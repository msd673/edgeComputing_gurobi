#include "gurobi_c++.h"
#include "bits/stdc++.h"
#include "math.h"
#include <queue>
using namespace std;

// EUs ESs
int n, k; 

// Bandwidth between the terminal and the edge server
double r_nk_e;

// Bandwidth between the terminal and the cloud
double r_nk_c;

// query executability vector
vector<vector<int>> e; 

// the amount of computation
// the result size
vector<int> c, w; 

// the computational capability
vector<int> F;

// Branch-and-Bound node
struct node
{
     double upper;
     double lower;
     int Nd_num;
     std::vector<std::vector<int>> d;
     std::vector<std::vector<int>> d_upper;
};

void processNode(queue<node>& Q, node& p, double& min_upper, vector<vector<int>>& best_D);
double calcProblem(node &p);
double calcTargetVal(node &p);
vector<vector<int>> readMatrixFromFile(const string& filename, int rows, int cols);
vector<int> readVectorFromFile(const string& filename, int size);
void initializeParameters();

int main(int argc,
         char *argv[])
{
     initializeParameters();

    // Output the values for verification
    cout << "Number of MUs (n): " << n << ", Number of ECs (k): " << k << endl;
    
    cout << "Matrix e: " << endl;
    for (const auto& row : e) {
        for (const auto& val : row) {
            cout << val << " ";
        }
        cout << endl;
    }
    
    cout << "Vector c: ";
    for (const auto& val : c) {
        cout << val << " ";
    }
    cout << endl;
    
    cout << "Vector w: ";
    for (const auto& val : w) {
        cout << val << " ";
    }
    cout << endl;
    
    cout << "Vector F: ";
    for (const auto& val : F) {
        cout << val << " ";
    }
    cout << endl;

     queue<node> Q;
     node p;
     p.Nd_num = 0;
     p.d = {};
     p.d_upper.assign(n, vector<int>(k, 0));
     calcTargetVal(p);
     double min_upper = p.upper;
     vector<vector<int>> best_D;
     cout << "min_upper:" << min_upper << endl;
     Q.push(p);
     while (!Q.empty())
     {
          p = Q.front();
          Q.pop();
          processNode(Q, p, min_upper, best_D);
          
     }

     cout << "bestD: " << endl;
     for (int i = 0; i < n; i++)
     {
          for (int j = 0; j < k; j++)
          {
               cout << best_D[i][j] << " ";
          }
          cout << endl;
     }
     vector<vector<double>> f(n, vector<double>(k));
     vector<double> fm(k);
     for (int i = 0; i < k; i++)
     {
          for (int j = 0; j < n; j++)
          {
               fm[i] += best_D[j][i] * e[j][i] * sqrt(c[j]);
          }
     }
     for (int i = 0; i < n; i++)
     {
          for (int j = 0; j < k; j++)
          {
               if (best_D[i][j] == 0)
               {
                    f[i][j] = 0;
                    cout << f[i][j] << " ";
                    continue;
               }
               f[i][j] = F[j] * sqrt(c[i]) / fm[j];
               cout << f[i][j] << " ";
          }
          cout << endl;
     }
     cout << "bestTarget: " << min_upper << endl;

     return 0;
}

void processNode(queue<node>& Q, node& p, double& min_upper, vector<vector<int>>& best_D) {
    int index = p.Nd_num;
    if (index >= n)
        return;

    vector<vector<int>> determinedD = p.d;
    node p_cloud;
    p_cloud.Nd_num = index + 1;
    determinedD.push_back({0, 0});
    p_cloud.d = determinedD;
    calcProblem(p_cloud);
    calcTargetVal(p_cloud);

    if (min_upper > p_cloud.upper) {
        min_upper = p_cloud.upper;
        best_D = p_cloud.d_upper;
    }

    vector<node> p_edge;
    for (int i = 0; i < k; i++) {
        vector<int> T(k, 0);
        if (e[index][i] == 1) {
            T[i] = 1;
            node p_e;
            p_e.Nd_num = index + 1;
            vector<vector<int>> determinedD_edge = p.d;
            determinedD_edge.push_back(T);
            p_e.d = determinedD_edge;
            calcProblem(p_e);
            calcTargetVal(p_e);

            if (min_upper > p_e.upper) {
                min_upper = p_e.upper;
                best_D = p_e.d_upper;
            }
            p_edge.push_back(p_e);
        }
    }

    if (p_cloud.lower <= min_upper)
        Q.push(p_cloud);
    
    for (int i = 0; i < p_edge.size(); i++) {
        if (p_edge[i].lower <= min_upper)
            Q.push(p_edge[i]);
    }
}

double calcTargetVal(node &p)
{
     std::vector<std::vector<int>> D_upper = p.d_upper;
     double obj = 0.0;
     for (int i = 0; i < k; i++)
     {
          double tmp = 0.0;
          for (int j = 0; j < n; j++)
          {
               tmp += D_upper[j][i] * e[j][i] * sqrt(c[j]);
          }
          obj += tmp * tmp / F[i];
     }
     for (int i = 0; i < k; i++)
     {
          for (int j = 0; j < n; j++)
          {
               obj += D_upper[j][i] * e[j][i] * w[j] / r_nk_e;
          }
     }
     for (int i = 0; i < n; i++)
     {
          int is_edge = 0.0;
          for (int j = 0; j < k; j++)
          {
               is_edge += D_upper[i][j] * e[i][j];
          }
          if (is_edge == 0)
               obj += w[i] / r_nk_c;
     }
     p.upper = obj;
     return obj;
}

double calcProblem(node &p)
{
     int Nd_num = p.Nd_num;
     std::vector<std::vector<int>> d = p.d;
     GRBEnv *env = 0;

     // Create an environment
     env = new GRBEnv();
     GRBModel model = GRBModel(*env);

     // Create variables
     GRBVar **D = 0;
     int dcnt = 0;

     try
     {
          D = new GRBVar *[n - Nd_num];
          for (int i = 0; i < n - Nd_num; i++)
          {
               D[i] = new GRBVar[k];
               dcnt++;
          }
          for (int i = 0; i < n - Nd_num; i++)
          {
               for (int j = 0; j < k; j++)
               {
                    int tmp = i + Nd_num;
                    ostringstream vname;
                    vname << "D_" << tmp << j;
                    D[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, vname.str());
               }
          }

          double obj1 = 0.0;
          for (int i = 0; i < k; i++)
          {
               double tmp = 0.0;
               for (int j = 0; j < Nd_num; j++)
               {
                    tmp += d[j][i] * e[j][i] * sqrt(c[j]);
               }
               obj1 += tmp * tmp / F[i];
          }
          for (int i = 0; i < k; i++)
          {
               for (int j = 0; j < Nd_num; j++)
               {
                    obj1 += d[j][i] * e[j][i] * w[j] / r_nk_e;
               }
          }
          for (int i = 0; i < Nd_num; i++)
          {
               double is_edge = 0.0;
               for (int j = 0; j < k; j++)
               {
                    is_edge += d[i][j] * e[i][j];
               }
               if (is_edge == 0)
                    obj1 += w[i] / r_nk_c;
          }

          GRBQuadExpr obj2 = 0.0;
          for (int i = 0; i < k; i++)
          {
               GRBLinExpr tmp = 0.0;
               for (int j = 0; j < n - Nd_num; j++)
               {
                    tmp += D[j][i] * e[j + Nd_num][i] * sqrt(c[j]);
               }
               obj2 += tmp * tmp / F[i];
          }
          GRBLinExpr obj3 = 0.0;
          for (int i = 0; i < k; i++)
          {
               for (int j = 0; j < n - Nd_num; j++)
               {
                    obj3 += D[j][i] * e[j + Nd_num][i] * w[j + Nd_num] / r_nk_e;
               }
          }
          for (int i = 0; i < n - Nd_num; i++)
          {
               GRBLinExpr is_edge = 0.0;
               for (int j = 0; j < k; j++)
               {
                    is_edge += D[i][j] * e[i + Nd_num][j];
               }
               obj3 += (1 - is_edge) * w[i + Nd_num] / r_nk_c;
          }
          if (n == Nd_num)
          {
               cout << "target::::::" << obj1 + obj2 + obj3 << endl;
          }

          model.setObjective(obj1 + obj2 + obj3);

          GRBLinExpr constr1 = 0.0;
          for (int i = 0; i < n - Nd_num; i++)
          {
               GRBLinExpr constr = 0.0;
               ostringstream cname;
               cname << "c" << i;
               for (int j = 0; j < k; j++)
               {
                    constr += D[i][j] * e[i + Nd_num][j];
               }
               model.addConstr(constr <= 1, cname.str());
          }

          // Optimize model
          model.optimize();
          std::vector<std::vector<int>> D_upper = p.d;

          for (int i = 0; i < n - Nd_num; i++)
          {
               std::vector<int> tmp;
               for (int j = 0; j < k; j++)
               {
                    cout << D[i][j].get(GRB_StringAttr_VarName) << " "
                         << D[i][j].get(GRB_DoubleAttr_X) << endl;
                    tmp.push_back((int)round(D[i][j].get(GRB_DoubleAttr_X)));
               }
               D_upper.push_back(tmp);
          }

          for (int i = 0; i < n; i++)
          {
               for (int j = 0; j < k; j++)
               {
                    cout << D_upper[i][j] << " ";
               }
               cout << endl;
          }

          cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
          p.lower = model.get(GRB_DoubleAttr_ObjVal);
          p.d_upper = D_upper;
     }
     catch (GRBException e)
     {
          cout << "Error code = " << e.getErrorCode() << endl;
          cout << e.getMessage() << endl;
     }
     catch (...)
     {
          cout << "Exception during optimization" << endl;
     }
     for (int i = 0; i < dcnt; ++i)
     {
          delete[] D[i];
     }
     delete[] D;
     cout << "----------------------------------" << endl;
     delete env;
     return model.get(GRB_DoubleAttr_ObjVal);
}

vector<vector<int>> readMatrixFromFile(const string& filename, int rows, int cols) {
    vector<vector<int>> matrix(rows, vector<int>(cols));
    ifstream file(filename);
    string line;
    int row = 0;
    
    while (getline(file, line) && row < rows) {
        stringstream ss(line);
        for (int col = 0; col < cols; ++col) {
            ss >> matrix[row][col];
        }
        ++row;
    }
    
    return matrix;
}

vector<int> readVectorFromFile(const string& filename, int size) {
    vector<int> vec(size);
    ifstream file(filename);
    string line;
    int idx = 0;
    
    while (getline(file, line) && idx < size) {
        stringstream ss(line);
        ss >> vec[idx];
        ++idx;
    }
    
    return vec;
}

void initializeParameters() {
    // Get user input for 'n' (EUs) and 'k' (ESs)
    cout << "Enter the number of EUs (n): ";
    cin >> n;
    cout << "Enter the number of ESs (k): ";
    cin >> k;

    // Get user input for 'r_nk_e' and 'r_nk_c'
    cout << "Enter bandwidth between the terminal and the edge server: ";
    cin >> r_nk_e;
    cout << "Enter Bandwidth between the terminal and the cloud: ";
    cin >> r_nk_c;
    
    string matrixFile;  // File containing matrix 'e'
    cout << "Enter the filename of query executability vector: ";
    cin >> matrixFile;

    string vectorCFile = "vector_c.txt"; // File containing vector 'c'
    cout << "Enter the filename of the amount of computation: ";
    cin >> vectorCFile;

    string vectorWFile = "vector_w.txt"; // File containing vector 'w'
    cout << "Enter the filename of the result size: ";
    cin >> vectorWFile;

    string vectorFFile = "vector_f.txt"; // File containing vector 'F'
    cout << "Enter the filename of the computational capability: ";
    cin >> vectorFFile;
    
    // Read the matrix 'e' and vectors 'c', 'w', and 'F' from respective files
    e = readMatrixFromFile(matrixFile, n, k);
    c = readVectorFromFile(vectorCFile, n);
    w = readVectorFromFile(vectorWFile, n);
    F = readVectorFromFile(vectorFFile, k);
}