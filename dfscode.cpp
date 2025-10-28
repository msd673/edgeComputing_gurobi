#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <functional>
#include <utility>
#include <set>
#include <regex>
#include <sstream>
#include <fstream> 
using namespace std;

struct Edge {
    int to;
    string elabel;
    Edge(int t, const string& e) : to(t), elabel(e) {}
};

struct Vertex {
    string label;
    vector<Edge> edges;
    Vertex(const string& l) : label(l) {} 
};

struct Graph {
    vector<Vertex> vertices;
    unordered_map<string, int> vertexMap;
    unordered_map<string, string> elementToParam;
    set<string> parameters;

    int addVertex(const string& label) {
        if (vertexMap.find(label) != vertexMap.end()) {
            return vertexMap[label]; 
        }
        int idx = vertices.size();
        vertices.emplace_back(label);
        vertexMap[label] = idx;
        return idx;
    }

    void addEdge(int from, int to, const string& predicate) {
        vertices[from].edges.emplace_back(to, predicate);
    }

    size_t size() const {
        return vertices.size();
    }

    int totalEdges() const {
        int count = 0;
        for (const auto& v : vertices) {
            count += v.edges.size();
        }
        return count;
    }

    bool isDirectedConnected() const {
        if (vertices.empty()) return true;
        vector<bool> visited(size(), false);
        function<void(int)> dfs = [&](int u) {
            visited[u] = true;
            for (const auto& e : vertices[u].edges) {
                if (!visited[e.to]) dfs(e.to);
            }
        };
        dfs(0);
        for (bool v : visited) if (!v) return false;
        return true;
    }

    void print() const {
        cout << "SPARQL convert graph structure:" << endl;
        for (size_t i = 0; i < vertices.size(); ++i) {
            cout << "vertex " << i << " (label: " << vertices[i].label << ") output: ";
            for (const auto& e : vertices[i].edges) {
                cout << i << "->" << e.to << "(" << e.elabel << ") ";
            }
            cout << endl;
        }

        if (!elementToParam.empty()) {
            cout << "\n element mapping to parameter:" << endl;
            for (const auto& [element, param] : elementToParam) {
                cout << element << " → " << param << endl;
            }
        }
    }
};

struct DFSEdge {
    int from;  
    int to;      
    string fromlabel;; 
    string elabel; 
    string tolabel;  
    DFSEdge(int f, int t, const string& fl, const string& el, const string tl)
        : from(f), to(t), fromlabel(fl), elabel(el), tolabel(tl) {}
        
    bool operator==(const DFSEdge& other) const {
        return from == other.from &&
               to == other.to &&
               fromlabel == other.fromlabel &&
               elabel == other.elabel &&
               tolabel == other.tolabel;
    }
};

struct DFSCode : public vector<DFSEdge> {
    void fromEdgeList(const Graph& g, const vector<pair<int, int>>& edge_list, set<int>& covered_vertices) {
        clear();
        covered_vertices.clear();
        for (const auto& [u, v] : edge_list) {
            string elabel;
            for (const auto& e : g.vertices[u].edges) {
                if (e.to == v) {
                    elabel = e.elabel;
                    break;
                }
            }
            push_back(DFSEdge(u, v, g.vertices[u].label, elabel, g.vertices[v].label));
            covered_vertices.insert(u);
            covered_vertices.insert(v);
        }
    }

    bool operator<(const DFSCode& other) const;
};

bool isLess(const DFSCode& a, const DFSCode& b) {
    int minSize = min(a.size(), b.size());
    for (int i = 0; i < minSize; ++i) {
        const DFSEdge& ea = a[i];
        const DFSEdge& eb = b[i];

        if (ea.fromlabel != eb.fromlabel) return ea.fromlabel < eb.fromlabel;
        if (ea.elabel != eb.elabel) return ea.elabel < eb.elabel;
        if (ea.tolabel != eb.tolabel) return ea.tolabel < eb.tolabel;
        if (ea.from != eb.from) return ea.from < eb.from;
        if (ea.to != eb.to) return ea.to < eb.to;
    }
    return a.size() < b.size();
}

DFSCode generateCodeFromRoot(const Graph& g, int root_from, int root_to, set<int>& covered_vertices) {
    vector<pair<int, int>> edge_list;
    vector<bool> visited(g.size(), false);

    function<void(int)> dfs = [&](int current) {
        visited[current] = true;
        for (const auto& e : g.vertices[current].edges) {
            int next = e.to;
            edge_list.emplace_back(current, next);
            if (!visited[next]) {
                dfs(next);
            }
        }
    };

    edge_list.emplace_back(root_from, root_to);
    visited[root_from] = true; 
    dfs(root_from); 

    DFSCode code;
    code.fromEdgeList(g, edge_list, covered_vertices);
    return code;
}

bool DFSCode::operator<(const DFSCode& other) const {
    return isLess(*this, other);
}

void printDFSCode(const DFSCode& code, const Graph& g) {
    if (code.empty()) {
        cout << "coding is null";
        return;
    }
    set<int> covered;
    for (const auto& e : code) {
        covered.insert(e.from);
        covered.insert(e.to);
    }
    for (size_t i = 0; i < code.size(); ++i) {
        const DFSEdge& e = code[i];
        if (i > 0) cout << " -> ";
        cout << "(" << e.from << "," << e.to << "," << e.fromlabel 
                  << ",\"" << e.elabel << "\"," << e.tolabel << ")";
    }
    cout << " [edge:" << code.size() << "/" << g.totalEdges() 
              << ", vertex" << covered.size() << "/" << g.size() << "]";
    if (covered.size() == g.size() && code.size() >= g.totalEdges()) {
        cout << "[valid]";
    }
    cout << endl;
}

pair<DFSCode, vector<DFSCode>> findMinAndPrunedDFSCode(const Graph& g) {
    vector<DFSCode> valid_codes; 
    DFSCode min_code;
    bool first_valid = true;
    int code_count = 0;  

    const int total_edges = g.totalEdges(); 
    if (total_edges == 0) {
        DFSCode single_node_code;
        single_node_code.push_back(DFSEdge(0, 0, g.vertices[0].label, "", g.vertices[0].label));
        valid_codes.push_back(single_node_code);
        cout << "generate coding:" << endl;
        cout << "coding 1:";
        printDFSCode(single_node_code, g);
        return {single_node_code, valid_codes};
    }

    if (!g.isDirectedConnected()) {
        cerr << "[Warning: the input graph is not a directed connected graph and may not have a valid code!]" << endl;
    }

    cout << "All generated codes (with pruning information):" << endl;

    for (int u = 0; u < (int)g.size(); ++u) {
        for (const auto& e : g.vertices[u].edges) {
            int v = e.to;
            set<int> covered_vertices;
            code_count++;  
            
            DFSCode current_code = generateCodeFromRoot(g, u, v, covered_vertices);
            
            cout << "coding " << code_count << ":";
            printDFSCode(current_code, g);
            
            if (current_code.size() < total_edges) {
                cout << " Reason for pruning: insufficient number of edges (" << current_code.size() << "<" << total_edges << ")" << endl;
                continue;
            }
            
            if (covered_vertices.size() != g.size()) {
                cout << " Reason for pruning: vertex not fully covered (" << covered_vertices.size() << "<" << g.size() << ")" << endl;
                continue;
            }

            valid_codes.push_back(current_code);
            cout << " Reserved as a valid code" << endl;
            
            if (first_valid || isLess(current_code, min_code)) {
                min_code = current_code;
                first_valid = false;
            }
        }
    }

    return {min_code, valid_codes};
}

string trim(const string& s) {
    auto start = s.begin();
    while (start != s.end() && isspace(*start)) {
        start++;
    }
    auto end = s.end();
    do {
        end--;
    } while (distance(start, end) > 0 && isspace(*end));
    return string(start, end + 1);
}

Graph sparqlToGraph(const string& sparqlQuery) {
    Graph graph;
    int paramCounter = 0;

    regex triplePattern(R"(\{([^}]+)\})");
    smatch match;
    string triplesStr;

    if (regex_search(sparqlQuery, match, triplePattern) && match.size() > 1) {
        triplesStr = match[1].str();
    } else {
        cerr << " No valid ternary patterns were found! " << endl;
        return graph;
    }

    vector<string> triples;
    stringstream ss(triplesStr);
    string triple;

    string buffer;
    char c;
    bool inQuotes = false;
    int braceCount = 0; 

    while (ss.get(c)) {
        if (c == '"') {
            inQuotes = !inQuotes;
        }
    
        if (c == '<') braceCount++;
        else if (c == '>') braceCount--;
    
        if (!inQuotes && braceCount == 0) {
            if (c == '.' || c == ';') {
                trim(buffer);
                if (!buffer.empty()) {
                    triples.push_back(buffer);
                }
                buffer.clear();
                continue;
            }
        }
    
        buffer += c;
    }
    trim(buffer);
    if (!buffer.empty()) {
        triples.push_back(buffer);
    }

    struct Triple {
        string s, p, o;
    };
    vector<Triple> parsedTriples;

    for (const auto& t : triples) {
        stringstream ss3(t);
        vector<string> parts;
        string part;
        while (ss3 >> part) {
            parts.push_back(part);
        }

        if (parts.size() != 3) {
            cerr << "Invalid ternary format:" << t << endl;
            continue;
        }

        parsedTriples.push_back({parts[0], parts[1], parts[2]});
    }

    sort(parsedTriples.begin(), parsedTriples.end(), [](const Triple& a, const Triple& b) {
        return a.p < b.p;
    });

    for (const auto& t : parsedTriples) {
        string s = t.s;
        string p = t.p;
        string o = t.o;

        if (graph.elementToParam.find(s) == graph.elementToParam.end()) {
            char letter = 'a' + (paramCounter % 26);
            string param = "?" + string(1, letter);
            if (paramCounter >= 26) {
                param = "?" + string(1, 'a' + (paramCounter / 26 - 1)) + string(1, letter);
            }
            graph.elementToParam[s] = param;
            graph.parameters.insert(param);
            paramCounter++;
        }
        s = graph.elementToParam[s];

        if (graph.elementToParam.find(o) == graph.elementToParam.end()) {
            char letter = 'a' + (paramCounter % 26);
            string param = "?" + string(1, letter);
            if (paramCounter >= 26) {
                param = "?" + string(1, 'a' + (paramCounter / 26 - 1)) + string(1, letter);
            }
            graph.elementToParam[o] = param;
            graph.parameters.insert(param);
            paramCounter++;
        }
        o = graph.elementToParam[o];

        int sIdx = graph.addVertex(s);
        int oIdx = graph.addVertex(o);
        graph.addEdge(sIdx, oIdx, p);
    }

    return graph;
}

vector<pair<string, string>> load_sparql_queries(const string& filename) {
    ifstream infile(filename);
    vector<pair<string, string>> queries;
    if (!infile.is_open()) {
        cerr << "Error: cannot open file " << filename << endl;
        return queries;
    }

    string line, title, query;
    ostringstream query_buffer;

    while (getline(infile, line)) {
        if (line.rfind("#", 0) == 0) { 
            if (!title.empty() && !query_buffer.str().empty()) {
                queries.emplace_back(title, query_buffer.str());
                query_buffer.str("");
                query_buffer.clear();
            }
            title = line.substr(1);
            if (!title.empty() && title[0] == ' ') title.erase(0, 1);
        } else {
            query_buffer << line << "\n";
        }
    }

    if (!title.empty() && !query_buffer.str().empty()) {
        queries.emplace_back(title, query_buffer.str());
    }

    return queries;
}


int main() {
    vector<pair<string, string>> sparql_queries = load_sparql_queries("queries.txt");

    map<DFSCode,string> all_min_dfscodes;

    for (const auto& [query_name, sparql] : sparql_queries) {
        cout << "\n====================================================" << endl;
        cout << "Processing quiries: " << query_name << endl;
        cout << "====================================================" << endl;

        Graph g = sparqlToGraph(sparql);

        cout << "--- graphical information ---" << endl;
        cout << "Vertex count:" << g.size() << ", the total number of sides:" << g.totalEdges() << endl;

        auto [min_code, valid_codes] = findMinAndPrunedDFSCode(g);
        
        all_min_dfscodes[min_code] = sparql;

        cout << "\n--- encoding results ---" << endl;
        cout << "Minimum valid DFS code:" << endl;
        printDFSCode(min_code, g);
    }

    cout << "\n====================================================" << endl;
    cout << "Minimum DFS code summary and Map take out test for all queries:" << endl;
    
    const string& test_sparql_key = sparql_queries[0].second;
    Graph g1 = sparqlToGraph(test_sparql_key);
    auto[test_key_code,code] = findMinAndPrunedDFSCode(g1); 
    if (all_min_dfscodes.count(test_key_code)) {
        cout << "The corresponding query is:" << all_min_dfscodes[test_key_code] << endl;
    } else{
        cout << "error!" << endl;
    }
    // ----------------------------------------------------
    cout << "====================================================" << endl;


    return 0;
}
