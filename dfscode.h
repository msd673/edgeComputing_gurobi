#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <climits>

struct Edge {
	int from;
	int to;
	int elabel;
	unsigned int id;
	Edge(): from(0), to(0), elabel(0), id(0) {};
};

struct Vertex {
    int label;
    std::vector<Edge> edges;
    Vertex(int l) : label(l) {}
};

struct Graph {
    std::vector<Vertex> vertices;
    bool directed;
    Graph(bool d) : directed(d) {}
    void addVertex(int label) {
        vertices.emplace_back(label);
    }
    void addEdge(int from, int to, int elabel) {
        vertices[from].edges.emplace_back(to, elabel);
        if (!directed) {
            vertices[to].edges.emplace_back(from, elabel);
        }
    }
};

struct DFSEdge {
    int from;
    int to;
    int fromlabel;
    int elabel;
    int tolabel;
    DFSEdge(int f, int t, int fl, int el, int tl)
        : from(f), to(t), fromlabel(fl), tolabel(tl), elabel(el) {}
};

struct DFSCode : public std::vector<DFSEdge> {
    void fromGraph(Graph &g, int start) {
        clear();
        std::vector<bool> visited(g.vertices.size(), false);
        std::function<void(int, int)> dfs = [&](int u, int parent) {
            visited[u] = true;
            for (const Edge &e : g.vertices[u].edges) {
                int v = e.to;
                if (v == parent && g.directed) continue;
                if (!g.directed && v < u && visited[v]) continue;
                push_back(DFSEdge(u, v, g.vertices[u].label, e.elabel, g.vertices[v].label));
                if (!visited[v]) {
                    dfs(v, u);
                }
            }
        };
        dfs(start, -1);
    }
};

bool isLess(const DFSCode &a, const DFSCode &b) {
    int minSize = std::min(a.size(), b.size());
    for (int i = 0; i < minSize; ++i) {
        const DFSEdge &ea = a[i];
        const DFSEdge &eb = b[i];

        if (ea.from != eb.from) return ea.from < eb.from;
        if (ea.to != eb.to) return ea.to < eb.to;
        if (ea.fromlabel != eb.fromlabel) return ea.fromlabel < eb.fromlabel;
        if (ea.elabel != eb.elabel) return ea.elabel < eb.elabel;
        if (ea.tolabel != eb.tolabel) return ea.tolabel < eb.tolabel;
    }
    return a.size() < b.size();
}

DFSCode getMinDFSCode(Graph &g) {
    DFSCode minCode;
    bool first = true;

    for (size_t i = 0; i < g.vertices.size(); ++i) {
        DFSCode code;
        code.fromGraph(g, i);

        if (first || isLess(code, minCode)) {
            minCode = code;
            first = false;
        }
    }

    return minCode;
}

void printDFSCode(const DFSCode &code) {
    for (const DFSEdge &e : code) {
        std::cout << "(" << e.from << "," << e.to << "," << e.fromlabel << "," 
                  << e.elabel << "," << e.tolabel << ") ";
    }
    std::cout << std::endl;
}

int main() {
    Graph g(false);
    g.addVertex(1);
    g.addVertex(2);
    g.addVertex(3);
    g.addVertex(15);
    g.addVertex(17);
    g.addEdge(0, 1, 5);
    g.addEdge(1, 2, 6);
    g.addEdge(0, 2, 7);

    DFSCode minCode = getMinDFSCode(g);
    std::cout << "minimal DFS coding: ";
    printDFSCode(minCode);

    return 0;
}