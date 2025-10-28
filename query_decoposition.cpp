#include <iostream>
#include <vector>
#include <string>
#include <regex>
#include <map>
#include <algorithm>
#include "gspan.h"

using namespace std;

vector<string> extractQueryPatterns(const string& sparql) {
    vector<string> patterns;
    map<string, string> varAliases; 
    vector<tuple<string, string, string>> triples;

    regex whereRegex(R"(WHERE\s*\{\s*([^}]+)\s*\})");
    smatch whereMatch;
    if (!regex_search(sparql, whereMatch, whereRegex) || whereMatch.size() < 2) {
        return patterns; 
    }
    string whereClause = whereMatch[1];

    regex tripleRegex(R"(([^.]+)\.)");
    sregex_iterator it(whereClause.begin(), whereClause.end(), tripleRegex);
    sregex_iterator end;

    for (; it != end; ++it) {
        string triple = it->str(1);
        triple.erase(remove_if(triple.begin(), triple.end(), ::isspace), triple.end());
        if (triple.empty()) continue;

        regex partsRegex(R"(([^<>"'\s]+)<([^>]+)>([^<>"'\s]+))");
        smatch partsMatch;
        if (regex_match(triple, partsMatch, partsRegex) && partsMatch.size() == 4) {
            string subj = partsMatch[1];
            string pred = partsMatch[2];
            string obj = partsMatch[3];
            triples.emplace_back(subj, pred, obj);
        }
    }

    vector<string> vars;
    for (const auto& t : triples) {
        if (get<0>(t).front() == '?') vars.push_back(get<0>(t));
        if (get<2>(t).front() == '?') vars.push_back(get<2>(t));
    }
    sort(vars.begin(), vars.end());
    vars.erase(unique(vars.begin(), vars.end()), vars.end());

    string aliasChars = "fpa";
    for (size_t i = 0; i < vars.size() && i < aliasChars.size(); ++i) {
        varAliases[vars[i]] = "?" + string(1, aliasChars[i]);
    }

    map<string, string> predMap = {
        {"rdfs:type", "type"},
        {"dbo:starring", "starring"}
    };

    for (const auto& t : triples) {
        string subj = get<0>(t);
        string pred = get<1>(t);
        string obj = get<2>(t);

        if (varAliases.count(subj)) subj = varAliases[subj];
        if (varAliases.count(obj)) obj = varAliases[obj];
        else if (obj.substr(0, 4) == "dbo:") { 
            obj = "?" + string(1, tolower(obj[4]));
        }

        string predShort = pred.substr(pred.find(':') + 1);
        if (predMap.count(pred)) predShort = predMap[pred];

        patterns.push_back(subj + "-" + predShort + "-" + obj);
    }

    return patterns;
}

int main() {
    string sparql;
    cin>>sparql;

    vector<string> patterns = extractQueryPatterns(sparql);
    
    cout << "Extracted query patterns:" << endl;
    for (const auto& p : patterns) {
        cout << p << endl;
    }

    return 0;
}
