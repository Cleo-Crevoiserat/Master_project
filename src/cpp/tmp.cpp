#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <exception>
#include <Eigen/Dense>
#include "gurobi_c++.h"
#include "Construction_Allocation.h"
#include "Construction_Part_1.h"

using namespace std;

// Trim whitespace from start and end of a string
string trim(const string& s) {
    size_t start = s.find_first_not_of(" \t\n\r");
    size_t end = s.find_last_not_of(" \t\n\r");
    return (start == string::npos) ? "" : s.substr(start, end - start + 1);
}

// Parse a vector<double> from "[x y z]" style string
vector<double> parseBracketVector(const string& input) {
    vector<double> result;
    size_t start = input.find('[');
    size_t end = input.find(']');
    if (start == string::npos || end == string::npos || end <= start)
        throw runtime_error("Invalid bracket vector format: " + input);

    string content = input.substr(start + 1, end - start - 1);
    stringstream ss(content);
    double val;
    while (ss >> val) {
        result.push_back(val);
    }
    return result;
}

// Parse a matrix from "[[1 2 3][4 5 6]]" style string
vector<vector<double>> parseBracketMatrix(const string& input) {
    vector<vector<double>> matrix;
    size_t pos = 0;
    while (true) {
        size_t start = input.find('[', pos);
        if (start == string::npos) break;
        size_t end = input.find(']', start);
        if (end == string::npos) throw runtime_error("Unmatched ']' in input");

        string rowStr = input.substr(start, end - start + 1);
        vector<double> row = parseBracketVector(rowStr);
        matrix.push_back(row);
        pos = end + 1;
    }
    return matrix;
}

int mainset(int argc, char* argv[]) {
    try {
        if (argc < 2) {
            cerr << "Usage:\n";
            cerr << "  " << argv[0] << " 0 m n \"[[W matrix]]\" nbQ \"[[Q matrix]]\"\n";
            cerr << "  " << argv[0] << " 1 nb_agent nb_object \"[N vector]\" \"[[Utility matrix]]\"\n";
            return 1;
        }

        GRBEnv env = GRBEnv(true);
        env.set(GRB_IntParam_InfUnbdInfo, 1);
        env.start();

        int prob = stoi(argv[1]);
        int argIndex = 2;

        if (prob == 0) {
            // Problem Wx <= b
            int m = stoi(argv[argIndex++]);
            int n = stoi(argv[argIndex++]);

            string wStr = argv[argIndex++];
            auto wVec = parseBracketMatrix(wStr);
            if ((int)wVec.size() != m)
                throw runtime_error("Row count mismatch in W matrix");
            if (!wVec.empty() && (int)wVec[0].size() != n)
                throw runtime_error("Column count mismatch in W matrix");

            Eigen::MatrixXd W(m, n);
            for (int i = 0; i < m; i++)
                for (int j = 0; j < n; j++)
                    W(i, j) = wVec[i][j];

            int nbQ = stoi(argv[argIndex++]);

            string qStr = argv[argIndex++];
            auto qVec = parseBracketMatrix(qStr);
            if ((int)qVec.size() != nbQ)
                throw runtime_error("Row count mismatch in Q matrix");
            for (const auto& row : qVec)
                if ((int)row.size() != m + 1)
                    throw runtime_error("Column count mismatch in Q matrix");

            vector<vector<double>> Q = qVec;

            Main(Q, W, n, m);
        }
        else if (prob == 1) {
            // Problem Fair Allocation
            int nb_agent = stoi(argv[argIndex++]);
            int nb_object = stoi(argv[argIndex++]);

            string nStr = argv[argIndex++];
            auto nVec = parseBracketVector(nStr);
            if ((int)nVec.size() != nb_object)
                throw runtime_error("Size mismatch in N vector");

            Eigen::VectorXd N(nb_object);
            for (int i = 0; i < nb_object; i++)
                N(i) = nVec[i];

            string utilStr = argv[argIndex++];
            auto utilVec = parseBracketMatrix(utilStr);
            if ((int)utilVec.size() != nb_agent)
                throw runtime_error("Agent count mismatch in Utility matrix");
            for (const auto& row : utilVec)
                if ((int)row.size() != nb_object)
                    throw runtime_error("Object count mismatch in Utility matrix rows");

            vector<vector<int>> Utility(nb_agent, vector<int>(nb_object));
            for (int i = 0; i < nb_agent; i++)
                for (int j = 0; j < nb_object; j++)
                    Utility[i][j] = static_cast<int>(utilVec[i][j]);

            Main_alloc(Utility, nb_object, nb_agent, N);
        }
        else {
            cerr << "Invalid problem type. Use 0 or 1.\n";
            return 1;
        }
    }
    catch (GRBException& e) {
        cerr << "Gurobi error: " << e.getMessage() << endl;
        return 1;
    }
    catch (exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    catch (...) {
        cerr << "Unknown error.\n";
        return 1;
    }

    return 0;
}

