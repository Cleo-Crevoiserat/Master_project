#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <Eigen/Dense>
#include "gurobi_c++.h"
#include "Construction_Allocation.h"
#include "Construction_Part_1.h"

using namespace std;

int main(int argc, char* argv[]) {
    try {
        if (argc < 2) {
            cerr << "Usage:\n";
            cerr << "  " << argv[0] << " 0 m n W_values... nbQ Q_values...\n";
            cerr << "  " << argv[0] << " 1 nb_agent nb_object N_values... Utility_values...\n";
            return 1;
        }

        GRBEnv env = GRBEnv(true);
        env.set(GRB_IntParam_InfUnbdInfo, 1);
        env.start();

        int prob = stoi(argv[1]);
        int argIndex = 2; // where we are in argv

        if (prob == 0) {
            // Wx <= b case
            int m = stoi(argv[argIndex++]);
            int n = stoi(argv[argIndex++]);

            Eigen::MatrixXd W(m, n);
            for (int i = 0; i < m; i++)
                for (int j = 0; j < n; j++)
                    W(i, j) = stod(argv[argIndex++]);

            int nbQ = stoi(argv[argIndex++]);
            vector<vector<double>> Q(nbQ, vector<double>(m + 1));
            for (int i = 0; i < nbQ; i++)
                for (int j = 0; j < m + 1; j++)
                    Q[i][j] = stod(argv[argIndex++]);

            Main(Q, W, n, m);
        }
        else if (prob == 1) {
            // Fair allocation case
            int nb_agent = stoi(argv[argIndex++]);
            int nb_object = stoi(argv[argIndex++]);

            Eigen::VectorXd N(nb_object);
            for (int i = 0; i < nb_object; i++)
                N(i) = stod(argv[argIndex++]);

            vector<vector<int>> Utility(nb_agent, vector<int>(nb_object));
            for (int i = 0; i < nb_agent; i++)
                for (int j = 0; j < nb_object; j++)
                    Utility[i][j] = stoi(argv[argIndex++]);

            Main_alloc(Utility, nb_object, nb_agent, N);
        }
        else {
            cerr << "Invalid problem type.\n";
            return 1;
        }
    }
    catch (GRBException& e) {
        cerr << "Gurobi error: " << e.getMessage() << endl;
    }
    catch (exception& e) {
        cerr << "Error: " << e.what() << endl;
    }
    catch (...) {
        cerr << "Unknown error.\n";
    }

    return 0;
}