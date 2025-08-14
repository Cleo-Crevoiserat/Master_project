#include <iostream>
#include <vector>
#include <algorithm>
#include <ciso646>
#include "gurobi_c++.h"
#include <chrono>
#include <cmath>
#include "Construction_Allocation.h"
#include "Useful_fct.h"
#include "Construction_Part_1.h"
#include "Solve_6.h"
using namespace std;

int maineeed() {
    try {
        GRBEnv env = GRBEnv(true);
        env.set(GRB_IntParam_InfUnbdInfo, 1);
        env.start();

        int prob;
        cout << "Chose the type of problem you want to solve (0 = Wx<=b, 1 = Fair Allocation) : ";
        cin >> prob;

        if (prob == 0) { // Case Wx <= b
            int m, n;
            cout << "Number of rows (m) : ";
            cin >> m;
            cout << "Number of columns (n) : ";
            cin >> n;

            Eigen::MatrixXd W(m, n);
            //cout << "Entrer les éléments de la matrice W (" << m << "x" << n << ") :\n";
            for (int i = 0; i < m; i++) {
                if (i == 0) {
                    cout << " Enter the coefficients of the first row of W : " << endl;
                }
                else if (i == 1) {
                    cout << " Enter the coefficients of the second row of W : " << endl;
                }
                else if (i == 2) {
                    cout << " Enter the coefficients of the third row of W : " << endl;
                }
                else {
                    cout << " Enter the coefficients of the " << i + 1 << "-th row of W : " << endl;
                }
                for (int j = 0; j < n; j++) {
                    cin >> W(i, j);
                }
            }

            int nbQ;
            cout << "Number of inequalities defining Q : ";
            cin >> nbQ;

            vector<vector<double>> Q(nbQ, vector<double>(n + m));
            //cout << "Enter les " << nbQ << " vecteurs Q (taille " <<  m+ 1 << " chacun) :\n";
            for (int i = 0; i < nbQ; i++) {
                if (i == 0) {
                    cout << " Enter the coefficients of the first inequality defining Q : " << endl;
                }
                else if (i == 1) {
                    cout << " Enter the coefficients of the second inequality defining Q : " << endl;
                }
                else if (i == 2) {
                    cout << " Enter the coefficients of the third inequality defining Q : " << endl;
                }
                else {
                    cout << " Enter the coefficients of the " << i + 1 << "-th inequality defining Q : "<<endl;
                }
                for (int j = 0; j < m+1; j++) {
                    cin >> Q[i][j];
                }
            }

            Main(Q, W, n, m);
        }

        else if (prob == 1) { // Cas Fair Allocation
            int nb_agent, nb_object;
            cout << "Number of agents : ";
            cin >> nb_agent;
            cout << "Number of type of objects : ";
            cin >> nb_object;

            Eigen::VectorXd N(nb_object);
            cout << "Number of items per type of objects: ";
            for (int i = 0; i < nb_object; i++) {
                cin >> N(i);
            }

            vector<vector<int>> Utility(nb_agent, vector<int>(nb_object));
            for (int i = 0; i < nb_agent; i++) {
                cout << "Enter the utility function of agent " << i+1<<" :\n";
                for (int j = 0; j < nb_object; j++) {
                    cin >> Utility[i][j];
                }
            }
            int nb_cell = 0;
            Main_alloc(Utility, nb_object, nb_agent, N);
        }

        else {
            cout << "Choix invalide.\n";
        }
    }
    catch (GRBException& e) {
        cout << "Erreur Gurobi : " << e.getMessage() << endl;
    }
    catch (...) {
        cout << "Erreur inconnue." << endl;
    }

    return 0;
}
