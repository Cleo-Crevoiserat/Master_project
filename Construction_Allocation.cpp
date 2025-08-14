#include <iostream>
#include <vector>
#include <ciso646>
#include <algorithm>
#include "gurobi_c++.h"
#include <chrono>
#include <fstream>
#include <cmath>
#include <Eigen/Dense>
#include "Usefull_fct.h"
#include "Solve_6.h"
#include "Construction_Part_1.h"
#include "Solve_6_4.h"
using namespace std;

vector<vector<double>> Q_create(const int& n, const int& m, Eigen::VectorXd& N, vector<vector<int>>& Utility) {
    // create the Polyhedron Q for envy free alocation as explained in the paper
    vector<vector<double>> Q;
    for (int i = 0; i < n; ++i) {
        vector<double> q_1(n + m * n + 2, 0);
        q_1[0] = 0;
        q_1[1 + i] = 1;
        vector<double> q_2(n + m * n + 2, 0);
        q_2[0] = 0;
        q_2[1 + i] = -1;
        for (int j = 0; j < m; ++j) {
            q_1[n + 1 + 1 + n * j + i] = Utility[i][j];
            q_2[n + 1 + 1 + n * j + i] = -Utility[i][j];
        }
        Q.push_back(q_1);
        Q.push_back(q_2);
    }
    vector<double> q_1(n + m * n + 2, 0);
    q_1[n + 1] = 1;
    vector<double> q_2(n + m * n + 2, 0);
    q_2[n + 1] = -1;
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            q_1[n + 2 + n * i + j] = Utility[j][i];
            q_2[n + 2 + n * i + j] = -Utility[j][i];
        }
    }
    q_1[0] = 1;
    q_2[0] = -1;
    Q.push_back(q_1);
    Q.push_back(q_2);
    for (int j = 0; j < m; ++j) {
        for (int i = 0; i < n; ++i) {
            vector<double> q_1(n + m * n + 2, 0);
            q_1[n + 2 + n * j + i] = 1;
            Q.push_back(q_1);
        }
    }
    for (int j = 0; j < m; ++j) {
        vector<double> q_1(n + m * n + 2, 0);
        q_1[0] = N[j];
        for (int i = 0; i < n; ++i) {
            q_1[n + 2 + n * j + i] = -1;
        }
        Q.push_back(q_1);
    }
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            vector<double> q_1(n + m * n + 2, 0);
            if (j != i) {
                for (int k = 0; k < m; ++k) {
                    q_1[n + 2 + k * n + i] = Utility[i][k];
                    q_1[n + 2 + k * n + j] = -Utility[i][k];
                }
                Q.push_back(q_1);
            }
        }
    }
    return Q;
}
vector<vector<double>> Q_create_case_2(const int& n, const int& m, Eigen::VectorXd& N, vector<vector<int>> Utility) {
    vector<vector<double>> Q;
    for (int i = 0; i < m; ++i) {
        vector<double> q_1(m + n + m * n + 2, 0);
        q_1[0] = N[i];
        q_1[i + 1] = -1;
        vector<double> q_2(m + n + m * n + 2, 0);
        q_2[0] = -N[i];
        q_2[i + 1] = 1;
        Q.push_back(q_1);
        Q.push_back(q_2);
    }
    //cout << "step 1 done" << endl;
    for (int i = 0; i < n; ++i) {
        vector<double> q_1(m + n + m * n + 2, 0);
        q_1[0] = 0;
        q_1[m + 1 + i] = 1;
        vector<double> q_2(m + n + m * n + 2, 0);
        q_2[0] = 0;
        q_2[m + 1 + i] = -1;
        for (int j = 0; j < m; ++j) {
            q_1[m + n + 1 + 1 + n * j + i] = Utility[i][j];
            q_2[m + n + 1 + 1 + n * j + i] = -Utility[i][j];
        }
        Q.push_back(q_1);
        Q.push_back(q_2);
    }
    //cout << "step 2 done" << endl;
    vector<double> q_1(m + n + m * n + 2, 0);
    q_1[m + n + 1] = 1;
    vector<double> q_2(m + n + m * n + 2, 0);
    q_2[m + n + 1] = -1;
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            q_1[m + n + 2 + n * i + j] = Utility[j][i];
            q_2[m + n + 2 + n * i + j] = -Utility[j][i];
        }
    }
    q_1[0] = 1;
    q_2[0] = -1;
    Q.push_back(q_1);
    Q.push_back(q_2);
    //cout << "step 3 done" << endl;
    for (int j = 0; j < m; ++j) {
        for (int i = 0; i < n; ++i) {
            vector<double> q_1(m + n + m * n + 2, 0);
            q_1[m + n + 2 + n * j + i] = 1;
            Q.push_back(q_1);
        }
    }
    //cout << "step 4 done" << endl;
    for (int j = 0; j < m; ++j) {
        vector<double> q_1(m + n + m * n + 2, 0);
        q_1[0] = N[j];
        for (int i = 0; i < n; ++i) {
            q_1[m + n + 2 + n * j + i] = -1;
        }
        Q.push_back(q_1);
    }
    //cout << "step 5 done" << endl;
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            vector<double> q_1(m + n + m * n + 2, 0);
            if (j != i) {
                for (int k = 0; k < m; ++k) {
                    q_1[m + n + 2 + k * n + i] = Utility[i][k];
                    q_1[m + n + 2 + k * n + j] = -Utility[i][k];
                }
                Q.push_back(q_1);
            }
        }
    }
    //cout << "step 6 done" << endl;
    return Q;
}
Eigen::MatrixXd Allocation_Matrix(const int& n, const int& m, Eigen::VectorXd& N, vector<vector<int>>& Utility) {
    Eigen::MatrixXd W = Eigen::MatrixXd::Zero(m + n + 1 + n * m, n * m);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            W(i, j + i * n) = 1;
        }
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            W(m + i, j * n + i) = -Utility[i][j];
        }
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            W(m + n, j * n + i) = -Utility[i][j];
        }
    }
    //cout << W;
    return W;
}
vector<Eigen::VectorXd> gen_C_Alloc(Eigen::MatrixXd& W, Eigen::VectorXd& N, int m, int n, int n_1, int m_1) {
    vector<Eigen::VectorXd> C;
    vector<Eigen::VectorXd> Coord;
    for (int i = 0; i < m_1; ++i) {
        Coord.push_back(Eigen::VectorXd::Ones(N(i)));
    }
    do {
        int pos = 0;
        Eigen::VectorXd coord = Eigen::VectorXd::Zero(n);
        for (int i = 0; i < m_1; ++i) {// Coord[i][k]=b means that the kth item of type  i is received by agent b
            for (int j = 0; j < N(i); ++j) {
                if (Coord[i][j] != 0) {
                    coord[pos + Coord[i][j]-1] += 1;
                }
            }
            pos += n_1;
        }
        Eigen::VectorXd c;
        c = W * coord;
        Eigen::VectorXd c_cut(m-m_1);
        for (int i =0; i < m-m_1; ++i) {
            c_cut(i)=c(i+m_1);
        }
        C.push_back(c_cut);
    } while (!(coord_Spe_C_Alloc(Coord, m_1 - 1, N[m_1 - 1] - 1, N, m_1, n_1)));
    return C;

}
Eigen::MatrixXd Allocation_expansion(Eigen::MatrixXd& W, int& m, int& n) {
    Eigen::MatrixXd test(m, n + m);
    for (size_t k = 0; k < n; ++k) {
        test.col(k) = W.col(k);
    };
    Eigen::MatrixXd E(m, m);
    E = Eigen::MatrixXd::Identity(m, m);
    for (size_t k = 0; k < m; ++k) {
        test.col(n + k) = E.col(k);
    }
    return test;
}

vector<Eigen::VectorXd> gen_C_Alloc_1(Eigen::MatrixXd& W, int bound, int m, int n,int mn, GRBEnv& env) {
    vector<Eigen::VectorXd> C_1;
    Eigen::VectorXd coord(m);
    int diff = m - mn;
    for (size_t i = 0; i < diff; ++i) {// the last drections will always be 0;
            coord[i] = -bound;
    }
    for (size_t i = 0; i < mn; ++i) {
        coord[i + diff] = 0;
    }
    do {
        vector<int> c;
        GRBModel model = GRBModel(env);
        model.set(GRB_IntParam_OutputFlag, 0);
        std::vector<GRBVar> vars(n-mn);
        for (size_t i = 0; i < n-mn; ++i) {
            vars[i] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_INTEGER, "x_" + std::to_string(i)); // Integer variables
        }
        int nb = 0;
        for (size_t i = 0; i < m-mn; ++i) {
            GRBLinExpr expr = 0;
            for (size_t j = 0; j < n-mn; ++j) {
                expr += W(i, j) * vars[j];
            }
            model.addConstr(expr == coord[i]);
        }
        model.optimize();
        int status0 = model.get(GRB_IntAttr_Status);
        if (model.get(GRB_IntAttr_SolCount) > 0) {
            C_1.push_back(coord);
        }
    } while (coord_base(coord, 0, bound, m-mn) == false);
    return C_1;
}

double Main_alloc(vector< vector<int>>& Utility, int m, int n, Eigen::VectorXd& N) {
    int nb_cell = 1;
    GRBEnv env = GRBEnv(true);
    env.set(GRB_IntParam_InfUnbdInfo, 1);
    env.start();
    auto start = std::chrono::high_resolution_clock::now();
    vector<vector<double>> Q;
    Q= Q_create(n, m, N, Utility);
    //Q = Q_create_case_2(n, m, N, Utility);
    //for (auto q : Q) {
    //    for (auto q_1 : q) {
    //        cout << q_1 << " ";
    //    }
    //    cout << endl;
    //}
    Eigen::MatrixXd W = Allocation_Matrix(n, m, N, Utility);
    int m_1 = m + n + 1 + n * m;
    int n_1 = n * m;
    Eigen::MatrixXd W2 = Allocation_expansion(W, m_1, n_1);
    //long long int bound = max_det_sub_mat_alloc(n, m, W2);
    long long int Delta = W.cwiseAbs().maxCoeff();// the max abs value of W
    long long int inf_norm = (W.cwiseAbs().rowwise().sum()).maxCoeff(); //the infinity norm of W
    long long int hadamard = sqrt(pow(2 * Delta, 2) + 1);// the hadamard inequality for Delta_m in order to be faster
    //cout << hadamard;
    long long int  bound_2 = inf_norm * m * pow(2 * m * Delta + 1, m);// previous bound
    for (int i = 0; i < n + m + 1; ++i) {

    }
    hadamard = pow(hadamard, n + m + 1);
    //cout << hadamard;
    hadamard = hadamard * inf_norm *(n_1);// *(n-m)=(n_1+n*m+1-(n*m+1))
    //cout << hadamard;
    cout << true << " " << (hadamard > 0);
    int bound=0;
    if ((hadamard < 0) and (bound_2 < 0)) {// the number is to big to be used
        bound = -1;
    }
    else if ((hadamard < 0)) {
        bound = bound_2;
    }
    else if ((bound_2 < 0)) {
        bound = hadamard;
    }
    else {
        if (hadamard < bound_2) {
            bound = hadamard;
        }
        else {
            bound = bound_2;
        }
    }
    int check = 1;
    int facto = factorial(n-1);
    for (int i = 0; i < m; ++i) {
        check = (check * factorial(n-1 + N(i)) / factorial(N(i))) / facto;// the number of point that have to be generated
    }
    cout << "finito ";
    cout << bound;
    if ((true) or (bound<0) or (check<=bound) or(check<=pow(bound,m+n+1))) {// we will have to do less operations
        cout << "on rentre";
        Q = Q_create(n, m, N, Utility);
        vector<Eigen::VectorXd> C;
        C = gen_C_Alloc(W, N, m_1, n_1, n, m);
        m_1 = m_1 - m;
        vector<int> coord(m_1, 0);
        std::vector<int> result(m_1, 0);
        GRBModel model = GRBModel(env);
        model.set(GRB_IntParam_OutputFlag, 0);
        std::vector<GRBVar> vars(m_1);
        for (size_t j = 0; j < m_1; ++j) {
            vars[j] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_INTEGER, "x_" + std::to_string(j));
        }
        model.setObjective(-vars[0], GRB_MAXIMIZE);
        int nb = 0;
        for (auto const& q : Q) {
            GRBLinExpr expr = q[0];
            for (size_t k = 0; k < m_1; ++k) {
                expr += (q[k + 1] * vars[k]);
            }
            model.addConstr(expr >= 0, "Q_constraint_" + std::to_string(nb));// here we add the constraint that defines Q
            nb += 1;
        }
        GRBLinExpr expr = 0;
        for (size_t k = n + 1; k < m_1; ++k) {
            expr += (vars[k] * 1);
        }
        model.addConstr(expr >= 1, "Q_constraint_" + std::to_string(nb));// we check if it exists at least one allocation which is ! [0,0,...,0] 
        model.write("modele_test.lp");
        model.optimize();
        int status1 = model.get(GRB_IntAttr_Status);
        if (status1 == GRB_OPTIMAL) {
            vector<vector<int>> fin = Solve6(C, Q, m_1, coord, env, result,nb_cell);
            cout << endl;
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::ratio<1>> duration = end - start;
            std::cout << endl << "Algorithm took " << duration.count() << " seconds.\n" << std::endl;
            if (fin[0][0] == 1) {// it means we founded a counter example
                vector<int> soluce(m_1 + m);
                for (int z = 0; z < m; ++z) {
                    soluce[z] = N(z);
                }
                for (int z = 0; z < m_1; ++z) {
                    soluce[z + m] = result[z];
                }
                cout << endl;
                GRBModel model1 = GRBModel(env);
                model1.set(GRB_IntParam_OutputFlag, 0);
                std::vector<GRBVar> vars1(n_1);
                for (size_t j = 0; j < n_1; ++j) {
                    vars1[j] = model1.addVar(0, GRB_INFINITY, 0.0, GRB_INTEGER, "x_" + std::to_string(j));
                }
                for (size_t j = 0; j < m_1 + m; ++j) {
                    GRBLinExpr expr = 0;
                    for (size_t k = 0; k < n_1; ++k) {
                        expr += W(j, k) * vars1[k];
                    }
                    model1.addConstr(expr <= soluce[j]);// here we add the constraint that makes Wx<=b
                    nb += 1;
                }
                model1.optimize();
                int status2 = model1.get(GRB_IntAttr_Status);
                if (status2 == GRB_INFEASIBLE) {

                    cout << " The allocation ";
                    for (int i = n+1; i < m_1;++i) {
                        cout << result[i] << " ";
                    }
                    cout << "is indeed fair";
                }
            }
            //return nb_cell;
            return duration.count();
        }
        else {
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::ratio<1>> duration = end - start;
            std::cout << "no solution and Algorithm took " << duration.count() << " seconds.\n" << std::endl;
            //return -12;
            return duration.count();
        }
    }
    else {
        Q = Q_create_case_2(n, m, N, Utility);
        vector<double> Time;
        auto start = std::chrono::high_resolution_clock::now();
        vector<int> result;
        vector<Eigen::VectorXd> C;
        n_1 = n_1 + m * n + 1+n+m;
        C = gen_C_Alloc_1(W2,bound, m_1, n_1, n * m, env);
        result = cone_W(n_1, m_1, env, Q, W2);
        if (result != vector<int>(0, m_1)) {
            cout << endl << "The counter-example is:" << endl;
            for (auto r : result) {
                cout << r << " ";
            }
            cout << endl;
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::ratio<1>> duration = end - start;
            std::cout << endl << "Algorithm took " << duration.count() << " seconds.\n" << std::endl;
            return duration.count();

        }
        Square_Matrices(W2, Q, C, m_1, n_1, env);
        cout << endl;
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::ratio<1>> duration = end - start;
        std::cout << endl << "Algorithm took " << duration.count() << " seconds.\n" << std::endl;
        return duration.count();

    }
}
