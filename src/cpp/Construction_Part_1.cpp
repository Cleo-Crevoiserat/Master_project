#include <iostream>
#include <vector>
#include <algorithm>
#include "gurobi_c++.h"
#include <chrono>
#include <ciso646>
#include <fstream>
#include <cmath>
#include <Eigen/Dense>
#include "Construction_Allocation.h"
#include "Useful_fct.h"
#include "Solve_6.h"
using namespace std;

Eigen::RowVectorXd Hyperplan_def(size_t const m, size_t const n, Eigen::MatrixXd const& W, vector<int> coord) {
    Eigen::MatrixXd W3(m - 1, m - 1);
    bool negative = false;
    Eigen::MatrixXd W2(m, m - 1);
    for (size_t k = 0; k < m - 1; ++k) {
        W2.col(k) = W.col(coord[k]);
    }
    Eigen::RowVectorXd rows(m + 1);
    rows(0) = 0;
    int test_pos = 0;// This is to check if perhaps the hyperplane we find is useless because there is elements of cone(W) on both sides
    int test_neg = 0;
    Eigen::VectorXd null_vect = Eigen::VectorXd::Zero(m + 1);
    for (size_t i = 0; i < m; ++i) {
        int val = 0;
        if (i == 0) {
            W3 << W2.block(1, 0, m - 1, m - 1);
            val = W3.determinant();
            rows(0, i + 1) = val;
            negative = true;
        }
        else if (i == m - 1) {
            if (negative) {
                W3 << W2.block(0, 0, m - 1, m - 1);
                val = -W3.determinant();
                rows(0, i + 1) = val;
            }
            else {
                W3 << W2.block(0, 0, m - 1, m - 1);
                val = W3.determinant();
                rows(0, i + 1) = val;
            }
        }
        else {
            if (negative) {
                W3 << W2.block(0, 0, i, m - 1), W2.block(i + 1, 0, m - 1 - i, m - 1);
                val = -W3.determinant();
                rows(0, i + 1) = val;
                negative = false;
            }
            else {
                W3 << W2.block(0, 0, i, m - 1), W2.block(i + 1, 0, m - 1 - i, m - 1);// 3 à remplacer par n-1 et 6 par m-1
                val = W3.determinant();
                rows(0, i + 1) = val;
                negative = true;
            }
        }
    }
    for (size_t j = 0; j < n; ++j) {// Now we check if the points e_i*W are all on the same side of the hyperplane 
        bool test_j = true;
        int tot = 0;
        for (size_t k = 0; k < m - 1; ++k) {
            if (j == coord[k]) {
                test_j = false;
                break;
            }
        }
        if (test_j) {
            for (int i = 0; i < m; ++i) {
                tot += W(i, j) * rows[i + 1];
            }
            if (tot > 0) {
                if (test_neg > 0) {// then we have elements on both side
                    return null_vect;
                }
                test_pos += 1;
            }
            else if (tot < 0) {
                if (test_pos > 0) {// then we have elements on both side
                    return null_vect;
                }
                test_neg += 1;
            }
        }
    }
    if (test_pos > 0) { // it means that the hyperplan is ax+by+...>=0
        return rows;
    }
    if (test_neg > 0) {// it means that the hyperplan is ax+by+...<=0
        return -rows; // so -(ax+by+...)>=0
    }
    else {
        return null_vect;
    }
}
vector<int> cone_W(const int n, const int m, const GRBEnv& env, vector<vector<double>>& Q, Eigen::MatrixXd W) {
    GRBModel model = GRBModel(env);
    std::vector<GRBVar> vars(m);
    int nb = 0;
    model.set(GRB_IntParam_OutputFlag, 0);
    for (size_t i = 0; i < m; ++i) {
        vars[i] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_INTEGER, "x_" + std::to_string(i)); // Integer variables
    }
    for (auto const& q : Q) {
        GRBLinExpr expr = q[0];
        for (size_t k = 0; k < m; ++k) {
            expr += (q[k + 1] * vars[k]);
        }
        model.addConstr(expr >= 0, "Q_constraint_" + std::to_string(nb));// here we add the constraint that defines Q
        nb += 1;
    }
    vector<Eigen::RowVectorXd> Ineq;
    vector<int> coord;
    for (int i = 0; i < m - 1; ++i) {
        coord.push_back(i);
    }
    do {
        Eigen::RowVectorXd ineq;
        ineq = Hyperplan_def(m, n, W, coord);
        if (ineq != Eigen::RowVectorXd::Zero(m + 1)) {
            bool InIneq = false;
            for (auto const in : Ineq) {
                if (in == ineq) {// We check that the inequality defining the hyperplan haven't already been added 
                    InIneq = true;
                    break;
                }
            }
            if (!(InIneq)) {
                Ineq.push_back(ineq);
                ineq = -ineq; 
                GRBLinExpr constraint = -1;
                //we know that the hyperplane is ax + by + ... >= 0 for a counter example outside we want ax+by+...<=-1
                // so -(ax+by+...)-1>=0
                for (size_t k = 0; k < m; ++k) {
                    constraint += (ineq[k + 1] * vars[k]);
                }
                GRBConstr c = model.addConstr(constraint >= 0, "new_constraint");
                model.optimize();
                int status0 = model.get(GRB_IntAttr_Status);
                if (status0 == GRB_INFEASIBLE) {//There is no integer point so we can continue
                    model.remove(c);
                }
                else if (model.get(GRB_IntAttr_SolCount) > 0) { // We have found a counter example
                    cout << "The counter example is: ";
                    vector<int> counter;
                    for (size_t i = 0; i < m; ++i) {
                        cout << vars[i].get(GRB_DoubleAttr_X) << " ";
                        counter.push_back(vars[i].get(GRB_DoubleAttr_X));
                    }
                    return counter;

                }
            }
        }
    } while (! (coord_spe_W(n, m - 1, coord)));
    return vector<int>(0, m);
}

Eigen::MatrixXd matrix_expansion(Eigen::MatrixXd& W, int& m, int& n) {
    Eigen::MatrixXd test(m, 2 *n + m);
    for (size_t k = 0; k < n; ++k) {
        test.col(k) = W.col(k);
        test.col(n + k) = -W.col(k);
    };
    Eigen::MatrixXd E(m, m);
    E = Eigen::MatrixXd::Identity(m, m);
    for (size_t k = 0; k < m; ++k) {
        test.col(2*n + k) = E.col(k);
    }
    return test;
}
vector<Eigen::VectorXd> find_P(vector<vector<double>>& Parall, int m) {
    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env);
    model.set(GRB_IntParam_OutputFlag, 0);
    std::vector<GRBVar> vars(m);
    for (size_t j = 0; j < m; ++j) {
        vars[j] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_INTEGER, "x_" + std::to_string(j)); // Integer variables
    }
    int nb = 0;
    for (auto const& q : Parall) {
        GRBLinExpr expr = q[0];
        for (size_t k = 0; k < m; ++k) {
            expr += (q[k + 1] * vars[k]);
        }
        model.addConstr(expr >= 0, "Q_constraint_" + std::to_string(nb));// here we add the constraint that defines Q
        nb += 1;
    }
    vector<Eigen::VectorXd> solutions;
    while (true) {
        model.optimize();
        if (model.get(GRB_IntAttr_Status) != GRB_OPTIMAL) {// there is no more solutions
            break;
        }
        Eigen::VectorXd solution;
        solution.resize(m);
        for (size_t j = 0; j < m; ++j) {
            solution(j) = vars[j].get(GRB_DoubleAttr_X);
        }
        solutions.push_back(solution);
        // We add a constraint to be sure not to obtain the same p again
        GRBQuadExpr expr1 = 0;
        for (size_t k = 0; k < m; ++k) {
            expr1 += (vars[k] - solution(k)) * (vars[k] - solution(k));
        }
        model.addQConstr(expr1 >= 1);
    }
    return solutions;
}

vector<Eigen::VectorXd> parallelogram(size_t const m, Eigen::MatrixXd const& W, vector<vector<double>>& Q) {
    vector<vector<double>> Parall(2 * m, vector<double>(m + 1, 0));//je connais déjà la taille
    Eigen::MatrixXd paral(2 * m, m + 1);
    Eigen::MatrixXd W3(m - 1, m - 1);
    bool negative = false;
    for (size_t j = 0; j < m; ++j) {
        Eigen::MatrixXd W2(m, m - 1);
        for (size_t k = 0; k < m; ++k) {
            if (k == j) {
            }
            else if (k < j) {
                W2.col(k) = W.col(k);
            }
            else {
                W2.col(k - 1) = W.col(k);
            }
        }
        int tot = 0;
        Eigen::MatrixXd rows(1, m + 1);
        rows(0, 0) = 0;
        for (size_t i = 0; i < m; ++i) {
            int val = 0;
            if (i == 0) {
                W3 << W2.block(1, 0, m - 1, m - 1);
                val = W3.determinant();
                tot += W(i, j) * val;
                rows(0, i + 1) = val;
                negative = true;
            }
            else if (i == m - 1) {
                if (negative) {
                    W3 << W2.block(0, 0, m - 1, m - 1);
                    val = -W3.determinant();
                    tot += W(i, j) * val;
                    rows(0, i + 1) = val;
                }
                else {
                    W3 << W2.block(0, 0, m - 1, m - 1);
                    val = W3.determinant();
                    tot += W(i, j) * val;
                    rows(0, i + 1) = val;
                }
            }
            else {
                if (negative) {
                    W3 << W2.block(0, 0, i, m - 1), W2.block(i + 1, 0, m - 1 - i, m - 1);
                    val = -W3.determinant();
                    tot += W(i, j) * val;
                    rows(0, i + 1) = val;
                    negative = false;
                }
                else {
                    W3 << W2.block(0, 0, i, m - 1), W2.block(i + 1, 0, m - 1 - i, m - 1);
                    val = W3.determinant();
                    tot += W(i, j) * val;
                    rows(0, i + 1) = val;
                    negative = true;
                }
            }
        }
        if (tot > 0) {//row(0,1)*x_1+...=tot 
            Parall[2 * j][0] = 0;
            Parall[2 * j + 1][0] = tot-1; // as we want to check row(0,1)*x_1+...<tot we will check row(0,1)*x_1+...<=tot-1
            for (size_t k = 1; k < m + 1; ++k) {
                Parall[2 * j][k] = rows(0, k);
                Parall[2 * j + 1][k] = -rows(0, k);
            }
            Q.push_back(Parall[2 * j]);// for Q to be in cone(W)
        }
        else if (tot < 0) {
            Parall[2 * j][0] = 0;
            Parall[2 * j + 1][0] = -(tot+1);
            for (size_t k = 1; k < m + 1; ++k) {
                Parall[2 * j][k] = -rows(0, k);// we want row(0,1)*x_1+...<0 so -row(0,1)*x_1-...>0
                Parall[2 * j + 1][k] = rows(0, k);
            }
            Q.push_back(Parall[2 * j]);// for Q to be in cone(W)
        }

    }
    vector<Eigen::VectorXd> P;
    P = find_P(Parall, m);
    return P;
}

vector<Eigen::VectorXd> solve_17(const vector<Eigen::VectorXd>& C, const Eigen::MatrixXd& W, const Eigen::MatrixXd& W_inv, const Eigen::VectorXd& p, const int& m) {
    vector<Eigen::VectorXd> C_1;
    for (const auto& c : C) {
        Eigen::VectorXd x = W_inv * (c - p);
        bool integer_point = true;// we need to check if x is an integer point
        for (int i = 0; i < m; ++i) {
            double roundedValue = round(x[i]);
            if (abs(x[i] - roundedValue) != 0) {
                integer_point = false;
                break;
            }
        }
        if (integer_point == true) {
            C_1.push_back(x);// we keep W^-1(C-P) because these are the element we want at the end
        }
    }
    return C_1;
}

vector<vector<double>> Image_Q(vector<vector<double>>& Q, Eigen::MatrixXd W, Eigen::VectorXd p, int m) {
    vector < vector<double>> Q_1;
    for (int i = 0; i < Q.size(); ++i) {// Q is as follows : Q[0][0] + x_1*Q[0][1] +x_2*Q[0][2]+... >=0
        vector<double> q = {};
        int t = 0;
        for (int k = 1; k < m + 1; ++k) {
            t += Q[i][k] * p[k - 1];// we adapt directly to Q-p
        }
        q.push_back(Q[i][0] + t);
        for (int j = 0; j < m; ++j) {// Q[0] is of size m+1 but the cst is not interesting
            int T = 0;
            for (int k = 1; k < m + 1; ++k) {
                T += Q[i][k] * W(k - 1, j);
            }
            q.push_back(T);
        }
        Q_1.push_back(q);
    }
    return Q_1;
}

int Square_Matrices(Eigen::MatrixXd& W, vector<vector<double>>& Q, vector<Eigen::VectorXd>& C,int m, int n,GRBEnv env) {
    Eigen::MatrixXd W_square(m, m);
    vector<int> coord(m, 0);
    for (int j = 0; j < m; ++j) {
        coord.push_back(j);
    }
    Eigen::VectorXd null_vect(m);
    null_vect.setZero();
    do {
        for (size_t i = 0; i < m; ++i) {
            W_square.col(i) = W.col(coord[i]);
        }
        if ((W_square.determinant() == 0) || (W_square.determinant() == 1) || (W_square.determinant() == -1)) {
        }
        else {
            vector< Eigen::VectorXd > new_C;
            vector<Eigen::VectorXd> P;
            Eigen::MatrixXd W_inv(m, m);
            W_inv = W_square.inverse();
            vector<vector<double>> Q1=Q;// we make a copie because we are going to modify Q1
            P = parallelogram(m,W_square,Q1);
            for (auto const p : P) {
                if (p != null_vect) {
                    new_C = solve_17(C, W_square, W_inv, p, m);
                    vector<vector<double>> new_Q;
                    new_Q = Image_Q(Q1, W_square, p, m);
                    GRBModel model = GRBModel(env);
                    std::vector<GRBVar> vars(m);
                    for (size_t i = 0; i < m; ++i) {
                        vars[i] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "x_" + std::to_string(i)); // continuous variables
                    }
                    model.set(GRB_IntParam_OutputFlag, 0);
                    int nb = 0;
                    for (auto const& q : new_Q) {
                        GRBLinExpr expr = q[0];
                        for (size_t k = 0; k < m; ++k) {
                            expr += (q[k + 1] * vars[k]);
                        }
                        model.addConstr(expr >= 0, "Q_constraint_" + std::to_string(nb));// here we add the constraint that defines Q
                        nb += 1;
                    }
                    model.optimize();
                    int status0 = model.get(GRB_IntAttr_Status);
                    if (status0 == GRB_OPTIMAL) {
                        vector<int> coord1(m, 0);
                        vector<int> result(m, 0);
                        vector<vector<int>> fin;
                        fin = Solve6(new_C, new_Q, m, coord1, env, result);
                        if (fin[0][0] == 1) {
                            Eigen::VectorXd result1(m);
                            for (int i = 0; i < m; ++i) {
                                result1(i) = result[i];
                            }
                            result1 = W_square * result1;
                            result1 += p; // in order to have the real b
                            GRBModel model1 = GRBModel(env);
                            model1.set(GRB_IntParam_OutputFlag, 0);
                            std::vector<GRBVar> vars1(n);
                            for (size_t j = 0; j < n; ++j) {
                                vars1[j] = model1.addVar(0, GRB_INFINITY, 0.0, GRB_INTEGER, "x_" + std::to_string(j));
                            }
                            for (size_t j = 0; j < m; ++j) {
                                GRBLinExpr expr = 0;
                                int nb = 0;
                                for (size_t k = 0; k < n; ++k) {
                                    expr += W(j, k) * vars1[k];
                                }
                                model1.addConstr(expr <= result[j]);
                                nb += 1;
                            }
                            model1.optimize();
                            int status2 = model1.get(GRB_IntAttr_Status);
                            if (status2 == GRB_INFEASIBLE) {
                                cout << "it is indeed a counter example";
                            }
                            else {
                                cout << "weird";
                            }
                        }
                        cout << endl;
                    }
                    else if (status0 == GRB_INFEASIBLE) {
                        // we just go to the next step
                    }
                    else {
                        cout << "Gurobi isn't working correctly on Q";
                    }
                }
            }
        }
    } while (! (coord_spe_W(n, m, coord)));
    cout << " There is no counterexample";
    return 0;
}

int Main(vector<vector<double>>& Q, Eigen::MatrixXd& W,int n, int m) {
    GRBEnv env = GRBEnv(true);
    env.set(GRB_IntParam_InfUnbdInfo, 1);
    env.start();
    vector<int> result;
    Eigen::MatrixXd W_2;
    W_2 = Allocation_expansion(W, m, n);
    int n_2 = n + m;
    W = matrix_expansion(W, m, n);
    n = 2 * n + m;
    result = cone_W(n, m, env,Q, W);
    if (result != vector<int>(0, m)) {
        cout << endl << "The counter example is" << endl;
        for (auto r : result) {
            cout << r << " ";
        }
        return 0;

    }
    vector<Eigen::VectorXd> C;
    long long int Delta = W.cwiseAbs().maxCoeff();// the max abs value of W
    long long int inf_norm= (W.cwiseAbs().rowwise().sum()).maxCoeff(); //the infinity norm of W
    long long int bound = max_det_sub_mat(n_2, m, W_2)* inf_norm*(n-m);//we do it on W_2 bcs the abs val of sub det of (W Id) is the same as the abs val of the sub det of (W -W Id)
    long long int  bound_2 = inf_norm*m*pow(2*m*Delta+1,m);
    if (bound > bound_2) {// we compare the two bounds
        bound = bound_2;
    }
    Eigen::VectorXd Max = W_2.cwiseAbs().rowwise().maxCoeff();// Max[i] contain the biggest max_j(W(i,j))
    C = gen_C_5(W, Max, bound, m, n);
    Square_Matrices(W, Q, C, m, n, env);
    return 0;
}