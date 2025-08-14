#ifndef ALLOCATION_H
#define ALLOCATION_H
#include <iostream>
#include <vector>
#include <algorithm>
#include "gurobi_c++.h"
#include <chrono>
#include <fstream>
#include <cmath>
#include <Eigen/Dense>
using namespace std;
bool coord_spe_W_4(int n, int m, vector<int>& coord);
vector<vector<int>> projection(const vector<vector<int>>& Q, const int& m, const int& n, GRBEnv& env);
vector<vector<int>> Q_create(const int& n, const int& m, Eigen::VectorXd N, vector<vector<int>> Utility);
vector<int> retour_alloc(const vector<vector<int>>& Q, const int& m, const int& n, GRBEnv& env, const vector<int>& b, vector<vector<int>> Utility);
Eigen::MatrixXd Allocation_Matrix(const int& n, const int& m, Eigen::VectorXd N, vector<vector<int>> Utility);
Eigen::RowVectorXd Hyperplan_def(size_t const m, size_t const n, Eigen::MatrixXd const& W, vector<int> coord);
vector<int> cone_W(const int n, const int m, const GRBEnv& env, vector<vector<int>> Q, Eigen::MatrixXd W);
Eigen::MatrixXd matrix_expansion(Eigen::MatrixXd& W, int& m, int& n);
void minimality_2(vector<Eigen::VectorXd>& C, size_t m);
bool coord_11(Eigen::VectorXd& coord, size_t direction, int born, size_t m);
vector<Eigen::VectorXd> find_P_3(vector<vector<int>> Parall, int m);
void keep_P_3(vector<Eigen::VectorXd>& P, Eigen::MatrixXd const& W, Eigen::MatrixXd const& W_inv, int m);
vector<Eigen::VectorXd> parallelograme_4(size_t const m, Eigen::MatrixXd const& W, Eigen::MatrixXd const& W_inv);
bool coord_Gen_C_poubelle(Eigen::VectorXd& coord, size_t direction, int born, int& max_0, size_t m);
vector<Eigen::VectorXd>  gen_C_5(Eigen::MatrixXd W, int bound, int m, int n);
vector<Eigen::VectorXd> solve_17_5(const vector<Eigen::VectorXd>& C, const Eigen::MatrixXd& W, const Eigen::MatrixXd& W_inv, const Eigen::VectorXd& p, const int& m);
void square_matrices(int m, int n, vector<int> coord, Eigen::MatrixXd W, const vector<Eigen::VectorXd>& C);
bool coord_Gen_C_100(Eigen::VectorXd& coord, size_t& direction, Eigen::VectorXd& Max, size_t n);
vector<vector<int>>gen_C_Alloc_v3(Eigen::MatrixXd W, Eigen::VectorXd N, int m, int n, int n_1, int m1);
int mainaaaaaaaaaaa();
#endif