#ifndef CONTRUCTION_PART_1_H
#define CONSTRUCTION_PART_1_H
#include <iostream>
#include <vector>
#include <algorithm>
#include "gurobi_c++.h"
#include <chrono>
#include <fstream>
#include <cmath>
#include <Eigen/Dense>
#include "Construction_Allocation.h"
#include "Useful_fct.h"
using namespace std;

Eigen::RowVectorXd Hyperplan_def(size_t const m, size_t const n, Eigen::MatrixXd const& W, vector<int> coord);
//return one of the hyperplane defining cone(W)
vector<int> cone_W(const int n, const int m, const GRBEnv& env, vector<vector<double>>& Q, Eigen::MatrixXd W);
// generates all the hyperplanes defining cone(W) and check wethere there is an integer point in Q outside of the cone(W)
Eigen::MatrixXd matrix_expansion(Eigen::MatrixXd& W, int& m, int& n);
// generates the new matrix going from the question Wx<=b to W'x'=b
vector<Eigen::VectorXd> find_P(vector<vector<double>>& Parall, int m);
// find all the integer point in the parallelograme of W
vector<Eigen::VectorXd> parallelogram(size_t const m, Eigen::MatrixXd const& W, vector<vector<double>>& Q);
// construct the paralellogram of (W) and find all integer point in it and add inequalities to Q in order for it to be in cone(W)
vector<Eigen::VectorXd> solve_17(const vector<Eigen::VectorXd>& C, const Eigen::MatrixXd& W, const Eigen::MatrixXd& W_inv, const Eigen::VectorXd& p, const int& m);
//keep only the elements in C such that c equiv p mod (the lattice of W)
vector<vector<double>> Image_Q(vector<vector<double>>& Q, Eigen::MatrixXd W, Eigen::VectorXd p, int m);
//goes from Q to Q'=W^-1(Q-p)
int Square_Matrices(Eigen::MatrixXd& W, vector<vector<double>>& Q, vector<Eigen::VectorXd>& C, int m, int n, GRBEnv env);
// we create the sub square matrices of W and go solve the end of the problem as in the paper
int Main(vector<vector<double>>& Q, Eigen::MatrixXd& W, int n, int m);
// the main for general matrix W of size mxn and polyhedron Q

#endif
#pragma once
