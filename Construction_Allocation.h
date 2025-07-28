#ifndef CONTRUCTION_ALLOCATION_H
#define CONSTRUCTION_ALLOCATION_H
#include <iostream>
#include <vector>
#include <algorithm>
#include "gurobi_c++.h"
#include <chrono>
#include <fstream>
#include <cmath>
#include <Eigen/Dense>
using namespace std;
vector<vector<double>> Q_create(const int& n, const int& m, Eigen::VectorXd& N, vector<vector<int>>& Utility);
// create the Polyhedron Q for envy free alocation as explained in the paper
Eigen::MatrixXd Allocation_Matrix(const int& n, const int& m, Eigen::VectorXd& N, vector<vector<int>>& Utility);
// create the W matrix that checks wether or not the allocation in Q is dominated (more detailed in the paper)
vector<Eigen::VectorXd> gen_C_Alloc(Eigen::MatrixXd& W, Eigen::VectorXd& N, int m, int n, int n_1, int m_1);
// Generate the set C depending on the multiplicity of each item
vector<Eigen::VectorXd> gen_C_Alloc_1(Eigen::MatrixXd& W, int bound, int m, int n, int mn, GRBEnv& env);
// generate all element in c in Cone (W) with infinity norm less than bound
Eigen::MatrixXd Allocation_expansion(Eigen::MatrixXd& W, int& m, int& n);
//// generates the new matrix going from the question Wx<=b to W'x'=b
int Main_alloc(vector< vector<int>>& Utility, int m, int n, Eigen::VectorXd& N);
//solve the fair allocation problem 
#endif
