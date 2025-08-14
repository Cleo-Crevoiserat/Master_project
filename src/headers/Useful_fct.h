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
using namespace std; 

bool coord_spe_W(int n, int m, vector<int>& coord);
// is used to create a loop choosing every set of m vectors in a set of n
vector< Eigen::VectorXd> minimality(vector<Eigen::VectorXd>& C, int m);
// keep only the minimal elements of the list C
bool coord_base(Eigen::VectorXd& coord, size_t direction, int born, size_t m);
// generate all integer points of dimension m with infinity norm less than born
bool coord_Spe_C_Alloc(vector<Eigen::VectorXd>& coord, int direction_1, int direction_2, Eigen::VectorXd N, int m, int n);
// this is the way to generate all elements in C depending on the N(i)
int max_det_sub_mat(int n, int m, Eigen::MatrixXd& W);
//compute the max sub det of W
int max_det_sub_mat_alloc(int n, int m, Eigen::MatrixXd& W);
// compute the max sub det of a mat for alloc
bool coord_spe_gen_C(Eigen::VectorXd& coord, size_t direction, size_t& face_nb, Eigen::VectorXd& Max, int born, size_t m);
//generate the coordonates of the different point that could be in C for a face which is bounded
bool Around_faces_C(Eigen::VectorXd& coord, size_t direction, size_t& face_nb, Eigen::VectorXd& Max, int born, size_t m);
//go through every faces 
vector<Eigen::VectorXd> gen_C_5(Eigen::MatrixXd W, Eigen::VectorXd& Max, int bound, int m, int n);
//Generate all elements in C
int factorial(int n);
// compute n!
#endif
#pragma once
