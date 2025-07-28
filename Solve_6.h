#ifndef SOLVE_6_H
#define SOLVE_6_H
#include <iostream>
#include <vector>
#include <algorithm>
#include "gurobi_c++.h"
#include <chrono>
#include <fstream>
#include <cmath>
#include <Eigen/Dense>
#include "Construction_Allocation.h"
#include "Usefull_fct.h"
using namespace std;

vector<bool> right(vector<int> actual_point, int min_y, const vector<vector<int>>& points, const vector<Eigen::VectorXd>& C, vector<int>& K, vector<bool>& max_K, vector<int> direction, const GRBEnv& env, vector<vector<double>>& Q, int m, vector<int>& result);
// check all the cells on the right of the one we are located at that point
vector<bool> bypass_2D_down(vector<int>& actual_point, int min_y, const vector<vector<int>>& points, const vector<Eigen::VectorXd>& C, vector<int>& K, vector<bool> max_K, vector<int> direction, const GRBEnv& env, vector<vector<double>>& Q, int m, vector<int>& result);
// will follow the vertices of polyhedra that are going "down" regarding the second axis
vector<bool> bypass_2D_up(vector<int>& actual_point, int min_y, const vector<vector<int>>& points, const vector<Eigen::VectorXd>& C, vector<int>& K, vector<bool> max_K, vector<int> direction, const GRBEnv& env, vector<vector<double>>& Q, int m, vector<int>& result);
// will follow the vertices of polyhedra that are going "up" regarding the second axis
bool bypass_2D(const vector<int>& starting_point, int min_y, const vector<vector<int>>& points, const vector<Eigen::VectorXd>& C, vector<int>& K, vector<bool> max_K, const GRBEnv& env, vector<vector<double>>& Q, int m, vector<int>& result);
//check all the cells that are needed to be checked when the direction> 2 are bounded
vector<vector<int>> check_cell(vector<int>& coord, const vector<vector<int>>& points, const vector<Eigen::VectorXd>& C, vector<int>& K, vector<bool> max_K, size_t m, int& n, const GRBEnv& env, vector<vector<double>>& Q, vector<int>& result,vector<int>& axes);
//check all the cells that are needed to be checked
vector<vector<int>> Solve6(vector<Eigen::VectorXd>& C, vector<vector<double>>& Q, size_t m, vector<int>& coord, const GRBEnv& env, vector<int>& result);
// split the space into cells
bool Coord(vector<int>& coord, size_t direction, vector<int>& K, vector<bool> max_K, size_t m, vector<int>& axes);
// gives the next cell we need "hyperplane " we need to check

#endif
#pragma once
