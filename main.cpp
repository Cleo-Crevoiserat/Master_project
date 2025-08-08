#include <iostream>
#include <vector>
#include <algorithm>
#include "gurobi_c++.h"
#include <chrono>
#include <cmath>
#include "Construction_Allocation.h"
#include "Usefull_fct.h"
#include "Construction_Part_1.h"
#include "Solve_6.h"
using namespace std;

int main(){
    GRBEnv env = GRBEnv(true);
    env.set(GRB_IntParam_InfUnbdInfo, 1);
    env.start();
    int prob = 1;// here chose either 0 or 1
    if (prob == 0) {// you want to solve Wx<=b 
        int m = 3;//you give the number of rows of your matrix
        int n = 2;//you give the number of columns of your matrix
        Eigen::MatrixXd W(3,2);//you define your matrix
        W << 1, 1,
            -2, 1,
            1, -2;
        vector<vector<double>> Q;// you define Q
        Q= {
        {2,0,-1,0},
        {2,-1,0,0},
        {2,0,0,-1},
        {1,1,0,0},
        {1,0,1,0},
        {1,0,0,1}
        };
        Main(Q, W, 2, 3);
    }
    if (prob == 1) {// you want to solve a fair allocation problem
        int nb_agent = 2;
        int nb_object = 4;
        Eigen::VectorXd N(nb_object);// define the nuber of item per object
        N << 1, 1, 1, 1;
        vector < vector<int>> Utility = { {279, 262, 233, 226}, {100, 484, 208, 208} };// you define your utility function
        Main_alloc(Utility, nb_object, nb_agent, N);
    }
    return 0;
}
