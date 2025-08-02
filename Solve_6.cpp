#include <iostream>
#include <vector>
#include <algorithm>
#include "gurobi_c++.h"
#include <chrono>
#include <cmath>
#include "Construction_Allocation.h"
#include "Usefull_fct.h"
#include "Construction_Part_1.h"
using namespace std;

vector<bool> right(vector<int> actual_point, int min_y, const vector<vector<int>>& points, const vector<Eigen::VectorXd>& C, vector<int>& K, vector<bool>& max_K, vector<int> direction, const GRBEnv& env, vector<vector<double>>& Q, int m,vector<int>& result,vector<int>& axes) {
    //
    vector<bool> T;
    bool equal = false;
    do {
        if (((actual_point[axes[0]] < K[axes[0]]) and (not max_K[axes[0]])) or ((actual_point[axes[0]] < K[axes[0]] - 1) and (max_K[axes[0]]))) {
            vector<int> point;
            actual_point[axes[0]] += 1;
            GRBModel model = GRBModel(env);
            model.set(GRB_IntParam_OutputFlag, 0);
            std::vector<GRBVar> vars(m);
            for (size_t i = 0; i < m; ++i) {
                vars[i] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "x_" + std::to_string(i)); // continuous variables
            }
            for (size_t k = 0; k < m; ++k) {
                if (actual_point[k] == K[k]) {// if in one direction we are at the last choice, we only need the last inequation
                    model.addConstr(vars[k] >= points[k][actual_point[k]] - 0.5, "constraint_" + std::to_string(k));
                }
                else if (actual_point[k] == -1) {// if we only need the first one
                    model.addConstr(vars[k] <= points[k][0] - 0.5, "constraint_" + std::to_string(k));
                }
                else {
                    model.addConstr(vars[k] >= points[k][actual_point[k]] - 0.5, "constraint_" + std::to_string(k));
                    model.addConstr(vars[k] <= points[k][actual_point[k] + 1] - 0.5, "constraint_" + std::to_string(k));
                }
            }
            int nb = 0;
            for (auto const& q : Q) {
                GRBLinExpr expr = q[0];
                for (size_t k = 0; k < m; ++k) {
                    expr += (q[k + 1] * vars[k]);
                }
                model.addConstr(expr >= 0, "Q_constraint_" + std::to_string(nb));// here we add the constraint that defines Q
                nb += 1;
            }
            model.optimize();
            int status0 = model.get(GRB_IntAttr_Status);
            if (status0 == GRB_UNBOUNDED) {
            }
            else if (status0 == GRB_INFEASIBLE) {// it means that we are outside of Q on the right
                T.push_back(false);
                T.push_back(false);// no need to continue on the right
                return T;
            }
            else if (model.get(GRB_IntAttr_SolCount) > 0) {
                vector<int> sol;
                for (size_t i = 0; i < m; ++i) {
                    sol.push_back(vars[i].get(GRB_DoubleAttr_X));
                }
                for (size_t i = 0; i < m; ++i) {// we create the element representing the hyperplan
                    if (actual_point[i] != -1) {
                        point.push_back(points[i][actual_point[i]]);
                    }
                }
                bool test = false;
                for (const auto& c : C) {
                    bool contained = true;
                    for (size_t j = 0; j < m; ++j) {
                        if (point[j] < c[j]) { // it means that the cell would not be in the hyperplane
                            contained = false;
                            break;
                        }
                    }
                    if (contained == true) {// the cell is contained in the union of the hyperplanes
                        test = true;
                        T.push_back(true);
                        T.push_back(true);// we continue
                        return T;
                    }
                }
                if (test == false) {  // we have an intersection with Q which is not in the union of hyperplans so we need to search for integer point
                    GRBModel model = GRBModel(env);
                    std::vector<GRBVar> vars(m);
                    model.set(GRB_IntParam_OutputFlag, 0);

                    for (size_t i = 0; i < m; ++i) {
                        vars[i] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_INTEGER, "x_" + std::to_string(i)); // Integer variables
                    }
                    for (size_t k = 0; k < m; ++k) {
                        if (actual_point[k] == K[k]) {// if in one direction we only need the last inequation
                            model.addConstr(vars[k] >= points[k][actual_point[k]] - 0.5, "constraint_" + std::to_string(k));
                        }
                        else if (actual_point[k] == -1) {// if we only need the first one
                            model.addConstr(vars[k] <= points[k][0] - 0.5, "constraint_" + std::to_string(k));
                        }
                        else {
                            model.addConstr(vars[k] >= points[k][actual_point[k]] - 0.5, "constraint_" + std::to_string(k));
                            model.addConstr(vars[k] <= points[k][actual_point[k] + 1] - 0.5, "constraint_" + std::to_string(k));
                        }
                    }
                    int nb = 0;
                    for (auto const& q : Q) {
                        GRBLinExpr expr = q[0];
                        for (size_t k = 0; k < m; ++k) {
                            expr += (q[k + 1] * vars[k]);
                        }
                        model.addConstr(expr >= 0, "Q_constraint_" + std::to_string(nb));// here we add the constraint that defines Q
                        nb += 1;
                    }
                    model.optimize();
                    int status0 = model.get(GRB_IntAttr_Status);
                    if (status0 == GRB_UNBOUNDED) {
                    }
                    else if (status0 == GRB_INFEASIBLE) {//the cell has an intersection with Q but there is no integer point so it's oke but we have to continue
                    }
                    else if (model.get(GRB_IntAttr_SolCount) > 0) {// on print la mauvaise solution si il y en a une
                        cout << std::endl;
                        cout << "The counter example is: ";
                        for (size_t i = 0; i < m; ++i) {
                            cout << vars[i].get(GRB_DoubleAttr_X) << " ";
                            result[i] = vars[i].get(GRB_DoubleAttr_X);
                        }
                        T.push_back(true);
                        T.push_back(false);// we stop because we found a counter example
                        return T;
                    }
                }
            }
        }
        else {
            T.push_back(true);
            T.push_back(true);
            return T;
        }
    } while (not(equal));
    return T;
}

vector<bool> bypass_2D_down(vector<int>& actual_point, int min_y, const vector<vector<int>>& points, const vector<Eigen::VectorXd>& C, vector<int>& K, vector<bool> max_K, vector<int> direction, const GRBEnv& env, vector<vector<double>>& Q, int m, vector<int>& result,vector<int>& axes ) {
    //T[0]= false: we stop T[1]: false -> we use 1,1 direct true -> 0,-1 T[2]=false if counter example found T[3]= we use bypass_up from here
    vector<bool> T;
    vector<int> point;
    actual_point[axes[0]] += direction[0];
    actual_point[axes[1]] += direction[1];
    GRBModel model = GRBModel(env);
    model.set(GRB_IntParam_OutputFlag, 0);
    std::vector<GRBVar> vars(m);
    //if (((direction[0] == 0 and direction[1] == 0) or (direction[0] == 1 and direction[1] == 1)) and (((actual_point[1] < K[1]) and (not max_K[1])) or ((actual_point[1] < K[1] - 1) and (max_K[1])))) {// if we are at the start or after going in [1,1] direction and that we can go up without being out of range
    //    if (not(up(actual_point, min_y, points, C, K, direction, env, Q, m, result))) {
    //        T.push_back(false);
    //        T.push_back(false);
    //        T.push_back(false);// il y eu un contre exemple au dessus
    //        return T;
    //    }
    //}// don't think we need it bcs right and bypass_up
    for (size_t i = 0; i < m; ++i) {
        vars[i] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "x_" + std::to_string(i)); // continuous variables
    }
    for (size_t k = 0; k < m; ++k) {
        if (actual_point[k] == K[k]) {
            model.addConstr(vars[k] >= points[k][actual_point[k]] - 0.5, "constraint_" + std::to_string(k));
        }
        else if (actual_point[k] == -1) {// if we only need the first one
            model.addConstr(vars[k] <= points[k][0] - 0.5, "constraint_" + std::to_string(k));
        }
        else {
            model.addConstr(vars[k] >= points[k][actual_point[k]] - 0.5, "constraint_" + std::to_string(k));
            model.addConstr(vars[k] <= points[k][actual_point[k] + 1] - 0.5, "constraint_" + std::to_string(k));
            
        }
    }
    int nb = 0;
    for (auto const& q : Q) {
        GRBLinExpr expr = q[0];
        for (size_t k = 0; k < m; ++k) {
            expr += (q[k + 1] * vars[k]);
        }
        model.addConstr(expr >= 0, "Q_constraint_" + std::to_string(nb));// here we add the constraint that defines Q
        nb += 1;
    }
    model.optimize();
    int status0 = model.get(GRB_IntAttr_Status);
    if (status0 == GRB_UNBOUNDED) {
    }
    else if (status0 == GRB_INFEASIBLE) {
        if (direction[1] == 1) {
            T.push_back(false);// we stop because it means that we are going up
            T.push_back(false);
            T.push_back(true);
            T.push_back(true);// we use bypass_2D_up from here
            return T;
        }
        else {
            T.push_back(true);// we continue 
            T.push_back(false);// but we move with {1,1} for direction
            T.push_back(true);
            return T;
        }
    }
    else if (model.get(GRB_IntAttr_SolCount) > 0) {
        if ((direction[0] == 1)) {// this cell has already been checked by right so we can just directly continue
            T.push_back(true);// we continue 
            T.push_back(true);
            T.push_back(true);
            return T;
        }
        vector<int> sol;
        for (size_t i = 0; i < m; ++i) {
            sol.push_back(vars[i].get(GRB_DoubleAttr_X));
        }
        for (size_t i = 0; i < m; ++i) {// creating the point representing the hyperplan
            if (actual_point[i] != -1) {
                point.push_back(points[i][actual_point[i]]);
            }
        }
        bool test = false;
        for (const auto& c : C) {
            bool equal = true;
            for (size_t j = 0; j < m; ++j) {
                if (point[j] < c[j]) { // it means that the cell would not be in the hyperplane
                    equal = false;
                    break;
                }
            }
            if (equal == true) {// the cell is contained in the union of the hyperplanes
                test = true;
                if (min_y >= point[1]) {// Q will not have smaller coordinated in the second axis and then everything bigger in the first direction will be contained in the hyperplane
                    T.push_back(false);// we stop
                    T.push_back(false);
                    T.push_back(true);
                    return T;
                }
                else {// we can still go down in the 2nd direction
                    T.push_back(true);// we continue
                    T.push_back(true);// in direction 0,-1
                    T.push_back(true);
                    return T;
                }

                break;
            }
        }
        if (test == false) { // so here we have an intersection with Q but which is not contained in the hyperplan so we need to check if there is an integer point
            GRBModel model = GRBModel(env);
            model.set(GRB_IntParam_OutputFlag, 0);
            std::vector<GRBVar> vars(m);

            for (size_t i = 0; i < m; ++i) {
                vars[i] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_INTEGER, "x_" + std::to_string(i)); // Integer variables
            }
            for (size_t k = 0; k < m; ++k) {
                if (actual_point[k] == K[k]) {// if in one direction we only need the last inequation
                    model.addConstr(vars[k] >= points[k][actual_point[k]] - 0.5, "constraint_" + std::to_string(k));
                }
                else if (actual_point[k] == -1) {// if we only need the first one
                    model.addConstr(vars[k] <= points[k][0] - 0.5, "constraint_" + std::to_string(k));
                }
                else {
                    model.addConstr(vars[k] >= points[k][actual_point[k]] - 0.5, "constraint_" + std::to_string(k));
                    model.addConstr(vars[k] <= points[k][actual_point[k] + 1] - 0.5, "constraint_" + std::to_string(k));
                }
            }
            int nb = 0;
            for (auto const& q : Q) {
                GRBLinExpr expr = q[0];
                for (size_t k = 0; k < m; ++k) {
                    expr += (q[k + 1] * vars[k]);
                }
                model.addConstr(expr >= 0, "Q_constraint_" + std::to_string(nb));// here we add the constraint that defines Q
                nb += 1;
            }
            model.optimize();
            int status0 = model.get(GRB_IntAttr_Status);
            if (status0 == GRB_UNBOUNDED) {
                int m = 0;
            }
            else if (status0 == GRB_INFEASIBLE) {//the cell has an intersection with Q but there is no integer point so it's oke
                vector<bool> L;// in that case we should check with bigger value on the first axis
                if (max_K[axes[0]]) {
                    if (actual_point[axes[0]] < K[axes[0]] - 1) {
                        L = right(actual_point, min_y, points, C, K, max_K, direction, env, Q, m,result,axes);// we are in the polyhedra so we can check on the right
                        if (L[0] == true) {
                            if (L[1] == false) {// we have found a counter example on the "right"
                                T.push_back(true);
                                T.push_back(true);
                                T.push_back(false);
                                return T;
                            }
                        }
                    }
                }
                else {
                    if (actual_point[axes[0]] < K[axes[0]]) {
                        L = right(actual_point, min_y, points, C, K, max_K, direction, env, Q, m,result,axes);// we are in the polyhedra so we can check on the right
                        if (L[0] == true) {
                            if (L[1] == false) { // we have found a counter example on the "right"
                                T.push_back(true);
                                T.push_back(true);
                                T.push_back(false);
                                return T;
                            }
                        }
                    }
                }
                T.push_back(true);
                T.push_back(true);// if we already at the limit of the second axis then we might only be able to go "down"
                T.push_back(true);
                return T;
            }
            else if (model.get(GRB_IntAttr_SolCount) > 0) {// there is a counter example
                cout << std::endl;
                cout << "The counter example is: ";
                for (size_t i = 0; i < m; ++i) {
                    cout << vars[i].get(GRB_DoubleAttr_X) << " ";
                    result[i] = vars[i].get(GRB_DoubleAttr_X);
                }
                T.push_back(true);
                T.push_back(true);
                T.push_back(false);// we have a counter example
                return T;
            }
        }
    }
    return T;
}

vector<bool> bypass_2D_up(vector<int>& actual_point, int min_y, const vector<vector<int>>& points, const vector<Eigen::VectorXd >& C, vector<int>& K, vector<bool> max_K, vector<int> direction, const GRBEnv& env, vector<vector<double>>& Q, int m,vector<int>& result,vector<int>& axes) {
    vector<bool> T;
    vector<int> point;
    actual_point[axes[0]] += direction[0];
    actual_point[axes[1]] += direction[1];
    GRBModel model = GRBModel(env);
    model.set(GRB_IntParam_OutputFlag, 0);
    std::vector<GRBVar> vars(m);
    //if ((direction[0] == 1 and direction[1] == -1) and (actual_point[1] > 0)) {// enlever ça mais add right
    //    if (not(down(actual_point, min_y, points, C, K, direction, env, Q, m,result))) {
    //        T.push_back(false);
    //        T.push_back(false);
    //        T.push_back(false);// il y eu un contre exemple en dessous
    //        return T;
    //    }
    //}
    for (size_t i = 0; i < m; ++i) {
        vars[i] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "x_" + std::to_string(i)); // continuous variables
    }
    for (size_t k = 0; k < m; ++k) {
        if (actual_point[k] == K[k]) {// if in one direction we are at the last choice, we only need the last inequation
            model.addConstr(vars[k] >= points[k][actual_point[k]] - 0.5, "constraint_" + std::to_string(k));
        }
        else if (actual_point[k] == -1) {// if we only need the first one
            model.addConstr(vars[k] <= points[k][0] - 0.5, "constraint_" + std::to_string(k));
        }
        else {
            model.addConstr(vars[k] >= points[k][actual_point[k]] - 0.5, "constraint_" + std::to_string(k));
            model.addConstr(vars[k] <= points[k][actual_point[k] + 1] - 0.5, "constraint_" + std::to_string(k));
        }
    }
    int nb = 0;
    for (auto const& q : Q) {
        GRBLinExpr expr = q[0];
        for (size_t k = 0; k < m; ++k) {
            expr += (q[k + 1] * vars[k]);
        }
        model.addConstr(expr >= 0, "Q_constraint_" + std::to_string(nb));// here we add the constraint that defines Q
        nb += 1;
    }
    model.optimize();
    int status0 = model.get(GRB_IntAttr_Status);
    if (status0 == GRB_UNBOUNDED) {
    }
    else if (status0 == GRB_INFEASIBLE) {
        if (direction[0] == 1) {//it means that Q is not going up anymore
            T.push_back(false);// we stop
            T.push_back(false);
            T.push_back(true);
            return T;
        }
        else {
            T.push_back(true);// we continue 
            T.push_back(false);// but we check by adding {1,-1} to the axis
            T.push_back(true);
            return T;
        }
    }
    else if (model.get(GRB_IntAttr_SolCount) > 0) {
        if (direction[1] == -1) {// these cell would have already been checked by right
            T.push_back(true);// we continue 
            T.push_back(true);
            T.push_back(true);
            return T;
        }
        vector<int> sol;
        for (size_t i = 0; i < m; ++i) {
            sol.push_back(vars[i].get(GRB_DoubleAttr_X));
        }
        for (size_t i = 0; i < m; ++i) {// we create the point representing the hyperplane
            if (actual_point[i] != -1) {
                point.push_back(points[i][actual_point[i]]);
            }
        }
        bool test = false;
        for (const auto& c : C) {
            bool equal = true;
            for (size_t j = 0; j < m; ++j) {
                if (point[j] < c[j]) { // it means that the cell would not be in the hyperplane
                    equal = false;
                    break;
                }
            }
            if (equal == true) {// the cell is contained in the union of the hyperplanes
                test = true;
                T.push_back(false);
                T.push_back(false);
                T.push_back(true);
                return T;
            }
        }
        if (test == false) {// so here we have an intersection with Q but which is not contained in the hyperplan so we need to check if there is an integer point
            GRBModel model = GRBModel(env);
            model.set(GRB_IntParam_OutputFlag, 0);
            std::vector<GRBVar> vars(m);

            for (size_t i = 0; i < m; ++i) {
                vars[i] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_INTEGER, "x_" + std::to_string(i)); // Integer variables
            }
            for (size_t k = 0; k < m; ++k) {
                if (actual_point[k] == K[k]) {// if in one direction we only need the last inequation
                    model.addConstr(vars[k] >= points[k][actual_point[k]] - 0.5, "constraint_" + std::to_string(k));
                }
                else if (actual_point[k] == -1) {// if we only need the first one
                    model.addConstr(vars[k] <= points[k][0] - 0.5, "constraint_" + std::to_string(k));
                }
                else {
                    model.addConstr(vars[k] >= points[k][actual_point[k]] - 0.5, "constraint_" + std::to_string(k));
                    model.addConstr(vars[k] <= points[k][actual_point[k] + 1] - 0.5, "constraint_" + std::to_string(k));
                }
            }
            int nb = 0;
            for (auto const& q : Q) {
                GRBLinExpr expr = q[0];
                for (size_t k = 0; k < m; ++k) {
                    expr += (q[k + 1] * vars[k]);
                }
                model.addConstr(expr >= 0, "Q_constraint_" + std::to_string(nb));// here we add the constraint that defines Q
                nb += 1;
            }
            model.optimize();
            int status0 = model.get(GRB_IntAttr_Status);
            if (status0 == GRB_UNBOUNDED) {
                int m = 0;
            }
            else if (status0 == GRB_INFEASIBLE) {//the cell has an intersection with Q but there is no integer point so it's oke but we have to check on the right
                vector<bool> L;
                bool remonte;
                if (max_K[axes[0]]) {
                    if (actual_point[axes[0]] < K[axes[0]] - 1) {
                        L = right(actual_point, min_y, points, C, K, max_K, direction, env, Q, m, result,axes);// we are in the polyhedra so we can check on the right
                        if (L[0] == true) {
                            if (L[1] == false) {// we have found a counter example on the right
                                T.push_back(true);
                                T.push_back(true);
                                T.push_back(false);
                                return T;
                            }
                        }
                    }
                }
                else {
                    if (actual_point[axes[0]] < K[axes[0]]) {
                        L = right(actual_point, min_y, points, C, K, max_K, direction, env, Q, m, result,axes);// we are in the polyhedra so we can check on the right
                        if (L[0] == true) {
                            if (L[1] == false) {// we have found a counter example on the right
                                T.push_back(true);
                                T.push_back(true);
                                T.push_back(false);
                                return T;
                            }
                        }
                    }
                }
                T.push_back(true);
                T.push_back(true);
                T.push_back(true);
                return T;
            }
            else if (model.get(GRB_IntAttr_SolCount) > 0) {
                cout << std::endl;
                cout << "The counter example is: ";
                for (size_t i = 0; i < m; ++i) {
                    cout << vars[i].get(GRB_DoubleAttr_X) << " ";
                    result[i] = vars[i].get(GRB_DoubleAttr_X);
                }
                T.push_back(true);
                T.push_back(true);
                T.push_back(false);
                return T;
            }
        }
    }
    return T;
}

bool bypass_2D(const vector<int>& starting_point, int min_y, const vector<vector<int>>& points, const vector<Eigen::VectorXd>& C, vector<int>& K, vector<bool> max_K, const GRBEnv& env, vector<vector<double>>& Q, int m,vector<int>& result,vector<int>& axes) {
    vector<int> direction = { 0,0 };
    vector<int> actual_point = starting_point;
    bool step;
    int nb_iter = 0;
    do {
        nb_iter += 1;
        vector<bool> T = bypass_2D_down(actual_point, min_y, points, C, K, max_K, direction, env, Q, m,result,axes);
        step = T[0];
        if (T[2] == false) {// we have found a counter exemple
            return false;
        }
        if (T[1] == true) {
            if ((actual_point[axes[1]] == 0) or (actual_point[axes[1]] <= min_y)) {// we can't go "down" anymore
                break; 
            }
            else {
                direction = { 0,-1 };
            }
        }
        if (T[1] == false) {
            if (max_K[axes[0]]) {
                if (actual_point[axes[0]]< K[axes[0]] - 1) {
                    if (((max_K[axes[1]]) and (actual_point[axes[1]] < K[axes[1]] - 1)) or ((not max_K[axes[1]]) and (actual_point[axes[1]] < K[axes[1]]))) {
                        direction = { 1,1 };
                    }
                    else {
                        direction = { 1,0 };
                    }
                }
                else {
                    break;
                }
            }
            else {
                if (actual_point[axes[0]] < K[axes[0]]) {
                    if (((max_K[axes[1]]) and (actual_point[axes[1]] < K[axes[1]] - 1)) or ((not max_K[axes[1]]) and (actual_point[axes[1]] < K[axes[1]]))) {
                        direction = { 1,1 };
                    }
                    else {
                        direction = { 1,0 };
                    }
                }
                else {
                    break;
                }
            }
        }
    } while (step);
    direction = { 0,1 };// we can go up one directly
    actual_point = starting_point;
    if (not (((actual_point[axes[1]] == K[axes[1]]) and (not max_K[axes[1]])) or ((actual_point[axes[1]] == K[axes[1]] - 1) and (max_K[axes[1]])))) {//if we are not already at the top
        do {
            nb_iter += 1;
            vector<bool> T = bypass_2D_up(actual_point, min_y, points, C, K, max_K, direction, env, Q, m, result,axes);
            step = T[0];
            if (T[2] == false) {// we have found a counter exemple
                return false;
            }
            if (T[1] == true) {
                if (((actual_point[axes[1]] == K[axes[1]]) and (not max_K[axes[1]])) or ((actual_point[axes[1]] == K[axes[1]] - 1) and (max_K[axes[1]]))) {// if we are at the "top" of the y axis then it's over
                    break;
                }
                else {
                    direction = { 0,1 };
                }
            }
            if (T[1] == false) {
                if (((actual_point[axes[0]]== K[axes[0]]) and (not max_K[axes[0]])) or ((actual_point[axes[0]] == K[axes[0]] - 1) and (max_K[axes[0]]))) {// if we are at the most on the right here we could only go down which is already check so it's over too
                    break;
                }
                else if (actual_point[axes[1]] == 0) {// we cannot go down anymore and the right has already been checked so it's over
                    break;
                }
                else {
                    direction = { 1,-1 };
                }
            }
        } while (step);
    }
    return true;
}
bool Coord(vector<int>& coord, size_t direction, vector<int>& K, vector<bool> max_K, size_t m, vector<int>& axes) {// in order to create a cartesian product
    if (m <= 2) {
        return true;// contourn_2D will do everything in one time
    }
    if (((axes[0] == m - 1) or (axes[1] == m - 1)) and ((axes[0] == m - 2) or (axes[1] == m - 2))) {
        if (direction == m - 3) {
            if (max_K[m - 3]) {
                if (coord[m - 3] == K[m - 3] - 1) {// As max_K is true if the values are smaller than the max
                    return true;
                }
                else {
                    coord[m - 3] += 1;
                    return false;
                }
            }
            else {
                if (coord[m - 3] == K[m - 3]) {
                    return true;
                }
                coord[m - 3] += 1;
                return false;
            }
        }
        else {
            if (max_K[direction]) {
                if (coord[direction] >= K[direction] - 1) {
                    coord[direction] = 0;
                    direction += 1;
                    return Coord(coord, direction, K, max_K, m, axes);
                }
                else {
                    coord[direction] += 1;
                    return false;
                }
            }
            else {
                if (coord[direction] >= K[direction]) {
                    coord[direction] = 0;
                    direction += 1;
                    return Coord(coord, direction, K, max_K, m, axes);
                }
                else {
                    coord[direction] += 1;
                    return false;
                }
            }
        }
    }
    else if ((axes[0] == m - 1) or (axes[1] == m - 1)) {
        if (direction == m - 2) {
            if (max_K[m - 2]) {
                if (coord[m - 2] == K[m - 2] - 1) {// As max_K is true if the values are smaller than the max
                    return true;
                }
                else {
                    coord[m - 2] += 1;
                    return false;
                }
            }
            else {
                if (coord[m - 2] == K[m - 2]) {
                    return true;
                }
                coord[m - 2] += 1;
                return false;
            }
        }
        else {
            if (((axes[0] == direction) or (axes[1] == direction))) {
                return Coord(coord, direction+1, K, max_K, m, axes);
            }
            if (max_K[direction]) {
                if (coord[direction] >= K[direction] - 1) {
                    coord[direction] = 0;
                    direction += 1;
                    return Coord(coord, direction, K, max_K, m, axes);
                }
                else {
                    coord[direction] += 1;
                    return false;
                }
            }
            else {
                if (coord[direction] >= K[direction]) {
                    coord[direction] = 0;
                    direction += 1;
                    return Coord(coord, direction, K, max_K, m, axes);
                }
                else {
                    coord[direction] += 1;
                    return false;
                }
            }
        }
    }
    else {
        if (direction == m - 1) {
            if (max_K[m - 1]) {
                if (coord[m - 1] == K[m - 1] - 1) {// As max_K is true if the values are smaller than the max
                    return true;
                }
                else {
                    coord[m - 1] += 1;
                    return false;
                }
            }
            else {
                if (coord[m - 1] == K[m - 1]) {
                    return true;
                }
                coord[m - 1] += 1;
                return false;
            }
        }
        else {
            if (((axes[0] == direction) or (axes[1] == direction))) {
                return Coord(coord, direction + 1, K, max_K, m, axes);
            }
            if (max_K[direction]) {
                if (coord[direction] >= K[direction] - 1) {
                    coord[direction] = 0;
                    direction += 1;
                    return Coord(coord, direction, K, max_K, m, axes);
                }
                else {
                    coord[direction] += 1;
                    return false;
                }
            }
            else {
                if (coord[direction] >= K[direction]) {
                    coord[direction] = 0;
                    direction += 1;
                    return Coord(coord, direction, K, max_K, m, axes);
                }
                else {
                    coord[direction] += 1;
                    return false;
                }
            }
        }
    }
}

vector<vector<int>> check_cell(vector<int>& coord, const vector<vector<int>>& points, const vector<Eigen::VectorXd>& C, vector<int>& K, vector<bool> max_K, size_t m, int& n, const GRBEnv& env, vector<vector<double>>& Q, vector<int>& result,vector<int>& axes) {
    size_t direction = 0;
    bool step = true;
    int min_y;
    do {
        vector<int> starting_point(m,0);
        GRBModel model = GRBModel(env);
        model.set(GRB_IntParam_OutputFlag, 0);
        std::vector<GRBVar> vars(m);
        for (size_t j = 0; j < m; ++j) {
            vars[j] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "x_" + std::to_string(j)); // continuous variables
        }
        model.setObjective(-vars[axes[0]], GRB_MAXIMIZE);// we searching the min for the first direction
        int nb = 0;
        if (m > 2) {
            for (size_t k = 0; k < m; ++k) {// we let x and y be free
                if ((k != axes[0]) and (k != axes[1])) {
                    if (coord[k] == K[k]) {// if in one direction we are at the last choice, we only need the last inequation
                        model.addConstr(vars[k] >= points[k][coord[k]] - 0.5, "constraint_" + std::to_string(k));
                    }
                    else if (coord[k] == -1) {// if we only need the first one
                        model.addConstr(vars[k] <= points[k][0] - 0.5, "constraint_" + std::to_string(k));
                    }
                    else {
                        model.addConstr(vars[k] >= points[k][coord[k]] - 0.5, "constraint_" + std::to_string(k));
                        model.addConstr(vars[k] <= points[k][coord[k] + 1] - 0.5, "constraint_" + std::to_string(k));
                    }
                }
            }
        }
        for (auto const& q : Q) {
            GRBLinExpr expr = q[0];
            for (size_t k = 0; k < m; ++k) {
                expr += (q[k + 1] * vars[k]);
            }
            model.addConstr(expr >= 0, "Q_constraint_" + std::to_string(nb));// here we add the constraint that defines Q
            nb += 1;
        }
        model.optimize();
        int status0 = model.get(GRB_IntAttr_Status);
        bool go_check = true;
        if (status0 == GRB_UNBOUNDED) {
        }
        else if (status0 == GRB_INFEASIBLE) {
            // std::cout << "The problem is infeasible in x !" << std::endl;
            go_check = false;
        }
        else if (status0 == GRB_OPTIMAL) {
            for (int i = 0; i < 2; i++) {
                vector<int> list0 = points[axes[i]];
                double l = round(list0.size() / 2);
                int l_step_supp = list0.size() - 1;
                int l_step_inf = 0;
                if (list0[0] - 0.5 < vars[axes[i]].get(GRB_DoubleAttr_X) and (l > 1)) {
                    while ((list0[l] - 0.5 < vars[axes[i]].get(GRB_DoubleAttr_X)) or (list0[l - 1] - 0.5 > vars[axes[i]].get(GRB_DoubleAttr_X))) {
                        if (l == list0.size() - 1) {// if we are at the end;
                            starting_point[axes[i]] = l;
                            break;
                        }
                        if (list0[l] - 0.5 < vars[axes[i]].get(GRB_DoubleAttr_X)) {
                            int step = l;
                            l = l + round(abs(l_step_supp - l) / 2);
                            l_step_inf = step;
                        }
                        if (list0[l - 1] - 0.5 > vars[axes[i]].get(GRB_DoubleAttr_X)) {
                            int step = l;
                            l = l - round(abs(l-l_step_inf) / 2);
                            l_step_supp = step;
                        }
                    }
                    if (l != list0.size()-1) {
                        starting_point[axes[i]] = l - 1; // l-1 because (list0[l-1]-0.5 < vars[i].get(GRB_DoubleAttr_X)<list0[l]
                    }
                }
                else {
                    starting_point[axes[i]] = 0;
                }
            }
        }
        else {
            go_check = false;
        }
        for (int i = 0; i < m; i++) {
            if ((i != axes[0]) and (i != axes[1])) {
                starting_point[i]= coord[i];
            }
        }
        GRBModel model1 = GRBModel(env);
        model1.set(GRB_IntParam_OutputFlag, 0);
        std::vector<GRBVar> vars1(m);
        for (size_t j = 0; j < m; ++j) {
            vars1[j] = model1.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "x_" + std::to_string(j)); // continuous variables
        }
        model1.setObjective(-vars1[axes[1]], GRB_MAXIMIZE);// on cherche le min pour y
        int nb1 = 0;
        if (m > 2) {
            for (size_t k = 0; k < m; ++k) {// we let x and y be free
                if ((k != axes[0]) and (k != axes[1])) {
                    if (coord[k] == K[k]) {// if in one direction we are at the last choice, we only need the last inequation
                        model1.addConstr(vars1[k] >= points[k][coord[k]] - 0.5, "constraint_" + std::to_string(k));
                    }
                    else if (coord[k] == -1) {// if we only need the first one
                        model1.addConstr(vars1[k] <= points[k][0] - 0.5, "constraint_" + std::to_string(k));
                    }
                    else {
                        model1.addConstr(vars1[k] >= points[k][coord[k]] - 0.5, "constraint_" + std::to_string(k));
                        model1.addConstr(vars1[k] <= points[k][coord[k] + 1] - 0.5, "constraint_" + std::to_string(k));
                    }
                }
            }
        }
        for (auto const& q : Q) {
            GRBLinExpr expr1 = q[0];
            for (size_t k = 0; k < m; ++k) {
                expr1 += (q[k + 1] * vars1[k]);
            }
            model1.addConstr(expr1 >= 0, "Q_constraint_" + std::to_string(nb1));// here we add the constraint that defines Q
            nb1 += 1;
        }
        model1.optimize();
        int status1 = model1.get(GRB_IntAttr_Status);
        if (status1 == GRB_UNBOUNDED) {
            std::cout << "The problem is unbounded!" << std::endl;
        }
        else if (status1 == GRB_INFEASIBLE) {// if it was infeasible it would have been already before
        }
        else {
            min_y = vars1[axes[1]].get(GRB_DoubleAttr_X);
        }
        if (go_check) {// we have an intersection with Q
            //cout << "starting point =";
            //for (auto s : starting_point) {
            //    cout << s << " ";
            //}
            //cout << endl ;
            step = bypass_2D(starting_point, min_y, points, C, K, max_K, env, Q, m,result,axes);
        }
        if (not step) {
            std::cout << "on a perdu !";
            return { {1},{0} };
        }
    } while (Coord (coord, direction, K, max_K, m,axes) == false);

    // Check if we have reached the end of coordinates

    return { {0},{0} };
}

vector<vector<int>> Solve6(vector<Eigen::VectorXd>& C, vector<vector<double>>& Q, size_t m, vector<int>& coord, const GRBEnv& env, vector<int>& result) {
    vector<vector<int>> List = {};
    vector<int> K = {};
    vector<bool> K_max = {};
    vector<int> minimal_elem = {};
    vector<int> axes;
    C = minimality(C, m);
    for (int i = 0; i < m; i++) {
        sort(C.begin(), C.end(), [i](const Eigen::VectorXd& a, const Eigen::VectorXd& b) {
            return a[i] < b[i]; // sort criteria is based on coordinate i
            });
        vector<int> list0;
        list0.push_back(C[0][i]);
        minimal_elem.push_back(C[0][i]);
        for (Eigen::VectorXd c : C) {
            if (list0.back() != c[i]) {
                list0.push_back(c[i]);
            }
        }
        GRBModel model = GRBModel(env);
        model.set(GRB_IntParam_OutputFlag, 0);
        std::vector<GRBVar> vars(m);

        // Add m integer variables to the model
        for (size_t j = 0; j < m; ++j) {
            vars[j] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "x_" + std::to_string(j)); // real variables
        }
        model.setObjective(-vars[i], GRB_MAXIMIZE); // we find the minimal value of direction i in Q
        int nb = 0;
        for (auto const& q : Q) {
            GRBLinExpr expr = q[0];
            for (size_t k = 0; k < m; ++k) {
                expr += (q[k + 1] * vars[k]);
            }
            model.addConstr(expr >= 0, "Q_constraint_" + std::to_string(nb));// here we add the constraint that defines Q
            nb += 1;
        }
        model.optimize();
        int status0 = model.get(GRB_IntAttr_Status);
        if (status0 == GRB_UNBOUNDED) {
            break;
        }
        else if (status0 == GRB_INFEASIBLE) {
            std::cout << "The problem is infeasible!, There is no points in Q" << std::endl;
        }
        else if (status0 == GRB_OPTIMAL) {
            double l = round(list0.size() / 2);
            int l_step_supp = list0.size();
            int l_step_inf = 0;
            if (l_step_supp > 1) {
                if (list0[0] - 0.5 <= vars[i].get(GRB_DoubleAttr_X) and (list0[l_step_supp - 1] - 0.5 >= vars[i].get(GRB_DoubleAttr_X))) {// if the min of Q in direction is contained between the min and max of the C[i] values
                    if (list0[0] - 0.5 == vars[i].get(GRB_DoubleAttr_X)) {

                    }
                    else if ((list0[l_step_supp - 1] - 0.5 == vars[i].get(GRB_DoubleAttr_X))) {
                        list0 = { list0[l_step_supp - 1] };
                    }
                    else {
                        while ((list0[l] - 0.5 < vars[i].get(GRB_DoubleAttr_X)) or (list0[l - 1] - 0.5 > vars[i].get(GRB_DoubleAttr_X))) {
                            if (list0[l] - 0.5 < vars[i].get(GRB_DoubleAttr_X)) {
                                int step = l;
                                l = l + round(abs(l_step_supp - l) / 2);
                                l_step_inf = step;
                            }
                            if (list0[l - 1] - 0.5 > vars[i].get(GRB_DoubleAttr_X)) {
                                int step = l;
                                l = l - round(abs(l-l_step_inf) / 2);
                                l_step_supp = step;
                            }
                            if (l == 0) {
                                break;
                            }

                        }
                        if (l == list0.size() - 1) {
                            list0.erase(list0.begin(), list0.begin() + l - 1);
                        }
                        else if (l - 2 >= 0) {
                            list0.erase(list0.begin(), list0.begin() + l - 1);
                        }// the only remaining case is if min is between the first and second value of the c[i] so then we just don't do anything
                    }
                }
                else if (list0[l_step_supp - 1] <= vars[i].get(GRB_DoubleAttr_X)) {// if the min is already bigger than everything
                    list0.erase(list0.begin(), list0.begin() + l_step_supp - 1);

                }
                else {// if the min is smaller than the range of the hyperplane 
                    // we then need to check if there is an integer vector with a value smaller in that direction
                    GRBModel model2 = GRBModel(env);
                    model2.set(GRB_IntParam_OutputFlag, 0);
                    std::vector<GRBVar> vars2(m);
                    int nb2 = 0;
                    for (size_t j = 0; j < m; ++j) {
                        vars2[j] = model2.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_INTEGER, "x_" + std::to_string(i)); // Integer variables
                    }
                    for (auto const& q : Q) {
                        GRBLinExpr expr2 = q[0];
                        for (size_t k = 0; k < m; ++k) {
                            expr2 += (q[k + 1] * vars2[k]);
                        }
                        model2.addConstr(expr2 >= 0, "Q_constraint_" + std::to_string(nb2));// here we add the constraint that defines Q
                        nb2 += 1;
                    }
                    model2.addConstr(list0[0] - vars2[i] - 1 >= 0);//vars2[i] <= list0[0]-1
                    GRBLinExpr linExpr2 = vars2[i];
                    model2.setObjective(linExpr2, GRB_MAXIMIZE);// we maximize the value but with setting up to be smaller than the min of C[i]
                    model2.optimize();
                    int status = model2.get(GRB_IntAttr_Status);
                    if (status == GRB_UNBOUNDED) {
                        //cout << "maximal in direction " << i << " = " << "infinity" << std::endl;
                    }
                    else if (status == GRB_INFEASIBLE) {
                        //std::cout << "The problem is infeasible!" << std::endl;
                    }
                    else if (status == GRB_OPTIMAL) {// if we found a counter example
                        cout << " counter example is" << std::endl;
                        for (size_t i = 0; i < m; ++i) {
                            cout << vars2[i].get(GRB_DoubleAttr_X) << ",";
                            result[i] = vars2[i].get(GRB_DoubleAttr_X);
                        }
                        return{ {1},{0} };
                    }
                }
            }
            else {// if size is 1
                if (list0[l_step_supp - 1] > vars[i].get(GRB_DoubleAttr_X)) {// if the min is smaller than the range of the hyperplane 
                    // we then need to check if there is an integer vector with a value smaller in that direction
                    GRBModel model2 = GRBModel(env);
                    model2.set(GRB_IntParam_OutputFlag, 0);
                    std::vector<GRBVar> vars2(m);
                    int nb2 = 0;
                    for (size_t j = 0; j < m; ++j) {
                        vars2[j] = model2.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_INTEGER, "x_" + std::to_string(i)); // Integer variables
                    }
                    for (auto const& q : Q) {
                        GRBLinExpr expr2 = q[0];
                        for (size_t k = 0; k < m; ++k) {
                            expr2 += (q[k + 1] * vars2[k]);
                        }
                        model2.addConstr(expr2 >= 0, "Q_constraint_" + std::to_string(nb2));// here we add the constraint that defines Q
                        nb2 += 1;
                    }
                    model2.addConstr(list0[0] - vars2[i] - 1 >= 0);//vars2[i] <= list0[0]-1
                    GRBLinExpr linExpr2 = vars2[i];
                    model2.setObjective(linExpr2, GRB_MAXIMIZE);// we maximize the value but with setting up to be smaller than the min of C[i]
                    model2.optimize();
                    int status = model2.get(GRB_IntAttr_Status);
                    if (status == GRB_UNBOUNDED) {
                    }
                    else if (status == GRB_INFEASIBLE) {
                    }
                    else if (status == GRB_OPTIMAL) {
                        cout << " counter example is" << std::endl;
                        for (size_t i = 0; i < m; ++i) {
                            cout << vars2[i].get(GRB_DoubleAttr_X) << ",";
                            result[i] = vars2[i].get(GRB_DoubleAttr_X);
                        }
                        return{ {1},{0} };
                    }
                }
            }
        }
        GRBModel model1 = GRBModel(env);
        model1.set(GRB_IntParam_OutputFlag, 0);
        std::vector<GRBVar> vars1(m);

        // Add m integer variables to the model
        for (size_t j = 0; j < m; ++j) {
            vars1[j] = model1.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "x_" + std::to_string(j)); // real variables
        }
        int nb1 = 0;
        for (auto const& q : Q) {
            GRBLinExpr expr1 = q[0];
            for (size_t k = 0; k < m; ++k) {
                expr1 += (q[k + 1] * vars1[k]);
            }
            model1.addConstr(expr1 >= 0, "Q_constraint_" + std::to_string(nb));// here we add the constraint that defines Q
            nb1 += 1;
        }
        GRBLinExpr linExpr1 = vars1[i];
        model1.setObjective(linExpr1, GRB_MAXIMIZE);// this time we maximize the value of direction i in Q
        model1.optimize();
        int status = model1.get(GRB_IntAttr_Status);
        if (status == GRB_UNBOUNDED) {
            K_max.push_back(false);// there is point of Q which are bigger than the max of the C[i]
        }
        else if (status == GRB_INFEASIBLE) {
            //std::cout << "The problem is infeasible!" << std::endl;
        }
        else if (status == GRB_OPTIMAL) {
            double l1 = round(list0.size() / 2);
            int l_step1_supp = list0.size();
            int l_step1_inf = 0;
            if (list0[l_step1_supp - 1] >= vars1[i].get(GRB_DoubleAttr_X)) {// if the max is contained in between the C[i]
                if (list0[l_step1_supp - 1] > vars1[i].get(GRB_DoubleAttr_X)) {
                    K_max.push_back(true);// K_max[i] is true if the max is contained in between the C[i]
                }
                else {
                    K_max.push_back(false);// if K_max[i]==C[i]_max then we have to check x[i]>=C[i] max-0.5
                }
                if (list0[l_step1_supp - 1] == vars1[i].get(GRB_DoubleAttr_X)) {
                }
                else {
                    while ((list0[l1] - 0.5 < vars1[i].get(GRB_DoubleAttr_X)) or (list0[l1 - 1] - 0.5 > vars1[i].get(GRB_DoubleAttr_X))) {
                        if (list0[l1] - 0.5 < vars1[i].get(GRB_DoubleAttr_X)) {
                            int step1 = l1;
                            l1 = l1 + round(abs(l_step1_supp - l1) / 2);
                            l_step1_inf = step1;
                        }
                        if (list0[l1 - 1] - 0.5 > vars1[i].get(GRB_DoubleAttr_X)) {
                            int step1 = l1;
                            l1 = l1 - round(abs(l1-l_step1_inf) / 2);
                            l_step1_supp = step1;
                        }
                    }
                    if (l1 == 0) {
                        list0.resize(1);
                    }
                    else if (l1 >= 0) {
                        list0.resize(l1 + 1);
                    }
                }
            }
            else {
                K_max.push_back(false);
            }
        }
        K.push_back(list0.size() - 1);
        if ((list0.size() > 2) and (axes.size() < 2)) {
            axes.push_back(i);
        }
        List.push_back(list0);// it contains for each direction what value would the intersection take
    }
    if (axes.size() < 2) {// almost all directions have only one element
        int j = 0;
        do {
            axes.push_back(j);
            j += 1;
        } while (axes.size() < 2);
    }
    int max_1 = 0;
    int max_2 = 0;
    for (int i = 0; i < m; ++i) {
        if (K[i] > max_1) {
            axes[1] = axes[0];
            max_2 = max_1;
            axes[0] = i;
            max_1 = K[i];
        }
        else if (K[i] > max_2) {
            axes[1] = i;
            max_2 = K[i];
        }
    }
    cout << "axes=";
    for (auto a : axes) {
        cout << a << " ";
    }
    cout << endl ;
    int n = 0;
    //cout << "prepa finito" << endl;
    for (auto l : List) {
        cout << l.size() << endl;
    }
    return check_cell(coord, List, C, K, K_max, m, n, env, Q, result, axes);
}
