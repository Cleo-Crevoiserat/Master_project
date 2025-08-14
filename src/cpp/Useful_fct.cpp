#include <iostream>
#include <vector>
#include <algorithm>
#include <ciso646>
#include "gurobi_c++.h"
#include <chrono>
#include <fstream>
#include <cmath>
#include <Eigen/Dense>
using namespace std;

bool coord_spe_W(int n, int m, vector<int>& coord) {// is used to create a loop choosing every set of m vectors in a set of n
    if (coord[m - 1] != n - 1) {
        coord[m - 1] += 1;
    }
    else {
        bool test = true;
        for (int i = m - 2; i > -1; --i) {
            if (coord[i] != n - m + i) {
                coord[i] += 1;
                test = false;
                for (size_t j = 1; j < m - i; ++j) {
                    coord[i + j] = coord[i] + j;
                }
                break;
            }
        }
        if (test) {
            return true;
        }
    }
    return false;
}

vector< Eigen::VectorXd> minimality(vector<Eigen::VectorXd>& C, int m) {
    vector< Eigen::VectorXd> new_C;
    new_C.push_back(C[0]);
    for (int i = 1; i < C.size();++i) {
        size_t j = 0;
        bool dominated_i = true;
        bool egal = true;
        while (j<new_C.size()) {
            bool dominated_j = true;
            dominated_i = true;
            egal = true;
            for (size_t k = 0; k < m; k++) {
                if (new_C[j][k] < C[i][k]) { //then c[j] is not dominated by c[i]
                    dominated_j = false;
                    egal= false;
                }
                else if (new_C[j][k] > C[i][k]) { // then c[i] is not dominated by c[j]
                    dominated_i = false;
                    egal =false;
                }
            }
            if (egal) {
                break;
            }
            else if (dominated_j) {// then every coordinate of c[j] is > c[i] so it's actually dominated
                new_C.erase(new_C.begin() + j);
            }
            else if (dominated_i) {
                break;
            }
            else {
                j++;
            }
        }
        if (egal or dominated_i) {// we do nothing

        }
        else {
            new_C.push_back(C[i]);
        }
    }
    return new_C;
}

bool coord_base(Eigen::VectorXd& coord, size_t direction, int born, size_t m) {// in order to create a cartesian product
    if (direction == m - 1) {
        if (coord[m - 1] == born) {
            return true;
        }
        coord[m - 1] += 1;
        return false;
    }
    else {
        if (coord[direction] == born) {
            coord[direction] = -born;
            direction += 1;
            return coord_base(coord, direction, born, m);
        }
        else {
            coord[direction] += 1;
            return false;
        }
    }
}

bool coord_Spe_C_Alloc(vector<Eigen::VectorXd>& coord, int direction_1, int direction_2, Eigen::VectorXd N, int m, int n) {
    if (direction_1 == 0) {
        if (coord[0][0] == n) {
            return true;
        }
    }
    if (coord[direction_1][direction_2] == n) {
        if (direction_2 == 0) {
            for (int i = 0; i < N[direction_1]; ++i) {
                coord[direction_1][i] = 1;
            }
            return coord_Spe_C_Alloc(coord, direction_1 - 1, N[direction_1-1]-1, N, m, n);
        }
        else {
            return coord_Spe_C_Alloc(coord, direction_1, direction_2-1, N, m, n);
        }
    }
    else {
        coord[direction_1][direction_2] += 1;
        int c = coord[direction_1][direction_2];
        for (int i = direction_2+1; i < N[direction_1]; ++i) {
            coord[direction_1][i] = c;
        }
        return false;
    }
}
int max_det_sub_mat(int n, int m, Eigen::MatrixXd & W) {
    Eigen::MatrixXd square_mat(m, m);
    vector<int> coord;
    int max=0;
    for (int j = 0; j < m; ++j) {
        coord.push_back(j);
    }
    do {
        for (size_t i = 0; i < m; ++i) {
            square_mat.col(i) = W.col(coord[i]);
        }
        int det = abs(square_mat.determinant());
        if (det > max) {
            max = det;
        }
    } while (not(coord_spe_W(n, m, coord)));
    return max;
}

int max_det_sub_mat_alloc(int n, int m, Eigen::MatrixXd& W) {
    Eigen::MatrixXd W2(m + n + 1, m * n + m + n + 1);
    W2 = W.block(0, 0, m + n + 1, m * n + m + n + 1);
    Eigen::MatrixXd square_mat(m + n + 1, m + n + 1);
    vector<int> coord;
    int max = 0;
    for (int j = 0; j < m+n+1; ++j) {
        coord.push_back(j);
    }
    do {
        for (size_t i = 0; i < m+n+1; ++i) {
            square_mat.col(i) = W2.col(coord[i]);
        }
        int det = abs(square_mat.determinant());
        if (det > max) {
            max = det;
        }
    } while (not(coord_spe_W(m * n + m + n + 1, m + n + 1, coord)));
    return max;
}

bool coord_spe_gen_C(Eigen::VectorXd& coord, size_t direction, size_t& face_nb, Eigen::VectorXd& Max, int born, size_t m) {// in order to create a cartesian product
    if (face_nb != m - 1) {// we split the two cases because it doesn't end at the same time
        if (direction == m - 1) {
            if (direction < face_nb) {
                if (coord[m - 1] == born - Max[m - 1]) {
                    return true;
                }
                coord[m - 1] += 1;
                return false;
            }
            else {
                if (coord[m - 1] == born ) {
                    return true;
                }
                coord[m - 1] += 1;
                return false;

            }
        }
        else {
            if (direction != face_nb) {
                if (direction < face_nb) {
                    if (coord[direction] == born - Max[direction]) {// for direction which are not the actual face there values goes between -born+Max to born-Max
                        coord[direction] = -born + Max[direction];
                        direction += 1;
                        return coord_spe_gen_C(coord, direction, face_nb, Max, born, m);
                    }
                    else {
                        coord[direction] += 1;
                        return false;
                    }
                }
                else {
                    if (coord[direction] == born) {
                        coord[direction] = -born ;
                        direction += 1;
                        return coord_spe_gen_C(coord, direction, face_nb, Max, born, m);
                    }
                    else {
                        coord[direction] += 1;
                        return false;
                    }
                }
            }
            else {//if face_nb=direction
                if (coord[direction] == -born + Max[direction] - 1) { // for direction ==face_nb we check in the values [-born;-born+Max]U[born-Max;born]
                    coord[direction] = born - Max[direction] + 1;
                    return false;
                }
                else if (coord[direction] == born) {
                    coord[direction] = -born;
                    direction += 1;
                    return coord_spe_gen_C(coord, direction, face_nb, Max, born, m);
                }
                else {
                    coord[direction] += 1;
                    return false;
                }

            }
        }
    }
    else {// now face==m-1
        if (direction == m - 1) {
            if (coord[m - 1] == born) {// it's the end
                return true;
            }
            else if (coord[m - 1] == -born + Max[m - 1] - 1) {// same arg as before
                coord[m - 1] = born - Max[m - 1] + 1;
                return false;
            }
            else {
                coord[m - 1] += 1;
                return false;
            }
        }
        else {
            if (direction != face_nb) {// we only have the case direction<face_nb as face_nb=m-1
                if (coord[direction] == born - Max[direction]) {
                    coord[direction] = -born + Max[direction];
                    direction += 1;
                    return coord_spe_gen_C(coord, direction, face_nb, Max, born, m);
                }
                else {
                    coord[direction] += 1;
                    return false;
                }
            }
        }

    }
}

bool Around_faces_C(Eigen::VectorXd& coord, size_t direction, size_t& face_nb, Eigen::VectorXd& Max, int born, size_t m) {
    if (coord_spe_gen_C(coord, direction, face_nb, Max, born, m)) {
        if (face_nb == m - 1) {
            return true;
        }
        else {
            face_nb += 1;// change the face
            for (size_t i = 0; i < m; ++i) {//generate the starting point
                if (i < face_nb) {
                    coord[i] = -born + Max[i];
                }
                else {
                    coord[i] = -born;
                }
            }
            return false;

        }
    }
    else {
        return false;
    }
}

vector<Eigen::VectorXd> gen_C(Eigen::MatrixXd W, Eigen::VectorXd& Max, int bound, int m, int n) {
    vector<Eigen::VectorXd> C_1;
    int w = 0;
    int tot = 0;
    GRBEnv env = GRBEnv();
    Eigen::VectorXd coord(m);
    size_t face_nb = 0;
    for (size_t i = 0; i < m; ++i) {
        if (i < face_nb) {
            coord[i] = -bound + Max[i];
        }
        else {
            coord[i] = -bound;
        }
    }
    GRBModel model = GRBModel(env);
    model.set(GRB_IntParam_OutputFlag, 0);
    std::vector<GRBVar> vars(n);
    for (size_t i = 0; i < n; ++i) {
        vars[i] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_INTEGER, "x_" + std::to_string(i)); // Integer variables from 0 to + inf
    }
    int nb = 0;
    vector<GRBConstr> Constr;
    for (size_t i = 0; i < m; ++i) {
        GRBLinExpr expr = 0;
        for (size_t j = 0; j < n; ++j) {
            expr += W(i, j) * vars[j];// we want W(i)*x=c(i)
        }
        GRBConstr c = model.addConstr(expr == coord[i], "matrix_equal_" + std::to_string(nb));
        Constr.push_back(c);
    }
    do {
        if (w == 0) {
            model.optimize();
        }
        else {
            for (int j = 0; j < m; ++j) {
                Constr[j].set(GRB_DoubleAttr_RHS, coord[j]);
            }
            model.optimize();
        }
        if (model.get(GRB_IntAttr_SolCount) > 0) {
            C_1.push_back(coord);

            if (face_nb != 0) {//as W contain the identity we know that all the points with the same coordinate except the first one which is bigger will also be in C
                for (int i = coord[0]; i < bound - Max[0]; ++i) {
                    coord[0] = i + 1;
                    C_1.push_back(coord);
                }
            }
        }
        w += 1;
    } while (Around_faces_C(coord, 0, face_nb, Max, bound, m) == false);
    return C_1;
}

int factorial(int n) {
    int result = 1;
    for (int i = 2; i <= n; ++i)
        result *= i;
    return result;
}
