#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 5.991464547107981;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with sympy 1.8                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_7406943224577723449) {
   out_7406943224577723449[0] = delta_x[0] + nom_x[0];
   out_7406943224577723449[1] = delta_x[1] + nom_x[1];
   out_7406943224577723449[2] = delta_x[2] + nom_x[2];
   out_7406943224577723449[3] = delta_x[3] + nom_x[3];
   out_7406943224577723449[4] = delta_x[4] + nom_x[4];
   out_7406943224577723449[5] = delta_x[5] + nom_x[5];
   out_7406943224577723449[6] = delta_x[6] + nom_x[6];
   out_7406943224577723449[7] = delta_x[7] + nom_x[7];
   out_7406943224577723449[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_6282749298357684830) {
   out_6282749298357684830[0] = -nom_x[0] + true_x[0];
   out_6282749298357684830[1] = -nom_x[1] + true_x[1];
   out_6282749298357684830[2] = -nom_x[2] + true_x[2];
   out_6282749298357684830[3] = -nom_x[3] + true_x[3];
   out_6282749298357684830[4] = -nom_x[4] + true_x[4];
   out_6282749298357684830[5] = -nom_x[5] + true_x[5];
   out_6282749298357684830[6] = -nom_x[6] + true_x[6];
   out_6282749298357684830[7] = -nom_x[7] + true_x[7];
   out_6282749298357684830[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_4025963510183789805) {
   out_4025963510183789805[0] = 1.0;
   out_4025963510183789805[1] = 0.0;
   out_4025963510183789805[2] = 0.0;
   out_4025963510183789805[3] = 0.0;
   out_4025963510183789805[4] = 0.0;
   out_4025963510183789805[5] = 0.0;
   out_4025963510183789805[6] = 0.0;
   out_4025963510183789805[7] = 0.0;
   out_4025963510183789805[8] = 0.0;
   out_4025963510183789805[9] = 0.0;
   out_4025963510183789805[10] = 1.0;
   out_4025963510183789805[11] = 0.0;
   out_4025963510183789805[12] = 0.0;
   out_4025963510183789805[13] = 0.0;
   out_4025963510183789805[14] = 0.0;
   out_4025963510183789805[15] = 0.0;
   out_4025963510183789805[16] = 0.0;
   out_4025963510183789805[17] = 0.0;
   out_4025963510183789805[18] = 0.0;
   out_4025963510183789805[19] = 0.0;
   out_4025963510183789805[20] = 1.0;
   out_4025963510183789805[21] = 0.0;
   out_4025963510183789805[22] = 0.0;
   out_4025963510183789805[23] = 0.0;
   out_4025963510183789805[24] = 0.0;
   out_4025963510183789805[25] = 0.0;
   out_4025963510183789805[26] = 0.0;
   out_4025963510183789805[27] = 0.0;
   out_4025963510183789805[28] = 0.0;
   out_4025963510183789805[29] = 0.0;
   out_4025963510183789805[30] = 1.0;
   out_4025963510183789805[31] = 0.0;
   out_4025963510183789805[32] = 0.0;
   out_4025963510183789805[33] = 0.0;
   out_4025963510183789805[34] = 0.0;
   out_4025963510183789805[35] = 0.0;
   out_4025963510183789805[36] = 0.0;
   out_4025963510183789805[37] = 0.0;
   out_4025963510183789805[38] = 0.0;
   out_4025963510183789805[39] = 0.0;
   out_4025963510183789805[40] = 1.0;
   out_4025963510183789805[41] = 0.0;
   out_4025963510183789805[42] = 0.0;
   out_4025963510183789805[43] = 0.0;
   out_4025963510183789805[44] = 0.0;
   out_4025963510183789805[45] = 0.0;
   out_4025963510183789805[46] = 0.0;
   out_4025963510183789805[47] = 0.0;
   out_4025963510183789805[48] = 0.0;
   out_4025963510183789805[49] = 0.0;
   out_4025963510183789805[50] = 1.0;
   out_4025963510183789805[51] = 0.0;
   out_4025963510183789805[52] = 0.0;
   out_4025963510183789805[53] = 0.0;
   out_4025963510183789805[54] = 0.0;
   out_4025963510183789805[55] = 0.0;
   out_4025963510183789805[56] = 0.0;
   out_4025963510183789805[57] = 0.0;
   out_4025963510183789805[58] = 0.0;
   out_4025963510183789805[59] = 0.0;
   out_4025963510183789805[60] = 1.0;
   out_4025963510183789805[61] = 0.0;
   out_4025963510183789805[62] = 0.0;
   out_4025963510183789805[63] = 0.0;
   out_4025963510183789805[64] = 0.0;
   out_4025963510183789805[65] = 0.0;
   out_4025963510183789805[66] = 0.0;
   out_4025963510183789805[67] = 0.0;
   out_4025963510183789805[68] = 0.0;
   out_4025963510183789805[69] = 0.0;
   out_4025963510183789805[70] = 1.0;
   out_4025963510183789805[71] = 0.0;
   out_4025963510183789805[72] = 0.0;
   out_4025963510183789805[73] = 0.0;
   out_4025963510183789805[74] = 0.0;
   out_4025963510183789805[75] = 0.0;
   out_4025963510183789805[76] = 0.0;
   out_4025963510183789805[77] = 0.0;
   out_4025963510183789805[78] = 0.0;
   out_4025963510183789805[79] = 0.0;
   out_4025963510183789805[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_2334097729121616737) {
   out_2334097729121616737[0] = state[0];
   out_2334097729121616737[1] = state[1];
   out_2334097729121616737[2] = state[2];
   out_2334097729121616737[3] = state[3];
   out_2334097729121616737[4] = state[4];
   out_2334097729121616737[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_2334097729121616737[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_2334097729121616737[7] = state[7];
   out_2334097729121616737[8] = state[8];
}
void F_fun(double *state, double dt, double *out_8837167548770486690) {
   out_8837167548770486690[0] = 1;
   out_8837167548770486690[1] = 0;
   out_8837167548770486690[2] = 0;
   out_8837167548770486690[3] = 0;
   out_8837167548770486690[4] = 0;
   out_8837167548770486690[5] = 0;
   out_8837167548770486690[6] = 0;
   out_8837167548770486690[7] = 0;
   out_8837167548770486690[8] = 0;
   out_8837167548770486690[9] = 0;
   out_8837167548770486690[10] = 1;
   out_8837167548770486690[11] = 0;
   out_8837167548770486690[12] = 0;
   out_8837167548770486690[13] = 0;
   out_8837167548770486690[14] = 0;
   out_8837167548770486690[15] = 0;
   out_8837167548770486690[16] = 0;
   out_8837167548770486690[17] = 0;
   out_8837167548770486690[18] = 0;
   out_8837167548770486690[19] = 0;
   out_8837167548770486690[20] = 1;
   out_8837167548770486690[21] = 0;
   out_8837167548770486690[22] = 0;
   out_8837167548770486690[23] = 0;
   out_8837167548770486690[24] = 0;
   out_8837167548770486690[25] = 0;
   out_8837167548770486690[26] = 0;
   out_8837167548770486690[27] = 0;
   out_8837167548770486690[28] = 0;
   out_8837167548770486690[29] = 0;
   out_8837167548770486690[30] = 1;
   out_8837167548770486690[31] = 0;
   out_8837167548770486690[32] = 0;
   out_8837167548770486690[33] = 0;
   out_8837167548770486690[34] = 0;
   out_8837167548770486690[35] = 0;
   out_8837167548770486690[36] = 0;
   out_8837167548770486690[37] = 0;
   out_8837167548770486690[38] = 0;
   out_8837167548770486690[39] = 0;
   out_8837167548770486690[40] = 1;
   out_8837167548770486690[41] = 0;
   out_8837167548770486690[42] = 0;
   out_8837167548770486690[43] = 0;
   out_8837167548770486690[44] = 0;
   out_8837167548770486690[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_8837167548770486690[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_8837167548770486690[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8837167548770486690[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8837167548770486690[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_8837167548770486690[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_8837167548770486690[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_8837167548770486690[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_8837167548770486690[53] = -9.8000000000000007*dt;
   out_8837167548770486690[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_8837167548770486690[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_8837167548770486690[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8837167548770486690[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8837167548770486690[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_8837167548770486690[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_8837167548770486690[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_8837167548770486690[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8837167548770486690[62] = 0;
   out_8837167548770486690[63] = 0;
   out_8837167548770486690[64] = 0;
   out_8837167548770486690[65] = 0;
   out_8837167548770486690[66] = 0;
   out_8837167548770486690[67] = 0;
   out_8837167548770486690[68] = 0;
   out_8837167548770486690[69] = 0;
   out_8837167548770486690[70] = 1;
   out_8837167548770486690[71] = 0;
   out_8837167548770486690[72] = 0;
   out_8837167548770486690[73] = 0;
   out_8837167548770486690[74] = 0;
   out_8837167548770486690[75] = 0;
   out_8837167548770486690[76] = 0;
   out_8837167548770486690[77] = 0;
   out_8837167548770486690[78] = 0;
   out_8837167548770486690[79] = 0;
   out_8837167548770486690[80] = 1;
}
void h_25(double *state, double *unused, double *out_3414124567549359716) {
   out_3414124567549359716[0] = state[6];
}
void H_25(double *state, double *unused, double *out_9047075177532746431) {
   out_9047075177532746431[0] = 0;
   out_9047075177532746431[1] = 0;
   out_9047075177532746431[2] = 0;
   out_9047075177532746431[3] = 0;
   out_9047075177532746431[4] = 0;
   out_9047075177532746431[5] = 0;
   out_9047075177532746431[6] = 1;
   out_9047075177532746431[7] = 0;
   out_9047075177532746431[8] = 0;
}
void h_24(double *state, double *unused, double *out_2828390139702683428) {
   out_2828390139702683428[0] = state[4];
   out_2828390139702683428[1] = state[5];
}
void H_24(double *state, double *unused, double *out_2824097089585287084) {
   out_2824097089585287084[0] = 0;
   out_2824097089585287084[1] = 0;
   out_2824097089585287084[2] = 0;
   out_2824097089585287084[3] = 0;
   out_2824097089585287084[4] = 1;
   out_2824097089585287084[5] = 0;
   out_2824097089585287084[6] = 0;
   out_2824097089585287084[7] = 0;
   out_2824097089585287084[8] = 0;
   out_2824097089585287084[9] = 0;
   out_2824097089585287084[10] = 0;
   out_2824097089585287084[11] = 0;
   out_2824097089585287084[12] = 0;
   out_2824097089585287084[13] = 0;
   out_2824097089585287084[14] = 1;
   out_2824097089585287084[15] = 0;
   out_2824097089585287084[16] = 0;
   out_2824097089585287084[17] = 0;
}
void h_30(double *state, double *unused, double *out_3268604743328704745) {
   out_3268604743328704745[0] = state[4];
}
void H_30(double *state, double *unused, double *out_8917736230389506361) {
   out_8917736230389506361[0] = 0;
   out_8917736230389506361[1] = 0;
   out_8917736230389506361[2] = 0;
   out_8917736230389506361[3] = 0;
   out_8917736230389506361[4] = 1;
   out_8917736230389506361[5] = 0;
   out_8917736230389506361[6] = 0;
   out_8917736230389506361[7] = 0;
   out_8917736230389506361[8] = 0;
}
void h_26(double *state, double *unused, double *out_8270278010989550390) {
   out_8270278010989550390[0] = state[7];
}
void H_26(double *state, double *unused, double *out_5305571858658690207) {
   out_5305571858658690207[0] = 0;
   out_5305571858658690207[1] = 0;
   out_5305571858658690207[2] = 0;
   out_5305571858658690207[3] = 0;
   out_5305571858658690207[4] = 0;
   out_5305571858658690207[5] = 0;
   out_5305571858658690207[6] = 0;
   out_5305571858658690207[7] = 1;
   out_5305571858658690207[8] = 0;
}
void h_27(double *state, double *unused, double *out_8152270617056173646) {
   out_8152270617056173646[0] = state[3];
}
void H_27(double *state, double *unused, double *out_6742972918589081450) {
   out_6742972918589081450[0] = 0;
   out_6742972918589081450[1] = 0;
   out_6742972918589081450[2] = 0;
   out_6742972918589081450[3] = 1;
   out_6742972918589081450[4] = 0;
   out_6742972918589081450[5] = 0;
   out_6742972918589081450[6] = 0;
   out_6742972918589081450[7] = 0;
   out_6742972918589081450[8] = 0;
}
void h_29(double *state, double *unused, double *out_9026589657322442701) {
   out_9026589657322442701[0] = state[1];
}
void H_29(double *state, double *unused, double *out_5029610191719530417) {
   out_5029610191719530417[0] = 0;
   out_5029610191719530417[1] = 1;
   out_5029610191719530417[2] = 0;
   out_5029610191719530417[3] = 0;
   out_5029610191719530417[4] = 0;
   out_5029610191719530417[5] = 0;
   out_5029610191719530417[6] = 0;
   out_5029610191719530417[7] = 0;
   out_5029610191719530417[8] = 0;
}
void h_28(double *state, double *unused, double *out_4606594481212752947) {
   out_4606594481212752947[0] = state[5];
   out_4606594481212752947[1] = state[6];
}
void H_28(double *state, double *unused, double *out_1354096657815228801) {
   out_1354096657815228801[0] = 0;
   out_1354096657815228801[1] = 0;
   out_1354096657815228801[2] = 0;
   out_1354096657815228801[3] = 0;
   out_1354096657815228801[4] = 0;
   out_1354096657815228801[5] = 1;
   out_1354096657815228801[6] = 0;
   out_1354096657815228801[7] = 0;
   out_1354096657815228801[8] = 0;
   out_1354096657815228801[9] = 0;
   out_1354096657815228801[10] = 0;
   out_1354096657815228801[11] = 0;
   out_1354096657815228801[12] = 0;
   out_1354096657815228801[13] = 0;
   out_1354096657815228801[14] = 0;
   out_1354096657815228801[15] = 1;
   out_1354096657815228801[16] = 0;
   out_1354096657815228801[17] = 0;
}
void h_31(double *state, double *unused, double *out_7051574095803382066) {
   out_7051574095803382066[0] = state[8];
}
void H_31(double *state, double *unused, double *out_9077721139409706859) {
   out_9077721139409706859[0] = 0;
   out_9077721139409706859[1] = 0;
   out_9077721139409706859[2] = 0;
   out_9077721139409706859[3] = 0;
   out_9077721139409706859[4] = 0;
   out_9077721139409706859[5] = 0;
   out_9077721139409706859[6] = 0;
   out_9077721139409706859[7] = 0;
   out_9077721139409706859[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_7406943224577723449) {
  err_fun(nom_x, delta_x, out_7406943224577723449);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_6282749298357684830) {
  inv_err_fun(nom_x, true_x, out_6282749298357684830);
}
void car_H_mod_fun(double *state, double *out_4025963510183789805) {
  H_mod_fun(state, out_4025963510183789805);
}
void car_f_fun(double *state, double dt, double *out_2334097729121616737) {
  f_fun(state,  dt, out_2334097729121616737);
}
void car_F_fun(double *state, double dt, double *out_8837167548770486690) {
  F_fun(state,  dt, out_8837167548770486690);
}
void car_h_25(double *state, double *unused, double *out_3414124567549359716) {
  h_25(state, unused, out_3414124567549359716);
}
void car_H_25(double *state, double *unused, double *out_9047075177532746431) {
  H_25(state, unused, out_9047075177532746431);
}
void car_h_24(double *state, double *unused, double *out_2828390139702683428) {
  h_24(state, unused, out_2828390139702683428);
}
void car_H_24(double *state, double *unused, double *out_2824097089585287084) {
  H_24(state, unused, out_2824097089585287084);
}
void car_h_30(double *state, double *unused, double *out_3268604743328704745) {
  h_30(state, unused, out_3268604743328704745);
}
void car_H_30(double *state, double *unused, double *out_8917736230389506361) {
  H_30(state, unused, out_8917736230389506361);
}
void car_h_26(double *state, double *unused, double *out_8270278010989550390) {
  h_26(state, unused, out_8270278010989550390);
}
void car_H_26(double *state, double *unused, double *out_5305571858658690207) {
  H_26(state, unused, out_5305571858658690207);
}
void car_h_27(double *state, double *unused, double *out_8152270617056173646) {
  h_27(state, unused, out_8152270617056173646);
}
void car_H_27(double *state, double *unused, double *out_6742972918589081450) {
  H_27(state, unused, out_6742972918589081450);
}
void car_h_29(double *state, double *unused, double *out_9026589657322442701) {
  h_29(state, unused, out_9026589657322442701);
}
void car_H_29(double *state, double *unused, double *out_5029610191719530417) {
  H_29(state, unused, out_5029610191719530417);
}
void car_h_28(double *state, double *unused, double *out_4606594481212752947) {
  h_28(state, unused, out_4606594481212752947);
}
void car_H_28(double *state, double *unused, double *out_1354096657815228801) {
  H_28(state, unused, out_1354096657815228801);
}
void car_h_31(double *state, double *unused, double *out_7051574095803382066) {
  h_31(state, unused, out_7051574095803382066);
}
void car_H_31(double *state, double *unused, double *out_9077721139409706859) {
  H_31(state, unused, out_9077721139409706859);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_init(car);
