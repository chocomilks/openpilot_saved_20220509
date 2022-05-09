#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_7406943224577723449);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_6282749298357684830);
void car_H_mod_fun(double *state, double *out_4025963510183789805);
void car_f_fun(double *state, double dt, double *out_2334097729121616737);
void car_F_fun(double *state, double dt, double *out_8837167548770486690);
void car_h_25(double *state, double *unused, double *out_3414124567549359716);
void car_H_25(double *state, double *unused, double *out_9047075177532746431);
void car_h_24(double *state, double *unused, double *out_2828390139702683428);
void car_H_24(double *state, double *unused, double *out_2824097089585287084);
void car_h_30(double *state, double *unused, double *out_3268604743328704745);
void car_H_30(double *state, double *unused, double *out_8917736230389506361);
void car_h_26(double *state, double *unused, double *out_8270278010989550390);
void car_H_26(double *state, double *unused, double *out_5305571858658690207);
void car_h_27(double *state, double *unused, double *out_8152270617056173646);
void car_H_27(double *state, double *unused, double *out_6742972918589081450);
void car_h_29(double *state, double *unused, double *out_9026589657322442701);
void car_H_29(double *state, double *unused, double *out_5029610191719530417);
void car_h_28(double *state, double *unused, double *out_4606594481212752947);
void car_H_28(double *state, double *unused, double *out_1354096657815228801);
void car_h_31(double *state, double *unused, double *out_7051574095803382066);
void car_H_31(double *state, double *unused, double *out_9077721139409706859);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}