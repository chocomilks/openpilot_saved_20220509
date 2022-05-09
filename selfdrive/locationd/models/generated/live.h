#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_3075168618539510648);
void live_err_fun(double *nom_x, double *delta_x, double *out_1090245281905906660);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_8851414583690127304);
void live_H_mod_fun(double *state, double *out_2951829955177005887);
void live_f_fun(double *state, double dt, double *out_3814539948165517167);
void live_F_fun(double *state, double dt, double *out_7875956373086618253);
void live_h_4(double *state, double *unused, double *out_308485439684569826);
void live_H_4(double *state, double *unused, double *out_448820613655623243);
void live_h_9(double *state, double *unused, double *out_7222384130018568816);
void live_H_9(double *state, double *unused, double *out_2440040938624456099);
void live_h_10(double *state, double *unused, double *out_7807615225781965650);
void live_H_10(double *state, double *unused, double *out_106603328486207745);
void live_h_12(double *state, double *unused, double *out_8240029641751488420);
void live_H_12(double *state, double *unused, double *out_7218307700026827249);
void live_h_31(double *state, double *unused, double *out_4739257077166324131);
void live_H_31(double *state, double *unused, double *out_8482873341357710658);
void live_h_32(double *state, double *unused, double *out_2545771817069536442);
void live_H_32(double *state, double *unused, double *out_8795495036245155686);
void live_h_13(double *state, double *unused, double *out_8364010069648976723);
void live_H_13(double *state, double *unused, double *out_2336756164620082167);
void live_h_14(double *state, double *unused, double *out_7222384130018568816);
void live_H_14(double *state, double *unused, double *out_2440040938624456099);
void live_h_33(double *state, double *unused, double *out_4663909138174759435);
void live_H_33(double *state, double *unused, double *out_5332316336718853054);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}