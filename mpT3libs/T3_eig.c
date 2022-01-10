/* Include Files */
#include "rt_nonfinite.h"
#include "T3_eig.h"
#include "sign.h"
#include "mpower.h"
#include "norm.h"
#include "sortElvs.h"
#include "T3_eig_rtwutil.h"

/* Function Definitions */

void T3_eig(float *H11_in, float *H12_in, float *H13_in,
            float *H22_in, float *H23_in, float *H33_in,
            float *El1, float *El2, float *El3,
            float *Ev11, float *Ev12, float *Ev13,
            float *Ev21, float *Ev22, float *Ev23,
            float *Ev31, float *Ev32, float *Ev33)
{
  float H11,H12,H13,H22,H23,H33;
  H11 = *H11_in;
  H12 = *H12_in;
  H13 = *H13_in;
  H22 = *H22_in;
  H23 = *H23_in;
  H33 = *H33_in;

  creal32_T T13;
  creal32_T T23;
  creal32_T T11;
  creal32_T fc0;
  float re;
  float im;
  float T13_re;
  float T13_im;
  float b_T13_re;
  float T12_re;
  float T12_im;
  float b_T12_re;
  float T23_re;
  float T23_im;
  float T11_re;
  float T11_im;
  float b_T11_re;
  float b_T12_im;
  float c_T11_re;
  float b_T11_im;
  float c_T11_im;
  float d_T11_re;
  float a_re;
  float a_im;
  float b_a_re;
  float c_a_re;
  creal32_T b_T11;
  float d_a_re;
  float b_a_im;
  float e_a_re;
  float T22_re;
  float T33_re;
  float f_a_re;
  float c_a_im;
  float b_T22_re;
  float T22_im;
  creal32_T a;
  creal32_T b_a;
  creal32_T c_T11;
  float e_T11_re;
  creal32_T fc1;
  float b_T13_im;
  float b_T23_re;
  float b_T23_im;
  float c_T22_re;
  float b_T22_im;
  creal32_T c_a;
  creal32_T d_a;
  float c_T13_re;
  float c_T12_re;
  creal32_T d_T11;
  float d_T11_im;
  float f_T11_re;
  creal32_T fc2;
  float V1_norm;
  float s;
  creal32_T e_T11;
  creal32_T f_T11;
  float v21_re;
  float g_a_re;
  creal32_T e_a;
  creal32_T f_a;
  creal32_T g_T11;
  creal32_T fc3;
  creal32_T g_a;
  creal32_T h_a;
  creal32_T h_T11;
  creal32_T fc4;
  float V11_re;
  float V11_im;
  creal32_T i_T11;
  creal32_T j_T11;
  creal32_T k_T11;
  float h_a_re;
  creal32_T l_T11;
  float i_a_re;
  float d_a_im;
  float j_a_re;
  float e_a_im;
  creal32_T i_a;
  creal32_T j_a;
  creal32_T m_T11;
  creal32_T fc5;
  creal32_T k_a;
  creal32_T l_a;
  creal32_T n_T11;
  creal32_T fc6;
  float b_re;
  creal32_T m_a;
  creal32_T n_a;
  creal32_T o_T11;
  float b_im;
  creal32_T fc7;
  creal32_T o_a;
  creal32_T p_a;
  creal32_T p_T11;
  creal32_T fc8;
  creal32_T q_T11;
  creal32_T r_T11;
  creal32_T s_T11;
  creal32_T t_T11;
  float V23_re;
  creal32_T q_a;
  creal32_T r_a;
  creal32_T u_T11;
  creal32_T fc9;
  creal32_T s_a;
  creal32_T t_a;
  creal32_T v_T11;
  creal32_T fc10;
  creal32_T u_a;
  creal32_T v_a;
  creal32_T w_T11;
  creal32_T fc11;
  creal32_T V23;
  creal32_T w_a;
  creal32_T x_T11;
  creal32_T fc12;
  float V12_re;
  float V12_im;
  creal32_T y_T11;
  creal32_T ab_T11;
  creal32_T bb_T11;
  creal32_T cb_T11;
  creal32_T x_a;
  creal32_T y_a;
  creal32_T db_T11;
  creal32_T fc13;
  creal32_T ab_a;
  creal32_T bb_a;
  creal32_T eb_T11;
  creal32_T fc14;
  creal32_T cb_a;
  creal32_T db_a;
  creal32_T fb_T11;
  creal32_T fc15;
  creal32_T eb_a;
  creal32_T fb_a;
  creal32_T gb_T11;
  creal32_T fc16;
  creal32_T hb_T11;
  creal32_T ib_T11;
  creal32_T jb_T11;
  creal32_T kb_T11;
  creal32_T gb_a;
  creal32_T hb_a;
  creal32_T lb_T11;
  creal32_T fc17;
  creal32_T ib_a;
  creal32_T jb_a;
  creal32_T mb_T11;
  creal32_T fc18;
  creal32_T kb_a;
  creal32_T lb_a;
  creal32_T nb_T11;
  creal32_T fc19;
  creal32_T b_V23;
  creal32_T mb_a;
  creal32_T ob_T11;
  creal32_T fc20;
  float V13_re;
  float V13_im;
  creal32_T pb_T11;
  creal32_T qb_T11;
  creal32_T nb_a;
  creal32_T ob_a;
  creal32_T rb_T11;
  creal32_T fc21;
  creal32_T pb_a;
  creal32_T qb_a;
  creal32_T sb_T11;
  creal32_T fc22;
  creal32_T tb_T11;
  creal32_T ub_T11;
  creal32_T rb_a;
  creal32_T sb_a;
  creal32_T vb_T11;
  creal32_T fc23;
  creal32_T tb_a;
  creal32_T ub_a;
  creal32_T wb_T11;
  creal32_T fc24;
  float V21_re;
  float V21_im;
  creal32_T xb_T11;
  creal32_T yb_T11;
  creal32_T ac_T11;
  creal32_T bc_T11;
  creal32_T vb_a;
  creal32_T wb_a;
  creal32_T cc_T11;
  creal32_T fc25;
  creal32_T xb_a;
  creal32_T yb_a;
  creal32_T dc_T11;
  creal32_T fc26;
  creal32_T ac_a;
  creal32_T bc_a;
  creal32_T ec_T11;
  creal32_T fc27;
  creal32_T cc_a;
  creal32_T dc_a;
  creal32_T fc_T11;
  creal32_T fc28;
  creal32_T gc_T11;
  creal32_T hc_T11;
  creal32_T ic_T11;
  creal32_T jc_T11;
  creal32_T ec_a;
  creal32_T fc_a;
  creal32_T kc_T11;
  creal32_T fc29;
  creal32_T gc_a;
  creal32_T hc_a;
  creal32_T lc_T11;
  creal32_T fc30;
  creal32_T ic_a;
  creal32_T jc_a;
  creal32_T mc_T11;
  creal32_T fc31;
  creal32_T c_V23;
  creal32_T kc_a;
  creal32_T nc_T11;
  creal32_T fc32;
  float V22_re;
  float V22_im;
  creal32_T oc_T11;
  creal32_T pc_T11;
  creal32_T qc_T11;
  creal32_T rc_T11;
  creal32_T lc_a;
  creal32_T mc_a;
  creal32_T sc_T11;
  creal32_T fc33;
  creal32_T nc_a;
  creal32_T oc_a;
  creal32_T tc_T11;
  creal32_T fc34;
  creal32_T pc_a;
  creal32_T qc_a;
  creal32_T uc_T11;
  creal32_T fc35;
  creal32_T rc_a;
  creal32_T sc_a;
  creal32_T vc_T11;
  creal32_T fc36;
  creal32_T wc_T11;
  creal32_T xc_T11;
  creal32_T yc_T11;
  creal32_T ad_T11;
  creal32_T tc_a;
  creal32_T uc_a;
  creal32_T bd_T11;
  creal32_T fc37;
  creal32_T vc_a;
  creal32_T wc_a;
  creal32_T cd_T11;
  creal32_T fc38;
  creal32_T xc_a;
  creal32_T yc_a;
  creal32_T dd_T11;
  creal32_T fc39;
  creal32_T d_V23;
  creal32_T ad_a;
  creal32_T ed_T11;
  creal32_T fc40;
  creal32_T V[9];
  int i0;
  float v11_re;
  float v11_im;
  float v12_re;
  float v12_im;
  float v13_re;
  float v22_re;
  float v22_im;
  float v23_re;
  creal32_T fd_T11;
  creal32_T gd_T11;
  creal32_T bd_a;
  creal32_T cd_a;
  creal32_T hd_T11;
  creal32_T fc41;
  creal32_T dd_a;
  creal32_T ed_a;
  creal32_T id_T11;
  creal32_T fc42;
  creal32_T jd_T11;
  creal32_T kd_T11;
  creal32_T ld_T11;
  creal32_T md_T11;
  creal32_T fd_a;
  creal32_T gd_a;
  creal32_T nd_T11;
  creal32_T fc43;
  creal32_T hd_a;
  creal32_T id_a;
  creal32_T od_T11;
  creal32_T fc44;
  creal32_T jd_a;
  creal32_T kd_a;
  creal32_T pd_T11;
  creal32_T fc45;
  creal32_T ld_a;
  creal32_T md_a;
  creal32_T qd_T11;
  creal32_T fc46;
  creal32_T rd_T11;
  creal32_T sd_T11;
  creal32_T td_T11;
  creal32_T ud_T11;
  creal32_T nd_a;
  creal32_T od_a;
  creal32_T vd_T11;
  creal32_T fc47;
  creal32_T pd_a;
  creal32_T qd_a;
  creal32_T wd_T11;
  creal32_T fc48;
  creal32_T rd_a;
  creal32_T sd_a;
  creal32_T xd_T11;
  creal32_T fc49;
  creal32_T td_a;
  creal32_T ud_a;
  creal32_T yd_T11;
  creal32_T fc50;
  float b_V23_re;
  float b_V12_re;
  float b_V22_re;
  float b_V13_re;
  float b_V11_re;
  float b_v21_re;
  float b_V21_re;
  float b_v22_re;
  float b_v11_re;
  float b_v23_re;
  float b_v12_re;
  float c_v11_re;
  float b_v13_re;
  float c_v12_re;
  float c_v21_re;
  float c_v13_re;
  float c_v22_re;
  float c_v23_re;
  T13.re = H13;
  T13.im = 0.0F;
  T23.re = H23;
  T23.im = 0.0F;
  T11.re = (H11 + H22) + H33;
  T11.im = 0.0F;
  fc0 = mpower(T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
           b_T12_im * 0.0F)) + c_T11_re;
  a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F +
           b_T12_im * H23)) + c_T11_im;
  b_a_re = (H11 + H22) + H33;
  c_a_re = (H11 + H22) + H33;
  b_T11.re = (H11 + H22) + H33;
  b_T11.im = 0.0F;
  fc0 = mpower(b_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  d_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  b_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  e_a_re = (H11 + H22) + H33;
  T11_re = H11 / 3.0F;
  T22_re = H22 / 3.0F;
  T33_re = H33 / 3.0F;
  f_a_re = b_a_re * b_a_re;
  c_a_im = b_a_re * 0.0F + 0.0F * b_a_re;
  if (c_a_im == 0.0F) {
    b_a_re = f_a_re / 9.0F;
    c_a_im = 0.0F;
  } else if (f_a_re == 0.0F) {
    b_a_re = 0.0F;
    c_a_im /= 9.0F;
  } else {
    b_a_re = f_a_re / 9.0F;
    c_a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  if (T13_im == 0.0F) {
    T13_re /= 3.0F;
    T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  a.re = (((((b_a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  a.im = (((((c_a_im + T12_im) + T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  fc0 = mpower(a);
  b_a.re = (a_re * a_re - a_im * a_im) - fc0.re;
  b_a.im = (a_re * a_im + a_im * a_re) - fc0.im;
  fc0 = b_mpower(b_a);
  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  b_T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    b_T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    b_T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  c_T11.re = (H11 + H22) + H33;
  c_T11.im = 0.0F;
  T11 = mpower(c_T11);
  if (T11.im == 0.0F) {
    re = T11.re / 27.0F;
    im = 0.0F;
  } else if (T11.re == 0.0F) {
    re = 0.0F;
    im = T11.im / 27.0F;
  } else {
    re = T11.re / 27.0F;
    im = T11.im / 27.0F;
  }

  c_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  d_T11_re = c_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = c_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  e_T11_re = d_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = d_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (e_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  fc1.re = ((((((fc0.re - T13_re) - T12_re) - b_T11_re) + re) + c_T11_re) +
            (b_T12_re * H23 - b_T12_im * 0.0F)) + d_T11_re;
  fc1.im = ((((((fc0.im - T13_im) - T12_im) - T11_im) + im) + b_T11_im) +
            (b_T12_re * 0.0F + b_T12_im * H23)) + c_T11_im;
  fc0 = c_mpower(fc1);
  a_re = c_a_re * c_a_re;
  a_im = c_a_re * 0.0F + 0.0F * c_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  if (T13_im == 0.0F) {
    T13_re /= 3.0F;
    T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  b_a_re = e_a_re * e_a_re;
  c_a_im = e_a_re * 0.0F + 0.0F * e_a_re;
  if (c_a_im == 0.0F) {
    b_a_re /= 9.0F;
    c_a_im = 0.0F;
  } else if (b_a_re == 0.0F) {
    b_a_re = 0.0F;
    c_a_im /= 9.0F;
  } else {
    b_a_re /= 9.0F;
    c_a_im /= 9.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  if (b_T12_im == 0.0F) {
    b_T12_re /= 3.0F;
    b_T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 3.0F;
  } else {
    b_T12_re /= 3.0F;
    b_T12_im /= 3.0F;
  }

  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  if (b_T13_im == 0.0F) {
    b_T13_re /= 3.0F;
    b_T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 3.0F;
  } else {
    b_T13_re /= 3.0F;
    b_T13_im /= 3.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  if (b_T23_im == 0.0F) {
    b_T23_re /= 3.0F;
    b_T23_im = 0.0F;
  } else if (b_T23_re == 0.0F) {
    b_T23_re = 0.0F;
    b_T23_im /= 3.0F;
  } else {
    b_T23_re /= 3.0F;
    b_T23_im /= 3.0F;
  }

  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 3.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 3.0F;
  } else {
    d_T11_re /= 3.0F;
    c_T11_im /= 3.0F;
  }

  e_T11_re = H11 * H33;
  c_a_re = H11 * 0.0F + 0.0F * H33;
  if (c_a_re == 0.0F) {
    e_T11_re /= 3.0F;
    c_a_re = 0.0F;
  } else if (e_T11_re == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 3.0F;
  } else {
    e_T11_re /= 3.0F;
    c_a_re /= 3.0F;
  }

  c_T22_re = H22 * H33;
  b_T22_im = H22 * 0.0F + 0.0F * H33;
  if (b_T22_im == 0.0F) {
    c_T22_re /= 3.0F;
    b_T22_im = 0.0F;
  } else if (c_T22_re == 0.0F) {
    c_T22_re = 0.0F;
    b_T22_im /= 3.0F;
  } else {
    c_T22_re /= 3.0F;
    b_T22_im /= 3.0F;
  }

  c_a.re = (((((b_a_re + b_T12_re) + b_T13_re) + b_T23_re) - d_T11_re) -
            e_T11_re) - c_T22_re;
  c_a.im = (((((c_a_im + b_T12_im) + b_T13_im) + b_T23_im) - c_T11_im) - c_a_re)
    - b_T22_im;
  T11 = mpower(c_a);
  d_a.re = (d_a_re * d_a_re - b_a_im * b_a_im) - T11.re;
  d_a.im = (d_a_re * b_a_im + b_a_im * d_a_re) - T11.im;
  T11 = b_mpower(d_a);
  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  c_T13_re = b_T13_re * H22 - b_T13_im * 0.0F;
  b_T13_im = b_T13_re * 0.0F + b_T13_im * H22;
  if (b_T13_im == 0.0F) {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im = 0.0F;
  } else if (c_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 2.0F;
  } else {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im /= 2.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  c_T12_re = b_T12_re * H33 - b_T12_im * 0.0F;
  b_T12_im = b_T12_re * 0.0F + b_T12_im * H33;
  if (b_T12_im == 0.0F) {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im = 0.0F;
  } else if (c_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 2.0F;
  } else {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im /= 2.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  d_T11_re = H11 * b_T23_re - 0.0F * b_T23_im;
  c_T11_im = H11 * b_T23_im + 0.0F * b_T23_re;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re /= 2.0F;
    c_T11_im /= 2.0F;
  }

  d_T11.re = (H11 + H22) + H33;
  d_T11.im = 0.0F;
  fc1 = mpower(d_T11);
  if (fc1.im == 0.0F) {
    re = fc1.re / 27.0F;
    im = 0.0F;
  } else if (fc1.re == 0.0F) {
    re = 0.0F;
    im = fc1.im / 27.0F;
  } else {
    re = fc1.re / 27.0F;
    im = fc1.im / 27.0F;
  }

  e_T11_re = (H11 + H22) + H33;
  c_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  c_T22_re = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  b_T22_im = e_T11_re * c_T12_re - 0.0F * c_T22_re;
  c_a_re = e_T11_re * c_T22_re + 0.0F * c_T12_re;
  if (c_a_re == 0.0F) {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re = 0.0F;
  } else if (b_T22_im == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 6.0F;
  } else {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re /= 6.0F;
  }

  c_T12_re = H12 * H13;
  c_T22_re = H12 * 0.0F + 0.0F * H13;
  b_T22_im = H11 * H22;
  d_T11_im = H11 * 0.0F + 0.0F * H22;
  f_T11_re = b_T22_im * H33 - d_T11_im * 0.0F;
  d_T11_im = b_T22_im * 0.0F + d_T11_im * H33;
  if (d_T11_im == 0.0F) {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im = 0.0F;
  } else if (f_T11_re == 0.0F) {
    b_T22_im = 0.0F;
    d_T11_im /= 2.0F;
  } else {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im /= 2.0F;
  }

  fc2.re = ((((((T11.re - b_T13_re) - b_T12_re) - d_T11_re) + re) + e_T11_re) +
            (c_T12_re * H23 - c_T22_re * 0.0F)) + b_T22_im;
  fc2.im = ((((((T11.im - b_T13_im) - b_T12_im) - c_T11_im) + im) + c_a_re) +
            (c_T12_re * 0.0F + c_T22_re * H23)) + d_T11_im;
  T11 = c_mpower(fc2);
  f_T11_re = (((((a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  b_T22_im = (((((a_im + T12_im) + T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  if (T11.im == 0.0F) {
    if (b_T22_im == 0.0F) {
      a_re = f_T11_re / T11.re;
      a_im = 0.0F;
    } else if (f_T11_re == 0.0F) {
      a_re = 0.0F;
      a_im = b_T22_im / T11.re;
    } else {
      a_re = f_T11_re / T11.re;
      a_im = b_T22_im / T11.re;
    }
  } else if (T11.re == 0.0F) {
    if (f_T11_re == 0.0F) {
      a_re = b_T22_im / T11.im;
      a_im = 0.0F;
    } else if (b_T22_im == 0.0F) {
      a_re = 0.0F;
      a_im = -(f_T11_re / T11.im);
    } else {
      a_re = b_T22_im / T11.im;
      a_im = -(f_T11_re / T11.im);
    }
  } else {
    c_T12_re = (float)fabs(T11.re);
    V1_norm = (float)fabs(T11.im);
    if (c_T12_re > V1_norm) {
      s = T11.im / T11.re;
      V1_norm = T11.re + s * T11.im;
      a_re = (f_T11_re + s * b_T22_im) / V1_norm;
      a_im = (b_T22_im - s * f_T11_re) / V1_norm;
    } else if (V1_norm == c_T12_re) {
      if (T11.re > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (T11.im > 0.0F) {
        V1_norm = 0.5F;
      } else {
        V1_norm = -0.5F;
      }

      a_re = (f_T11_re * s + b_T22_im * V1_norm) / c_T12_re;
      a_im = (b_T22_im * s - f_T11_re * V1_norm) / c_T12_re;
    } else {
      s = T11.re / T11.im;
      V1_norm = T11.im + s * T11.re;
      a_re = (s * f_T11_re + b_T22_im) / V1_norm;
      a_im = (s * b_T22_im - f_T11_re) / V1_norm;
    }
  }

  a_re += ((T11_re + T22_re) + T33_re) + fc0.re;
  a_im += fc0.im;
  e_T11.re = (H11 + H22) + H33;
  e_T11.im = 0.0F;
  fc0 = mpower(e_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  b_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  b_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  c_a_re = (H11 + H22) + H33;
  d_a_re = (H11 + H22) + H33;
  f_T11.re = (H11 + H22) + H33;
  f_T11.im = 0.0F;
  fc0 = mpower(f_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  e_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  c_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  f_a_re = (H11 + H22) + H33;
  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  T12_re = H12 * H33;
  T12_im = H12 * 0.0F + 0.0F * H33;
  fc0 = mpower(T23);
  T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * H13;
  T11_im = H11 * 0.0F + 0.0F * H13;
  c_T13_re = H13 * H22;
  v21_re = H13 * 0.0F + 0.0F * H22;
  f_T11_re = (((T13_re * H23 - T13_im * 0.0F) - (T12_re * H13 - T12_im * 0.0F))
              + fc0.re) - (T22_re * H23 - T22_im * 0.0F);
  b_T22_im = (((T13_re * 0.0F + T13_im * H23) - (T12_re * 0.0F + T12_im * H13))
              + fc0.im) - (T22_re * 0.0F + T22_im * H23);
  d_T11_im = (((H12 * b_T13_re - 0.0F * b_T13_im) - (H12 * T23_re - 0.0F *
    T23_im)) - (T11_re * H23 - T11_im * 0.0F)) + (c_T13_re * H23 - v21_re * 0.0F);
  c_T22_re = (((H12 * b_T13_im + 0.0F * b_T13_re) - (H12 * T23_im + 0.0F *
    T23_re)) - (T11_re * 0.0F + T11_im * H23)) + (c_T13_re * 0.0F + v21_re * H23);
  if (c_T22_re == 0.0F) {
    if (b_T22_im == 0.0F) {
      T13_re = f_T11_re / d_T11_im;
      T13_im = 0.0F;
    } else if (f_T11_re == 0.0F) {
      T13_re = 0.0F;
      T13_im = b_T22_im / d_T11_im;
    } else {
      T13_re = f_T11_re / d_T11_im;
      T13_im = b_T22_im / d_T11_im;
    }
  } else if (d_T11_im == 0.0F) {
    if (f_T11_re == 0.0F) {
      T13_re = b_T22_im / c_T22_re;
      T13_im = 0.0F;
    } else if (b_T22_im == 0.0F) {
      T13_re = 0.0F;
      T13_im = -(f_T11_re / c_T22_re);
    } else {
      T13_re = b_T22_im / c_T22_re;
      T13_im = -(f_T11_re / c_T22_re);
    }
  } else {
    c_T12_re = (float)fabs(d_T11_im);
    if (c_T12_re > c_T22_re) {
      s = c_T22_re / d_T11_im;
      V1_norm = d_T11_im + s * c_T22_re;
      T13_re = (f_T11_re + s * b_T22_im) / V1_norm;
      T13_im = (b_T22_im - s * f_T11_re) / V1_norm;
    } else if (c_T22_re == c_T12_re) {
      if (d_T11_im > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (f_T11_re * s + b_T22_im * -0.5F > 0.0F) {
        T13_re = ((real32_T)rtInf);
      } else if (f_T11_re * s + b_T22_im * -0.5F < 0.0F) {
        T13_re = ((real32_T)rtMinusInf);
      } else {
        T13_re = ((real32_T)rtNaN);
      }

      if (b_T22_im * s - f_T11_re * -0.5F > 0.0F) {
        T13_im = ((real32_T)rtInf);
      } else if (b_T22_im * s - f_T11_re * -0.5F < 0.0F) {
        T13_im = ((real32_T)rtMinusInf);
      } else {
        T13_im = ((real32_T)rtNaN);
      }
    } else {
      s = d_T11_im / c_T22_re;
      V1_norm = c_T22_re + s * d_T11_im;
      T13_re = (s * f_T11_re + b_T22_im) / V1_norm;
      T13_im = (s * b_T22_im - f_T11_re) / V1_norm;
    }
  }

  g_a_re = a_re * a_re - a_im * a_im;
  a_im = a_re * a_im + a_im * a_re;
  T23_re = H23 * g_a_re - 0.0F * a_im;
  T23_im = H23 * a_im + 0.0F * g_a_re;
  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * H13;
  T11_im = H11 * 0.0F + 0.0F * H13;
  c_T13_re = H13 * H22;
  v21_re = H13 * 0.0F + 0.0F * H22;
  d_T11_im = (((H12 * b_T13_re - 0.0F * b_T13_im) - (H12 * b_T23_re - 0.0F *
    b_T23_im)) - (T11_re * H23 - T11_im * 0.0F)) + (c_T13_re * H23 - v21_re *
    0.0F);
  c_T22_re = (((H12 * b_T13_im + 0.0F * b_T13_re) - (H12 * b_T23_im + 0.0F *
    b_T23_re)) - (T11_re * 0.0F + T11_im * H23)) + (c_T13_re * 0.0F + v21_re *
    H23);
  if (c_T22_re == 0.0F) {
    if (T23_im == 0.0F) {
      b_T23_re = T23_re / d_T11_im;
      T23_im = 0.0F;
    } else if (T23_re == 0.0F) {
      b_T23_re = 0.0F;
      T23_im /= d_T11_im;
    } else {
      b_T23_re = T23_re / d_T11_im;
      T23_im /= d_T11_im;
    }
  } else if (d_T11_im == 0.0F) {
    if (T23_re == 0.0F) {
      b_T23_re = T23_im / c_T22_re;
      T23_im = 0.0F;
    } else if (T23_im == 0.0F) {
      b_T23_re = 0.0F;
      T23_im = -(T23_re / c_T22_re);
    } else {
      b_T23_re = T23_im / c_T22_re;
      T23_im = -(T23_re / c_T22_re);
    }
  } else {
    c_T12_re = (float)fabs(d_T11_im);
    if (c_T12_re > c_T22_re) {
      s = c_T22_re / d_T11_im;
      V1_norm = d_T11_im + s * c_T22_re;
      b_T23_re = (T23_re + s * T23_im) / V1_norm;
      T23_im = (T23_im - s * T23_re) / V1_norm;
    } else if (c_T22_re == c_T12_re) {
      if (d_T11_im > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (T23_re * s + T23_im * -0.5F > 0.0F) {
        b_T23_re = ((real32_T)rtInf);
      } else if (T23_re * s + T23_im * -0.5F < 0.0F) {
        b_T23_re = ((real32_T)rtMinusInf);
      } else {
        b_T23_re = ((real32_T)rtNaN);
      }

      if (T23_im * s - T23_re * -0.5F > 0.0F) {
        T23_im = ((real32_T)rtInf);
      } else if (T23_im * s - T23_re * -0.5F < 0.0F) {
        T23_im = ((real32_T)rtMinusInf);
      } else {
        T23_im = ((real32_T)rtNaN);
      }
    } else {
      s = d_T11_im / c_T22_re;
      V1_norm = c_T22_re + s * d_T11_im;
      b_T23_re = (s * T23_re + T23_im) / V1_norm;
      T23_im = (s * T23_im - T23_re) / V1_norm;
    }
  }

  T11_re = H11 / 3.0F;
  T22_re = H22 / 3.0F;
  T33_re = H33 / 3.0F;
  a_re = c_a_re * c_a_re;
  a_im = c_a_re * 0.0F + 0.0F * c_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  if (b_T13_im == 0.0F) {
    b_T13_re /= 3.0F;
    b_T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 3.0F;
  } else {
    b_T13_re /= 3.0F;
    b_T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  if (b_T23_im == 0.0F) {
    T23_re /= 3.0F;
    b_T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    b_T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    b_T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  e_a.re = (((((a_re + T12_re) + b_T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  e_a.im = (((((a_im + T12_im) + b_T13_im) + b_T23_im) - T11_im) - b_T11_im) -
    T22_im;
  fc0 = mpower(e_a);
  f_a.re = (b_a_re * b_a_re - b_a_im * b_a_im) - fc0.re;
  f_a.im = (b_a_re * b_a_im + b_a_im * b_a_re) - fc0.im;
  fc0 = b_mpower(f_a);
  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  c_T13_re = b_T13_re * H22 - b_T13_im * 0.0F;
  b_T13_im = b_T13_re * 0.0F + b_T13_im * H22;
  if (b_T13_im == 0.0F) {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im = 0.0F;
  } else if (c_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 2.0F;
  } else {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  b_T11_re = H11 * T23_re - 0.0F * b_T23_im;
  T11_im = H11 * b_T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    b_T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    b_T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  g_T11.re = (H11 + H22) + H33;
  g_T11.im = 0.0F;
  T11 = mpower(g_T11);
  if (T11.im == 0.0F) {
    re = T11.re / 27.0F;
    im = 0.0F;
  } else if (T11.re == 0.0F) {
    re = 0.0F;
    im = T11.im / 27.0F;
  } else {
    re = T11.re / 27.0F;
    im = T11.im / 27.0F;
  }

  c_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  d_T11_re = c_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = c_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  e_T11_re = d_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = d_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (e_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  fc3.re = ((((((fc0.re - b_T13_re) - T12_re) - b_T11_re) + re) + c_T11_re) +
            (b_T12_re * H23 - b_T12_im * 0.0F)) + d_T11_re;
  fc3.im = ((((((fc0.im - b_T13_im) - T12_im) - T11_im) + im) + b_T11_im) +
            (b_T12_re * 0.0F + b_T12_im * H23)) + c_T11_im;
  fc0 = c_mpower(fc3);
  a_re = d_a_re * d_a_re;
  a_im = d_a_re * 0.0F + 0.0F * d_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  if (b_T13_im == 0.0F) {
    b_T13_re /= 3.0F;
    b_T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 3.0F;
  } else {
    b_T13_re /= 3.0F;
    b_T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  if (b_T23_im == 0.0F) {
    T23_re /= 3.0F;
    b_T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    b_T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    b_T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  b_a_re = f_a_re * f_a_re;
  b_a_im = f_a_re * 0.0F + 0.0F * f_a_re;
  if (b_a_im == 0.0F) {
    b_a_re /= 9.0F;
    b_a_im = 0.0F;
  } else if (b_a_re == 0.0F) {
    b_a_re = 0.0F;
    b_a_im /= 9.0F;
  } else {
    b_a_re /= 9.0F;
    b_a_im /= 9.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  if (b_T12_im == 0.0F) {
    b_T12_re /= 3.0F;
    b_T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 3.0F;
  } else {
    b_T12_re /= 3.0F;
    b_T12_im /= 3.0F;
  }

  c_T13_re = H13 * H13;
  v21_re = H13 * 0.0F + 0.0F * H13;
  if (v21_re == 0.0F) {
    c_T13_re /= 3.0F;
    v21_re = 0.0F;
  } else if (c_T13_re == 0.0F) {
    c_T13_re = 0.0F;
    v21_re /= 3.0F;
  } else {
    c_T13_re /= 3.0F;
    v21_re /= 3.0F;
  }

  V1_norm = H23 * H23;
  s = H23 * 0.0F + 0.0F * H23;
  if (s == 0.0F) {
    V1_norm /= 3.0F;
    s = 0.0F;
  } else if (V1_norm == 0.0F) {
    V1_norm = 0.0F;
    s /= 3.0F;
  } else {
    V1_norm /= 3.0F;
    s /= 3.0F;
  }

  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 3.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 3.0F;
  } else {
    d_T11_re /= 3.0F;
    c_T11_im /= 3.0F;
  }

  e_T11_re = H11 * H33;
  c_a_re = H11 * 0.0F + 0.0F * H33;
  if (c_a_re == 0.0F) {
    e_T11_re /= 3.0F;
    c_a_re = 0.0F;
  } else if (e_T11_re == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 3.0F;
  } else {
    e_T11_re /= 3.0F;
    c_a_re /= 3.0F;
  }

  c_T22_re = H22 * H33;
  b_T22_im = H22 * 0.0F + 0.0F * H33;
  if (b_T22_im == 0.0F) {
    c_T22_re /= 3.0F;
    b_T22_im = 0.0F;
  } else if (c_T22_re == 0.0F) {
    c_T22_re = 0.0F;
    b_T22_im /= 3.0F;
  } else {
    c_T22_re /= 3.0F;
    b_T22_im /= 3.0F;
  }

  g_a.re = (((((b_a_re + b_T12_re) + c_T13_re) + V1_norm) - d_T11_re) - e_T11_re)
    - c_T22_re;
  g_a.im = (((((b_a_im + b_T12_im) + v21_re) + s) - c_T11_im) - c_a_re) -
    b_T22_im;
  T11 = mpower(g_a);
  h_a.re = (e_a_re * e_a_re - c_a_im * c_a_im) - T11.re;
  h_a.im = (e_a_re * c_a_im + c_a_im * e_a_re) - T11.im;
  T11 = b_mpower(h_a);
  c_T13_re = H13 * H13;
  v21_re = H13 * 0.0F + 0.0F * H13;
  V1_norm = c_T13_re * H22 - v21_re * 0.0F;
  v21_re = c_T13_re * 0.0F + v21_re * H22;
  if (v21_re == 0.0F) {
    c_T13_re = V1_norm / 2.0F;
    v21_re = 0.0F;
  } else if (V1_norm == 0.0F) {
    c_T13_re = 0.0F;
    v21_re /= 2.0F;
  } else {
    c_T13_re = V1_norm / 2.0F;
    v21_re /= 2.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  c_T12_re = b_T12_re * H33 - b_T12_im * 0.0F;
  b_T12_im = b_T12_re * 0.0F + b_T12_im * H33;
  if (b_T12_im == 0.0F) {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im = 0.0F;
  } else if (c_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 2.0F;
  } else {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im /= 2.0F;
  }

  V1_norm = H23 * H23;
  s = H23 * 0.0F + 0.0F * H23;
  d_T11_re = H11 * V1_norm - 0.0F * s;
  c_T11_im = H11 * s + 0.0F * V1_norm;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re /= 2.0F;
    c_T11_im /= 2.0F;
  }

  h_T11.re = (H11 + H22) + H33;
  h_T11.im = 0.0F;
  fc1 = mpower(h_T11);
  if (fc1.im == 0.0F) {
    re = fc1.re / 27.0F;
    im = 0.0F;
  } else if (fc1.re == 0.0F) {
    re = 0.0F;
    im = fc1.im / 27.0F;
  } else {
    re = fc1.re / 27.0F;
    im = fc1.im / 27.0F;
  }

  e_T11_re = (H11 + H22) + H33;
  c_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  c_T22_re = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  b_T22_im = e_T11_re * c_T12_re - 0.0F * c_T22_re;
  c_a_re = e_T11_re * c_T22_re + 0.0F * c_T12_re;
  if (c_a_re == 0.0F) {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re = 0.0F;
  } else if (b_T22_im == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 6.0F;
  } else {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re /= 6.0F;
  }

  c_T12_re = H12 * H13;
  c_T22_re = H12 * 0.0F + 0.0F * H13;
  b_T22_im = H11 * H22;
  d_T11_im = H11 * 0.0F + 0.0F * H22;
  f_T11_re = b_T22_im * H33 - d_T11_im * 0.0F;
  d_T11_im = b_T22_im * 0.0F + d_T11_im * H33;
  if (d_T11_im == 0.0F) {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im = 0.0F;
  } else if (f_T11_re == 0.0F) {
    b_T22_im = 0.0F;
    d_T11_im /= 2.0F;
  } else {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im /= 2.0F;
  }

  fc4.re = ((((((T11.re - c_T13_re) - b_T12_re) - d_T11_re) + re) + e_T11_re) +
            (c_T12_re * H23 - c_T22_re * 0.0F)) + b_T22_im;
  fc4.im = ((((((T11.im - v21_re) - b_T12_im) - c_T11_im) + im) + c_a_re) +
            (c_T12_re * 0.0F + c_T22_re * H23)) + d_T11_im;
  T11 = c_mpower(fc4);
  f_T11_re = (((((a_re + T12_re) + b_T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  b_T22_im = (((((a_im + T12_im) + b_T13_im) + b_T23_im) - T11_im) - b_T11_im) -
    T22_im;
  if (T11.im == 0.0F) {
    if (b_T22_im == 0.0F) {
      a_re = f_T11_re / T11.re;
      a_im = 0.0F;
    } else if (f_T11_re == 0.0F) {
      a_re = 0.0F;
      a_im = b_T22_im / T11.re;
    } else {
      a_re = f_T11_re / T11.re;
      a_im = b_T22_im / T11.re;
    }
  } else if (T11.re == 0.0F) {
    if (f_T11_re == 0.0F) {
      a_re = b_T22_im / T11.im;
      a_im = 0.0F;
    } else if (b_T22_im == 0.0F) {
      a_re = 0.0F;
      a_im = -(f_T11_re / T11.im);
    } else {
      a_re = b_T22_im / T11.im;
      a_im = -(f_T11_re / T11.im);
    }
  } else {
    c_T12_re = (float)fabs(T11.re);
    V1_norm = (float)fabs(T11.im);
    if (c_T12_re > V1_norm) {
      s = T11.im / T11.re;
      V1_norm = T11.re + s * T11.im;
      a_re = (f_T11_re + s * b_T22_im) / V1_norm;
      a_im = (b_T22_im - s * f_T11_re) / V1_norm;
    } else if (V1_norm == c_T12_re) {
      if (T11.re > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (T11.im > 0.0F) {
        V1_norm = 0.5F;
      } else {
        V1_norm = -0.5F;
      }

      a_re = (f_T11_re * s + b_T22_im * V1_norm) / c_T12_re;
      a_im = (b_T22_im * s - f_T11_re * V1_norm) / c_T12_re;
    } else {
      s = T11.re / T11.im;
      V1_norm = T11.im + s * T11.re;
      a_re = (s * f_T11_re + b_T22_im) / V1_norm;
      a_im = (s * b_T22_im - f_T11_re) / V1_norm;
    }
  }

  T12_re = (H12 * H13 + H22 * H23) + H23 * H33;
  T12_im = ((H12 * 0.0F + 0.0F * H13) + (H22 * 0.0F + 0.0F * H23)) + (H23 * 0.0F
    + 0.0F * H33);
  T11_re = (((T11_re + T22_re) + T33_re) + fc0.re) + a_re;
  T11_im = fc0.im + a_im;
  b_T12_re = T12_re * T11_re - T12_im * T11_im;
  T12_im = T12_re * T11_im + T12_im * T11_re;
  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * H13;
  T11_im = H11 * 0.0F + 0.0F * H13;
  c_T13_re = H13 * H22;
  v21_re = H13 * 0.0F + 0.0F * H22;
  d_T11_im = (((H12 * b_T13_re - 0.0F * b_T13_im) - (H12 * T23_re - 0.0F *
    b_T23_im)) - (T11_re * H23 - T11_im * 0.0F)) + (c_T13_re * H23 - v21_re *
    0.0F);
  c_T22_re = (((H12 * b_T13_im + 0.0F * b_T13_re) - (H12 * b_T23_im + 0.0F *
    T23_re)) - (T11_re * 0.0F + T11_im * H23)) + (c_T13_re * 0.0F + v21_re * H23);
  if (c_T22_re == 0.0F) {
    if (T12_im == 0.0F) {
      T12_re = b_T12_re / d_T11_im;
      T12_im = 0.0F;
    } else if (b_T12_re == 0.0F) {
      T12_re = 0.0F;
      T12_im /= d_T11_im;
    } else {
      T12_re = b_T12_re / d_T11_im;
      T12_im /= d_T11_im;
    }
  } else if (d_T11_im == 0.0F) {
    if (b_T12_re == 0.0F) {
      T12_re = T12_im / c_T22_re;
      T12_im = 0.0F;
    } else if (T12_im == 0.0F) {
      T12_re = 0.0F;
      T12_im = -(b_T12_re / c_T22_re);
    } else {
      T12_re = T12_im / c_T22_re;
      T12_im = -(b_T12_re / c_T22_re);
    }
  } else {
    c_T12_re = (float)fabs(d_T11_im);
    if (c_T12_re > c_T22_re) {
      s = c_T22_re / d_T11_im;
      V1_norm = d_T11_im + s * c_T22_re;
      T12_re = (b_T12_re + s * T12_im) / V1_norm;
      T12_im = (T12_im - s * b_T12_re) / V1_norm;
    } else if (c_T22_re == c_T12_re) {
      if (d_T11_im > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (b_T12_re * s + T12_im * -0.5F > 0.0F) {
        T12_re = ((real32_T)rtInf);
      } else if (b_T12_re * s + T12_im * -0.5F < 0.0F) {
        T12_re = ((real32_T)rtMinusInf);
      } else {
        T12_re = ((real32_T)rtNaN);
      }

      if (T12_im * s - b_T12_re * -0.5F > 0.0F) {
        T12_im = ((real32_T)rtInf);
      } else if (T12_im * s - b_T12_re * -0.5F < 0.0F) {
        T12_im = ((real32_T)rtMinusInf);
      } else {
        T12_im = ((real32_T)rtNaN);
      }
    } else {
      s = d_T11_im / c_T22_re;
      V1_norm = c_T22_re + s * d_T11_im;
      T12_re = (s * b_T12_re + T12_im) / V1_norm;
      T12_im = (s * T12_im - b_T12_re) / V1_norm;
    }
  }

  V11_re = (T13_re - b_T23_re) + T12_re;
  V11_im = (T13_im - T23_im) + T12_im;
  i_T11.re = (H11 + H22) + H33;
  i_T11.im = 0.0F;
  fc0 = mpower(i_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
           b_T12_im * 0.0F)) + c_T11_re;
  a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F +
           b_T12_im * H23)) + c_T11_im;
  b_a_re = (H11 + H22) + H33;
  c_a_re = (H11 + H22) + H33;
  j_T11.re = (H11 + H22) + H33;
  j_T11.im = 0.0F;
  fc0 = mpower(j_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  d_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  b_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  e_a_re = (H11 + H22) + H33;
  k_T11.re = (H11 + H22) + H33;
  k_T11.im = 0.0F;
  fc0 = mpower(k_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  f_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  c_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  g_a_re = (H11 + H22) + H33;
  h_a_re = (H11 + H22) + H33;
  l_T11.re = (H11 + H22) + H33;
  l_T11.im = 0.0F;
  fc0 = mpower(l_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  i_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  d_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  j_a_re = (H11 + H22) + H33;
  T11_re = H11 / 3.0F;
  T22_re = H22 / 3.0F;
  T33_re = H33 / 3.0F;
  V1_norm = b_a_re * b_a_re;
  e_a_im = b_a_re * 0.0F + 0.0F * b_a_re;
  if (e_a_im == 0.0F) {
    b_a_re = V1_norm / 9.0F;
    e_a_im = 0.0F;
  } else if (V1_norm == 0.0F) {
    b_a_re = 0.0F;
    e_a_im /= 9.0F;
  } else {
    b_a_re = V1_norm / 9.0F;
    e_a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  if (T13_im == 0.0F) {
    T13_re /= 3.0F;
    T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  i_a.re = (((((b_a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  i_a.im = (((((e_a_im + T12_im) + T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  fc0 = mpower(i_a);
  j_a.re = (a_re * a_re - a_im * a_im) - fc0.re;
  j_a.im = (a_re * a_im + a_im * a_re) - fc0.im;
  fc0 = b_mpower(j_a);
  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  b_T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    b_T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    b_T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  m_T11.re = (H11 + H22) + H33;
  m_T11.im = 0.0F;
  T11 = mpower(m_T11);
  if (T11.im == 0.0F) {
    re = T11.re / 27.0F;
    im = 0.0F;
  } else if (T11.re == 0.0F) {
    re = 0.0F;
    im = T11.im / 27.0F;
  } else {
    re = T11.re / 27.0F;
    im = T11.im / 27.0F;
  }

  c_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  d_T11_re = c_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = c_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  e_T11_re = d_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = d_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (e_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  fc5.re = ((((((fc0.re - T13_re) - T12_re) - b_T11_re) + re) + c_T11_re) +
            (b_T12_re * H23 - b_T12_im * 0.0F)) + d_T11_re;
  fc5.im = ((((((fc0.im - T13_im) - T12_im) - T11_im) + im) + b_T11_im) +
            (b_T12_re * 0.0F + b_T12_im * H23)) + c_T11_im;
  fc0 = c_mpower(fc5);
  a_re = c_a_re * c_a_re;
  a_im = c_a_re * 0.0F + 0.0F * c_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  if (T13_im == 0.0F) {
    T13_re /= 3.0F;
    T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  b_a_re = e_a_re * e_a_re;
  e_a_im = e_a_re * 0.0F + 0.0F * e_a_re;
  if (e_a_im == 0.0F) {
    b_a_re /= 9.0F;
    e_a_im = 0.0F;
  } else if (b_a_re == 0.0F) {
    b_a_re = 0.0F;
    e_a_im /= 9.0F;
  } else {
    b_a_re /= 9.0F;
    e_a_im /= 9.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  if (b_T12_im == 0.0F) {
    b_T12_re /= 3.0F;
    b_T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 3.0F;
  } else {
    b_T12_re /= 3.0F;
    b_T12_im /= 3.0F;
  }

  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  if (b_T13_im == 0.0F) {
    b_T13_re /= 3.0F;
    b_T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 3.0F;
  } else {
    b_T13_re /= 3.0F;
    b_T13_im /= 3.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  if (b_T23_im == 0.0F) {
    b_T23_re /= 3.0F;
    b_T23_im = 0.0F;
  } else if (b_T23_re == 0.0F) {
    b_T23_re = 0.0F;
    b_T23_im /= 3.0F;
  } else {
    b_T23_re /= 3.0F;
    b_T23_im /= 3.0F;
  }

  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 3.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 3.0F;
  } else {
    d_T11_re /= 3.0F;
    c_T11_im /= 3.0F;
  }

  e_T11_re = H11 * H33;
  c_a_re = H11 * 0.0F + 0.0F * H33;
  if (c_a_re == 0.0F) {
    e_T11_re /= 3.0F;
    c_a_re = 0.0F;
  } else if (e_T11_re == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 3.0F;
  } else {
    e_T11_re /= 3.0F;
    c_a_re /= 3.0F;
  }

  c_T22_re = H22 * H33;
  b_T22_im = H22 * 0.0F + 0.0F * H33;
  if (b_T22_im == 0.0F) {
    c_T22_re /= 3.0F;
    b_T22_im = 0.0F;
  } else if (c_T22_re == 0.0F) {
    c_T22_re = 0.0F;
    b_T22_im /= 3.0F;
  } else {
    c_T22_re /= 3.0F;
    b_T22_im /= 3.0F;
  }

  k_a.re = (((((b_a_re + b_T12_re) + b_T13_re) + b_T23_re) - d_T11_re) -
            e_T11_re) - c_T22_re;
  k_a.im = (((((e_a_im + b_T12_im) + b_T13_im) + b_T23_im) - c_T11_im) - c_a_re)
    - b_T22_im;
  T11 = mpower(k_a);
  l_a.re = (d_a_re * d_a_re - b_a_im * b_a_im) - T11.re;
  l_a.im = (d_a_re * b_a_im + b_a_im * d_a_re) - T11.im;
  T11 = b_mpower(l_a);
  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  c_T13_re = b_T13_re * H22 - b_T13_im * 0.0F;
  b_T13_im = b_T13_re * 0.0F + b_T13_im * H22;
  if (b_T13_im == 0.0F) {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im = 0.0F;
  } else if (c_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 2.0F;
  } else {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im /= 2.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  c_T12_re = b_T12_re * H33 - b_T12_im * 0.0F;
  b_T12_im = b_T12_re * 0.0F + b_T12_im * H33;
  if (b_T12_im == 0.0F) {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im = 0.0F;
  } else if (c_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 2.0F;
  } else {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im /= 2.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  d_T11_re = H11 * b_T23_re - 0.0F * b_T23_im;
  c_T11_im = H11 * b_T23_im + 0.0F * b_T23_re;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re /= 2.0F;
    c_T11_im /= 2.0F;
  }

  n_T11.re = (H11 + H22) + H33;
  n_T11.im = 0.0F;
  fc1 = mpower(n_T11);
  if (fc1.im == 0.0F) {
    re = fc1.re / 27.0F;
    im = 0.0F;
  } else if (fc1.re == 0.0F) {
    re = 0.0F;
    im = fc1.im / 27.0F;
  } else {
    re = fc1.re / 27.0F;
    im = fc1.im / 27.0F;
  }

  e_T11_re = (H11 + H22) + H33;
  c_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  c_T22_re = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  b_T22_im = e_T11_re * c_T12_re - 0.0F * c_T22_re;
  c_a_re = e_T11_re * c_T22_re + 0.0F * c_T12_re;
  if (c_a_re == 0.0F) {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re = 0.0F;
  } else if (b_T22_im == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 6.0F;
  } else {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re /= 6.0F;
  }

  c_T12_re = H12 * H13;
  c_T22_re = H12 * 0.0F + 0.0F * H13;
  b_T22_im = H11 * H22;
  d_T11_im = H11 * 0.0F + 0.0F * H22;
  f_T11_re = b_T22_im * H33 - d_T11_im * 0.0F;
  d_T11_im = b_T22_im * 0.0F + d_T11_im * H33;
  if (d_T11_im == 0.0F) {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im = 0.0F;
  } else if (f_T11_re == 0.0F) {
    b_T22_im = 0.0F;
    d_T11_im /= 2.0F;
  } else {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im /= 2.0F;
  }

  fc6.re = ((((((T11.re - b_T13_re) - b_T12_re) - d_T11_re) + re) + e_T11_re) +
            (c_T12_re * H23 - c_T22_re * 0.0F)) + b_T22_im;
  fc6.im = ((((((T11.im - b_T13_im) - b_T12_im) - c_T11_im) + im) + c_a_re) +
            (c_T12_re * 0.0F + c_T22_re * H23)) + d_T11_im;
  T11 = c_mpower(fc6);
  f_T11_re = (((((a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  b_T22_im = (((((a_im + T12_im) + T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  if (T11.im == 0.0F) {
    if (b_T22_im == 0.0F) {
      a_re = f_T11_re / T11.re;
      a_im = 0.0F;
    } else if (f_T11_re == 0.0F) {
      a_re = 0.0F;
      a_im = b_T22_im / T11.re;
    } else {
      a_re = f_T11_re / T11.re;
      a_im = b_T22_im / T11.re;
    }
  } else if (T11.re == 0.0F) {
    if (f_T11_re == 0.0F) {
      a_re = b_T22_im / T11.im;
      a_im = 0.0F;
    } else if (b_T22_im == 0.0F) {
      a_re = 0.0F;
      a_im = -(f_T11_re / T11.im);
    } else {
      a_re = b_T22_im / T11.im;
      a_im = -(f_T11_re / T11.im);
    }
  } else {
    c_T12_re = (float)fabs(T11.re);
    V1_norm = (float)fabs(T11.im);
    if (c_T12_re > V1_norm) {
      s = T11.im / T11.re;
      V1_norm = T11.re + s * T11.im;
      a_re = (f_T11_re + s * b_T22_im) / V1_norm;
      a_im = (b_T22_im - s * f_T11_re) / V1_norm;
    } else if (V1_norm == c_T12_re) {
      if (T11.re > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (T11.im > 0.0F) {
        V1_norm = 0.5F;
      } else {
        V1_norm = -0.5F;
      }

      a_re = (f_T11_re * s + b_T22_im * V1_norm) / c_T12_re;
      a_im = (b_T22_im * s - f_T11_re * V1_norm) / c_T12_re;
    } else {
      s = T11.re / T11.im;
      V1_norm = T11.im + s * T11.re;
      a_re = (s * f_T11_re + b_T22_im) / V1_norm;
      a_im = (s * b_T22_im - f_T11_re) / V1_norm;
    }
  }

  re = 1.73205078F * (fc0.re - a_re);
  im = 1.73205078F * (fc0.im - a_im);
  b_re = re * 0.0F - im;
  im = re + im * 0.0F;
  if (im == 0.0F) {
    re = b_re / 2.0F;
    im = 0.0F;
  } else if (b_re == 0.0F) {
    re = 0.0F;
    im /= 2.0F;
  } else {
    re = b_re / 2.0F;
    im /= 2.0F;
  }

  a_re = g_a_re * g_a_re;
  a_im = g_a_re * 0.0F + 0.0F * g_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  if (T13_im == 0.0F) {
    T13_re /= 3.0F;
    T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  m_a.re = (((((a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  m_a.im = (((((a_im + T12_im) + T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  fc0 = mpower(m_a);
  n_a.re = (f_a_re * f_a_re - c_a_im * c_a_im) - fc0.re;
  n_a.im = (f_a_re * c_a_im + c_a_im * f_a_re) - fc0.im;
  fc0 = b_mpower(n_a);
  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  b_T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    b_T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    b_T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  o_T11.re = (H11 + H22) + H33;
  o_T11.im = 0.0F;
  T11 = mpower(o_T11);
  if (T11.im == 0.0F) {
    b_re = T11.re / 27.0F;
    b_im = 0.0F;
  } else if (T11.re == 0.0F) {
    b_re = 0.0F;
    b_im = T11.im / 27.0F;
  } else {
    b_re = T11.re / 27.0F;
    b_im = T11.im / 27.0F;
  }

  c_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  d_T11_re = c_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = c_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  e_T11_re = d_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = d_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (e_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  fc7.re = ((((((fc0.re - T13_re) - T12_re) - b_T11_re) + b_re) + c_T11_re) +
            (b_T12_re * H23 - b_T12_im * 0.0F)) + d_T11_re;
  fc7.im = ((((((fc0.im - T13_im) - T12_im) - T11_im) + b_im) + b_T11_im) +
            (b_T12_re * 0.0F + b_T12_im * H23)) + c_T11_im;
  fc0 = c_mpower(fc7);
  if (fc0.im == 0.0F) {
    b_re = fc0.re / 2.0F;
    b_im = 0.0F;
  } else if (fc0.re == 0.0F) {
    b_re = 0.0F;
    b_im = fc0.im / 2.0F;
  } else {
    b_re = fc0.re / 2.0F;
    b_im = fc0.im / 2.0F;
  }

  a_re = h_a_re * h_a_re;
  a_im = h_a_re * 0.0F + 0.0F * h_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  if (T13_im == 0.0F) {
    T13_re /= 3.0F;
    T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  b_a_re = j_a_re * j_a_re;
  b_a_im = j_a_re * 0.0F + 0.0F * j_a_re;
  if (b_a_im == 0.0F) {
    b_a_re /= 9.0F;
    b_a_im = 0.0F;
  } else if (b_a_re == 0.0F) {
    b_a_re = 0.0F;
    b_a_im /= 9.0F;
  } else {
    b_a_re /= 9.0F;
    b_a_im /= 9.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  if (b_T12_im == 0.0F) {
    b_T12_re /= 3.0F;
    b_T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 3.0F;
  } else {
    b_T12_re /= 3.0F;
    b_T12_im /= 3.0F;
  }

  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  if (b_T13_im == 0.0F) {
    b_T13_re /= 3.0F;
    b_T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 3.0F;
  } else {
    b_T13_re /= 3.0F;
    b_T13_im /= 3.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  if (b_T23_im == 0.0F) {
    b_T23_re /= 3.0F;
    b_T23_im = 0.0F;
  } else if (b_T23_re == 0.0F) {
    b_T23_re = 0.0F;
    b_T23_im /= 3.0F;
  } else {
    b_T23_re /= 3.0F;
    b_T23_im /= 3.0F;
  }

  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 3.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 3.0F;
  } else {
    d_T11_re /= 3.0F;
    c_T11_im /= 3.0F;
  }

  e_T11_re = H11 * H33;
  c_a_re = H11 * 0.0F + 0.0F * H33;
  if (c_a_re == 0.0F) {
    e_T11_re /= 3.0F;
    c_a_re = 0.0F;
  } else if (e_T11_re == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 3.0F;
  } else {
    e_T11_re /= 3.0F;
    c_a_re /= 3.0F;
  }

  c_T22_re = H22 * H33;
  b_T22_im = H22 * 0.0F + 0.0F * H33;
  if (b_T22_im == 0.0F) {
    c_T22_re /= 3.0F;
    b_T22_im = 0.0F;
  } else if (c_T22_re == 0.0F) {
    c_T22_re = 0.0F;
    b_T22_im /= 3.0F;
  } else {
    c_T22_re /= 3.0F;
    b_T22_im /= 3.0F;
  }

  o_a.re = (((((b_a_re + b_T12_re) + b_T13_re) + b_T23_re) - d_T11_re) -
            e_T11_re) - c_T22_re;
  o_a.im = (((((b_a_im + b_T12_im) + b_T13_im) + b_T23_im) - c_T11_im) - c_a_re)
    - b_T22_im;
  fc0 = mpower(o_a);
  p_a.re = (i_a_re * i_a_re - d_a_im * d_a_im) - fc0.re;
  p_a.im = (i_a_re * d_a_im + d_a_im * i_a_re) - fc0.im;
  fc0 = b_mpower(p_a);
  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  c_T13_re = b_T13_re * H22 - b_T13_im * 0.0F;
  b_T13_im = b_T13_re * 0.0F + b_T13_im * H22;
  if (b_T13_im == 0.0F) {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im = 0.0F;
  } else if (c_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 2.0F;
  } else {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im /= 2.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  c_T12_re = b_T12_re * H33 - b_T12_im * 0.0F;
  b_T12_im = b_T12_re * 0.0F + b_T12_im * H33;
  if (b_T12_im == 0.0F) {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im = 0.0F;
  } else if (c_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 2.0F;
  } else {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im /= 2.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  d_T11_re = H11 * b_T23_re - 0.0F * b_T23_im;
  c_T11_im = H11 * b_T23_im + 0.0F * b_T23_re;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re /= 2.0F;
    c_T11_im /= 2.0F;
  }

  p_T11.re = (H11 + H22) + H33;
  p_T11.im = 0.0F;
  T11 = mpower(p_T11);
  if (T11.im == 0.0F) {
    b_a_re = T11.re / 27.0F;
    d_a_re = 0.0F;
  } else if (T11.re == 0.0F) {
    b_a_re = 0.0F;
    d_a_re = T11.im / 27.0F;
  } else {
    b_a_re = T11.re / 27.0F;
    d_a_re = T11.im / 27.0F;
  }

  e_T11_re = (H11 + H22) + H33;
  c_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  c_T22_re = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  b_T22_im = e_T11_re * c_T12_re - 0.0F * c_T22_re;
  c_a_re = e_T11_re * c_T22_re + 0.0F * c_T12_re;
  if (c_a_re == 0.0F) {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re = 0.0F;
  } else if (b_T22_im == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 6.0F;
  } else {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re /= 6.0F;
  }

  c_T12_re = H12 * H13;
  c_T22_re = H12 * 0.0F + 0.0F * H13;
  b_T22_im = H11 * H22;
  d_T11_im = H11 * 0.0F + 0.0F * H22;
  f_T11_re = b_T22_im * H33 - d_T11_im * 0.0F;
  d_T11_im = b_T22_im * 0.0F + d_T11_im * H33;
  if (d_T11_im == 0.0F) {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im = 0.0F;
  } else if (f_T11_re == 0.0F) {
    b_T22_im = 0.0F;
    d_T11_im /= 2.0F;
  } else {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im /= 2.0F;
  }

  fc8.re = ((((((fc0.re - b_T13_re) - b_T12_re) - d_T11_re) + b_a_re) + e_T11_re)
            + (c_T12_re * H23 - c_T22_re * 0.0F)) + b_T22_im;
  fc8.im = ((((((fc0.im - b_T13_im) - b_T12_im) - c_T11_im) + d_a_re) + c_a_re)
            + (c_T12_re * 0.0F + c_T22_re * H23)) + d_T11_im;
  fc0 = c_mpower(fc8);
  f_T11_re = (((((a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  b_T22_im = (((((a_im + T12_im) + T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  d_T11_im = 2.0F * fc0.re;
  c_T22_re = 2.0F * fc0.im;
  if (c_T22_re == 0.0F) {
    if (b_T22_im == 0.0F) {
      a_re = f_T11_re / d_T11_im;
      a_im = 0.0F;
    } else if (f_T11_re == 0.0F) {
      a_re = 0.0F;
      a_im = b_T22_im / d_T11_im;
    } else {
      a_re = f_T11_re / d_T11_im;
      a_im = b_T22_im / d_T11_im;
    }
  } else if (d_T11_im == 0.0F) {
    if (f_T11_re == 0.0F) {
      a_re = b_T22_im / c_T22_re;
      a_im = 0.0F;
    } else if (b_T22_im == 0.0F) {
      a_re = 0.0F;
      a_im = -(f_T11_re / c_T22_re);
    } else {
      a_re = b_T22_im / c_T22_re;
      a_im = -(f_T11_re / c_T22_re);
    }
  } else {
    c_T12_re = (float)fabs(d_T11_im);
    V1_norm = (float)fabs(c_T22_re);
    if (c_T12_re > V1_norm) {
      s = c_T22_re / d_T11_im;
      V1_norm = d_T11_im + s * c_T22_re;
      a_re = (f_T11_re + s * b_T22_im) / V1_norm;
      a_im = (b_T22_im - s * f_T11_re) / V1_norm;
    } else if (V1_norm == c_T12_re) {
      if (d_T11_im > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (c_T22_re > 0.0F) {
        V1_norm = 0.5F;
      } else {
        V1_norm = -0.5F;
      }

      a_re = (f_T11_re * s + b_T22_im * V1_norm) / c_T12_re;
      a_im = (b_T22_im * s - f_T11_re * V1_norm) / c_T12_re;
    } else {
      s = d_T11_im / c_T22_re;
      V1_norm = c_T22_re + s * d_T11_im;
      a_re = (s * f_T11_re + b_T22_im) / V1_norm;
      a_im = (s * b_T22_im - f_T11_re) / V1_norm;
    }
  }

  a_re = ((((T11_re + T22_re) + T33_re) - re) - b_re) - a_re;
  a_im = ((0.0F - im) - b_im) - a_im;
  q_T11.re = (H11 + H22) + H33;
  q_T11.im = 0.0F;
  fc0 = mpower(q_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  b_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  b_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  c_a_re = (H11 + H22) + H33;
  d_a_re = (H11 + H22) + H33;
  r_T11.re = (H11 + H22) + H33;
  r_T11.im = 0.0F;
  fc0 = mpower(r_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  e_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  c_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  f_a_re = (H11 + H22) + H33;
  s_T11.re = (H11 + H22) + H33;
  s_T11.im = 0.0F;
  fc0 = mpower(s_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  g_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  d_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  h_a_re = (H11 + H22) + H33;
  i_a_re = (H11 + H22) + H33;
  t_T11.re = (H11 + H22) + H33;
  t_T11.im = 0.0F;
  fc0 = mpower(t_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  j_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  e_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  V23_re = (H11 + H22) + H33;
  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  T12_re = H12 * H33;
  T12_im = H12 * 0.0F + 0.0F * H33;
  fc0 = mpower(T23);
  T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * H13;
  T11_im = H11 * 0.0F + 0.0F * H13;
  c_T13_re = H13 * H22;
  v21_re = H13 * 0.0F + 0.0F * H22;
  f_T11_re = (((T13_re * H23 - T13_im * 0.0F) - (T12_re * H13 - T12_im * 0.0F))
              + fc0.re) - (T22_re * H23 - T22_im * 0.0F);
  b_T22_im = (((T13_re * 0.0F + T13_im * H23) - (T12_re * 0.0F + T12_im * H13))
              + fc0.im) - (T22_re * 0.0F + T22_im * H23);
  d_T11_im = (((H12 * b_T13_re - 0.0F * b_T13_im) - (H12 * T23_re - 0.0F *
    T23_im)) - (T11_re * H23 - T11_im * 0.0F)) + (c_T13_re * H23 - v21_re * 0.0F);
  c_T22_re = (((H12 * b_T13_im + 0.0F * b_T13_re) - (H12 * T23_im + 0.0F *
    T23_re)) - (T11_re * 0.0F + T11_im * H23)) + (c_T13_re * 0.0F + v21_re * H23);
  if (c_T22_re == 0.0F) {
    if (b_T22_im == 0.0F) {
      T13_re = f_T11_re / d_T11_im;
      T13_im = 0.0F;
    } else if (f_T11_re == 0.0F) {
      T13_re = 0.0F;
      T13_im = b_T22_im / d_T11_im;
    } else {
      T13_re = f_T11_re / d_T11_im;
      T13_im = b_T22_im / d_T11_im;
    }
  } else if (d_T11_im == 0.0F) {
    if (f_T11_re == 0.0F) {
      T13_re = b_T22_im / c_T22_re;
      T13_im = 0.0F;
    } else if (b_T22_im == 0.0F) {
      T13_re = 0.0F;
      T13_im = -(f_T11_re / c_T22_re);
    } else {
      T13_re = b_T22_im / c_T22_re;
      T13_im = -(f_T11_re / c_T22_re);
    }
  } else {
    c_T12_re = (float)fabs(d_T11_im);
    if (c_T12_re > c_T22_re) {
      s = c_T22_re / d_T11_im;
      V1_norm = d_T11_im + s * c_T22_re;
      T13_re = (f_T11_re + s * b_T22_im) / V1_norm;
      T13_im = (b_T22_im - s * f_T11_re) / V1_norm;
    } else if (c_T22_re == c_T12_re) {
      if (d_T11_im > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (f_T11_re * s + b_T22_im * -0.5F > 0.0F) {
        T13_re = ((real32_T)rtInf);
      } else if (f_T11_re * s + b_T22_im * -0.5F < 0.0F) {
        T13_re = ((real32_T)rtMinusInf);
      } else {
        T13_re = ((real32_T)rtNaN);
      }

      if (b_T22_im * s - f_T11_re * -0.5F > 0.0F) {
        T13_im = ((real32_T)rtInf);
      } else if (b_T22_im * s - f_T11_re * -0.5F < 0.0F) {
        T13_im = ((real32_T)rtMinusInf);
      } else {
        T13_im = ((real32_T)rtNaN);
      }
    } else {
      s = d_T11_im / c_T22_re;
      V1_norm = c_T22_re + s * d_T11_im;
      T13_re = (s * f_T11_re + b_T22_im) / V1_norm;
      T13_im = (s * b_T22_im - f_T11_re) / V1_norm;
    }
  }

  V1_norm = a_re * a_re - a_im * a_im;
  a_im = a_re * a_im + a_im * a_re;
  T23_re = H23 * V1_norm - 0.0F * a_im;
  T23_im = H23 * a_im + 0.0F * V1_norm;
  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * H13;
  T11_im = H11 * 0.0F + 0.0F * H13;
  c_T13_re = H13 * H22;
  v21_re = H13 * 0.0F + 0.0F * H22;
  d_T11_im = (((H12 * b_T13_re - 0.0F * b_T13_im) - (H12 * b_T23_re - 0.0F *
    b_T23_im)) - (T11_re * H23 - T11_im * 0.0F)) + (c_T13_re * H23 - v21_re *
    0.0F);
  c_T22_re = (((H12 * b_T13_im + 0.0F * b_T13_re) - (H12 * b_T23_im + 0.0F *
    b_T23_re)) - (T11_re * 0.0F + T11_im * H23)) + (c_T13_re * 0.0F + v21_re *
    H23);
  if (c_T22_re == 0.0F) {
    if (T23_im == 0.0F) {
      b_T23_re = T23_re / d_T11_im;
      T23_im = 0.0F;
    } else if (T23_re == 0.0F) {
      b_T23_re = 0.0F;
      T23_im /= d_T11_im;
    } else {
      b_T23_re = T23_re / d_T11_im;
      T23_im /= d_T11_im;
    }
  } else if (d_T11_im == 0.0F) {
    if (T23_re == 0.0F) {
      b_T23_re = T23_im / c_T22_re;
      T23_im = 0.0F;
    } else if (T23_im == 0.0F) {
      b_T23_re = 0.0F;
      T23_im = -(T23_re / c_T22_re);
    } else {
      b_T23_re = T23_im / c_T22_re;
      T23_im = -(T23_re / c_T22_re);
    }
  } else {
    c_T12_re = (float)fabs(d_T11_im);
    if (c_T12_re > c_T22_re) {
      s = c_T22_re / d_T11_im;
      V1_norm = d_T11_im + s * c_T22_re;
      b_T23_re = (T23_re + s * T23_im) / V1_norm;
      T23_im = (T23_im - s * T23_re) / V1_norm;
    } else if (c_T22_re == c_T12_re) {
      if (d_T11_im > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (T23_re * s + T23_im * -0.5F > 0.0F) {
        b_T23_re = ((real32_T)rtInf);
      } else if (T23_re * s + T23_im * -0.5F < 0.0F) {
        b_T23_re = ((real32_T)rtMinusInf);
      } else {
        b_T23_re = ((real32_T)rtNaN);
      }

      if (T23_im * s - T23_re * -0.5F > 0.0F) {
        T23_im = ((real32_T)rtInf);
      } else if (T23_im * s - T23_re * -0.5F < 0.0F) {
        T23_im = ((real32_T)rtMinusInf);
      } else {
        T23_im = ((real32_T)rtNaN);
      }
    } else {
      s = d_T11_im / c_T22_re;
      V1_norm = c_T22_re + s * d_T11_im;
      b_T23_re = (s * T23_re + T23_im) / V1_norm;
      T23_im = (s * T23_im - T23_re) / V1_norm;
    }
  }

  T11_re = H11 / 3.0F;
  T22_re = H22 / 3.0F;
  T33_re = H33 / 3.0F;
  a_re = c_a_re * c_a_re;
  a_im = c_a_re * 0.0F + 0.0F * c_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  if (b_T13_im == 0.0F) {
    b_T13_re /= 3.0F;
    b_T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 3.0F;
  } else {
    b_T13_re /= 3.0F;
    b_T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  if (b_T23_im == 0.0F) {
    T23_re /= 3.0F;
    b_T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    b_T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    b_T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  q_a.re = (((((a_re + T12_re) + b_T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  q_a.im = (((((a_im + T12_im) + b_T13_im) + b_T23_im) - T11_im) - b_T11_im) -
    T22_im;
  fc0 = mpower(q_a);
  r_a.re = (b_a_re * b_a_re - b_a_im * b_a_im) - fc0.re;
  r_a.im = (b_a_re * b_a_im + b_a_im * b_a_re) - fc0.im;
  fc0 = b_mpower(r_a);
  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  c_T13_re = b_T13_re * H22 - b_T13_im * 0.0F;
  b_T13_im = b_T13_re * 0.0F + b_T13_im * H22;
  if (b_T13_im == 0.0F) {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im = 0.0F;
  } else if (c_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 2.0F;
  } else {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  b_T11_re = H11 * T23_re - 0.0F * b_T23_im;
  T11_im = H11 * b_T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    b_T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    b_T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  u_T11.re = (H11 + H22) + H33;
  u_T11.im = 0.0F;
  T11 = mpower(u_T11);
  if (T11.im == 0.0F) {
    re = T11.re / 27.0F;
    im = 0.0F;
  } else if (T11.re == 0.0F) {
    re = 0.0F;
    im = T11.im / 27.0F;
  } else {
    re = T11.re / 27.0F;
    im = T11.im / 27.0F;
  }

  c_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  d_T11_re = c_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = c_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  e_T11_re = d_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = d_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (e_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  fc9.re = ((((((fc0.re - b_T13_re) - T12_re) - b_T11_re) + re) + c_T11_re) +
            (b_T12_re * H23 - b_T12_im * 0.0F)) + d_T11_re;
  fc9.im = ((((((fc0.im - b_T13_im) - T12_im) - T11_im) + im) + b_T11_im) +
            (b_T12_re * 0.0F + b_T12_im * H23)) + c_T11_im;
  fc0 = c_mpower(fc9);
  a_re = d_a_re * d_a_re;
  a_im = d_a_re * 0.0F + 0.0F * d_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  if (b_T13_im == 0.0F) {
    b_T13_re /= 3.0F;
    b_T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 3.0F;
  } else {
    b_T13_re /= 3.0F;
    b_T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  if (b_T23_im == 0.0F) {
    T23_re /= 3.0F;
    b_T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    b_T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    b_T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  b_a_re = f_a_re * f_a_re;
  b_a_im = f_a_re * 0.0F + 0.0F * f_a_re;
  if (b_a_im == 0.0F) {
    b_a_re /= 9.0F;
    b_a_im = 0.0F;
  } else if (b_a_re == 0.0F) {
    b_a_re = 0.0F;
    b_a_im /= 9.0F;
  } else {
    b_a_re /= 9.0F;
    b_a_im /= 9.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  if (b_T12_im == 0.0F) {
    b_T12_re /= 3.0F;
    b_T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 3.0F;
  } else {
    b_T12_re /= 3.0F;
    b_T12_im /= 3.0F;
  }

  c_T13_re = H13 * H13;
  v21_re = H13 * 0.0F + 0.0F * H13;
  if (v21_re == 0.0F) {
    c_T13_re /= 3.0F;
    v21_re = 0.0F;
  } else if (c_T13_re == 0.0F) {
    c_T13_re = 0.0F;
    v21_re /= 3.0F;
  } else {
    c_T13_re /= 3.0F;
    v21_re /= 3.0F;
  }

  V1_norm = H23 * H23;
  s = H23 * 0.0F + 0.0F * H23;
  if (s == 0.0F) {
    V1_norm /= 3.0F;
    s = 0.0F;
  } else if (V1_norm == 0.0F) {
    V1_norm = 0.0F;
    s /= 3.0F;
  } else {
    V1_norm /= 3.0F;
    s /= 3.0F;
  }

  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 3.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 3.0F;
  } else {
    d_T11_re /= 3.0F;
    c_T11_im /= 3.0F;
  }

  e_T11_re = H11 * H33;
  c_a_re = H11 * 0.0F + 0.0F * H33;
  if (c_a_re == 0.0F) {
    e_T11_re /= 3.0F;
    c_a_re = 0.0F;
  } else if (e_T11_re == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 3.0F;
  } else {
    e_T11_re /= 3.0F;
    c_a_re /= 3.0F;
  }

  c_T22_re = H22 * H33;
  b_T22_im = H22 * 0.0F + 0.0F * H33;
  if (b_T22_im == 0.0F) {
    c_T22_re /= 3.0F;
    b_T22_im = 0.0F;
  } else if (c_T22_re == 0.0F) {
    c_T22_re = 0.0F;
    b_T22_im /= 3.0F;
  } else {
    c_T22_re /= 3.0F;
    b_T22_im /= 3.0F;
  }

  s_a.re = (((((b_a_re + b_T12_re) + c_T13_re) + V1_norm) - d_T11_re) - e_T11_re)
    - c_T22_re;
  s_a.im = (((((b_a_im + b_T12_im) + v21_re) + s) - c_T11_im) - c_a_re) -
    b_T22_im;
  T11 = mpower(s_a);
  t_a.re = (e_a_re * e_a_re - c_a_im * c_a_im) - T11.re;
  t_a.im = (e_a_re * c_a_im + c_a_im * e_a_re) - T11.im;
  T11 = b_mpower(t_a);
  c_T13_re = H13 * H13;
  v21_re = H13 * 0.0F + 0.0F * H13;
  V1_norm = c_T13_re * H22 - v21_re * 0.0F;
  v21_re = c_T13_re * 0.0F + v21_re * H22;
  if (v21_re == 0.0F) {
    c_T13_re = V1_norm / 2.0F;
    v21_re = 0.0F;
  } else if (V1_norm == 0.0F) {
    c_T13_re = 0.0F;
    v21_re /= 2.0F;
  } else {
    c_T13_re = V1_norm / 2.0F;
    v21_re /= 2.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  c_T12_re = b_T12_re * H33 - b_T12_im * 0.0F;
  b_T12_im = b_T12_re * 0.0F + b_T12_im * H33;
  if (b_T12_im == 0.0F) {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im = 0.0F;
  } else if (c_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 2.0F;
  } else {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im /= 2.0F;
  }

  V1_norm = H23 * H23;
  s = H23 * 0.0F + 0.0F * H23;
  d_T11_re = H11 * V1_norm - 0.0F * s;
  c_T11_im = H11 * s + 0.0F * V1_norm;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re /= 2.0F;
    c_T11_im /= 2.0F;
  }

  v_T11.re = (H11 + H22) + H33;
  v_T11.im = 0.0F;
  fc1 = mpower(v_T11);
  if (fc1.im == 0.0F) {
    re = fc1.re / 27.0F;
    im = 0.0F;
  } else if (fc1.re == 0.0F) {
    re = 0.0F;
    im = fc1.im / 27.0F;
  } else {
    re = fc1.re / 27.0F;
    im = fc1.im / 27.0F;
  }

  e_T11_re = (H11 + H22) + H33;
  c_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  c_T22_re = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  b_T22_im = e_T11_re * c_T12_re - 0.0F * c_T22_re;
  c_a_re = e_T11_re * c_T22_re + 0.0F * c_T12_re;
  if (c_a_re == 0.0F) {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re = 0.0F;
  } else if (b_T22_im == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 6.0F;
  } else {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re /= 6.0F;
  }

  c_T12_re = H12 * H13;
  c_T22_re = H12 * 0.0F + 0.0F * H13;
  b_T22_im = H11 * H22;
  d_T11_im = H11 * 0.0F + 0.0F * H22;
  f_T11_re = b_T22_im * H33 - d_T11_im * 0.0F;
  d_T11_im = b_T22_im * 0.0F + d_T11_im * H33;
  if (d_T11_im == 0.0F) {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im = 0.0F;
  } else if (f_T11_re == 0.0F) {
    b_T22_im = 0.0F;
    d_T11_im /= 2.0F;
  } else {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im /= 2.0F;
  }

  fc10.re = ((((((T11.re - c_T13_re) - b_T12_re) - d_T11_re) + re) + e_T11_re) +
             (c_T12_re * H23 - c_T22_re * 0.0F)) + b_T22_im;
  fc10.im = ((((((T11.im - v21_re) - b_T12_im) - c_T11_im) + im) + c_a_re) +
             (c_T12_re * 0.0F + c_T22_re * H23)) + d_T11_im;
  T11 = c_mpower(fc10);
  f_T11_re = (((((a_re + T12_re) + b_T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  b_T22_im = (((((a_im + T12_im) + b_T13_im) + b_T23_im) - T11_im) - b_T11_im) -
    T22_im;
  if (T11.im == 0.0F) {
    if (b_T22_im == 0.0F) {
      a_re = f_T11_re / T11.re;
      a_im = 0.0F;
    } else if (f_T11_re == 0.0F) {
      a_re = 0.0F;
      a_im = b_T22_im / T11.re;
    } else {
      a_re = f_T11_re / T11.re;
      a_im = b_T22_im / T11.re;
    }
  } else if (T11.re == 0.0F) {
    if (f_T11_re == 0.0F) {
      a_re = b_T22_im / T11.im;
      a_im = 0.0F;
    } else if (b_T22_im == 0.0F) {
      a_re = 0.0F;
      a_im = -(f_T11_re / T11.im);
    } else {
      a_re = b_T22_im / T11.im;
      a_im = -(f_T11_re / T11.im);
    }
  } else {
    c_T12_re = (float)fabs(T11.re);
    V1_norm = (float)fabs(T11.im);
    if (c_T12_re > V1_norm) {
      s = T11.im / T11.re;
      V1_norm = T11.re + s * T11.im;
      a_re = (f_T11_re + s * b_T22_im) / V1_norm;
      a_im = (b_T22_im - s * f_T11_re) / V1_norm;
    } else if (V1_norm == c_T12_re) {
      if (T11.re > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (T11.im > 0.0F) {
        V1_norm = 0.5F;
      } else {
        V1_norm = -0.5F;
      }

      a_re = (f_T11_re * s + b_T22_im * V1_norm) / c_T12_re;
      a_im = (b_T22_im * s - f_T11_re * V1_norm) / c_T12_re;
    } else {
      s = T11.re / T11.im;
      V1_norm = T11.im + s * T11.re;
      a_re = (s * f_T11_re + b_T22_im) / V1_norm;
      a_im = (s * b_T22_im - f_T11_re) / V1_norm;
    }
  }

  re = 1.73205078F * (fc0.re - a_re);
  im = 1.73205078F * (fc0.im - a_im);
  b_re = re * 0.0F - im;
  im = re + im * 0.0F;
  if (im == 0.0F) {
    re = b_re / 2.0F;
    im = 0.0F;
  } else if (b_re == 0.0F) {
    re = 0.0F;
    im /= 2.0F;
  } else {
    re = b_re / 2.0F;
    im /= 2.0F;
  }

  a_re = h_a_re * h_a_re;
  a_im = h_a_re * 0.0F + 0.0F * h_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  if (b_T13_im == 0.0F) {
    b_T13_re /= 3.0F;
    b_T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 3.0F;
  } else {
    b_T13_re /= 3.0F;
    b_T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  if (b_T23_im == 0.0F) {
    T23_re /= 3.0F;
    b_T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    b_T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    b_T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  u_a.re = (((((a_re + T12_re) + b_T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  u_a.im = (((((a_im + T12_im) + b_T13_im) + b_T23_im) - T11_im) - b_T11_im) -
    T22_im;
  fc0 = mpower(u_a);
  v_a.re = (g_a_re * g_a_re - d_a_im * d_a_im) - fc0.re;
  v_a.im = (g_a_re * d_a_im + d_a_im * g_a_re) - fc0.im;
  fc0 = b_mpower(v_a);
  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  c_T13_re = b_T13_re * H22 - b_T13_im * 0.0F;
  b_T13_im = b_T13_re * 0.0F + b_T13_im * H22;
  if (b_T13_im == 0.0F) {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im = 0.0F;
  } else if (c_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 2.0F;
  } else {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  b_T11_re = H11 * T23_re - 0.0F * b_T23_im;
  T11_im = H11 * b_T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    b_T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    b_T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  w_T11.re = (H11 + H22) + H33;
  w_T11.im = 0.0F;
  T11 = mpower(w_T11);
  if (T11.im == 0.0F) {
    b_re = T11.re / 27.0F;
    b_im = 0.0F;
  } else if (T11.re == 0.0F) {
    b_re = 0.0F;
    b_im = T11.im / 27.0F;
  } else {
    b_re = T11.re / 27.0F;
    b_im = T11.im / 27.0F;
  }

  c_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  d_T11_re = c_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = c_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  e_T11_re = d_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = d_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (e_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  fc11.re = ((((((fc0.re - b_T13_re) - T12_re) - b_T11_re) + b_re) + c_T11_re) +
             (b_T12_re * H23 - b_T12_im * 0.0F)) + d_T11_re;
  fc11.im = ((((((fc0.im - b_T13_im) - T12_im) - T11_im) + b_im) + b_T11_im) +
             (b_T12_re * 0.0F + b_T12_im * H23)) + c_T11_im;
  fc0 = c_mpower(fc11);
  if (fc0.im == 0.0F) {
    b_re = fc0.re / 2.0F;
    b_im = 0.0F;
  } else if (fc0.re == 0.0F) {
    b_re = 0.0F;
    b_im = fc0.im / 2.0F;
  } else {
    b_re = fc0.re / 2.0F;
    b_im = fc0.im / 2.0F;
  }

  a_re = i_a_re * i_a_re;
  a_im = i_a_re * 0.0F + 0.0F * i_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  if (b_T13_im == 0.0F) {
    b_T13_re /= 3.0F;
    b_T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 3.0F;
  } else {
    b_T13_re /= 3.0F;
    b_T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  if (b_T23_im == 0.0F) {
    T23_re /= 3.0F;
    b_T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    b_T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    b_T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  V1_norm = V23_re * V23_re;
  c_T12_re = V23_re * 0.0F + 0.0F * V23_re;
  if (c_T12_re == 0.0F) {
    V23_re = V1_norm / 9.0F;
    c_T12_re = 0.0F;
  } else if (V1_norm == 0.0F) {
    V23_re = 0.0F;
    c_T12_re /= 9.0F;
  } else {
    V23_re = V1_norm / 9.0F;
    c_T12_re /= 9.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  if (b_T12_im == 0.0F) {
    b_T12_re /= 3.0F;
    b_T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 3.0F;
  } else {
    b_T12_re /= 3.0F;
    b_T12_im /= 3.0F;
  }

  c_T13_re = H13 * H13;
  v21_re = H13 * 0.0F + 0.0F * H13;
  if (v21_re == 0.0F) {
    c_T13_re /= 3.0F;
    v21_re = 0.0F;
  } else if (c_T13_re == 0.0F) {
    c_T13_re = 0.0F;
    v21_re /= 3.0F;
  } else {
    c_T13_re /= 3.0F;
    v21_re /= 3.0F;
  }

  V1_norm = H23 * H23;
  s = H23 * 0.0F + 0.0F * H23;
  if (s == 0.0F) {
    V1_norm /= 3.0F;
    s = 0.0F;
  } else if (V1_norm == 0.0F) {
    V1_norm = 0.0F;
    s /= 3.0F;
  } else {
    V1_norm /= 3.0F;
    s /= 3.0F;
  }

  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 3.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 3.0F;
  } else {
    d_T11_re /= 3.0F;
    c_T11_im /= 3.0F;
  }

  e_T11_re = H11 * H33;
  c_a_re = H11 * 0.0F + 0.0F * H33;
  if (c_a_re == 0.0F) {
    e_T11_re /= 3.0F;
    c_a_re = 0.0F;
  } else if (e_T11_re == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 3.0F;
  } else {
    e_T11_re /= 3.0F;
    c_a_re /= 3.0F;
  }

  c_T22_re = H22 * H33;
  b_T22_im = H22 * 0.0F + 0.0F * H33;
  if (b_T22_im == 0.0F) {
    c_T22_re /= 3.0F;
    b_T22_im = 0.0F;
  } else if (c_T22_re == 0.0F) {
    c_T22_re = 0.0F;
    b_T22_im /= 3.0F;
  } else {
    c_T22_re /= 3.0F;
    b_T22_im /= 3.0F;
  }

  V23.re = (((((V23_re + b_T12_re) + c_T13_re) + V1_norm) - d_T11_re) - e_T11_re)
    - c_T22_re;
  V23.im = (((((c_T12_re + b_T12_im) + v21_re) + s) - c_T11_im) - c_a_re) -
    b_T22_im;
  fc0 = mpower(V23);
  w_a.re = (j_a_re * j_a_re - e_a_im * e_a_im) - fc0.re;
  w_a.im = (j_a_re * e_a_im + e_a_im * j_a_re) - fc0.im;
  fc0 = b_mpower(w_a);
  c_T13_re = H13 * H13;
  v21_re = H13 * 0.0F + 0.0F * H13;
  V1_norm = c_T13_re * H22 - v21_re * 0.0F;
  v21_re = c_T13_re * 0.0F + v21_re * H22;
  if (v21_re == 0.0F) {
    c_T13_re = V1_norm / 2.0F;
    v21_re = 0.0F;
  } else if (V1_norm == 0.0F) {
    c_T13_re = 0.0F;
    v21_re /= 2.0F;
  } else {
    c_T13_re = V1_norm / 2.0F;
    v21_re /= 2.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  c_T12_re = b_T12_re * H33 - b_T12_im * 0.0F;
  b_T12_im = b_T12_re * 0.0F + b_T12_im * H33;
  if (b_T12_im == 0.0F) {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im = 0.0F;
  } else if (c_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 2.0F;
  } else {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im /= 2.0F;
  }

  V1_norm = H23 * H23;
  s = H23 * 0.0F + 0.0F * H23;
  d_T11_re = H11 * V1_norm - 0.0F * s;
  c_T11_im = H11 * s + 0.0F * V1_norm;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re /= 2.0F;
    c_T11_im /= 2.0F;
  }

  x_T11.re = (H11 + H22) + H33;
  x_T11.im = 0.0F;
  T11 = mpower(x_T11);
  if (T11.im == 0.0F) {
    b_a_re = T11.re / 27.0F;
    d_a_re = 0.0F;
  } else if (T11.re == 0.0F) {
    b_a_re = 0.0F;
    d_a_re = T11.im / 27.0F;
  } else {
    b_a_re = T11.re / 27.0F;
    d_a_re = T11.im / 27.0F;
  }

  e_T11_re = (H11 + H22) + H33;
  c_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  c_T22_re = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  b_T22_im = e_T11_re * c_T12_re - 0.0F * c_T22_re;
  c_a_re = e_T11_re * c_T22_re + 0.0F * c_T12_re;
  if (c_a_re == 0.0F) {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re = 0.0F;
  } else if (b_T22_im == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 6.0F;
  } else {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re /= 6.0F;
  }

  c_T12_re = H12 * H13;
  c_T22_re = H12 * 0.0F + 0.0F * H13;
  b_T22_im = H11 * H22;
  d_T11_im = H11 * 0.0F + 0.0F * H22;
  f_T11_re = b_T22_im * H33 - d_T11_im * 0.0F;
  d_T11_im = b_T22_im * 0.0F + d_T11_im * H33;
  if (d_T11_im == 0.0F) {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im = 0.0F;
  } else if (f_T11_re == 0.0F) {
    b_T22_im = 0.0F;
    d_T11_im /= 2.0F;
  } else {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im /= 2.0F;
  }

  fc12.re = ((((((fc0.re - c_T13_re) - b_T12_re) - d_T11_re) + b_a_re) +
              e_T11_re) + (c_T12_re * H23 - c_T22_re * 0.0F)) + b_T22_im;
  fc12.im = ((((((fc0.im - v21_re) - b_T12_im) - c_T11_im) + d_a_re) + c_a_re) +
             (c_T12_re * 0.0F + c_T22_re * H23)) + d_T11_im;
  fc0 = c_mpower(fc12);
  f_T11_re = (((((a_re + T12_re) + b_T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  b_T22_im = (((((a_im + T12_im) + b_T13_im) + b_T23_im) - T11_im) - b_T11_im) -
    T22_im;
  d_T11_im = 2.0F * fc0.re;
  c_T22_re = 2.0F * fc0.im;
  if (c_T22_re == 0.0F) {
    if (b_T22_im == 0.0F) {
      a_re = f_T11_re / d_T11_im;
      a_im = 0.0F;
    } else if (f_T11_re == 0.0F) {
      a_re = 0.0F;
      a_im = b_T22_im / d_T11_im;
    } else {
      a_re = f_T11_re / d_T11_im;
      a_im = b_T22_im / d_T11_im;
    }
  } else if (d_T11_im == 0.0F) {
    if (f_T11_re == 0.0F) {
      a_re = b_T22_im / c_T22_re;
      a_im = 0.0F;
    } else if (b_T22_im == 0.0F) {
      a_re = 0.0F;
      a_im = -(f_T11_re / c_T22_re);
    } else {
      a_re = b_T22_im / c_T22_re;
      a_im = -(f_T11_re / c_T22_re);
    }
  } else {
    c_T12_re = (float)fabs(d_T11_im);
    V1_norm = (float)fabs(c_T22_re);
    if (c_T12_re > V1_norm) {
      s = c_T22_re / d_T11_im;
      V1_norm = d_T11_im + s * c_T22_re;
      a_re = (f_T11_re + s * b_T22_im) / V1_norm;
      a_im = (b_T22_im - s * f_T11_re) / V1_norm;
    } else if (V1_norm == c_T12_re) {
      if (d_T11_im > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (c_T22_re > 0.0F) {
        V1_norm = 0.5F;
      } else {
        V1_norm = -0.5F;
      }

      a_re = (f_T11_re * s + b_T22_im * V1_norm) / c_T12_re;
      a_im = (b_T22_im * s - f_T11_re * V1_norm) / c_T12_re;
    } else {
      s = d_T11_im / c_T22_re;
      V1_norm = c_T22_re + s * d_T11_im;
      a_re = (s * f_T11_re + b_T22_im) / V1_norm;
      a_im = (s * b_T22_im - f_T11_re) / V1_norm;
    }
  }

  T12_re = (H12 * H13 + H22 * H23) + H23 * H33;
  T12_im = ((H12 * 0.0F + 0.0F * H13) + (H22 * 0.0F + 0.0F * H23)) + (H23 * 0.0F
    + 0.0F * H33);
  T11_re = ((((T11_re + T22_re) + T33_re) - re) - b_re) - a_re;
  T11_im = ((0.0F - im) - b_im) - a_im;
  b_T12_re = T12_re * T11_re - T12_im * T11_im;
  T12_im = T12_re * T11_im + T12_im * T11_re;
  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * H13;
  T11_im = H11 * 0.0F + 0.0F * H13;
  c_T13_re = H13 * H22;
  v21_re = H13 * 0.0F + 0.0F * H22;
  d_T11_im = (((H12 * b_T13_re - 0.0F * b_T13_im) - (H12 * T23_re - 0.0F *
    b_T23_im)) - (T11_re * H23 - T11_im * 0.0F)) + (c_T13_re * H23 - v21_re *
    0.0F);
  c_T22_re = (((H12 * b_T13_im + 0.0F * b_T13_re) - (H12 * b_T23_im + 0.0F *
    T23_re)) - (T11_re * 0.0F + T11_im * H23)) + (c_T13_re * 0.0F + v21_re * H23);
  if (c_T22_re == 0.0F) {
    if (T12_im == 0.0F) {
      T12_re = b_T12_re / d_T11_im;
      T12_im = 0.0F;
    } else if (b_T12_re == 0.0F) {
      T12_re = 0.0F;
      T12_im /= d_T11_im;
    } else {
      T12_re = b_T12_re / d_T11_im;
      T12_im /= d_T11_im;
    }
  } else if (d_T11_im == 0.0F) {
    if (b_T12_re == 0.0F) {
      T12_re = T12_im / c_T22_re;
      T12_im = 0.0F;
    } else if (T12_im == 0.0F) {
      T12_re = 0.0F;
      T12_im = -(b_T12_re / c_T22_re);
    } else {
      T12_re = T12_im / c_T22_re;
      T12_im = -(b_T12_re / c_T22_re);
    }
  } else {
    c_T12_re = (float)fabs(d_T11_im);
    if (c_T12_re > c_T22_re) {
      s = c_T22_re / d_T11_im;
      V1_norm = d_T11_im + s * c_T22_re;
      T12_re = (b_T12_re + s * T12_im) / V1_norm;
      T12_im = (T12_im - s * b_T12_re) / V1_norm;
    } else if (c_T22_re == c_T12_re) {
      if (d_T11_im > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (b_T12_re * s + T12_im * -0.5F > 0.0F) {
        T12_re = ((real32_T)rtInf);
      } else if (b_T12_re * s + T12_im * -0.5F < 0.0F) {
        T12_re = ((real32_T)rtMinusInf);
      } else {
        T12_re = ((real32_T)rtNaN);
      }

      if (T12_im * s - b_T12_re * -0.5F > 0.0F) {
        T12_im = ((real32_T)rtInf);
      } else if (T12_im * s - b_T12_re * -0.5F < 0.0F) {
        T12_im = ((real32_T)rtMinusInf);
      } else {
        T12_im = ((real32_T)rtNaN);
      }
    } else {
      s = d_T11_im / c_T22_re;
      V1_norm = c_T22_re + s * d_T11_im;
      T12_re = (s * b_T12_re + T12_im) / V1_norm;
      T12_im = (s * T12_im - b_T12_re) / V1_norm;
    }
  }

  V12_re = (T13_re - b_T23_re) + T12_re;
  V12_im = (T13_im - T23_im) + T12_im;
  y_T11.re = (H11 + H22) + H33;
  y_T11.im = 0.0F;
  fc0 = mpower(y_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
           b_T12_im * 0.0F)) + c_T11_re;
  a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F +
           b_T12_im * H23)) + c_T11_im;
  b_a_re = (H11 + H22) + H33;
  c_a_re = (H11 + H22) + H33;
  ab_T11.re = (H11 + H22) + H33;
  ab_T11.im = 0.0F;
  fc0 = mpower(ab_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  d_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  b_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  e_a_re = (H11 + H22) + H33;
  bb_T11.re = (H11 + H22) + H33;
  bb_T11.im = 0.0F;
  fc0 = mpower(bb_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  f_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  c_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  g_a_re = (H11 + H22) + H33;
  h_a_re = (H11 + H22) + H33;
  cb_T11.re = (H11 + H22) + H33;
  cb_T11.im = 0.0F;
  fc0 = mpower(cb_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  i_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  d_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  j_a_re = (H11 + H22) + H33;
  T11_re = H11 / 3.0F;
  T22_re = H22 / 3.0F;
  T33_re = H33 / 3.0F;
  V1_norm = b_a_re * b_a_re;
  e_a_im = b_a_re * 0.0F + 0.0F * b_a_re;
  if (e_a_im == 0.0F) {
    b_a_re = V1_norm / 9.0F;
    e_a_im = 0.0F;
  } else if (V1_norm == 0.0F) {
    b_a_re = 0.0F;
    e_a_im /= 9.0F;
  } else {
    b_a_re = V1_norm / 9.0F;
    e_a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  if (T13_im == 0.0F) {
    T13_re /= 3.0F;
    T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  x_a.re = (((((b_a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  x_a.im = (((((e_a_im + T12_im) + T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  fc0 = mpower(x_a);
  y_a.re = (a_re * a_re - a_im * a_im) - fc0.re;
  y_a.im = (a_re * a_im + a_im * a_re) - fc0.im;
  fc0 = b_mpower(y_a);
  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  b_T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    b_T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    b_T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  db_T11.re = (H11 + H22) + H33;
  db_T11.im = 0.0F;
  T11 = mpower(db_T11);
  if (T11.im == 0.0F) {
    re = T11.re / 27.0F;
    im = 0.0F;
  } else if (T11.re == 0.0F) {
    re = 0.0F;
    im = T11.im / 27.0F;
  } else {
    re = T11.re / 27.0F;
    im = T11.im / 27.0F;
  }

  c_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  d_T11_re = c_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = c_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  e_T11_re = d_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = d_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (e_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  fc13.re = ((((((fc0.re - T13_re) - T12_re) - b_T11_re) + re) + c_T11_re) +
             (b_T12_re * H23 - b_T12_im * 0.0F)) + d_T11_re;
  fc13.im = ((((((fc0.im - T13_im) - T12_im) - T11_im) + im) + b_T11_im) +
             (b_T12_re * 0.0F + b_T12_im * H23)) + c_T11_im;
  fc0 = c_mpower(fc13);
  a_re = c_a_re * c_a_re;
  a_im = c_a_re * 0.0F + 0.0F * c_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  if (T13_im == 0.0F) {
    T13_re /= 3.0F;
    T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  b_a_re = e_a_re * e_a_re;
  e_a_im = e_a_re * 0.0F + 0.0F * e_a_re;
  if (e_a_im == 0.0F) {
    b_a_re /= 9.0F;
    e_a_im = 0.0F;
  } else if (b_a_re == 0.0F) {
    b_a_re = 0.0F;
    e_a_im /= 9.0F;
  } else {
    b_a_re /= 9.0F;
    e_a_im /= 9.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  if (b_T12_im == 0.0F) {
    b_T12_re /= 3.0F;
    b_T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 3.0F;
  } else {
    b_T12_re /= 3.0F;
    b_T12_im /= 3.0F;
  }

  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  if (b_T13_im == 0.0F) {
    b_T13_re /= 3.0F;
    b_T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 3.0F;
  } else {
    b_T13_re /= 3.0F;
    b_T13_im /= 3.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  if (b_T23_im == 0.0F) {
    b_T23_re /= 3.0F;
    b_T23_im = 0.0F;
  } else if (b_T23_re == 0.0F) {
    b_T23_re = 0.0F;
    b_T23_im /= 3.0F;
  } else {
    b_T23_re /= 3.0F;
    b_T23_im /= 3.0F;
  }

  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 3.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 3.0F;
  } else {
    d_T11_re /= 3.0F;
    c_T11_im /= 3.0F;
  }

  e_T11_re = H11 * H33;
  c_a_re = H11 * 0.0F + 0.0F * H33;
  if (c_a_re == 0.0F) {
    e_T11_re /= 3.0F;
    c_a_re = 0.0F;
  } else if (e_T11_re == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 3.0F;
  } else {
    e_T11_re /= 3.0F;
    c_a_re /= 3.0F;
  }

  c_T22_re = H22 * H33;
  b_T22_im = H22 * 0.0F + 0.0F * H33;
  if (b_T22_im == 0.0F) {
    c_T22_re /= 3.0F;
    b_T22_im = 0.0F;
  } else if (c_T22_re == 0.0F) {
    c_T22_re = 0.0F;
    b_T22_im /= 3.0F;
  } else {
    c_T22_re /= 3.0F;
    b_T22_im /= 3.0F;
  }

  ab_a.re = (((((b_a_re + b_T12_re) + b_T13_re) + b_T23_re) - d_T11_re) -
             e_T11_re) - c_T22_re;
  ab_a.im = (((((e_a_im + b_T12_im) + b_T13_im) + b_T23_im) - c_T11_im) - c_a_re)
    - b_T22_im;
  T11 = mpower(ab_a);
  bb_a.re = (d_a_re * d_a_re - b_a_im * b_a_im) - T11.re;
  bb_a.im = (d_a_re * b_a_im + b_a_im * d_a_re) - T11.im;
  T11 = b_mpower(bb_a);
  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  c_T13_re = b_T13_re * H22 - b_T13_im * 0.0F;
  b_T13_im = b_T13_re * 0.0F + b_T13_im * H22;
  if (b_T13_im == 0.0F) {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im = 0.0F;
  } else if (c_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 2.0F;
  } else {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im /= 2.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  c_T12_re = b_T12_re * H33 - b_T12_im * 0.0F;
  b_T12_im = b_T12_re * 0.0F + b_T12_im * H33;
  if (b_T12_im == 0.0F) {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im = 0.0F;
  } else if (c_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 2.0F;
  } else {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im /= 2.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  d_T11_re = H11 * b_T23_re - 0.0F * b_T23_im;
  c_T11_im = H11 * b_T23_im + 0.0F * b_T23_re;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re /= 2.0F;
    c_T11_im /= 2.0F;
  }

  eb_T11.re = (H11 + H22) + H33;
  eb_T11.im = 0.0F;
  fc1 = mpower(eb_T11);
  if (fc1.im == 0.0F) {
    re = fc1.re / 27.0F;
    im = 0.0F;
  } else if (fc1.re == 0.0F) {
    re = 0.0F;
    im = fc1.im / 27.0F;
  } else {
    re = fc1.re / 27.0F;
    im = fc1.im / 27.0F;
  }

  e_T11_re = (H11 + H22) + H33;
  c_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  c_T22_re = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  b_T22_im = e_T11_re * c_T12_re - 0.0F * c_T22_re;
  c_a_re = e_T11_re * c_T22_re + 0.0F * c_T12_re;
  if (c_a_re == 0.0F) {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re = 0.0F;
  } else if (b_T22_im == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 6.0F;
  } else {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re /= 6.0F;
  }

  c_T12_re = H12 * H13;
  c_T22_re = H12 * 0.0F + 0.0F * H13;
  b_T22_im = H11 * H22;
  d_T11_im = H11 * 0.0F + 0.0F * H22;
  f_T11_re = b_T22_im * H33 - d_T11_im * 0.0F;
  d_T11_im = b_T22_im * 0.0F + d_T11_im * H33;
  if (d_T11_im == 0.0F) {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im = 0.0F;
  } else if (f_T11_re == 0.0F) {
    b_T22_im = 0.0F;
    d_T11_im /= 2.0F;
  } else {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im /= 2.0F;
  }

  fc14.re = ((((((T11.re - b_T13_re) - b_T12_re) - d_T11_re) + re) + e_T11_re) +
             (c_T12_re * H23 - c_T22_re * 0.0F)) + b_T22_im;
  fc14.im = ((((((T11.im - b_T13_im) - b_T12_im) - c_T11_im) + im) + c_a_re) +
             (c_T12_re * 0.0F + c_T22_re * H23)) + d_T11_im;
  T11 = c_mpower(fc14);
  f_T11_re = (((((a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  b_T22_im = (((((a_im + T12_im) + T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  if (T11.im == 0.0F) {
    if (b_T22_im == 0.0F) {
      a_re = f_T11_re / T11.re;
      a_im = 0.0F;
    } else if (f_T11_re == 0.0F) {
      a_re = 0.0F;
      a_im = b_T22_im / T11.re;
    } else {
      a_re = f_T11_re / T11.re;
      a_im = b_T22_im / T11.re;
    }
  } else if (T11.re == 0.0F) {
    if (f_T11_re == 0.0F) {
      a_re = b_T22_im / T11.im;
      a_im = 0.0F;
    } else if (b_T22_im == 0.0F) {
      a_re = 0.0F;
      a_im = -(f_T11_re / T11.im);
    } else {
      a_re = b_T22_im / T11.im;
      a_im = -(f_T11_re / T11.im);
    }
  } else {
    c_T12_re = (float)fabs(T11.re);
    V1_norm = (float)fabs(T11.im);
    if (c_T12_re > V1_norm) {
      s = T11.im / T11.re;
      V1_norm = T11.re + s * T11.im;
      a_re = (f_T11_re + s * b_T22_im) / V1_norm;
      a_im = (b_T22_im - s * f_T11_re) / V1_norm;
    } else if (V1_norm == c_T12_re) {
      if (T11.re > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (T11.im > 0.0F) {
        V1_norm = 0.5F;
      } else {
        V1_norm = -0.5F;
      }

      a_re = (f_T11_re * s + b_T22_im * V1_norm) / c_T12_re;
      a_im = (b_T22_im * s - f_T11_re * V1_norm) / c_T12_re;
    } else {
      s = T11.re / T11.im;
      V1_norm = T11.im + s * T11.re;
      a_re = (s * f_T11_re + b_T22_im) / V1_norm;
      a_im = (s * b_T22_im - f_T11_re) / V1_norm;
    }
  }

  re = 1.73205078F * (fc0.re - a_re);
  im = 1.73205078F * (fc0.im - a_im);
  b_re = re * 0.0F - im;
  im = re + im * 0.0F;
  if (im == 0.0F) {
    re = b_re / 2.0F;
    im = 0.0F;
  } else if (b_re == 0.0F) {
    re = 0.0F;
    im /= 2.0F;
  } else {
    re = b_re / 2.0F;
    im /= 2.0F;
  }

  a_re = g_a_re * g_a_re;
  a_im = g_a_re * 0.0F + 0.0F * g_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  if (T13_im == 0.0F) {
    T13_re /= 3.0F;
    T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  cb_a.re = (((((a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  cb_a.im = (((((a_im + T12_im) + T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  fc0 = mpower(cb_a);
  db_a.re = (f_a_re * f_a_re - c_a_im * c_a_im) - fc0.re;
  db_a.im = (f_a_re * c_a_im + c_a_im * f_a_re) - fc0.im;
  fc0 = b_mpower(db_a);
  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  b_T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    b_T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    b_T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  fb_T11.re = (H11 + H22) + H33;
  fb_T11.im = 0.0F;
  T11 = mpower(fb_T11);
  if (T11.im == 0.0F) {
    b_re = T11.re / 27.0F;
    b_im = 0.0F;
  } else if (T11.re == 0.0F) {
    b_re = 0.0F;
    b_im = T11.im / 27.0F;
  } else {
    b_re = T11.re / 27.0F;
    b_im = T11.im / 27.0F;
  }

  c_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  d_T11_re = c_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = c_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  e_T11_re = d_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = d_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (e_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  fc15.re = ((((((fc0.re - T13_re) - T12_re) - b_T11_re) + b_re) + c_T11_re) +
             (b_T12_re * H23 - b_T12_im * 0.0F)) + d_T11_re;
  fc15.im = ((((((fc0.im - T13_im) - T12_im) - T11_im) + b_im) + b_T11_im) +
             (b_T12_re * 0.0F + b_T12_im * H23)) + c_T11_im;
  fc0 = c_mpower(fc15);
  if (fc0.im == 0.0F) {
    b_re = fc0.re / 2.0F;
    b_im = 0.0F;
  } else if (fc0.re == 0.0F) {
    b_re = 0.0F;
    b_im = fc0.im / 2.0F;
  } else {
    b_re = fc0.re / 2.0F;
    b_im = fc0.im / 2.0F;
  }

  a_re = h_a_re * h_a_re;
  a_im = h_a_re * 0.0F + 0.0F * h_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  if (T13_im == 0.0F) {
    T13_re /= 3.0F;
    T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  b_a_re = j_a_re * j_a_re;
  b_a_im = j_a_re * 0.0F + 0.0F * j_a_re;
  if (b_a_im == 0.0F) {
    b_a_re /= 9.0F;
    b_a_im = 0.0F;
  } else if (b_a_re == 0.0F) {
    b_a_re = 0.0F;
    b_a_im /= 9.0F;
  } else {
    b_a_re /= 9.0F;
    b_a_im /= 9.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  if (b_T12_im == 0.0F) {
    b_T12_re /= 3.0F;
    b_T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 3.0F;
  } else {
    b_T12_re /= 3.0F;
    b_T12_im /= 3.0F;
  }

  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  if (b_T13_im == 0.0F) {
    b_T13_re /= 3.0F;
    b_T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 3.0F;
  } else {
    b_T13_re /= 3.0F;
    b_T13_im /= 3.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  if (b_T23_im == 0.0F) {
    b_T23_re /= 3.0F;
    b_T23_im = 0.0F;
  } else if (b_T23_re == 0.0F) {
    b_T23_re = 0.0F;
    b_T23_im /= 3.0F;
  } else {
    b_T23_re /= 3.0F;
    b_T23_im /= 3.0F;
  }

  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 3.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 3.0F;
  } else {
    d_T11_re /= 3.0F;
    c_T11_im /= 3.0F;
  }

  e_T11_re = H11 * H33;
  c_a_re = H11 * 0.0F + 0.0F * H33;
  if (c_a_re == 0.0F) {
    e_T11_re /= 3.0F;
    c_a_re = 0.0F;
  } else if (e_T11_re == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 3.0F;
  } else {
    e_T11_re /= 3.0F;
    c_a_re /= 3.0F;
  }

  c_T22_re = H22 * H33;
  b_T22_im = H22 * 0.0F + 0.0F * H33;
  if (b_T22_im == 0.0F) {
    c_T22_re /= 3.0F;
    b_T22_im = 0.0F;
  } else if (c_T22_re == 0.0F) {
    c_T22_re = 0.0F;
    b_T22_im /= 3.0F;
  } else {
    c_T22_re /= 3.0F;
    b_T22_im /= 3.0F;
  }

  eb_a.re = (((((b_a_re + b_T12_re) + b_T13_re) + b_T23_re) - d_T11_re) -
             e_T11_re) - c_T22_re;
  eb_a.im = (((((b_a_im + b_T12_im) + b_T13_im) + b_T23_im) - c_T11_im) - c_a_re)
    - b_T22_im;
  fc0 = mpower(eb_a);
  fb_a.re = (i_a_re * i_a_re - d_a_im * d_a_im) - fc0.re;
  fb_a.im = (i_a_re * d_a_im + d_a_im * i_a_re) - fc0.im;
  fc0 = b_mpower(fb_a);
  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  c_T13_re = b_T13_re * H22 - b_T13_im * 0.0F;
  b_T13_im = b_T13_re * 0.0F + b_T13_im * H22;
  if (b_T13_im == 0.0F) {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im = 0.0F;
  } else if (c_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 2.0F;
  } else {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im /= 2.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  c_T12_re = b_T12_re * H33 - b_T12_im * 0.0F;
  b_T12_im = b_T12_re * 0.0F + b_T12_im * H33;
  if (b_T12_im == 0.0F) {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im = 0.0F;
  } else if (c_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 2.0F;
  } else {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im /= 2.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  d_T11_re = H11 * b_T23_re - 0.0F * b_T23_im;
  c_T11_im = H11 * b_T23_im + 0.0F * b_T23_re;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re /= 2.0F;
    c_T11_im /= 2.0F;
  }

  gb_T11.re = (H11 + H22) + H33;
  gb_T11.im = 0.0F;
  T11 = mpower(gb_T11);
  if (T11.im == 0.0F) {
    b_a_re = T11.re / 27.0F;
    d_a_re = 0.0F;
  } else if (T11.re == 0.0F) {
    b_a_re = 0.0F;
    d_a_re = T11.im / 27.0F;
  } else {
    b_a_re = T11.re / 27.0F;
    d_a_re = T11.im / 27.0F;
  }

  e_T11_re = (H11 + H22) + H33;
  c_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  c_T22_re = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  b_T22_im = e_T11_re * c_T12_re - 0.0F * c_T22_re;
  c_a_re = e_T11_re * c_T22_re + 0.0F * c_T12_re;
  if (c_a_re == 0.0F) {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re = 0.0F;
  } else if (b_T22_im == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 6.0F;
  } else {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re /= 6.0F;
  }

  c_T12_re = H12 * H13;
  c_T22_re = H12 * 0.0F + 0.0F * H13;
  b_T22_im = H11 * H22;
  d_T11_im = H11 * 0.0F + 0.0F * H22;
  f_T11_re = b_T22_im * H33 - d_T11_im * 0.0F;
  d_T11_im = b_T22_im * 0.0F + d_T11_im * H33;
  if (d_T11_im == 0.0F) {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im = 0.0F;
  } else if (f_T11_re == 0.0F) {
    b_T22_im = 0.0F;
    d_T11_im /= 2.0F;
  } else {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im /= 2.0F;
  }

  fc16.re = ((((((fc0.re - b_T13_re) - b_T12_re) - d_T11_re) + b_a_re) +
              e_T11_re) + (c_T12_re * H23 - c_T22_re * 0.0F)) + b_T22_im;
  fc16.im = ((((((fc0.im - b_T13_im) - b_T12_im) - c_T11_im) + d_a_re) + c_a_re)
             + (c_T12_re * 0.0F + c_T22_re * H23)) + d_T11_im;
  fc0 = c_mpower(fc16);
  f_T11_re = (((((a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  b_T22_im = (((((a_im + T12_im) + T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  d_T11_im = 2.0F * fc0.re;
  c_T22_re = 2.0F * fc0.im;
  if (c_T22_re == 0.0F) {
    if (b_T22_im == 0.0F) {
      a_re = f_T11_re / d_T11_im;
      a_im = 0.0F;
    } else if (f_T11_re == 0.0F) {
      a_re = 0.0F;
      a_im = b_T22_im / d_T11_im;
    } else {
      a_re = f_T11_re / d_T11_im;
      a_im = b_T22_im / d_T11_im;
    }
  } else if (d_T11_im == 0.0F) {
    if (f_T11_re == 0.0F) {
      a_re = b_T22_im / c_T22_re;
      a_im = 0.0F;
    } else if (b_T22_im == 0.0F) {
      a_re = 0.0F;
      a_im = -(f_T11_re / c_T22_re);
    } else {
      a_re = b_T22_im / c_T22_re;
      a_im = -(f_T11_re / c_T22_re);
    }
  } else {
    c_T12_re = (float)fabs(d_T11_im);
    V1_norm = (float)fabs(c_T22_re);
    if (c_T12_re > V1_norm) {
      s = c_T22_re / d_T11_im;
      V1_norm = d_T11_im + s * c_T22_re;
      a_re = (f_T11_re + s * b_T22_im) / V1_norm;
      a_im = (b_T22_im - s * f_T11_re) / V1_norm;
    } else if (V1_norm == c_T12_re) {
      if (d_T11_im > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (c_T22_re > 0.0F) {
        V1_norm = 0.5F;
      } else {
        V1_norm = -0.5F;
      }

      a_re = (f_T11_re * s + b_T22_im * V1_norm) / c_T12_re;
      a_im = (b_T22_im * s - f_T11_re * V1_norm) / c_T12_re;
    } else {
      s = d_T11_im / c_T22_re;
      V1_norm = c_T22_re + s * d_T11_im;
      a_re = (s * f_T11_re + b_T22_im) / V1_norm;
      a_im = (s * b_T22_im - f_T11_re) / V1_norm;
    }
  }

  a_re = ((((T11_re + T22_re) + T33_re) + re) - b_re) - a_re;
  a_im = (im - b_im) - a_im;
  hb_T11.re = (H11 + H22) + H33;
  hb_T11.im = 0.0F;
  fc0 = mpower(hb_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  b_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  b_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  c_a_re = (H11 + H22) + H33;
  d_a_re = (H11 + H22) + H33;
  ib_T11.re = (H11 + H22) + H33;
  ib_T11.im = 0.0F;
  fc0 = mpower(ib_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  e_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  c_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  f_a_re = (H11 + H22) + H33;
  jb_T11.re = (H11 + H22) + H33;
  jb_T11.im = 0.0F;
  fc0 = mpower(jb_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  g_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  d_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  h_a_re = (H11 + H22) + H33;
  i_a_re = (H11 + H22) + H33;
  kb_T11.re = (H11 + H22) + H33;
  kb_T11.im = 0.0F;
  fc0 = mpower(kb_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  j_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  e_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  V23_re = (H11 + H22) + H33;
  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  T12_re = H12 * H33;
  T12_im = H12 * 0.0F + 0.0F * H33;
  fc0 = mpower(T23);
  T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * H13;
  T11_im = H11 * 0.0F + 0.0F * H13;
  c_T13_re = H13 * H22;
  v21_re = H13 * 0.0F + 0.0F * H22;
  f_T11_re = (((T13_re * H23 - T13_im * 0.0F) - (T12_re * H13 - T12_im * 0.0F))
              + fc0.re) - (T22_re * H23 - T22_im * 0.0F);
  b_T22_im = (((T13_re * 0.0F + T13_im * H23) - (T12_re * 0.0F + T12_im * H13))
              + fc0.im) - (T22_re * 0.0F + T22_im * H23);
  d_T11_im = (((H12 * b_T13_re - 0.0F * b_T13_im) - (H12 * T23_re - 0.0F *
    T23_im)) - (T11_re * H23 - T11_im * 0.0F)) + (c_T13_re * H23 - v21_re * 0.0F);
  c_T22_re = (((H12 * b_T13_im + 0.0F * b_T13_re) - (H12 * T23_im + 0.0F *
    T23_re)) - (T11_re * 0.0F + T11_im * H23)) + (c_T13_re * 0.0F + v21_re * H23);
  if (c_T22_re == 0.0F) {
    if (b_T22_im == 0.0F) {
      T13_re = f_T11_re / d_T11_im;
      T13_im = 0.0F;
    } else if (f_T11_re == 0.0F) {
      T13_re = 0.0F;
      T13_im = b_T22_im / d_T11_im;
    } else {
      T13_re = f_T11_re / d_T11_im;
      T13_im = b_T22_im / d_T11_im;
    }
  } else if (d_T11_im == 0.0F) {
    if (f_T11_re == 0.0F) {
      T13_re = b_T22_im / c_T22_re;
      T13_im = 0.0F;
    } else if (b_T22_im == 0.0F) {
      T13_re = 0.0F;
      T13_im = -(f_T11_re / c_T22_re);
    } else {
      T13_re = b_T22_im / c_T22_re;
      T13_im = -(f_T11_re / c_T22_re);
    }
  } else {
    c_T12_re = (float)fabs(d_T11_im);
    if (c_T12_re > c_T22_re) {
      s = c_T22_re / d_T11_im;
      V1_norm = d_T11_im + s * c_T22_re;
      T13_re = (f_T11_re + s * b_T22_im) / V1_norm;
      T13_im = (b_T22_im - s * f_T11_re) / V1_norm;
    } else if (c_T22_re == c_T12_re) {
      if (d_T11_im > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (f_T11_re * s + b_T22_im * -0.5F > 0.0F) {
        T13_re = ((real32_T)rtInf);
      } else if (f_T11_re * s + b_T22_im * -0.5F < 0.0F) {
        T13_re = ((real32_T)rtMinusInf);
      } else {
        T13_re = ((real32_T)rtNaN);
      }

      if (b_T22_im * s - f_T11_re * -0.5F > 0.0F) {
        T13_im = ((real32_T)rtInf);
      } else if (b_T22_im * s - f_T11_re * -0.5F < 0.0F) {
        T13_im = ((real32_T)rtMinusInf);
      } else {
        T13_im = ((real32_T)rtNaN);
      }
    } else {
      s = d_T11_im / c_T22_re;
      V1_norm = c_T22_re + s * d_T11_im;
      T13_re = (s * f_T11_re + b_T22_im) / V1_norm;
      T13_im = (s * b_T22_im - f_T11_re) / V1_norm;
    }
  }

  V1_norm = a_re * a_re - a_im * a_im;
  a_im = a_re * a_im + a_im * a_re;
  T23_re = H23 * V1_norm - 0.0F * a_im;
  T23_im = H23 * a_im + 0.0F * V1_norm;
  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * H13;
  T11_im = H11 * 0.0F + 0.0F * H13;
  c_T13_re = H13 * H22;
  v21_re = H13 * 0.0F + 0.0F * H22;
  d_T11_im = (((H12 * b_T13_re - 0.0F * b_T13_im) - (H12 * b_T23_re - 0.0F *
    b_T23_im)) - (T11_re * H23 - T11_im * 0.0F)) + (c_T13_re * H23 - v21_re *
    0.0F);
  c_T22_re = (((H12 * b_T13_im + 0.0F * b_T13_re) - (H12 * b_T23_im + 0.0F *
    b_T23_re)) - (T11_re * 0.0F + T11_im * H23)) + (c_T13_re * 0.0F + v21_re *
    H23);
  if (c_T22_re == 0.0F) {
    if (T23_im == 0.0F) {
      b_T23_re = T23_re / d_T11_im;
      T23_im = 0.0F;
    } else if (T23_re == 0.0F) {
      b_T23_re = 0.0F;
      T23_im /= d_T11_im;
    } else {
      b_T23_re = T23_re / d_T11_im;
      T23_im /= d_T11_im;
    }
  } else if (d_T11_im == 0.0F) {
    if (T23_re == 0.0F) {
      b_T23_re = T23_im / c_T22_re;
      T23_im = 0.0F;
    } else if (T23_im == 0.0F) {
      b_T23_re = 0.0F;
      T23_im = -(T23_re / c_T22_re);
    } else {
      b_T23_re = T23_im / c_T22_re;
      T23_im = -(T23_re / c_T22_re);
    }
  } else {
    c_T12_re = (float)fabs(d_T11_im);
    if (c_T12_re > c_T22_re) {
      s = c_T22_re / d_T11_im;
      V1_norm = d_T11_im + s * c_T22_re;
      b_T23_re = (T23_re + s * T23_im) / V1_norm;
      T23_im = (T23_im - s * T23_re) / V1_norm;
    } else if (c_T22_re == c_T12_re) {
      if (d_T11_im > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (T23_re * s + T23_im * -0.5F > 0.0F) {
        b_T23_re = ((real32_T)rtInf);
      } else if (T23_re * s + T23_im * -0.5F < 0.0F) {
        b_T23_re = ((real32_T)rtMinusInf);
      } else {
        b_T23_re = ((real32_T)rtNaN);
      }

      if (T23_im * s - T23_re * -0.5F > 0.0F) {
        T23_im = ((real32_T)rtInf);
      } else if (T23_im * s - T23_re * -0.5F < 0.0F) {
        T23_im = ((real32_T)rtMinusInf);
      } else {
        T23_im = ((real32_T)rtNaN);
      }
    } else {
      s = d_T11_im / c_T22_re;
      V1_norm = c_T22_re + s * d_T11_im;
      b_T23_re = (s * T23_re + T23_im) / V1_norm;
      T23_im = (s * T23_im - T23_re) / V1_norm;
    }
  }

  T11_re = H11 / 3.0F;
  T22_re = H22 / 3.0F;
  T33_re = H33 / 3.0F;
  a_re = c_a_re * c_a_re;
  a_im = c_a_re * 0.0F + 0.0F * c_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  if (b_T13_im == 0.0F) {
    b_T13_re /= 3.0F;
    b_T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 3.0F;
  } else {
    b_T13_re /= 3.0F;
    b_T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  if (b_T23_im == 0.0F) {
    T23_re /= 3.0F;
    b_T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    b_T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    b_T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  gb_a.re = (((((a_re + T12_re) + b_T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  gb_a.im = (((((a_im + T12_im) + b_T13_im) + b_T23_im) - T11_im) - b_T11_im) -
    T22_im;
  fc0 = mpower(gb_a);
  hb_a.re = (b_a_re * b_a_re - b_a_im * b_a_im) - fc0.re;
  hb_a.im = (b_a_re * b_a_im + b_a_im * b_a_re) - fc0.im;
  fc0 = b_mpower(hb_a);
  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  c_T13_re = b_T13_re * H22 - b_T13_im * 0.0F;
  b_T13_im = b_T13_re * 0.0F + b_T13_im * H22;
  if (b_T13_im == 0.0F) {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im = 0.0F;
  } else if (c_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 2.0F;
  } else {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  b_T11_re = H11 * T23_re - 0.0F * b_T23_im;
  T11_im = H11 * b_T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    b_T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    b_T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  lb_T11.re = (H11 + H22) + H33;
  lb_T11.im = 0.0F;
  T11 = mpower(lb_T11);
  if (T11.im == 0.0F) {
    re = T11.re / 27.0F;
    im = 0.0F;
  } else if (T11.re == 0.0F) {
    re = 0.0F;
    im = T11.im / 27.0F;
  } else {
    re = T11.re / 27.0F;
    im = T11.im / 27.0F;
  }

  c_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  d_T11_re = c_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = c_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  e_T11_re = d_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = d_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (e_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  fc17.re = ((((((fc0.re - b_T13_re) - T12_re) - b_T11_re) + re) + c_T11_re) +
             (b_T12_re * H23 - b_T12_im * 0.0F)) + d_T11_re;
  fc17.im = ((((((fc0.im - b_T13_im) - T12_im) - T11_im) + im) + b_T11_im) +
             (b_T12_re * 0.0F + b_T12_im * H23)) + c_T11_im;
  fc0 = c_mpower(fc17);
  a_re = d_a_re * d_a_re;
  a_im = d_a_re * 0.0F + 0.0F * d_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  if (b_T13_im == 0.0F) {
    b_T13_re /= 3.0F;
    b_T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 3.0F;
  } else {
    b_T13_re /= 3.0F;
    b_T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  if (b_T23_im == 0.0F) {
    T23_re /= 3.0F;
    b_T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    b_T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    b_T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  b_a_re = f_a_re * f_a_re;
  b_a_im = f_a_re * 0.0F + 0.0F * f_a_re;
  if (b_a_im == 0.0F) {
    b_a_re /= 9.0F;
    b_a_im = 0.0F;
  } else if (b_a_re == 0.0F) {
    b_a_re = 0.0F;
    b_a_im /= 9.0F;
  } else {
    b_a_re /= 9.0F;
    b_a_im /= 9.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  if (b_T12_im == 0.0F) {
    b_T12_re /= 3.0F;
    b_T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 3.0F;
  } else {
    b_T12_re /= 3.0F;
    b_T12_im /= 3.0F;
  }

  c_T13_re = H13 * H13;
  v21_re = H13 * 0.0F + 0.0F * H13;
  if (v21_re == 0.0F) {
    c_T13_re /= 3.0F;
    v21_re = 0.0F;
  } else if (c_T13_re == 0.0F) {
    c_T13_re = 0.0F;
    v21_re /= 3.0F;
  } else {
    c_T13_re /= 3.0F;
    v21_re /= 3.0F;
  }

  V1_norm = H23 * H23;
  s = H23 * 0.0F + 0.0F * H23;
  if (s == 0.0F) {
    V1_norm /= 3.0F;
    s = 0.0F;
  } else if (V1_norm == 0.0F) {
    V1_norm = 0.0F;
    s /= 3.0F;
  } else {
    V1_norm /= 3.0F;
    s /= 3.0F;
  }

  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 3.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 3.0F;
  } else {
    d_T11_re /= 3.0F;
    c_T11_im /= 3.0F;
  }

  e_T11_re = H11 * H33;
  c_a_re = H11 * 0.0F + 0.0F * H33;
  if (c_a_re == 0.0F) {
    e_T11_re /= 3.0F;
    c_a_re = 0.0F;
  } else if (e_T11_re == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 3.0F;
  } else {
    e_T11_re /= 3.0F;
    c_a_re /= 3.0F;
  }

  c_T22_re = H22 * H33;
  b_T22_im = H22 * 0.0F + 0.0F * H33;
  if (b_T22_im == 0.0F) {
    c_T22_re /= 3.0F;
    b_T22_im = 0.0F;
  } else if (c_T22_re == 0.0F) {
    c_T22_re = 0.0F;
    b_T22_im /= 3.0F;
  } else {
    c_T22_re /= 3.0F;
    b_T22_im /= 3.0F;
  }

  ib_a.re = (((((b_a_re + b_T12_re) + c_T13_re) + V1_norm) - d_T11_re) -
             e_T11_re) - c_T22_re;
  ib_a.im = (((((b_a_im + b_T12_im) + v21_re) + s) - c_T11_im) - c_a_re) -
    b_T22_im;
  T11 = mpower(ib_a);
  jb_a.re = (e_a_re * e_a_re - c_a_im * c_a_im) - T11.re;
  jb_a.im = (e_a_re * c_a_im + c_a_im * e_a_re) - T11.im;
  T11 = b_mpower(jb_a);
  c_T13_re = H13 * H13;
  v21_re = H13 * 0.0F + 0.0F * H13;
  V1_norm = c_T13_re * H22 - v21_re * 0.0F;
  v21_re = c_T13_re * 0.0F + v21_re * H22;
  if (v21_re == 0.0F) {
    c_T13_re = V1_norm / 2.0F;
    v21_re = 0.0F;
  } else if (V1_norm == 0.0F) {
    c_T13_re = 0.0F;
    v21_re /= 2.0F;
  } else {
    c_T13_re = V1_norm / 2.0F;
    v21_re /= 2.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  c_T12_re = b_T12_re * H33 - b_T12_im * 0.0F;
  b_T12_im = b_T12_re * 0.0F + b_T12_im * H33;
  if (b_T12_im == 0.0F) {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im = 0.0F;
  } else if (c_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 2.0F;
  } else {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im /= 2.0F;
  }

  V1_norm = H23 * H23;
  s = H23 * 0.0F + 0.0F * H23;
  d_T11_re = H11 * V1_norm - 0.0F * s;
  c_T11_im = H11 * s + 0.0F * V1_norm;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re /= 2.0F;
    c_T11_im /= 2.0F;
  }

  mb_T11.re = (H11 + H22) + H33;
  mb_T11.im = 0.0F;
  fc1 = mpower(mb_T11);
  if (fc1.im == 0.0F) {
    re = fc1.re / 27.0F;
    im = 0.0F;
  } else if (fc1.re == 0.0F) {
    re = 0.0F;
    im = fc1.im / 27.0F;
  } else {
    re = fc1.re / 27.0F;
    im = fc1.im / 27.0F;
  }

  e_T11_re = (H11 + H22) + H33;
  c_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  c_T22_re = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  b_T22_im = e_T11_re * c_T12_re - 0.0F * c_T22_re;
  c_a_re = e_T11_re * c_T22_re + 0.0F * c_T12_re;
  if (c_a_re == 0.0F) {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re = 0.0F;
  } else if (b_T22_im == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 6.0F;
  } else {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re /= 6.0F;
  }

  c_T12_re = H12 * H13;
  c_T22_re = H12 * 0.0F + 0.0F * H13;
  b_T22_im = H11 * H22;
  d_T11_im = H11 * 0.0F + 0.0F * H22;
  f_T11_re = b_T22_im * H33 - d_T11_im * 0.0F;
  d_T11_im = b_T22_im * 0.0F + d_T11_im * H33;
  if (d_T11_im == 0.0F) {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im = 0.0F;
  } else if (f_T11_re == 0.0F) {
    b_T22_im = 0.0F;
    d_T11_im /= 2.0F;
  } else {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im /= 2.0F;
  }

  fc18.re = ((((((T11.re - c_T13_re) - b_T12_re) - d_T11_re) + re) + e_T11_re) +
             (c_T12_re * H23 - c_T22_re * 0.0F)) + b_T22_im;
  fc18.im = ((((((T11.im - v21_re) - b_T12_im) - c_T11_im) + im) + c_a_re) +
             (c_T12_re * 0.0F + c_T22_re * H23)) + d_T11_im;
  T11 = c_mpower(fc18);
  f_T11_re = (((((a_re + T12_re) + b_T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  b_T22_im = (((((a_im + T12_im) + b_T13_im) + b_T23_im) - T11_im) - b_T11_im) -
    T22_im;
  if (T11.im == 0.0F) {
    if (b_T22_im == 0.0F) {
      a_re = f_T11_re / T11.re;
      a_im = 0.0F;
    } else if (f_T11_re == 0.0F) {
      a_re = 0.0F;
      a_im = b_T22_im / T11.re;
    } else {
      a_re = f_T11_re / T11.re;
      a_im = b_T22_im / T11.re;
    }
  } else if (T11.re == 0.0F) {
    if (f_T11_re == 0.0F) {
      a_re = b_T22_im / T11.im;
      a_im = 0.0F;
    } else if (b_T22_im == 0.0F) {
      a_re = 0.0F;
      a_im = -(f_T11_re / T11.im);
    } else {
      a_re = b_T22_im / T11.im;
      a_im = -(f_T11_re / T11.im);
    }
  } else {
    c_T12_re = (float)fabs(T11.re);
    V1_norm = (float)fabs(T11.im);
    if (c_T12_re > V1_norm) {
      s = T11.im / T11.re;
      V1_norm = T11.re + s * T11.im;
      a_re = (f_T11_re + s * b_T22_im) / V1_norm;
      a_im = (b_T22_im - s * f_T11_re) / V1_norm;
    } else if (V1_norm == c_T12_re) {
      if (T11.re > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (T11.im > 0.0F) {
        V1_norm = 0.5F;
      } else {
        V1_norm = -0.5F;
      }

      a_re = (f_T11_re * s + b_T22_im * V1_norm) / c_T12_re;
      a_im = (b_T22_im * s - f_T11_re * V1_norm) / c_T12_re;
    } else {
      s = T11.re / T11.im;
      V1_norm = T11.im + s * T11.re;
      a_re = (s * f_T11_re + b_T22_im) / V1_norm;
      a_im = (s * b_T22_im - f_T11_re) / V1_norm;
    }
  }

  re = 1.73205078F * (fc0.re - a_re);
  im = 1.73205078F * (fc0.im - a_im);
  b_re = re * 0.0F - im;
  im = re + im * 0.0F;
  if (im == 0.0F) {
    re = b_re / 2.0F;
    im = 0.0F;
  } else if (b_re == 0.0F) {
    re = 0.0F;
    im /= 2.0F;
  } else {
    re = b_re / 2.0F;
    im /= 2.0F;
  }

  a_re = h_a_re * h_a_re;
  a_im = h_a_re * 0.0F + 0.0F * h_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  if (b_T13_im == 0.0F) {
    b_T13_re /= 3.0F;
    b_T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 3.0F;
  } else {
    b_T13_re /= 3.0F;
    b_T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  if (b_T23_im == 0.0F) {
    T23_re /= 3.0F;
    b_T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    b_T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    b_T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  kb_a.re = (((((a_re + T12_re) + b_T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  kb_a.im = (((((a_im + T12_im) + b_T13_im) + b_T23_im) - T11_im) - b_T11_im) -
    T22_im;
  fc0 = mpower(kb_a);
  lb_a.re = (g_a_re * g_a_re - d_a_im * d_a_im) - fc0.re;
  lb_a.im = (g_a_re * d_a_im + d_a_im * g_a_re) - fc0.im;
  fc0 = b_mpower(lb_a);
  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  c_T13_re = b_T13_re * H22 - b_T13_im * 0.0F;
  b_T13_im = b_T13_re * 0.0F + b_T13_im * H22;
  if (b_T13_im == 0.0F) {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im = 0.0F;
  } else if (c_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 2.0F;
  } else {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  b_T11_re = H11 * T23_re - 0.0F * b_T23_im;
  T11_im = H11 * b_T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    b_T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    b_T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  nb_T11.re = (H11 + H22) + H33;
  nb_T11.im = 0.0F;
  T11 = mpower(nb_T11);
  if (T11.im == 0.0F) {
    b_re = T11.re / 27.0F;
    b_im = 0.0F;
  } else if (T11.re == 0.0F) {
    b_re = 0.0F;
    b_im = T11.im / 27.0F;
  } else {
    b_re = T11.re / 27.0F;
    b_im = T11.im / 27.0F;
  }

  c_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  d_T11_re = c_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = c_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  e_T11_re = d_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = d_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (e_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  fc19.re = ((((((fc0.re - b_T13_re) - T12_re) - b_T11_re) + b_re) + c_T11_re) +
             (b_T12_re * H23 - b_T12_im * 0.0F)) + d_T11_re;
  fc19.im = ((((((fc0.im - b_T13_im) - T12_im) - T11_im) + b_im) + b_T11_im) +
             (b_T12_re * 0.0F + b_T12_im * H23)) + c_T11_im;
  fc0 = c_mpower(fc19);
  if (fc0.im == 0.0F) {
    b_re = fc0.re / 2.0F;
    b_im = 0.0F;
  } else if (fc0.re == 0.0F) {
    b_re = 0.0F;
    b_im = fc0.im / 2.0F;
  } else {
    b_re = fc0.re / 2.0F;
    b_im = fc0.im / 2.0F;
  }

  a_re = i_a_re * i_a_re;
  a_im = i_a_re * 0.0F + 0.0F * i_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  if (b_T13_im == 0.0F) {
    b_T13_re /= 3.0F;
    b_T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 3.0F;
  } else {
    b_T13_re /= 3.0F;
    b_T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  if (b_T23_im == 0.0F) {
    T23_re /= 3.0F;
    b_T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    b_T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    b_T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  V1_norm = V23_re * V23_re;
  c_T12_re = V23_re * 0.0F + 0.0F * V23_re;
  if (c_T12_re == 0.0F) {
    V23_re = V1_norm / 9.0F;
    c_T12_re = 0.0F;
  } else if (V1_norm == 0.0F) {
    V23_re = 0.0F;
    c_T12_re /= 9.0F;
  } else {
    V23_re = V1_norm / 9.0F;
    c_T12_re /= 9.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  if (b_T12_im == 0.0F) {
    b_T12_re /= 3.0F;
    b_T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 3.0F;
  } else {
    b_T12_re /= 3.0F;
    b_T12_im /= 3.0F;
  }

  c_T13_re = H13 * H13;
  v21_re = H13 * 0.0F + 0.0F * H13;
  if (v21_re == 0.0F) {
    c_T13_re /= 3.0F;
    v21_re = 0.0F;
  } else if (c_T13_re == 0.0F) {
    c_T13_re = 0.0F;
    v21_re /= 3.0F;
  } else {
    c_T13_re /= 3.0F;
    v21_re /= 3.0F;
  }

  V1_norm = H23 * H23;
  s = H23 * 0.0F + 0.0F * H23;
  if (s == 0.0F) {
    V1_norm /= 3.0F;
    s = 0.0F;
  } else if (V1_norm == 0.0F) {
    V1_norm = 0.0F;
    s /= 3.0F;
  } else {
    V1_norm /= 3.0F;
    s /= 3.0F;
  }

  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 3.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 3.0F;
  } else {
    d_T11_re /= 3.0F;
    c_T11_im /= 3.0F;
  }

  e_T11_re = H11 * H33;
  c_a_re = H11 * 0.0F + 0.0F * H33;
  if (c_a_re == 0.0F) {
    e_T11_re /= 3.0F;
    c_a_re = 0.0F;
  } else if (e_T11_re == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 3.0F;
  } else {
    e_T11_re /= 3.0F;
    c_a_re /= 3.0F;
  }

  c_T22_re = H22 * H33;
  b_T22_im = H22 * 0.0F + 0.0F * H33;
  if (b_T22_im == 0.0F) {
    c_T22_re /= 3.0F;
    b_T22_im = 0.0F;
  } else if (c_T22_re == 0.0F) {
    c_T22_re = 0.0F;
    b_T22_im /= 3.0F;
  } else {
    c_T22_re /= 3.0F;
    b_T22_im /= 3.0F;
  }

  b_V23.re = (((((V23_re + b_T12_re) + c_T13_re) + V1_norm) - d_T11_re) -
              e_T11_re) - c_T22_re;
  b_V23.im = (((((c_T12_re + b_T12_im) + v21_re) + s) - c_T11_im) - c_a_re) -
    b_T22_im;
  fc0 = mpower(b_V23);
  mb_a.re = (j_a_re * j_a_re - e_a_im * e_a_im) - fc0.re;
  mb_a.im = (j_a_re * e_a_im + e_a_im * j_a_re) - fc0.im;
  fc0 = b_mpower(mb_a);
  c_T13_re = H13 * H13;
  v21_re = H13 * 0.0F + 0.0F * H13;
  V1_norm = c_T13_re * H22 - v21_re * 0.0F;
  v21_re = c_T13_re * 0.0F + v21_re * H22;
  if (v21_re == 0.0F) {
    c_T13_re = V1_norm / 2.0F;
    v21_re = 0.0F;
  } else if (V1_norm == 0.0F) {
    c_T13_re = 0.0F;
    v21_re /= 2.0F;
  } else {
    c_T13_re = V1_norm / 2.0F;
    v21_re /= 2.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  c_T12_re = b_T12_re * H33 - b_T12_im * 0.0F;
  b_T12_im = b_T12_re * 0.0F + b_T12_im * H33;
  if (b_T12_im == 0.0F) {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im = 0.0F;
  } else if (c_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 2.0F;
  } else {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im /= 2.0F;
  }

  V1_norm = H23 * H23;
  s = H23 * 0.0F + 0.0F * H23;
  d_T11_re = H11 * V1_norm - 0.0F * s;
  c_T11_im = H11 * s + 0.0F * V1_norm;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re /= 2.0F;
    c_T11_im /= 2.0F;
  }

  ob_T11.re = (H11 + H22) + H33;
  ob_T11.im = 0.0F;
  T11 = mpower(ob_T11);
  if (T11.im == 0.0F) {
    b_a_re = T11.re / 27.0F;
    d_a_re = 0.0F;
  } else if (T11.re == 0.0F) {
    b_a_re = 0.0F;
    d_a_re = T11.im / 27.0F;
  } else {
    b_a_re = T11.re / 27.0F;
    d_a_re = T11.im / 27.0F;
  }

  e_T11_re = (H11 + H22) + H33;
  c_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  c_T22_re = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  b_T22_im = e_T11_re * c_T12_re - 0.0F * c_T22_re;
  c_a_re = e_T11_re * c_T22_re + 0.0F * c_T12_re;
  if (c_a_re == 0.0F) {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re = 0.0F;
  } else if (b_T22_im == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 6.0F;
  } else {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re /= 6.0F;
  }

  c_T12_re = H12 * H13;
  c_T22_re = H12 * 0.0F + 0.0F * H13;
  b_T22_im = H11 * H22;
  d_T11_im = H11 * 0.0F + 0.0F * H22;
  f_T11_re = b_T22_im * H33 - d_T11_im * 0.0F;
  d_T11_im = b_T22_im * 0.0F + d_T11_im * H33;
  if (d_T11_im == 0.0F) {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im = 0.0F;
  } else if (f_T11_re == 0.0F) {
    b_T22_im = 0.0F;
    d_T11_im /= 2.0F;
  } else {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im /= 2.0F;
  }

  fc20.re = ((((((fc0.re - c_T13_re) - b_T12_re) - d_T11_re) + b_a_re) +
              e_T11_re) + (c_T12_re * H23 - c_T22_re * 0.0F)) + b_T22_im;
  fc20.im = ((((((fc0.im - v21_re) - b_T12_im) - c_T11_im) + d_a_re) + c_a_re) +
             (c_T12_re * 0.0F + c_T22_re * H23)) + d_T11_im;
  fc0 = c_mpower(fc20);
  f_T11_re = (((((a_re + T12_re) + b_T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  b_T22_im = (((((a_im + T12_im) + b_T13_im) + b_T23_im) - T11_im) - b_T11_im) -
    T22_im;
  d_T11_im = 2.0F * fc0.re;
  c_T22_re = 2.0F * fc0.im;
  if (c_T22_re == 0.0F) {
    if (b_T22_im == 0.0F) {
      a_re = f_T11_re / d_T11_im;
      a_im = 0.0F;
    } else if (f_T11_re == 0.0F) {
      a_re = 0.0F;
      a_im = b_T22_im / d_T11_im;
    } else {
      a_re = f_T11_re / d_T11_im;
      a_im = b_T22_im / d_T11_im;
    }
  } else if (d_T11_im == 0.0F) {
    if (f_T11_re == 0.0F) {
      a_re = b_T22_im / c_T22_re;
      a_im = 0.0F;
    } else if (b_T22_im == 0.0F) {
      a_re = 0.0F;
      a_im = -(f_T11_re / c_T22_re);
    } else {
      a_re = b_T22_im / c_T22_re;
      a_im = -(f_T11_re / c_T22_re);
    }
  } else {
    c_T12_re = (float)fabs(d_T11_im);
    V1_norm = (float)fabs(c_T22_re);
    if (c_T12_re > V1_norm) {
      s = c_T22_re / d_T11_im;
      V1_norm = d_T11_im + s * c_T22_re;
      a_re = (f_T11_re + s * b_T22_im) / V1_norm;
      a_im = (b_T22_im - s * f_T11_re) / V1_norm;
    } else if (V1_norm == c_T12_re) {
      if (d_T11_im > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (c_T22_re > 0.0F) {
        V1_norm = 0.5F;
      } else {
        V1_norm = -0.5F;
      }

      a_re = (f_T11_re * s + b_T22_im * V1_norm) / c_T12_re;
      a_im = (b_T22_im * s - f_T11_re * V1_norm) / c_T12_re;
    } else {
      s = d_T11_im / c_T22_re;
      V1_norm = c_T22_re + s * d_T11_im;
      a_re = (s * f_T11_re + b_T22_im) / V1_norm;
      a_im = (s * b_T22_im - f_T11_re) / V1_norm;
    }
  }

  T12_re = (H12 * H13 + H22 * H23) + H23 * H33;
  T12_im = ((H12 * 0.0F + 0.0F * H13) + (H22 * 0.0F + 0.0F * H23)) + (H23 * 0.0F
    + 0.0F * H33);
  T11_re = ((((T11_re + T22_re) + T33_re) + re) - b_re) - a_re;
  T11_im = (im - b_im) - a_im;
  b_T12_re = T12_re * T11_re - T12_im * T11_im;
  T12_im = T12_re * T11_im + T12_im * T11_re;
  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * H13;
  T11_im = H11 * 0.0F + 0.0F * H13;
  c_T13_re = H13 * H22;
  v21_re = H13 * 0.0F + 0.0F * H22;
  d_T11_im = (((H12 * b_T13_re - 0.0F * b_T13_im) - (H12 * T23_re - 0.0F *
    b_T23_im)) - (T11_re * H23 - T11_im * 0.0F)) + (c_T13_re * H23 - v21_re *
    0.0F);
  c_T22_re = (((H12 * b_T13_im + 0.0F * b_T13_re) - (H12 * b_T23_im + 0.0F *
    T23_re)) - (T11_re * 0.0F + T11_im * H23)) + (c_T13_re * 0.0F + v21_re * H23);
  if (c_T22_re == 0.0F) {
    if (T12_im == 0.0F) {
      T12_re = b_T12_re / d_T11_im;
      T12_im = 0.0F;
    } else if (b_T12_re == 0.0F) {
      T12_re = 0.0F;
      T12_im /= d_T11_im;
    } else {
      T12_re = b_T12_re / d_T11_im;
      T12_im /= d_T11_im;
    }
  } else if (d_T11_im == 0.0F) {
    if (b_T12_re == 0.0F) {
      T12_re = T12_im / c_T22_re;
      T12_im = 0.0F;
    } else if (T12_im == 0.0F) {
      T12_re = 0.0F;
      T12_im = -(b_T12_re / c_T22_re);
    } else {
      T12_re = T12_im / c_T22_re;
      T12_im = -(b_T12_re / c_T22_re);
    }
  } else {
    c_T12_re = (float)fabs(d_T11_im);
    if (c_T12_re > c_T22_re) {
      s = c_T22_re / d_T11_im;
      V1_norm = d_T11_im + s * c_T22_re;
      T12_re = (b_T12_re + s * T12_im) / V1_norm;
      T12_im = (T12_im - s * b_T12_re) / V1_norm;
    } else if (c_T22_re == c_T12_re) {
      if (d_T11_im > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (b_T12_re * s + T12_im * -0.5F > 0.0F) {
        T12_re = ((real32_T)rtInf);
      } else if (b_T12_re * s + T12_im * -0.5F < 0.0F) {
        T12_re = ((real32_T)rtMinusInf);
      } else {
        T12_re = ((real32_T)rtNaN);
      }

      if (T12_im * s - b_T12_re * -0.5F > 0.0F) {
        T12_im = ((real32_T)rtInf);
      } else if (T12_im * s - b_T12_re * -0.5F < 0.0F) {
        T12_im = ((real32_T)rtMinusInf);
      } else {
        T12_im = ((real32_T)rtNaN);
      }
    } else {
      s = d_T11_im / c_T22_re;
      V1_norm = c_T22_re + s * d_T11_im;
      T12_re = (s * b_T12_re + T12_im) / V1_norm;
      T12_im = (s * T12_im - b_T12_re) / V1_norm;
    }
  }

  V13_re = (T13_re - b_T23_re) + T12_re;
  V13_im = (T13_im - T23_im) + T12_im;
  pb_T11.re = (H11 + H22) + H33;
  pb_T11.im = 0.0F;
  fc0 = mpower(pb_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
           b_T12_im * 0.0F)) + c_T11_re;
  a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F +
           b_T12_im * H23)) + c_T11_im;
  b_a_re = (H11 + H22) + H33;
  c_a_re = (H11 + H22) + H33;
  qb_T11.re = (H11 + H22) + H33;
  qb_T11.im = 0.0F;
  fc0 = mpower(qb_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  d_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  b_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  e_a_re = (H11 + H22) + H33;
  T11_re = H11 / 3.0F;
  T22_re = H22 / 3.0F;
  T33_re = H33 / 3.0F;
  f_a_re = b_a_re * b_a_re;
  c_a_im = b_a_re * 0.0F + 0.0F * b_a_re;
  if (c_a_im == 0.0F) {
    b_a_re = f_a_re / 9.0F;
    c_a_im = 0.0F;
  } else if (f_a_re == 0.0F) {
    b_a_re = 0.0F;
    c_a_im /= 9.0F;
  } else {
    b_a_re = f_a_re / 9.0F;
    c_a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  if (T13_im == 0.0F) {
    T13_re /= 3.0F;
    T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  nb_a.re = (((((b_a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  nb_a.im = (((((c_a_im + T12_im) + T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  fc0 = mpower(nb_a);
  ob_a.re = (a_re * a_re - a_im * a_im) - fc0.re;
  ob_a.im = (a_re * a_im + a_im * a_re) - fc0.im;
  fc0 = b_mpower(ob_a);
  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  b_T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    b_T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    b_T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  rb_T11.re = (H11 + H22) + H33;
  rb_T11.im = 0.0F;
  T11 = mpower(rb_T11);
  if (T11.im == 0.0F) {
    re = T11.re / 27.0F;
    im = 0.0F;
  } else if (T11.re == 0.0F) {
    re = 0.0F;
    im = T11.im / 27.0F;
  } else {
    re = T11.re / 27.0F;
    im = T11.im / 27.0F;
  }

  c_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  d_T11_re = c_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = c_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  e_T11_re = d_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = d_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (e_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  fc21.re = ((((((fc0.re - T13_re) - T12_re) - b_T11_re) + re) + c_T11_re) +
             (b_T12_re * H23 - b_T12_im * 0.0F)) + d_T11_re;
  fc21.im = ((((((fc0.im - T13_im) - T12_im) - T11_im) + im) + b_T11_im) +
             (b_T12_re * 0.0F + b_T12_im * H23)) + c_T11_im;
  fc0 = c_mpower(fc21);
  a_re = c_a_re * c_a_re;
  a_im = c_a_re * 0.0F + 0.0F * c_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  if (T13_im == 0.0F) {
    T13_re /= 3.0F;
    T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  b_a_re = e_a_re * e_a_re;
  c_a_im = e_a_re * 0.0F + 0.0F * e_a_re;
  if (c_a_im == 0.0F) {
    b_a_re /= 9.0F;
    c_a_im = 0.0F;
  } else if (b_a_re == 0.0F) {
    b_a_re = 0.0F;
    c_a_im /= 9.0F;
  } else {
    b_a_re /= 9.0F;
    c_a_im /= 9.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  if (b_T12_im == 0.0F) {
    b_T12_re /= 3.0F;
    b_T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 3.0F;
  } else {
    b_T12_re /= 3.0F;
    b_T12_im /= 3.0F;
  }

  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  if (b_T13_im == 0.0F) {
    b_T13_re /= 3.0F;
    b_T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 3.0F;
  } else {
    b_T13_re /= 3.0F;
    b_T13_im /= 3.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  if (b_T23_im == 0.0F) {
    b_T23_re /= 3.0F;
    b_T23_im = 0.0F;
  } else if (b_T23_re == 0.0F) {
    b_T23_re = 0.0F;
    b_T23_im /= 3.0F;
  } else {
    b_T23_re /= 3.0F;
    b_T23_im /= 3.0F;
  }

  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 3.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 3.0F;
  } else {
    d_T11_re /= 3.0F;
    c_T11_im /= 3.0F;
  }

  e_T11_re = H11 * H33;
  c_a_re = H11 * 0.0F + 0.0F * H33;
  if (c_a_re == 0.0F) {
    e_T11_re /= 3.0F;
    c_a_re = 0.0F;
  } else if (e_T11_re == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 3.0F;
  } else {
    e_T11_re /= 3.0F;
    c_a_re /= 3.0F;
  }

  c_T22_re = H22 * H33;
  b_T22_im = H22 * 0.0F + 0.0F * H33;
  if (b_T22_im == 0.0F) {
    c_T22_re /= 3.0F;
    b_T22_im = 0.0F;
  } else if (c_T22_re == 0.0F) {
    c_T22_re = 0.0F;
    b_T22_im /= 3.0F;
  } else {
    c_T22_re /= 3.0F;
    b_T22_im /= 3.0F;
  }

  pb_a.re = (((((b_a_re + b_T12_re) + b_T13_re) + b_T23_re) - d_T11_re) -
             e_T11_re) - c_T22_re;
  pb_a.im = (((((c_a_im + b_T12_im) + b_T13_im) + b_T23_im) - c_T11_im) - c_a_re)
    - b_T22_im;
  T11 = mpower(pb_a);
  qb_a.re = (d_a_re * d_a_re - b_a_im * b_a_im) - T11.re;
  qb_a.im = (d_a_re * b_a_im + b_a_im * d_a_re) - T11.im;
  T11 = b_mpower(qb_a);
  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  c_T13_re = b_T13_re * H22 - b_T13_im * 0.0F;
  b_T13_im = b_T13_re * 0.0F + b_T13_im * H22;
  if (b_T13_im == 0.0F) {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im = 0.0F;
  } else if (c_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 2.0F;
  } else {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im /= 2.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  c_T12_re = b_T12_re * H33 - b_T12_im * 0.0F;
  b_T12_im = b_T12_re * 0.0F + b_T12_im * H33;
  if (b_T12_im == 0.0F) {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im = 0.0F;
  } else if (c_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 2.0F;
  } else {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im /= 2.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  d_T11_re = H11 * b_T23_re - 0.0F * b_T23_im;
  c_T11_im = H11 * b_T23_im + 0.0F * b_T23_re;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re /= 2.0F;
    c_T11_im /= 2.0F;
  }

  sb_T11.re = (H11 + H22) + H33;
  sb_T11.im = 0.0F;
  fc1 = mpower(sb_T11);
  if (fc1.im == 0.0F) {
    re = fc1.re / 27.0F;
    im = 0.0F;
  } else if (fc1.re == 0.0F) {
    re = 0.0F;
    im = fc1.im / 27.0F;
  } else {
    re = fc1.re / 27.0F;
    im = fc1.im / 27.0F;
  }

  e_T11_re = (H11 + H22) + H33;
  c_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  c_T22_re = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  b_T22_im = e_T11_re * c_T12_re - 0.0F * c_T22_re;
  c_a_re = e_T11_re * c_T22_re + 0.0F * c_T12_re;
  if (c_a_re == 0.0F) {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re = 0.0F;
  } else if (b_T22_im == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 6.0F;
  } else {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re /= 6.0F;
  }

  c_T12_re = H12 * H13;
  c_T22_re = H12 * 0.0F + 0.0F * H13;
  b_T22_im = H11 * H22;
  d_T11_im = H11 * 0.0F + 0.0F * H22;
  f_T11_re = b_T22_im * H33 - d_T11_im * 0.0F;
  d_T11_im = b_T22_im * 0.0F + d_T11_im * H33;
  if (d_T11_im == 0.0F) {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im = 0.0F;
  } else if (f_T11_re == 0.0F) {
    b_T22_im = 0.0F;
    d_T11_im /= 2.0F;
  } else {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im /= 2.0F;
  }

  fc22.re = ((((((T11.re - b_T13_re) - b_T12_re) - d_T11_re) + re) + e_T11_re) +
             (c_T12_re * H23 - c_T22_re * 0.0F)) + b_T22_im;
  fc22.im = ((((((T11.im - b_T13_im) - b_T12_im) - c_T11_im) + im) + c_a_re) +
             (c_T12_re * 0.0F + c_T22_re * H23)) + d_T11_im;
  T11 = c_mpower(fc22);
  f_T11_re = (((((a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  b_T22_im = (((((a_im + T12_im) + T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  if (T11.im == 0.0F) {
    if (b_T22_im == 0.0F) {
      a_re = f_T11_re / T11.re;
      a_im = 0.0F;
    } else if (f_T11_re == 0.0F) {
      a_re = 0.0F;
      a_im = b_T22_im / T11.re;
    } else {
      a_re = f_T11_re / T11.re;
      a_im = b_T22_im / T11.re;
    }
  } else if (T11.re == 0.0F) {
    if (f_T11_re == 0.0F) {
      a_re = b_T22_im / T11.im;
      a_im = 0.0F;
    } else if (b_T22_im == 0.0F) {
      a_re = 0.0F;
      a_im = -(f_T11_re / T11.im);
    } else {
      a_re = b_T22_im / T11.im;
      a_im = -(f_T11_re / T11.im);
    }
  } else {
    c_T12_re = (float)fabs(T11.re);
    V1_norm = (float)fabs(T11.im);
    if (c_T12_re > V1_norm) {
      s = T11.im / T11.re;
      V1_norm = T11.re + s * T11.im;
      a_re = (f_T11_re + s * b_T22_im) / V1_norm;
      a_im = (b_T22_im - s * f_T11_re) / V1_norm;
    } else if (V1_norm == c_T12_re) {
      if (T11.re > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (T11.im > 0.0F) {
        V1_norm = 0.5F;
      } else {
        V1_norm = -0.5F;
      }

      a_re = (f_T11_re * s + b_T22_im * V1_norm) / c_T12_re;
      a_im = (b_T22_im * s - f_T11_re * V1_norm) / c_T12_re;
    } else {
      s = T11.re / T11.im;
      V1_norm = T11.im + s * T11.re;
      a_re = (s * f_T11_re + b_T22_im) / V1_norm;
      a_im = (s * b_T22_im - f_T11_re) / V1_norm;
    }
  }

  a_re += ((T11_re + T22_re) + T33_re) + fc0.re;
  a_im += fc0.im;
  tb_T11.re = (H11 + H22) + H33;
  tb_T11.im = 0.0F;
  fc0 = mpower(tb_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  b_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  b_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  c_a_re = (H11 + H22) + H33;
  d_a_re = (H11 + H22) + H33;
  ub_T11.re = (H11 + H22) + H33;
  ub_T11.im = 0.0F;
  fc0 = mpower(ub_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  e_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  c_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  f_a_re = (H11 + H22) + H33;
  g_a_re = a_re * a_re - a_im * a_im;
  a_im = a_re * a_im + a_im * a_re;
  T13_re = H13 * g_a_re - 0.0F * a_im;
  T13_im = H13 * a_im + 0.0F * g_a_re;
  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * H13;
  T11_im = H11 * 0.0F + 0.0F * H13;
  c_T13_re = H13 * H22;
  v21_re = H13 * 0.0F + 0.0F * H22;
  d_T11_im = (((H12 * b_T13_re - 0.0F * b_T13_im) - (H12 * T23_re - 0.0F *
    T23_im)) - (T11_re * H23 - T11_im * 0.0F)) + (c_T13_re * H23 - v21_re * 0.0F);
  c_T22_re = (((H12 * b_T13_im + 0.0F * b_T13_re) - (H12 * T23_im + 0.0F *
    T23_re)) - (T11_re * 0.0F + T11_im * H23)) + (c_T13_re * 0.0F + v21_re * H23);
  if (c_T22_re == 0.0F) {
    if (T13_im == 0.0F) {
      b_T13_re = T13_re / d_T11_im;
      T13_im = 0.0F;
    } else if (T13_re == 0.0F) {
      b_T13_re = 0.0F;
      T13_im /= d_T11_im;
    } else {
      b_T13_re = T13_re / d_T11_im;
      T13_im /= d_T11_im;
    }
  } else if (d_T11_im == 0.0F) {
    if (T13_re == 0.0F) {
      b_T13_re = T13_im / c_T22_re;
      T13_im = 0.0F;
    } else if (T13_im == 0.0F) {
      b_T13_re = 0.0F;
      T13_im = -(T13_re / c_T22_re);
    } else {
      b_T13_re = T13_im / c_T22_re;
      T13_im = -(T13_re / c_T22_re);
    }
  } else {
    c_T12_re = (float)fabs(d_T11_im);
    if (c_T12_re > c_T22_re) {
      s = c_T22_re / d_T11_im;
      V1_norm = d_T11_im + s * c_T22_re;
      b_T13_re = (T13_re + s * T13_im) / V1_norm;
      T13_im = (T13_im - s * T13_re) / V1_norm;
    } else if (c_T22_re == c_T12_re) {
      if (d_T11_im > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (T13_re * s + T13_im * -0.5F > 0.0F) {
        b_T13_re = ((real32_T)rtInf);
      } else if (T13_re * s + T13_im * -0.5F < 0.0F) {
        b_T13_re = ((real32_T)rtMinusInf);
      } else {
        b_T13_re = ((real32_T)rtNaN);
      }

      if (T13_im * s - T13_re * -0.5F > 0.0F) {
        T13_im = ((real32_T)rtInf);
      } else if (T13_im * s - T13_re * -0.5F < 0.0F) {
        T13_im = ((real32_T)rtMinusInf);
      } else {
        T13_im = ((real32_T)rtNaN);
      }
    } else {
      s = d_T11_im / c_T22_re;
      V1_norm = c_T22_re + s * d_T11_im;
      b_T13_re = (s * T13_re + T13_im) / V1_norm;
      T13_im = (s * T13_im - T13_re) / V1_norm;
    }
  }

  fc0 = mpower(T13);
  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * H33;
  T11_im = H11 * 0.0F + 0.0F * H33;
  T12_re = H12 * H33;
  T12_im = H12 * 0.0F + 0.0F * H33;
  T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  b_T11_re = H11 * H13;
  b_T11_im = H11 * 0.0F + 0.0F * H13;
  c_T13_re = H13 * H22;
  v21_re = H13 * 0.0F + 0.0F * H22;
  f_T11_re = ((fc0.re + (H13 * T23_re - 0.0F * T23_im)) - (T11_re * H13 - T11_im
    * 0.0F)) - (T12_re * H23 - T12_im * 0.0F);
  b_T22_im = ((fc0.im + (H13 * T23_im + 0.0F * T23_re)) - (T11_re * 0.0F +
    T11_im * H13)) - (T12_re * 0.0F + T12_im * H23);
  d_T11_im = (((H12 * T13_re - 0.0F * b_T13_im) - (H12 * b_T23_re - 0.0F *
    b_T23_im)) - (b_T11_re * H23 - b_T11_im * 0.0F)) + (c_T13_re * H23 - v21_re *
    0.0F);
  c_T22_re = (((H12 * b_T13_im + 0.0F * T13_re) - (H12 * b_T23_im + 0.0F *
    b_T23_re)) - (b_T11_re * 0.0F + b_T11_im * H23)) + (c_T13_re * 0.0F + v21_re
    * H23);
  if (c_T22_re == 0.0F) {
    if (b_T22_im == 0.0F) {
      re = f_T11_re / d_T11_im;
      im = 0.0F;
    } else if (f_T11_re == 0.0F) {
      re = 0.0F;
      im = b_T22_im / d_T11_im;
    } else {
      re = f_T11_re / d_T11_im;
      im = b_T22_im / d_T11_im;
    }
  } else if (d_T11_im == 0.0F) {
    if (f_T11_re == 0.0F) {
      re = b_T22_im / c_T22_re;
      im = 0.0F;
    } else if (b_T22_im == 0.0F) {
      re = 0.0F;
      im = -(f_T11_re / c_T22_re);
    } else {
      re = b_T22_im / c_T22_re;
      im = -(f_T11_re / c_T22_re);
    }
  } else {
    c_T12_re = (float)fabs(d_T11_im);
    if (c_T12_re > c_T22_re) {
      s = c_T22_re / d_T11_im;
      V1_norm = d_T11_im + s * c_T22_re;
      re = (f_T11_re + s * b_T22_im) / V1_norm;
      im = (b_T22_im - s * f_T11_re) / V1_norm;
    } else if (c_T22_re == c_T12_re) {
      if (d_T11_im > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (f_T11_re * s + b_T22_im * -0.5F > 0.0F) {
        re = ((real32_T)rtInf);
      } else if (f_T11_re * s + b_T22_im * -0.5F < 0.0F) {
        re = ((real32_T)rtMinusInf);
      } else {
        re = ((real32_T)rtNaN);
      }

      if (b_T22_im * s - f_T11_re * -0.5F > 0.0F) {
        im = ((real32_T)rtInf);
      } else if (b_T22_im * s - f_T11_re * -0.5F < 0.0F) {
        im = ((real32_T)rtMinusInf);
      } else {
        im = ((real32_T)rtNaN);
      }
    } else {
      s = d_T11_im / c_T22_re;
      V1_norm = c_T22_re + s * d_T11_im;
      re = (s * f_T11_re + b_T22_im) / V1_norm;
      im = (s * b_T22_im - f_T11_re) / V1_norm;
    }
  }

  T11_re = H11 / 3.0F;
  T22_re = H22 / 3.0F;
  T33_re = H33 / 3.0F;
  a_re = c_a_re * c_a_re;
  a_im = c_a_re * 0.0F + 0.0F * c_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  if (b_T13_im == 0.0F) {
    T13_re /= 3.0F;
    b_T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    b_T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    b_T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  rb_a.re = (((((a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  rb_a.im = (((((a_im + T12_im) + b_T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  fc0 = mpower(rb_a);
  sb_a.re = (b_a_re * b_a_re - b_a_im * b_a_im) - fc0.re;
  sb_a.im = (b_a_re * b_a_im + b_a_im * b_a_re) - fc0.im;
  fc0 = b_mpower(sb_a);
  T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  c_T13_re = T13_re * H22 - b_T13_im * 0.0F;
  b_T13_im = T13_re * 0.0F + b_T13_im * H22;
  if (b_T13_im == 0.0F) {
    T13_re = c_T13_re / 2.0F;
    b_T13_im = 0.0F;
  } else if (c_T13_re == 0.0F) {
    T13_re = 0.0F;
    b_T13_im /= 2.0F;
  } else {
    T13_re = c_T13_re / 2.0F;
    b_T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  b_T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    b_T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    b_T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  vb_T11.re = (H11 + H22) + H33;
  vb_T11.im = 0.0F;
  T11 = mpower(vb_T11);
  if (T11.im == 0.0F) {
    b_re = T11.re / 27.0F;
    b_im = 0.0F;
  } else if (T11.re == 0.0F) {
    b_re = 0.0F;
    b_im = T11.im / 27.0F;
  } else {
    b_re = T11.re / 27.0F;
    b_im = T11.im / 27.0F;
  }

  c_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  d_T11_re = c_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = c_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  e_T11_re = d_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = d_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (e_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  fc23.re = ((((((fc0.re - T13_re) - T12_re) - b_T11_re) + b_re) + c_T11_re) +
             (b_T12_re * H23 - b_T12_im * 0.0F)) + d_T11_re;
  fc23.im = ((((((fc0.im - b_T13_im) - T12_im) - T11_im) + b_im) + b_T11_im) +
             (b_T12_re * 0.0F + b_T12_im * H23)) + c_T11_im;
  fc0 = c_mpower(fc23);
  a_re = d_a_re * d_a_re;
  a_im = d_a_re * 0.0F + 0.0F * d_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  if (b_T13_im == 0.0F) {
    T13_re /= 3.0F;
    b_T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    b_T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    b_T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  b_a_re = f_a_re * f_a_re;
  b_a_im = f_a_re * 0.0F + 0.0F * f_a_re;
  if (b_a_im == 0.0F) {
    b_a_re /= 9.0F;
    b_a_im = 0.0F;
  } else if (b_a_re == 0.0F) {
    b_a_re = 0.0F;
    b_a_im /= 9.0F;
  } else {
    b_a_re /= 9.0F;
    b_a_im /= 9.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  if (b_T12_im == 0.0F) {
    b_T12_re /= 3.0F;
    b_T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 3.0F;
  } else {
    b_T12_re /= 3.0F;
    b_T12_im /= 3.0F;
  }

  c_T13_re = H13 * H13;
  v21_re = H13 * 0.0F + 0.0F * H13;
  if (v21_re == 0.0F) {
    c_T13_re /= 3.0F;
    v21_re = 0.0F;
  } else if (c_T13_re == 0.0F) {
    c_T13_re = 0.0F;
    v21_re /= 3.0F;
  } else {
    c_T13_re /= 3.0F;
    v21_re /= 3.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  if (b_T23_im == 0.0F) {
    b_T23_re /= 3.0F;
    b_T23_im = 0.0F;
  } else if (b_T23_re == 0.0F) {
    b_T23_re = 0.0F;
    b_T23_im /= 3.0F;
  } else {
    b_T23_re /= 3.0F;
    b_T23_im /= 3.0F;
  }

  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 3.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 3.0F;
  } else {
    d_T11_re /= 3.0F;
    c_T11_im /= 3.0F;
  }

  e_T11_re = H11 * H33;
  c_a_re = H11 * 0.0F + 0.0F * H33;
  if (c_a_re == 0.0F) {
    e_T11_re /= 3.0F;
    c_a_re = 0.0F;
  } else if (e_T11_re == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 3.0F;
  } else {
    e_T11_re /= 3.0F;
    c_a_re /= 3.0F;
  }

  c_T22_re = H22 * H33;
  b_T22_im = H22 * 0.0F + 0.0F * H33;
  if (b_T22_im == 0.0F) {
    c_T22_re /= 3.0F;
    b_T22_im = 0.0F;
  } else if (c_T22_re == 0.0F) {
    c_T22_re = 0.0F;
    b_T22_im /= 3.0F;
  } else {
    c_T22_re /= 3.0F;
    b_T22_im /= 3.0F;
  }

  tb_a.re = (((((b_a_re + b_T12_re) + c_T13_re) + b_T23_re) - d_T11_re) -
             e_T11_re) - c_T22_re;
  tb_a.im = (((((b_a_im + b_T12_im) + v21_re) + b_T23_im) - c_T11_im) - c_a_re)
    - b_T22_im;
  T11 = mpower(tb_a);
  ub_a.re = (e_a_re * e_a_re - c_a_im * c_a_im) - T11.re;
  ub_a.im = (e_a_re * c_a_im + c_a_im * e_a_re) - T11.im;
  T11 = b_mpower(ub_a);
  c_T13_re = H13 * H13;
  v21_re = H13 * 0.0F + 0.0F * H13;
  V1_norm = c_T13_re * H22 - v21_re * 0.0F;
  v21_re = c_T13_re * 0.0F + v21_re * H22;
  if (v21_re == 0.0F) {
    c_T13_re = V1_norm / 2.0F;
    v21_re = 0.0F;
  } else if (V1_norm == 0.0F) {
    c_T13_re = 0.0F;
    v21_re /= 2.0F;
  } else {
    c_T13_re = V1_norm / 2.0F;
    v21_re /= 2.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  c_T12_re = b_T12_re * H33 - b_T12_im * 0.0F;
  b_T12_im = b_T12_re * 0.0F + b_T12_im * H33;
  if (b_T12_im == 0.0F) {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im = 0.0F;
  } else if (c_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 2.0F;
  } else {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im /= 2.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  d_T11_re = H11 * b_T23_re - 0.0F * b_T23_im;
  c_T11_im = H11 * b_T23_im + 0.0F * b_T23_re;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re /= 2.0F;
    c_T11_im /= 2.0F;
  }

  wb_T11.re = (H11 + H22) + H33;
  wb_T11.im = 0.0F;
  fc1 = mpower(wb_T11);
  if (fc1.im == 0.0F) {
    b_re = fc1.re / 27.0F;
    b_im = 0.0F;
  } else if (fc1.re == 0.0F) {
    b_re = 0.0F;
    b_im = fc1.im / 27.0F;
  } else {
    b_re = fc1.re / 27.0F;
    b_im = fc1.im / 27.0F;
  }

  e_T11_re = (H11 + H22) + H33;
  c_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  c_T22_re = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  b_T22_im = e_T11_re * c_T12_re - 0.0F * c_T22_re;
  c_a_re = e_T11_re * c_T22_re + 0.0F * c_T12_re;
  if (c_a_re == 0.0F) {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re = 0.0F;
  } else if (b_T22_im == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 6.0F;
  } else {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re /= 6.0F;
  }

  c_T12_re = H12 * H13;
  c_T22_re = H12 * 0.0F + 0.0F * H13;
  b_T22_im = H11 * H22;
  d_T11_im = H11 * 0.0F + 0.0F * H22;
  f_T11_re = b_T22_im * H33 - d_T11_im * 0.0F;
  d_T11_im = b_T22_im * 0.0F + d_T11_im * H33;
  if (d_T11_im == 0.0F) {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im = 0.0F;
  } else if (f_T11_re == 0.0F) {
    b_T22_im = 0.0F;
    d_T11_im /= 2.0F;
  } else {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im /= 2.0F;
  }

  fc24.re = ((((((T11.re - c_T13_re) - b_T12_re) - d_T11_re) + b_re) + e_T11_re)
             + (c_T12_re * H23 - c_T22_re * 0.0F)) + b_T22_im;
  fc24.im = ((((((T11.im - v21_re) - b_T12_im) - c_T11_im) + b_im) + c_a_re) +
             (c_T12_re * 0.0F + c_T22_re * H23)) + d_T11_im;
  T11 = c_mpower(fc24);
  f_T11_re = (((((a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  b_T22_im = (((((a_im + T12_im) + b_T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  if (T11.im == 0.0F) {
    if (b_T22_im == 0.0F) {
      a_re = f_T11_re / T11.re;
      a_im = 0.0F;
    } else if (f_T11_re == 0.0F) {
      a_re = 0.0F;
      a_im = b_T22_im / T11.re;
    } else {
      a_re = f_T11_re / T11.re;
      a_im = b_T22_im / T11.re;
    }
  } else if (T11.re == 0.0F) {
    if (f_T11_re == 0.0F) {
      a_re = b_T22_im / T11.im;
      a_im = 0.0F;
    } else if (b_T22_im == 0.0F) {
      a_re = 0.0F;
      a_im = -(f_T11_re / T11.im);
    } else {
      a_re = b_T22_im / T11.im;
      a_im = -(f_T11_re / T11.im);
    }
  } else {
    c_T12_re = (float)fabs(T11.re);
    V1_norm = (float)fabs(T11.im);
    if (c_T12_re > V1_norm) {
      s = T11.im / T11.re;
      V1_norm = T11.re + s * T11.im;
      a_re = (f_T11_re + s * b_T22_im) / V1_norm;
      a_im = (b_T22_im - s * f_T11_re) / V1_norm;
    } else if (V1_norm == c_T12_re) {
      if (T11.re > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (T11.im > 0.0F) {
        V1_norm = 0.5F;
      } else {
        V1_norm = -0.5F;
      }

      a_re = (f_T11_re * s + b_T22_im * V1_norm) / c_T12_re;
      a_im = (b_T22_im * s - f_T11_re * V1_norm) / c_T12_re;
    } else {
      s = T11.re / T11.im;
      V1_norm = T11.im + s * T11.re;
      a_re = (s * f_T11_re + b_T22_im) / V1_norm;
      a_im = (s * b_T22_im - f_T11_re) / V1_norm;
    }
  }

  b_T11_re = (H11 * H13 + H12 * H23) + H13 * H33;
  T11_im = ((H11 * 0.0F + 0.0F * H13) + (H12 * 0.0F + 0.0F * H23)) + (H13 * 0.0F
    + 0.0F * H33);
  T11_re = (((T11_re + T22_re) + T33_re) + fc0.re) + a_re;
  b_T11_im = fc0.im + a_im;
  c_T11_re = b_T11_re * T11_re - T11_im * b_T11_im;
  T11_im = b_T11_re * b_T11_im + T11_im * T11_re;
  T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * H13;
  b_T11_im = H11 * 0.0F + 0.0F * H13;
  c_T13_re = H13 * H22;
  v21_re = H13 * 0.0F + 0.0F * H22;
  d_T11_im = (((H12 * T13_re - 0.0F * b_T13_im) - (H12 * T23_re - 0.0F * T23_im))
              - (T11_re * H23 - b_T11_im * 0.0F)) + (c_T13_re * H23 - v21_re *
    0.0F);
  c_T22_re = (((H12 * b_T13_im + 0.0F * T13_re) - (H12 * T23_im + 0.0F * T23_re))
              - (T11_re * 0.0F + b_T11_im * H23)) + (c_T13_re * 0.0F + v21_re *
    H23);
  if (c_T22_re == 0.0F) {
    if (T11_im == 0.0F) {
      T11_re = c_T11_re / d_T11_im;
      T11_im = 0.0F;
    } else if (c_T11_re == 0.0F) {
      T11_re = 0.0F;
      T11_im /= d_T11_im;
    } else {
      T11_re = c_T11_re / d_T11_im;
      T11_im /= d_T11_im;
    }
  } else if (d_T11_im == 0.0F) {
    if (c_T11_re == 0.0F) {
      T11_re = T11_im / c_T22_re;
      T11_im = 0.0F;
    } else if (T11_im == 0.0F) {
      T11_re = 0.0F;
      T11_im = -(c_T11_re / c_T22_re);
    } else {
      T11_re = T11_im / c_T22_re;
      T11_im = -(c_T11_re / c_T22_re);
    }
  } else {
    c_T12_re = (float)fabs(d_T11_im);
    if (c_T12_re > c_T22_re) {
      s = c_T22_re / d_T11_im;
      V1_norm = d_T11_im + s * c_T22_re;
      T11_re = (c_T11_re + s * T11_im) / V1_norm;
      T11_im = (T11_im - s * c_T11_re) / V1_norm;
    } else if (c_T22_re == c_T12_re) {
      if (d_T11_im > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (c_T11_re * s + T11_im * -0.5F > 0.0F) {
        T11_re = ((real32_T)rtInf);
      } else if (c_T11_re * s + T11_im * -0.5F < 0.0F) {
        T11_re = ((real32_T)rtMinusInf);
      } else {
        T11_re = ((real32_T)rtNaN);
      }

      if (T11_im * s - c_T11_re * -0.5F > 0.0F) {
        T11_im = ((real32_T)rtInf);
      } else if (T11_im * s - c_T11_re * -0.5F < 0.0F) {
        T11_im = ((real32_T)rtMinusInf);
      } else {
        T11_im = ((real32_T)rtNaN);
      }
    } else {
      s = d_T11_im / c_T22_re;
      V1_norm = c_T22_re + s * d_T11_im;
      T11_re = (s * c_T11_re + T11_im) / V1_norm;
      T11_im = (s * T11_im - c_T11_re) / V1_norm;
    }
  }

  V21_re = (b_T13_re - re) - T11_re;
  V21_im = (T13_im - im) - T11_im;
  xb_T11.re = (H11 + H22) + H33;
  xb_T11.im = 0.0F;
  fc0 = mpower(xb_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
           b_T12_im * 0.0F)) + c_T11_re;
  a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F +
           b_T12_im * H23)) + c_T11_im;
  b_a_re = (H11 + H22) + H33;
  c_a_re = (H11 + H22) + H33;
  yb_T11.re = (H11 + H22) + H33;
  yb_T11.im = 0.0F;
  fc0 = mpower(yb_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  d_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  b_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  e_a_re = (H11 + H22) + H33;
  ac_T11.re = (H11 + H22) + H33;
  ac_T11.im = 0.0F;
  fc0 = mpower(ac_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  f_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  c_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  g_a_re = (H11 + H22) + H33;
  h_a_re = (H11 + H22) + H33;
  bc_T11.re = (H11 + H22) + H33;
  bc_T11.im = 0.0F;
  fc0 = mpower(bc_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  i_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  d_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  j_a_re = (H11 + H22) + H33;
  T11_re = H11 / 3.0F;
  T22_re = H22 / 3.0F;
  T33_re = H33 / 3.0F;
  V1_norm = b_a_re * b_a_re;
  e_a_im = b_a_re * 0.0F + 0.0F * b_a_re;
  if (e_a_im == 0.0F) {
    b_a_re = V1_norm / 9.0F;
    e_a_im = 0.0F;
  } else if (V1_norm == 0.0F) {
    b_a_re = 0.0F;
    e_a_im /= 9.0F;
  } else {
    b_a_re = V1_norm / 9.0F;
    e_a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  if (T13_im == 0.0F) {
    T13_re /= 3.0F;
    T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  vb_a.re = (((((b_a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  vb_a.im = (((((e_a_im + T12_im) + T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  fc0 = mpower(vb_a);
  wb_a.re = (a_re * a_re - a_im * a_im) - fc0.re;
  wb_a.im = (a_re * a_im + a_im * a_re) - fc0.im;
  fc0 = b_mpower(wb_a);
  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  b_T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    b_T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    b_T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  cc_T11.re = (H11 + H22) + H33;
  cc_T11.im = 0.0F;
  T11 = mpower(cc_T11);
  if (T11.im == 0.0F) {
    re = T11.re / 27.0F;
    im = 0.0F;
  } else if (T11.re == 0.0F) {
    re = 0.0F;
    im = T11.im / 27.0F;
  } else {
    re = T11.re / 27.0F;
    im = T11.im / 27.0F;
  }

  c_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  d_T11_re = c_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = c_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  e_T11_re = d_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = d_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (e_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  fc25.re = ((((((fc0.re - T13_re) - T12_re) - b_T11_re) + re) + c_T11_re) +
             (b_T12_re * H23 - b_T12_im * 0.0F)) + d_T11_re;
  fc25.im = ((((((fc0.im - T13_im) - T12_im) - T11_im) + im) + b_T11_im) +
             (b_T12_re * 0.0F + b_T12_im * H23)) + c_T11_im;
  fc0 = c_mpower(fc25);
  a_re = c_a_re * c_a_re;
  a_im = c_a_re * 0.0F + 0.0F * c_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  if (T13_im == 0.0F) {
    T13_re /= 3.0F;
    T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  b_a_re = e_a_re * e_a_re;
  e_a_im = e_a_re * 0.0F + 0.0F * e_a_re;
  if (e_a_im == 0.0F) {
    b_a_re /= 9.0F;
    e_a_im = 0.0F;
  } else if (b_a_re == 0.0F) {
    b_a_re = 0.0F;
    e_a_im /= 9.0F;
  } else {
    b_a_re /= 9.0F;
    e_a_im /= 9.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  if (b_T12_im == 0.0F) {
    b_T12_re /= 3.0F;
    b_T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 3.0F;
  } else {
    b_T12_re /= 3.0F;
    b_T12_im /= 3.0F;
  }

  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  if (b_T13_im == 0.0F) {
    b_T13_re /= 3.0F;
    b_T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 3.0F;
  } else {
    b_T13_re /= 3.0F;
    b_T13_im /= 3.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  if (b_T23_im == 0.0F) {
    b_T23_re /= 3.0F;
    b_T23_im = 0.0F;
  } else if (b_T23_re == 0.0F) {
    b_T23_re = 0.0F;
    b_T23_im /= 3.0F;
  } else {
    b_T23_re /= 3.0F;
    b_T23_im /= 3.0F;
  }

  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 3.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 3.0F;
  } else {
    d_T11_re /= 3.0F;
    c_T11_im /= 3.0F;
  }

  e_T11_re = H11 * H33;
  c_a_re = H11 * 0.0F + 0.0F * H33;
  if (c_a_re == 0.0F) {
    e_T11_re /= 3.0F;
    c_a_re = 0.0F;
  } else if (e_T11_re == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 3.0F;
  } else {
    e_T11_re /= 3.0F;
    c_a_re /= 3.0F;
  }

  c_T22_re = H22 * H33;
  b_T22_im = H22 * 0.0F + 0.0F * H33;
  if (b_T22_im == 0.0F) {
    c_T22_re /= 3.0F;
    b_T22_im = 0.0F;
  } else if (c_T22_re == 0.0F) {
    c_T22_re = 0.0F;
    b_T22_im /= 3.0F;
  } else {
    c_T22_re /= 3.0F;
    b_T22_im /= 3.0F;
  }

  xb_a.re = (((((b_a_re + b_T12_re) + b_T13_re) + b_T23_re) - d_T11_re) -
             e_T11_re) - c_T22_re;
  xb_a.im = (((((e_a_im + b_T12_im) + b_T13_im) + b_T23_im) - c_T11_im) - c_a_re)
    - b_T22_im;
  T11 = mpower(xb_a);
  yb_a.re = (d_a_re * d_a_re - b_a_im * b_a_im) - T11.re;
  yb_a.im = (d_a_re * b_a_im + b_a_im * d_a_re) - T11.im;
  T11 = b_mpower(yb_a);
  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  c_T13_re = b_T13_re * H22 - b_T13_im * 0.0F;
  b_T13_im = b_T13_re * 0.0F + b_T13_im * H22;
  if (b_T13_im == 0.0F) {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im = 0.0F;
  } else if (c_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 2.0F;
  } else {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im /= 2.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  c_T12_re = b_T12_re * H33 - b_T12_im * 0.0F;
  b_T12_im = b_T12_re * 0.0F + b_T12_im * H33;
  if (b_T12_im == 0.0F) {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im = 0.0F;
  } else if (c_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 2.0F;
  } else {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im /= 2.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  d_T11_re = H11 * b_T23_re - 0.0F * b_T23_im;
  c_T11_im = H11 * b_T23_im + 0.0F * b_T23_re;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re /= 2.0F;
    c_T11_im /= 2.0F;
  }

  dc_T11.re = (H11 + H22) + H33;
  dc_T11.im = 0.0F;
  fc1 = mpower(dc_T11);
  if (fc1.im == 0.0F) {
    re = fc1.re / 27.0F;
    im = 0.0F;
  } else if (fc1.re == 0.0F) {
    re = 0.0F;
    im = fc1.im / 27.0F;
  } else {
    re = fc1.re / 27.0F;
    im = fc1.im / 27.0F;
  }

  e_T11_re = (H11 + H22) + H33;
  c_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  c_T22_re = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  b_T22_im = e_T11_re * c_T12_re - 0.0F * c_T22_re;
  c_a_re = e_T11_re * c_T22_re + 0.0F * c_T12_re;
  if (c_a_re == 0.0F) {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re = 0.0F;
  } else if (b_T22_im == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 6.0F;
  } else {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re /= 6.0F;
  }

  c_T12_re = H12 * H13;
  c_T22_re = H12 * 0.0F + 0.0F * H13;
  b_T22_im = H11 * H22;
  d_T11_im = H11 * 0.0F + 0.0F * H22;
  f_T11_re = b_T22_im * H33 - d_T11_im * 0.0F;
  d_T11_im = b_T22_im * 0.0F + d_T11_im * H33;
  if (d_T11_im == 0.0F) {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im = 0.0F;
  } else if (f_T11_re == 0.0F) {
    b_T22_im = 0.0F;
    d_T11_im /= 2.0F;
  } else {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im /= 2.0F;
  }

  fc26.re = ((((((T11.re - b_T13_re) - b_T12_re) - d_T11_re) + re) + e_T11_re) +
             (c_T12_re * H23 - c_T22_re * 0.0F)) + b_T22_im;
  fc26.im = ((((((T11.im - b_T13_im) - b_T12_im) - c_T11_im) + im) + c_a_re) +
             (c_T12_re * 0.0F + c_T22_re * H23)) + d_T11_im;
  T11 = c_mpower(fc26);
  f_T11_re = (((((a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  b_T22_im = (((((a_im + T12_im) + T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  if (T11.im == 0.0F) {
    if (b_T22_im == 0.0F) {
      a_re = f_T11_re / T11.re;
      a_im = 0.0F;
    } else if (f_T11_re == 0.0F) {
      a_re = 0.0F;
      a_im = b_T22_im / T11.re;
    } else {
      a_re = f_T11_re / T11.re;
      a_im = b_T22_im / T11.re;
    }
  } else if (T11.re == 0.0F) {
    if (f_T11_re == 0.0F) {
      a_re = b_T22_im / T11.im;
      a_im = 0.0F;
    } else if (b_T22_im == 0.0F) {
      a_re = 0.0F;
      a_im = -(f_T11_re / T11.im);
    } else {
      a_re = b_T22_im / T11.im;
      a_im = -(f_T11_re / T11.im);
    }
  } else {
    c_T12_re = (float)fabs(T11.re);
    V1_norm = (float)fabs(T11.im);
    if (c_T12_re > V1_norm) {
      s = T11.im / T11.re;
      V1_norm = T11.re + s * T11.im;
      a_re = (f_T11_re + s * b_T22_im) / V1_norm;
      a_im = (b_T22_im - s * f_T11_re) / V1_norm;
    } else if (V1_norm == c_T12_re) {
      if (T11.re > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (T11.im > 0.0F) {
        V1_norm = 0.5F;
      } else {
        V1_norm = -0.5F;
      }

      a_re = (f_T11_re * s + b_T22_im * V1_norm) / c_T12_re;
      a_im = (b_T22_im * s - f_T11_re * V1_norm) / c_T12_re;
    } else {
      s = T11.re / T11.im;
      V1_norm = T11.im + s * T11.re;
      a_re = (s * f_T11_re + b_T22_im) / V1_norm;
      a_im = (s * b_T22_im - f_T11_re) / V1_norm;
    }
  }

  re = 1.73205078F * (fc0.re - a_re);
  im = 1.73205078F * (fc0.im - a_im);
  b_re = re * 0.0F - im;
  im = re + im * 0.0F;
  if (im == 0.0F) {
    re = b_re / 2.0F;
    im = 0.0F;
  } else if (b_re == 0.0F) {
    re = 0.0F;
    im /= 2.0F;
  } else {
    re = b_re / 2.0F;
    im /= 2.0F;
  }

  a_re = g_a_re * g_a_re;
  a_im = g_a_re * 0.0F + 0.0F * g_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  if (T13_im == 0.0F) {
    T13_re /= 3.0F;
    T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  ac_a.re = (((((a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  ac_a.im = (((((a_im + T12_im) + T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  fc0 = mpower(ac_a);
  bc_a.re = (f_a_re * f_a_re - c_a_im * c_a_im) - fc0.re;
  bc_a.im = (f_a_re * c_a_im + c_a_im * f_a_re) - fc0.im;
  fc0 = b_mpower(bc_a);
  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  b_T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    b_T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    b_T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  ec_T11.re = (H11 + H22) + H33;
  ec_T11.im = 0.0F;
  T11 = mpower(ec_T11);
  if (T11.im == 0.0F) {
    b_re = T11.re / 27.0F;
    b_im = 0.0F;
  } else if (T11.re == 0.0F) {
    b_re = 0.0F;
    b_im = T11.im / 27.0F;
  } else {
    b_re = T11.re / 27.0F;
    b_im = T11.im / 27.0F;
  }

  c_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  d_T11_re = c_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = c_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  e_T11_re = d_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = d_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (e_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  fc27.re = ((((((fc0.re - T13_re) - T12_re) - b_T11_re) + b_re) + c_T11_re) +
             (b_T12_re * H23 - b_T12_im * 0.0F)) + d_T11_re;
  fc27.im = ((((((fc0.im - T13_im) - T12_im) - T11_im) + b_im) + b_T11_im) +
             (b_T12_re * 0.0F + b_T12_im * H23)) + c_T11_im;
  fc0 = c_mpower(fc27);
  if (fc0.im == 0.0F) {
    b_re = fc0.re / 2.0F;
    b_im = 0.0F;
  } else if (fc0.re == 0.0F) {
    b_re = 0.0F;
    b_im = fc0.im / 2.0F;
  } else {
    b_re = fc0.re / 2.0F;
    b_im = fc0.im / 2.0F;
  }

  a_re = h_a_re * h_a_re;
  a_im = h_a_re * 0.0F + 0.0F * h_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  if (T13_im == 0.0F) {
    T13_re /= 3.0F;
    T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  b_a_re = j_a_re * j_a_re;
  b_a_im = j_a_re * 0.0F + 0.0F * j_a_re;
  if (b_a_im == 0.0F) {
    b_a_re /= 9.0F;
    b_a_im = 0.0F;
  } else if (b_a_re == 0.0F) {
    b_a_re = 0.0F;
    b_a_im /= 9.0F;
  } else {
    b_a_re /= 9.0F;
    b_a_im /= 9.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  if (b_T12_im == 0.0F) {
    b_T12_re /= 3.0F;
    b_T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 3.0F;
  } else {
    b_T12_re /= 3.0F;
    b_T12_im /= 3.0F;
  }

  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  if (b_T13_im == 0.0F) {
    b_T13_re /= 3.0F;
    b_T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 3.0F;
  } else {
    b_T13_re /= 3.0F;
    b_T13_im /= 3.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  if (b_T23_im == 0.0F) {
    b_T23_re /= 3.0F;
    b_T23_im = 0.0F;
  } else if (b_T23_re == 0.0F) {
    b_T23_re = 0.0F;
    b_T23_im /= 3.0F;
  } else {
    b_T23_re /= 3.0F;
    b_T23_im /= 3.0F;
  }

  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 3.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 3.0F;
  } else {
    d_T11_re /= 3.0F;
    c_T11_im /= 3.0F;
  }

  e_T11_re = H11 * H33;
  c_a_re = H11 * 0.0F + 0.0F * H33;
  if (c_a_re == 0.0F) {
    e_T11_re /= 3.0F;
    c_a_re = 0.0F;
  } else if (e_T11_re == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 3.0F;
  } else {
    e_T11_re /= 3.0F;
    c_a_re /= 3.0F;
  }

  c_T22_re = H22 * H33;
  b_T22_im = H22 * 0.0F + 0.0F * H33;
  if (b_T22_im == 0.0F) {
    c_T22_re /= 3.0F;
    b_T22_im = 0.0F;
  } else if (c_T22_re == 0.0F) {
    c_T22_re = 0.0F;
    b_T22_im /= 3.0F;
  } else {
    c_T22_re /= 3.0F;
    b_T22_im /= 3.0F;
  }

  cc_a.re = (((((b_a_re + b_T12_re) + b_T13_re) + b_T23_re) - d_T11_re) -
             e_T11_re) - c_T22_re;
  cc_a.im = (((((b_a_im + b_T12_im) + b_T13_im) + b_T23_im) - c_T11_im) - c_a_re)
    - b_T22_im;
  fc0 = mpower(cc_a);
  dc_a.re = (i_a_re * i_a_re - d_a_im * d_a_im) - fc0.re;
  dc_a.im = (i_a_re * d_a_im + d_a_im * i_a_re) - fc0.im;
  fc0 = b_mpower(dc_a);
  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  c_T13_re = b_T13_re * H22 - b_T13_im * 0.0F;
  b_T13_im = b_T13_re * 0.0F + b_T13_im * H22;
  if (b_T13_im == 0.0F) {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im = 0.0F;
  } else if (c_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 2.0F;
  } else {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im /= 2.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  c_T12_re = b_T12_re * H33 - b_T12_im * 0.0F;
  b_T12_im = b_T12_re * 0.0F + b_T12_im * H33;
  if (b_T12_im == 0.0F) {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im = 0.0F;
  } else if (c_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 2.0F;
  } else {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im /= 2.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  d_T11_re = H11 * b_T23_re - 0.0F * b_T23_im;
  c_T11_im = H11 * b_T23_im + 0.0F * b_T23_re;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re /= 2.0F;
    c_T11_im /= 2.0F;
  }

  fc_T11.re = (H11 + H22) + H33;
  fc_T11.im = 0.0F;
  T11 = mpower(fc_T11);
  if (T11.im == 0.0F) {
    b_a_re = T11.re / 27.0F;
    d_a_re = 0.0F;
  } else if (T11.re == 0.0F) {
    b_a_re = 0.0F;
    d_a_re = T11.im / 27.0F;
  } else {
    b_a_re = T11.re / 27.0F;
    d_a_re = T11.im / 27.0F;
  }

  e_T11_re = (H11 + H22) + H33;
  c_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  c_T22_re = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  b_T22_im = e_T11_re * c_T12_re - 0.0F * c_T22_re;
  c_a_re = e_T11_re * c_T22_re + 0.0F * c_T12_re;
  if (c_a_re == 0.0F) {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re = 0.0F;
  } else if (b_T22_im == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 6.0F;
  } else {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re /= 6.0F;
  }

  c_T12_re = H12 * H13;
  c_T22_re = H12 * 0.0F + 0.0F * H13;
  b_T22_im = H11 * H22;
  d_T11_im = H11 * 0.0F + 0.0F * H22;
  f_T11_re = b_T22_im * H33 - d_T11_im * 0.0F;
  d_T11_im = b_T22_im * 0.0F + d_T11_im * H33;
  if (d_T11_im == 0.0F) {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im = 0.0F;
  } else if (f_T11_re == 0.0F) {
    b_T22_im = 0.0F;
    d_T11_im /= 2.0F;
  } else {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im /= 2.0F;
  }

  fc28.re = ((((((fc0.re - b_T13_re) - b_T12_re) - d_T11_re) + b_a_re) +
              e_T11_re) + (c_T12_re * H23 - c_T22_re * 0.0F)) + b_T22_im;
  fc28.im = ((((((fc0.im - b_T13_im) - b_T12_im) - c_T11_im) + d_a_re) + c_a_re)
             + (c_T12_re * 0.0F + c_T22_re * H23)) + d_T11_im;
  fc0 = c_mpower(fc28);
  f_T11_re = (((((a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  b_T22_im = (((((a_im + T12_im) + T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  d_T11_im = 2.0F * fc0.re;
  c_T22_re = 2.0F * fc0.im;
  if (c_T22_re == 0.0F) {
    if (b_T22_im == 0.0F) {
      a_re = f_T11_re / d_T11_im;
      a_im = 0.0F;
    } else if (f_T11_re == 0.0F) {
      a_re = 0.0F;
      a_im = b_T22_im / d_T11_im;
    } else {
      a_re = f_T11_re / d_T11_im;
      a_im = b_T22_im / d_T11_im;
    }
  } else if (d_T11_im == 0.0F) {
    if (f_T11_re == 0.0F) {
      a_re = b_T22_im / c_T22_re;
      a_im = 0.0F;
    } else if (b_T22_im == 0.0F) {
      a_re = 0.0F;
      a_im = -(f_T11_re / c_T22_re);
    } else {
      a_re = b_T22_im / c_T22_re;
      a_im = -(f_T11_re / c_T22_re);
    }
  } else {
    c_T12_re = (float)fabs(d_T11_im);
    V1_norm = (float)fabs(c_T22_re);
    if (c_T12_re > V1_norm) {
      s = c_T22_re / d_T11_im;
      V1_norm = d_T11_im + s * c_T22_re;
      a_re = (f_T11_re + s * b_T22_im) / V1_norm;
      a_im = (b_T22_im - s * f_T11_re) / V1_norm;
    } else if (V1_norm == c_T12_re) {
      if (d_T11_im > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (c_T22_re > 0.0F) {
        V1_norm = 0.5F;
      } else {
        V1_norm = -0.5F;
      }

      a_re = (f_T11_re * s + b_T22_im * V1_norm) / c_T12_re;
      a_im = (b_T22_im * s - f_T11_re * V1_norm) / c_T12_re;
    } else {
      s = d_T11_im / c_T22_re;
      V1_norm = c_T22_re + s * d_T11_im;
      a_re = (s * f_T11_re + b_T22_im) / V1_norm;
      a_im = (s * b_T22_im - f_T11_re) / V1_norm;
    }
  }

  a_re = ((((T11_re + T22_re) + T33_re) - re) - b_re) - a_re;
  a_im = ((0.0F - im) - b_im) - a_im;
  gc_T11.re = (H11 + H22) + H33;
  gc_T11.im = 0.0F;
  fc0 = mpower(gc_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  b_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  b_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  c_a_re = (H11 + H22) + H33;
  d_a_re = (H11 + H22) + H33;
  hc_T11.re = (H11 + H22) + H33;
  hc_T11.im = 0.0F;
  fc0 = mpower(hc_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  e_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  c_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  f_a_re = (H11 + H22) + H33;
  ic_T11.re = (H11 + H22) + H33;
  ic_T11.im = 0.0F;
  fc0 = mpower(ic_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  g_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  d_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  h_a_re = (H11 + H22) + H33;
  i_a_re = (H11 + H22) + H33;
  jc_T11.re = (H11 + H22) + H33;
  jc_T11.im = 0.0F;
  fc0 = mpower(jc_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  j_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  e_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  V23_re = (H11 + H22) + H33;
  V1_norm = a_re * a_re - a_im * a_im;
  a_im = a_re * a_im + a_im * a_re;
  T13_re = H13 * V1_norm - 0.0F * a_im;
  T13_im = H13 * a_im + 0.0F * V1_norm;
  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * H13;
  T11_im = H11 * 0.0F + 0.0F * H13;
  c_T13_re = H13 * H22;
  v21_re = H13 * 0.0F + 0.0F * H22;
  d_T11_im = (((H12 * b_T13_re - 0.0F * b_T13_im) - (H12 * T23_re - 0.0F *
    T23_im)) - (T11_re * H23 - T11_im * 0.0F)) + (c_T13_re * H23 - v21_re * 0.0F);
  c_T22_re = (((H12 * b_T13_im + 0.0F * b_T13_re) - (H12 * T23_im + 0.0F *
    T23_re)) - (T11_re * 0.0F + T11_im * H23)) + (c_T13_re * 0.0F + v21_re * H23);
  if (c_T22_re == 0.0F) {
    if (T13_im == 0.0F) {
      b_T13_re = T13_re / d_T11_im;
      T13_im = 0.0F;
    } else if (T13_re == 0.0F) {
      b_T13_re = 0.0F;
      T13_im /= d_T11_im;
    } else {
      b_T13_re = T13_re / d_T11_im;
      T13_im /= d_T11_im;
    }
  } else if (d_T11_im == 0.0F) {
    if (T13_re == 0.0F) {
      b_T13_re = T13_im / c_T22_re;
      T13_im = 0.0F;
    } else if (T13_im == 0.0F) {
      b_T13_re = 0.0F;
      T13_im = -(T13_re / c_T22_re);
    } else {
      b_T13_re = T13_im / c_T22_re;
      T13_im = -(T13_re / c_T22_re);
    }
  } else {
    c_T12_re = (float)fabs(d_T11_im);
    if (c_T12_re > c_T22_re) {
      s = c_T22_re / d_T11_im;
      V1_norm = d_T11_im + s * c_T22_re;
      b_T13_re = (T13_re + s * T13_im) / V1_norm;
      T13_im = (T13_im - s * T13_re) / V1_norm;
    } else if (c_T22_re == c_T12_re) {
      if (d_T11_im > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (T13_re * s + T13_im * -0.5F > 0.0F) {
        b_T13_re = ((real32_T)rtInf);
      } else if (T13_re * s + T13_im * -0.5F < 0.0F) {
        b_T13_re = ((real32_T)rtMinusInf);
      } else {
        b_T13_re = ((real32_T)rtNaN);
      }

      if (T13_im * s - T13_re * -0.5F > 0.0F) {
        T13_im = ((real32_T)rtInf);
      } else if (T13_im * s - T13_re * -0.5F < 0.0F) {
        T13_im = ((real32_T)rtMinusInf);
      } else {
        T13_im = ((real32_T)rtNaN);
      }
    } else {
      s = d_T11_im / c_T22_re;
      V1_norm = c_T22_re + s * d_T11_im;
      b_T13_re = (s * T13_re + T13_im) / V1_norm;
      T13_im = (s * T13_im - T13_re) / V1_norm;
    }
  }

  fc0 = mpower(T13);
  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * H33;
  T11_im = H11 * 0.0F + 0.0F * H33;
  T12_re = H12 * H33;
  T12_im = H12 * 0.0F + 0.0F * H33;
  T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  b_T11_re = H11 * H13;
  b_T11_im = H11 * 0.0F + 0.0F * H13;
  c_T13_re = H13 * H22;
  v21_re = H13 * 0.0F + 0.0F * H22;
  f_T11_re = ((fc0.re + (H13 * T23_re - 0.0F * T23_im)) - (T11_re * H13 - T11_im
    * 0.0F)) - (T12_re * H23 - T12_im * 0.0F);
  b_T22_im = ((fc0.im + (H13 * T23_im + 0.0F * T23_re)) - (T11_re * 0.0F +
    T11_im * H13)) - (T12_re * 0.0F + T12_im * H23);
  d_T11_im = (((H12 * T13_re - 0.0F * b_T13_im) - (H12 * b_T23_re - 0.0F *
    b_T23_im)) - (b_T11_re * H23 - b_T11_im * 0.0F)) + (c_T13_re * H23 - v21_re *
    0.0F);
  c_T22_re = (((H12 * b_T13_im + 0.0F * T13_re) - (H12 * b_T23_im + 0.0F *
    b_T23_re)) - (b_T11_re * 0.0F + b_T11_im * H23)) + (c_T13_re * 0.0F + v21_re
    * H23);
  if (c_T22_re == 0.0F) {
    if (b_T22_im == 0.0F) {
      re = f_T11_re / d_T11_im;
      im = 0.0F;
    } else if (f_T11_re == 0.0F) {
      re = 0.0F;
      im = b_T22_im / d_T11_im;
    } else {
      re = f_T11_re / d_T11_im;
      im = b_T22_im / d_T11_im;
    }
  } else if (d_T11_im == 0.0F) {
    if (f_T11_re == 0.0F) {
      re = b_T22_im / c_T22_re;
      im = 0.0F;
    } else if (b_T22_im == 0.0F) {
      re = 0.0F;
      im = -(f_T11_re / c_T22_re);
    } else {
      re = b_T22_im / c_T22_re;
      im = -(f_T11_re / c_T22_re);
    }
  } else {
    c_T12_re = (float)fabs(d_T11_im);
    if (c_T12_re > c_T22_re) {
      s = c_T22_re / d_T11_im;
      V1_norm = d_T11_im + s * c_T22_re;
      re = (f_T11_re + s * b_T22_im) / V1_norm;
      im = (b_T22_im - s * f_T11_re) / V1_norm;
    } else if (c_T22_re == c_T12_re) {
      if (d_T11_im > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (f_T11_re * s + b_T22_im * -0.5F > 0.0F) {
        re = ((real32_T)rtInf);
      } else if (f_T11_re * s + b_T22_im * -0.5F < 0.0F) {
        re = ((real32_T)rtMinusInf);
      } else {
        re = ((real32_T)rtNaN);
      }

      if (b_T22_im * s - f_T11_re * -0.5F > 0.0F) {
        im = ((real32_T)rtInf);
      } else if (b_T22_im * s - f_T11_re * -0.5F < 0.0F) {
        im = ((real32_T)rtMinusInf);
      } else {
        im = ((real32_T)rtNaN);
      }
    } else {
      s = d_T11_im / c_T22_re;
      V1_norm = c_T22_re + s * d_T11_im;
      re = (s * f_T11_re + b_T22_im) / V1_norm;
      im = (s * b_T22_im - f_T11_re) / V1_norm;
    }
  }

  T11_re = H11 / 3.0F;
  T22_re = H22 / 3.0F;
  T33_re = H33 / 3.0F;
  a_re = c_a_re * c_a_re;
  a_im = c_a_re * 0.0F + 0.0F * c_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  if (b_T13_im == 0.0F) {
    T13_re /= 3.0F;
    b_T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    b_T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    b_T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  ec_a.re = (((((a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  ec_a.im = (((((a_im + T12_im) + b_T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  fc0 = mpower(ec_a);
  fc_a.re = (b_a_re * b_a_re - b_a_im * b_a_im) - fc0.re;
  fc_a.im = (b_a_re * b_a_im + b_a_im * b_a_re) - fc0.im;
  fc0 = b_mpower(fc_a);
  T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  c_T13_re = T13_re * H22 - b_T13_im * 0.0F;
  b_T13_im = T13_re * 0.0F + b_T13_im * H22;
  if (b_T13_im == 0.0F) {
    T13_re = c_T13_re / 2.0F;
    b_T13_im = 0.0F;
  } else if (c_T13_re == 0.0F) {
    T13_re = 0.0F;
    b_T13_im /= 2.0F;
  } else {
    T13_re = c_T13_re / 2.0F;
    b_T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  b_T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    b_T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    b_T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  kc_T11.re = (H11 + H22) + H33;
  kc_T11.im = 0.0F;
  T11 = mpower(kc_T11);
  if (T11.im == 0.0F) {
    b_re = T11.re / 27.0F;
    b_im = 0.0F;
  } else if (T11.re == 0.0F) {
    b_re = 0.0F;
    b_im = T11.im / 27.0F;
  } else {
    b_re = T11.re / 27.0F;
    b_im = T11.im / 27.0F;
  }

  c_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  d_T11_re = c_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = c_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  e_T11_re = d_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = d_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (e_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  fc29.re = ((((((fc0.re - T13_re) - T12_re) - b_T11_re) + b_re) + c_T11_re) +
             (b_T12_re * H23 - b_T12_im * 0.0F)) + d_T11_re;
  fc29.im = ((((((fc0.im - b_T13_im) - T12_im) - T11_im) + b_im) + b_T11_im) +
             (b_T12_re * 0.0F + b_T12_im * H23)) + c_T11_im;
  fc0 = c_mpower(fc29);
  a_re = d_a_re * d_a_re;
  a_im = d_a_re * 0.0F + 0.0F * d_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  if (b_T13_im == 0.0F) {
    T13_re /= 3.0F;
    b_T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    b_T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    b_T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  b_a_re = f_a_re * f_a_re;
  b_a_im = f_a_re * 0.0F + 0.0F * f_a_re;
  if (b_a_im == 0.0F) {
    b_a_re /= 9.0F;
    b_a_im = 0.0F;
  } else if (b_a_re == 0.0F) {
    b_a_re = 0.0F;
    b_a_im /= 9.0F;
  } else {
    b_a_re /= 9.0F;
    b_a_im /= 9.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  if (b_T12_im == 0.0F) {
    b_T12_re /= 3.0F;
    b_T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 3.0F;
  } else {
    b_T12_re /= 3.0F;
    b_T12_im /= 3.0F;
  }

  c_T13_re = H13 * H13;
  v21_re = H13 * 0.0F + 0.0F * H13;
  if (v21_re == 0.0F) {
    c_T13_re /= 3.0F;
    v21_re = 0.0F;
  } else if (c_T13_re == 0.0F) {
    c_T13_re = 0.0F;
    v21_re /= 3.0F;
  } else {
    c_T13_re /= 3.0F;
    v21_re /= 3.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  if (b_T23_im == 0.0F) {
    b_T23_re /= 3.0F;
    b_T23_im = 0.0F;
  } else if (b_T23_re == 0.0F) {
    b_T23_re = 0.0F;
    b_T23_im /= 3.0F;
  } else {
    b_T23_re /= 3.0F;
    b_T23_im /= 3.0F;
  }

  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 3.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 3.0F;
  } else {
    d_T11_re /= 3.0F;
    c_T11_im /= 3.0F;
  }

  e_T11_re = H11 * H33;
  c_a_re = H11 * 0.0F + 0.0F * H33;
  if (c_a_re == 0.0F) {
    e_T11_re /= 3.0F;
    c_a_re = 0.0F;
  } else if (e_T11_re == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 3.0F;
  } else {
    e_T11_re /= 3.0F;
    c_a_re /= 3.0F;
  }

  c_T22_re = H22 * H33;
  b_T22_im = H22 * 0.0F + 0.0F * H33;
  if (b_T22_im == 0.0F) {
    c_T22_re /= 3.0F;
    b_T22_im = 0.0F;
  } else if (c_T22_re == 0.0F) {
    c_T22_re = 0.0F;
    b_T22_im /= 3.0F;
  } else {
    c_T22_re /= 3.0F;
    b_T22_im /= 3.0F;
  }

  gc_a.re = (((((b_a_re + b_T12_re) + c_T13_re) + b_T23_re) - d_T11_re) -
             e_T11_re) - c_T22_re;
  gc_a.im = (((((b_a_im + b_T12_im) + v21_re) + b_T23_im) - c_T11_im) - c_a_re)
    - b_T22_im;
  T11 = mpower(gc_a);
  hc_a.re = (e_a_re * e_a_re - c_a_im * c_a_im) - T11.re;
  hc_a.im = (e_a_re * c_a_im + c_a_im * e_a_re) - T11.im;
  T11 = b_mpower(hc_a);
  c_T13_re = H13 * H13;
  v21_re = H13 * 0.0F + 0.0F * H13;
  V1_norm = c_T13_re * H22 - v21_re * 0.0F;
  v21_re = c_T13_re * 0.0F + v21_re * H22;
  if (v21_re == 0.0F) {
    c_T13_re = V1_norm / 2.0F;
    v21_re = 0.0F;
  } else if (V1_norm == 0.0F) {
    c_T13_re = 0.0F;
    v21_re /= 2.0F;
  } else {
    c_T13_re = V1_norm / 2.0F;
    v21_re /= 2.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  c_T12_re = b_T12_re * H33 - b_T12_im * 0.0F;
  b_T12_im = b_T12_re * 0.0F + b_T12_im * H33;
  if (b_T12_im == 0.0F) {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im = 0.0F;
  } else if (c_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 2.0F;
  } else {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im /= 2.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  d_T11_re = H11 * b_T23_re - 0.0F * b_T23_im;
  c_T11_im = H11 * b_T23_im + 0.0F * b_T23_re;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re /= 2.0F;
    c_T11_im /= 2.0F;
  }

  lc_T11.re = (H11 + H22) + H33;
  lc_T11.im = 0.0F;
  fc1 = mpower(lc_T11);
  if (fc1.im == 0.0F) {
    b_re = fc1.re / 27.0F;
    b_im = 0.0F;
  } else if (fc1.re == 0.0F) {
    b_re = 0.0F;
    b_im = fc1.im / 27.0F;
  } else {
    b_re = fc1.re / 27.0F;
    b_im = fc1.im / 27.0F;
  }

  e_T11_re = (H11 + H22) + H33;
  c_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  c_T22_re = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  b_T22_im = e_T11_re * c_T12_re - 0.0F * c_T22_re;
  c_a_re = e_T11_re * c_T22_re + 0.0F * c_T12_re;
  if (c_a_re == 0.0F) {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re = 0.0F;
  } else if (b_T22_im == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 6.0F;
  } else {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re /= 6.0F;
  }

  c_T12_re = H12 * H13;
  c_T22_re = H12 * 0.0F + 0.0F * H13;
  b_T22_im = H11 * H22;
  d_T11_im = H11 * 0.0F + 0.0F * H22;
  f_T11_re = b_T22_im * H33 - d_T11_im * 0.0F;
  d_T11_im = b_T22_im * 0.0F + d_T11_im * H33;
  if (d_T11_im == 0.0F) {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im = 0.0F;
  } else if (f_T11_re == 0.0F) {
    b_T22_im = 0.0F;
    d_T11_im /= 2.0F;
  } else {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im /= 2.0F;
  }

  fc30.re = ((((((T11.re - c_T13_re) - b_T12_re) - d_T11_re) + b_re) + e_T11_re)
             + (c_T12_re * H23 - c_T22_re * 0.0F)) + b_T22_im;
  fc30.im = ((((((T11.im - v21_re) - b_T12_im) - c_T11_im) + b_im) + c_a_re) +
             (c_T12_re * 0.0F + c_T22_re * H23)) + d_T11_im;
  T11 = c_mpower(fc30);
  f_T11_re = (((((a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  b_T22_im = (((((a_im + T12_im) + b_T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  if (T11.im == 0.0F) {
    if (b_T22_im == 0.0F) {
      a_re = f_T11_re / T11.re;
      a_im = 0.0F;
    } else if (f_T11_re == 0.0F) {
      a_re = 0.0F;
      a_im = b_T22_im / T11.re;
    } else {
      a_re = f_T11_re / T11.re;
      a_im = b_T22_im / T11.re;
    }
  } else if (T11.re == 0.0F) {
    if (f_T11_re == 0.0F) {
      a_re = b_T22_im / T11.im;
      a_im = 0.0F;
    } else if (b_T22_im == 0.0F) {
      a_re = 0.0F;
      a_im = -(f_T11_re / T11.im);
    } else {
      a_re = b_T22_im / T11.im;
      a_im = -(f_T11_re / T11.im);
    }
  } else {
    c_T12_re = (float)fabs(T11.re);
    V1_norm = (float)fabs(T11.im);
    if (c_T12_re > V1_norm) {
      s = T11.im / T11.re;
      V1_norm = T11.re + s * T11.im;
      a_re = (f_T11_re + s * b_T22_im) / V1_norm;
      a_im = (b_T22_im - s * f_T11_re) / V1_norm;
    } else if (V1_norm == c_T12_re) {
      if (T11.re > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (T11.im > 0.0F) {
        V1_norm = 0.5F;
      } else {
        V1_norm = -0.5F;
      }

      a_re = (f_T11_re * s + b_T22_im * V1_norm) / c_T12_re;
      a_im = (b_T22_im * s - f_T11_re * V1_norm) / c_T12_re;
    } else {
      s = T11.re / T11.im;
      V1_norm = T11.im + s * T11.re;
      a_re = (s * f_T11_re + b_T22_im) / V1_norm;
      a_im = (s * b_T22_im - f_T11_re) / V1_norm;
    }
  }

  b_re = 1.73205078F * (fc0.re - a_re);
  b_im = 1.73205078F * (fc0.im - a_im);
  b_a_re = b_re * 0.0F - b_im;
  b_im = b_re + b_im * 0.0F;
  if (b_im == 0.0F) {
    b_re = b_a_re / 2.0F;
    b_im = 0.0F;
  } else if (b_a_re == 0.0F) {
    b_re = 0.0F;
    b_im /= 2.0F;
  } else {
    b_re = b_a_re / 2.0F;
    b_im /= 2.0F;
  }

  a_re = h_a_re * h_a_re;
  a_im = h_a_re * 0.0F + 0.0F * h_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  if (b_T13_im == 0.0F) {
    T13_re /= 3.0F;
    b_T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    b_T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    b_T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  ic_a.re = (((((a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  ic_a.im = (((((a_im + T12_im) + b_T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  fc0 = mpower(ic_a);
  jc_a.re = (g_a_re * g_a_re - d_a_im * d_a_im) - fc0.re;
  jc_a.im = (g_a_re * d_a_im + d_a_im * g_a_re) - fc0.im;
  fc0 = b_mpower(jc_a);
  T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  c_T13_re = T13_re * H22 - b_T13_im * 0.0F;
  b_T13_im = T13_re * 0.0F + b_T13_im * H22;
  if (b_T13_im == 0.0F) {
    T13_re = c_T13_re / 2.0F;
    b_T13_im = 0.0F;
  } else if (c_T13_re == 0.0F) {
    T13_re = 0.0F;
    b_T13_im /= 2.0F;
  } else {
    T13_re = c_T13_re / 2.0F;
    b_T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  b_T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    b_T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    b_T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  mc_T11.re = (H11 + H22) + H33;
  mc_T11.im = 0.0F;
  T11 = mpower(mc_T11);
  if (T11.im == 0.0F) {
    b_a_re = T11.re / 27.0F;
    d_a_re = 0.0F;
  } else if (T11.re == 0.0F) {
    b_a_re = 0.0F;
    d_a_re = T11.im / 27.0F;
  } else {
    b_a_re = T11.re / 27.0F;
    d_a_re = T11.im / 27.0F;
  }

  c_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  d_T11_re = c_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = c_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  e_T11_re = d_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = d_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (e_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  fc31.re = ((((((fc0.re - T13_re) - T12_re) - b_T11_re) + b_a_re) + c_T11_re) +
             (b_T12_re * H23 - b_T12_im * 0.0F)) + d_T11_re;
  fc31.im = ((((((fc0.im - b_T13_im) - T12_im) - T11_im) + d_a_re) + b_T11_im) +
             (b_T12_re * 0.0F + b_T12_im * H23)) + c_T11_im;
  fc0 = c_mpower(fc31);
  if (fc0.im == 0.0F) {
    b_a_re = fc0.re / 2.0F;
    d_a_re = 0.0F;
  } else if (fc0.re == 0.0F) {
    b_a_re = 0.0F;
    d_a_re = fc0.im / 2.0F;
  } else {
    b_a_re = fc0.re / 2.0F;
    d_a_re = fc0.im / 2.0F;
  }

  a_re = i_a_re * i_a_re;
  a_im = i_a_re * 0.0F + 0.0F * i_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  if (b_T13_im == 0.0F) {
    T13_re /= 3.0F;
    b_T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    b_T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    b_T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  V1_norm = V23_re * V23_re;
  c_T12_re = V23_re * 0.0F + 0.0F * V23_re;
  if (c_T12_re == 0.0F) {
    V23_re = V1_norm / 9.0F;
    c_T12_re = 0.0F;
  } else if (V1_norm == 0.0F) {
    V23_re = 0.0F;
    c_T12_re /= 9.0F;
  } else {
    V23_re = V1_norm / 9.0F;
    c_T12_re /= 9.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  if (b_T12_im == 0.0F) {
    b_T12_re /= 3.0F;
    b_T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 3.0F;
  } else {
    b_T12_re /= 3.0F;
    b_T12_im /= 3.0F;
  }

  c_T13_re = H13 * H13;
  v21_re = H13 * 0.0F + 0.0F * H13;
  if (v21_re == 0.0F) {
    c_T13_re /= 3.0F;
    v21_re = 0.0F;
  } else if (c_T13_re == 0.0F) {
    c_T13_re = 0.0F;
    v21_re /= 3.0F;
  } else {
    c_T13_re /= 3.0F;
    v21_re /= 3.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  if (b_T23_im == 0.0F) {
    b_T23_re /= 3.0F;
    b_T23_im = 0.0F;
  } else if (b_T23_re == 0.0F) {
    b_T23_re = 0.0F;
    b_T23_im /= 3.0F;
  } else {
    b_T23_re /= 3.0F;
    b_T23_im /= 3.0F;
  }

  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 3.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 3.0F;
  } else {
    d_T11_re /= 3.0F;
    c_T11_im /= 3.0F;
  }

  e_T11_re = H11 * H33;
  c_a_re = H11 * 0.0F + 0.0F * H33;
  if (c_a_re == 0.0F) {
    e_T11_re /= 3.0F;
    c_a_re = 0.0F;
  } else if (e_T11_re == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 3.0F;
  } else {
    e_T11_re /= 3.0F;
    c_a_re /= 3.0F;
  }

  c_T22_re = H22 * H33;
  b_T22_im = H22 * 0.0F + 0.0F * H33;
  if (b_T22_im == 0.0F) {
    c_T22_re /= 3.0F;
    b_T22_im = 0.0F;
  } else if (c_T22_re == 0.0F) {
    c_T22_re = 0.0F;
    b_T22_im /= 3.0F;
  } else {
    c_T22_re /= 3.0F;
    b_T22_im /= 3.0F;
  }

  c_V23.re = (((((V23_re + b_T12_re) + c_T13_re) + b_T23_re) - d_T11_re) -
              e_T11_re) - c_T22_re;
  c_V23.im = (((((c_T12_re + b_T12_im) + v21_re) + b_T23_im) - c_T11_im) -
              c_a_re) - b_T22_im;
  fc0 = mpower(c_V23);
  kc_a.re = (j_a_re * j_a_re - e_a_im * e_a_im) - fc0.re;
  kc_a.im = (j_a_re * e_a_im + e_a_im * j_a_re) - fc0.im;
  fc0 = b_mpower(kc_a);
  c_T13_re = H13 * H13;
  v21_re = H13 * 0.0F + 0.0F * H13;
  V1_norm = c_T13_re * H22 - v21_re * 0.0F;
  v21_re = c_T13_re * 0.0F + v21_re * H22;
  if (v21_re == 0.0F) {
    c_T13_re = V1_norm / 2.0F;
    v21_re = 0.0F;
  } else if (V1_norm == 0.0F) {
    c_T13_re = 0.0F;
    v21_re /= 2.0F;
  } else {
    c_T13_re = V1_norm / 2.0F;
    v21_re /= 2.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  c_T12_re = b_T12_re * H33 - b_T12_im * 0.0F;
  b_T12_im = b_T12_re * 0.0F + b_T12_im * H33;
  if (b_T12_im == 0.0F) {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im = 0.0F;
  } else if (c_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 2.0F;
  } else {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im /= 2.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  d_T11_re = H11 * b_T23_re - 0.0F * b_T23_im;
  c_T11_im = H11 * b_T23_im + 0.0F * b_T23_re;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re /= 2.0F;
    c_T11_im /= 2.0F;
  }

  nc_T11.re = (H11 + H22) + H33;
  nc_T11.im = 0.0F;
  T11 = mpower(nc_T11);
  if (T11.im == 0.0F) {
    V1_norm = T11.re / 27.0F;
    s = 0.0F;
  } else if (T11.re == 0.0F) {
    V1_norm = 0.0F;
    s = T11.im / 27.0F;
  } else {
    V1_norm = T11.re / 27.0F;
    s = T11.im / 27.0F;
  }

  e_T11_re = (H11 + H22) + H33;
  c_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  c_T22_re = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  b_T22_im = e_T11_re * c_T12_re - 0.0F * c_T22_re;
  c_a_re = e_T11_re * c_T22_re + 0.0F * c_T12_re;
  if (c_a_re == 0.0F) {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re = 0.0F;
  } else if (b_T22_im == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 6.0F;
  } else {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re /= 6.0F;
  }

  c_T12_re = H12 * H13;
  c_T22_re = H12 * 0.0F + 0.0F * H13;
  b_T22_im = H11 * H22;
  d_T11_im = H11 * 0.0F + 0.0F * H22;
  f_T11_re = b_T22_im * H33 - d_T11_im * 0.0F;
  d_T11_im = b_T22_im * 0.0F + d_T11_im * H33;
  if (d_T11_im == 0.0F) {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im = 0.0F;
  } else if (f_T11_re == 0.0F) {
    b_T22_im = 0.0F;
    d_T11_im /= 2.0F;
  } else {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im /= 2.0F;
  }

  fc32.re = ((((((fc0.re - c_T13_re) - b_T12_re) - d_T11_re) + V1_norm) +
              e_T11_re) + (c_T12_re * H23 - c_T22_re * 0.0F)) + b_T22_im;
  fc32.im = ((((((fc0.im - v21_re) - b_T12_im) - c_T11_im) + s) + c_a_re) +
             (c_T12_re * 0.0F + c_T22_re * H23)) + d_T11_im;
  fc0 = c_mpower(fc32);
  f_T11_re = (((((a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  b_T22_im = (((((a_im + T12_im) + b_T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  d_T11_im = 2.0F * fc0.re;
  c_T22_re = 2.0F * fc0.im;
  if (c_T22_re == 0.0F) {
    if (b_T22_im == 0.0F) {
      a_re = f_T11_re / d_T11_im;
      a_im = 0.0F;
    } else if (f_T11_re == 0.0F) {
      a_re = 0.0F;
      a_im = b_T22_im / d_T11_im;
    } else {
      a_re = f_T11_re / d_T11_im;
      a_im = b_T22_im / d_T11_im;
    }
  } else if (d_T11_im == 0.0F) {
    if (f_T11_re == 0.0F) {
      a_re = b_T22_im / c_T22_re;
      a_im = 0.0F;
    } else if (b_T22_im == 0.0F) {
      a_re = 0.0F;
      a_im = -(f_T11_re / c_T22_re);
    } else {
      a_re = b_T22_im / c_T22_re;
      a_im = -(f_T11_re / c_T22_re);
    }
  } else {
    c_T12_re = (float)fabs(d_T11_im);
    V1_norm = (float)fabs(c_T22_re);
    if (c_T12_re > V1_norm) {
      s = c_T22_re / d_T11_im;
      V1_norm = d_T11_im + s * c_T22_re;
      a_re = (f_T11_re + s * b_T22_im) / V1_norm;
      a_im = (b_T22_im - s * f_T11_re) / V1_norm;
    } else if (V1_norm == c_T12_re) {
      if (d_T11_im > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (c_T22_re > 0.0F) {
        V1_norm = 0.5F;
      } else {
        V1_norm = -0.5F;
      }

      a_re = (f_T11_re * s + b_T22_im * V1_norm) / c_T12_re;
      a_im = (b_T22_im * s - f_T11_re * V1_norm) / c_T12_re;
    } else {
      s = d_T11_im / c_T22_re;
      V1_norm = c_T22_re + s * d_T11_im;
      a_re = (s * f_T11_re + b_T22_im) / V1_norm;
      a_im = (s * b_T22_im - f_T11_re) / V1_norm;
    }
  }

  b_T11_re = (H11 * H13 + H12 * H23) + H13 * H33;
  T11_im = ((H11 * 0.0F + 0.0F * H13) + (H12 * 0.0F + 0.0F * H23)) + (H13 * 0.0F
    + 0.0F * H33);
  T11_re = ((((T11_re + T22_re) + T33_re) - b_re) - b_a_re) - a_re;
  b_T11_im = ((0.0F - b_im) - d_a_re) - a_im;
  c_T11_re = b_T11_re * T11_re - T11_im * b_T11_im;
  T11_im = b_T11_re * b_T11_im + T11_im * T11_re;
  T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * H13;
  b_T11_im = H11 * 0.0F + 0.0F * H13;
  c_T13_re = H13 * H22;
  v21_re = H13 * 0.0F + 0.0F * H22;
  d_T11_im = (((H12 * T13_re - 0.0F * b_T13_im) - (H12 * T23_re - 0.0F * T23_im))
              - (T11_re * H23 - b_T11_im * 0.0F)) + (c_T13_re * H23 - v21_re *
    0.0F);
  c_T22_re = (((H12 * b_T13_im + 0.0F * T13_re) - (H12 * T23_im + 0.0F * T23_re))
              - (T11_re * 0.0F + b_T11_im * H23)) + (c_T13_re * 0.0F + v21_re *
    H23);
  if (c_T22_re == 0.0F) {
    if (T11_im == 0.0F) {
      T11_re = c_T11_re / d_T11_im;
      T11_im = 0.0F;
    } else if (c_T11_re == 0.0F) {
      T11_re = 0.0F;
      T11_im /= d_T11_im;
    } else {
      T11_re = c_T11_re / d_T11_im;
      T11_im /= d_T11_im;
    }
  } else if (d_T11_im == 0.0F) {
    if (c_T11_re == 0.0F) {
      T11_re = T11_im / c_T22_re;
      T11_im = 0.0F;
    } else if (T11_im == 0.0F) {
      T11_re = 0.0F;
      T11_im = -(c_T11_re / c_T22_re);
    } else {
      T11_re = T11_im / c_T22_re;
      T11_im = -(c_T11_re / c_T22_re);
    }
  } else {
    c_T12_re = (float)fabs(d_T11_im);
    if (c_T12_re > c_T22_re) {
      s = c_T22_re / d_T11_im;
      V1_norm = d_T11_im + s * c_T22_re;
      T11_re = (c_T11_re + s * T11_im) / V1_norm;
      T11_im = (T11_im - s * c_T11_re) / V1_norm;
    } else if (c_T22_re == c_T12_re) {
      if (d_T11_im > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (c_T11_re * s + T11_im * -0.5F > 0.0F) {
        T11_re = ((real32_T)rtInf);
      } else if (c_T11_re * s + T11_im * -0.5F < 0.0F) {
        T11_re = ((real32_T)rtMinusInf);
      } else {
        T11_re = ((real32_T)rtNaN);
      }

      if (T11_im * s - c_T11_re * -0.5F > 0.0F) {
        T11_im = ((real32_T)rtInf);
      } else if (T11_im * s - c_T11_re * -0.5F < 0.0F) {
        T11_im = ((real32_T)rtMinusInf);
      } else {
        T11_im = ((real32_T)rtNaN);
      }
    } else {
      s = d_T11_im / c_T22_re;
      V1_norm = c_T22_re + s * d_T11_im;
      T11_re = (s * c_T11_re + T11_im) / V1_norm;
      T11_im = (s * T11_im - c_T11_re) / V1_norm;
    }
  }

  V22_re = (b_T13_re - re) - T11_re;
  V22_im = (T13_im - im) - T11_im;
  oc_T11.re = (H11 + H22) + H33;
  oc_T11.im = 0.0F;
  fc0 = mpower(oc_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
           b_T12_im * 0.0F)) + c_T11_re;
  a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F +
           b_T12_im * H23)) + c_T11_im;
  b_a_re = (H11 + H22) + H33;
  c_a_re = (H11 + H22) + H33;
  pc_T11.re = (H11 + H22) + H33;
  pc_T11.im = 0.0F;
  fc0 = mpower(pc_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  d_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  b_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  e_a_re = (H11 + H22) + H33;
  qc_T11.re = (H11 + H22) + H33;
  qc_T11.im = 0.0F;
  fc0 = mpower(qc_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  f_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  c_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  g_a_re = (H11 + H22) + H33;
  h_a_re = (H11 + H22) + H33;
  rc_T11.re = (H11 + H22) + H33;
  rc_T11.im = 0.0F;
  fc0 = mpower(rc_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  i_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  d_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  j_a_re = (H11 + H22) + H33;
  T11_re = H11 / 3.0F;
  T22_re = H22 / 3.0F;
  T33_re = H33 / 3.0F;
  V1_norm = b_a_re * b_a_re;
  e_a_im = b_a_re * 0.0F + 0.0F * b_a_re;
  if (e_a_im == 0.0F) {
    b_a_re = V1_norm / 9.0F;
    e_a_im = 0.0F;
  } else if (V1_norm == 0.0F) {
    b_a_re = 0.0F;
    e_a_im /= 9.0F;
  } else {
    b_a_re = V1_norm / 9.0F;
    e_a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  if (T13_im == 0.0F) {
    T13_re /= 3.0F;
    T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  lc_a.re = (((((b_a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  lc_a.im = (((((e_a_im + T12_im) + T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  fc0 = mpower(lc_a);
  mc_a.re = (a_re * a_re - a_im * a_im) - fc0.re;
  mc_a.im = (a_re * a_im + a_im * a_re) - fc0.im;
  fc0 = b_mpower(mc_a);
  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  b_T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    b_T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    b_T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  sc_T11.re = (H11 + H22) + H33;
  sc_T11.im = 0.0F;
  T11 = mpower(sc_T11);
  if (T11.im == 0.0F) {
    re = T11.re / 27.0F;
    im = 0.0F;
  } else if (T11.re == 0.0F) {
    re = 0.0F;
    im = T11.im / 27.0F;
  } else {
    re = T11.re / 27.0F;
    im = T11.im / 27.0F;
  }

  c_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  d_T11_re = c_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = c_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  e_T11_re = d_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = d_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (e_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  fc33.re = ((((((fc0.re - T13_re) - T12_re) - b_T11_re) + re) + c_T11_re) +
             (b_T12_re * H23 - b_T12_im * 0.0F)) + d_T11_re;
  fc33.im = ((((((fc0.im - T13_im) - T12_im) - T11_im) + im) + b_T11_im) +
             (b_T12_re * 0.0F + b_T12_im * H23)) + c_T11_im;
  fc0 = c_mpower(fc33);
  a_re = c_a_re * c_a_re;
  a_im = c_a_re * 0.0F + 0.0F * c_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  if (T13_im == 0.0F) {
    T13_re /= 3.0F;
    T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  b_a_re = e_a_re * e_a_re;
  e_a_im = e_a_re * 0.0F + 0.0F * e_a_re;
  if (e_a_im == 0.0F) {
    b_a_re /= 9.0F;
    e_a_im = 0.0F;
  } else if (b_a_re == 0.0F) {
    b_a_re = 0.0F;
    e_a_im /= 9.0F;
  } else {
    b_a_re /= 9.0F;
    e_a_im /= 9.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  if (b_T12_im == 0.0F) {
    b_T12_re /= 3.0F;
    b_T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 3.0F;
  } else {
    b_T12_re /= 3.0F;
    b_T12_im /= 3.0F;
  }

  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  if (b_T13_im == 0.0F) {
    b_T13_re /= 3.0F;
    b_T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 3.0F;
  } else {
    b_T13_re /= 3.0F;
    b_T13_im /= 3.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  if (b_T23_im == 0.0F) {
    b_T23_re /= 3.0F;
    b_T23_im = 0.0F;
  } else if (b_T23_re == 0.0F) {
    b_T23_re = 0.0F;
    b_T23_im /= 3.0F;
  } else {
    b_T23_re /= 3.0F;
    b_T23_im /= 3.0F;
  }

  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 3.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 3.0F;
  } else {
    d_T11_re /= 3.0F;
    c_T11_im /= 3.0F;
  }

  e_T11_re = H11 * H33;
  c_a_re = H11 * 0.0F + 0.0F * H33;
  if (c_a_re == 0.0F) {
    e_T11_re /= 3.0F;
    c_a_re = 0.0F;
  } else if (e_T11_re == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 3.0F;
  } else {
    e_T11_re /= 3.0F;
    c_a_re /= 3.0F;
  }

  c_T22_re = H22 * H33;
  b_T22_im = H22 * 0.0F + 0.0F * H33;
  if (b_T22_im == 0.0F) {
    c_T22_re /= 3.0F;
    b_T22_im = 0.0F;
  } else if (c_T22_re == 0.0F) {
    c_T22_re = 0.0F;
    b_T22_im /= 3.0F;
  } else {
    c_T22_re /= 3.0F;
    b_T22_im /= 3.0F;
  }

  nc_a.re = (((((b_a_re + b_T12_re) + b_T13_re) + b_T23_re) - d_T11_re) -
             e_T11_re) - c_T22_re;
  nc_a.im = (((((e_a_im + b_T12_im) + b_T13_im) + b_T23_im) - c_T11_im) - c_a_re)
    - b_T22_im;
  T11 = mpower(nc_a);
  oc_a.re = (d_a_re * d_a_re - b_a_im * b_a_im) - T11.re;
  oc_a.im = (d_a_re * b_a_im + b_a_im * d_a_re) - T11.im;
  T11 = b_mpower(oc_a);
  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  c_T13_re = b_T13_re * H22 - b_T13_im * 0.0F;
  b_T13_im = b_T13_re * 0.0F + b_T13_im * H22;
  if (b_T13_im == 0.0F) {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im = 0.0F;
  } else if (c_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 2.0F;
  } else {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im /= 2.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  c_T12_re = b_T12_re * H33 - b_T12_im * 0.0F;
  b_T12_im = b_T12_re * 0.0F + b_T12_im * H33;
  if (b_T12_im == 0.0F) {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im = 0.0F;
  } else if (c_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 2.0F;
  } else {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im /= 2.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  d_T11_re = H11 * b_T23_re - 0.0F * b_T23_im;
  c_T11_im = H11 * b_T23_im + 0.0F * b_T23_re;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re /= 2.0F;
    c_T11_im /= 2.0F;
  }

  tc_T11.re = (H11 + H22) + H33;
  tc_T11.im = 0.0F;
  fc1 = mpower(tc_T11);
  if (fc1.im == 0.0F) {
    re = fc1.re / 27.0F;
    im = 0.0F;
  } else if (fc1.re == 0.0F) {
    re = 0.0F;
    im = fc1.im / 27.0F;
  } else {
    re = fc1.re / 27.0F;
    im = fc1.im / 27.0F;
  }

  e_T11_re = (H11 + H22) + H33;
  c_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  c_T22_re = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  b_T22_im = e_T11_re * c_T12_re - 0.0F * c_T22_re;
  c_a_re = e_T11_re * c_T22_re + 0.0F * c_T12_re;
  if (c_a_re == 0.0F) {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re = 0.0F;
  } else if (b_T22_im == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 6.0F;
  } else {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re /= 6.0F;
  }

  c_T12_re = H12 * H13;
  c_T22_re = H12 * 0.0F + 0.0F * H13;
  b_T22_im = H11 * H22;
  d_T11_im = H11 * 0.0F + 0.0F * H22;
  f_T11_re = b_T22_im * H33 - d_T11_im * 0.0F;
  d_T11_im = b_T22_im * 0.0F + d_T11_im * H33;
  if (d_T11_im == 0.0F) {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im = 0.0F;
  } else if (f_T11_re == 0.0F) {
    b_T22_im = 0.0F;
    d_T11_im /= 2.0F;
  } else {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im /= 2.0F;
  }

  fc34.re = ((((((T11.re - b_T13_re) - b_T12_re) - d_T11_re) + re) + e_T11_re) +
             (c_T12_re * H23 - c_T22_re * 0.0F)) + b_T22_im;
  fc34.im = ((((((T11.im - b_T13_im) - b_T12_im) - c_T11_im) + im) + c_a_re) +
             (c_T12_re * 0.0F + c_T22_re * H23)) + d_T11_im;
  T11 = c_mpower(fc34);
  f_T11_re = (((((a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  b_T22_im = (((((a_im + T12_im) + T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  if (T11.im == 0.0F) {
    if (b_T22_im == 0.0F) {
      a_re = f_T11_re / T11.re;
      a_im = 0.0F;
    } else if (f_T11_re == 0.0F) {
      a_re = 0.0F;
      a_im = b_T22_im / T11.re;
    } else {
      a_re = f_T11_re / T11.re;
      a_im = b_T22_im / T11.re;
    }
  } else if (T11.re == 0.0F) {
    if (f_T11_re == 0.0F) {
      a_re = b_T22_im / T11.im;
      a_im = 0.0F;
    } else if (b_T22_im == 0.0F) {
      a_re = 0.0F;
      a_im = -(f_T11_re / T11.im);
    } else {
      a_re = b_T22_im / T11.im;
      a_im = -(f_T11_re / T11.im);
    }
  } else {
    c_T12_re = (float)fabs(T11.re);
    V1_norm = (float)fabs(T11.im);
    if (c_T12_re > V1_norm) {
      s = T11.im / T11.re;
      V1_norm = T11.re + s * T11.im;
      a_re = (f_T11_re + s * b_T22_im) / V1_norm;
      a_im = (b_T22_im - s * f_T11_re) / V1_norm;
    } else if (V1_norm == c_T12_re) {
      if (T11.re > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (T11.im > 0.0F) {
        V1_norm = 0.5F;
      } else {
        V1_norm = -0.5F;
      }

      a_re = (f_T11_re * s + b_T22_im * V1_norm) / c_T12_re;
      a_im = (b_T22_im * s - f_T11_re * V1_norm) / c_T12_re;
    } else {
      s = T11.re / T11.im;
      V1_norm = T11.im + s * T11.re;
      a_re = (s * f_T11_re + b_T22_im) / V1_norm;
      a_im = (s * b_T22_im - f_T11_re) / V1_norm;
    }
  }

  re = 1.73205078F * (fc0.re - a_re);
  im = 1.73205078F * (fc0.im - a_im);
  b_re = re * 0.0F - im;
  im = re + im * 0.0F;
  if (im == 0.0F) {
    re = b_re / 2.0F;
    im = 0.0F;
  } else if (b_re == 0.0F) {
    re = 0.0F;
    im /= 2.0F;
  } else {
    re = b_re / 2.0F;
    im /= 2.0F;
  }

  a_re = g_a_re * g_a_re;
  a_im = g_a_re * 0.0F + 0.0F * g_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  if (T13_im == 0.0F) {
    T13_re /= 3.0F;
    T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  pc_a.re = (((((a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  pc_a.im = (((((a_im + T12_im) + T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  fc0 = mpower(pc_a);
  qc_a.re = (f_a_re * f_a_re - c_a_im * c_a_im) - fc0.re;
  qc_a.im = (f_a_re * c_a_im + c_a_im * f_a_re) - fc0.im;
  fc0 = b_mpower(qc_a);
  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  b_T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    b_T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    b_T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  uc_T11.re = (H11 + H22) + H33;
  uc_T11.im = 0.0F;
  T11 = mpower(uc_T11);
  if (T11.im == 0.0F) {
    b_re = T11.re / 27.0F;
    b_im = 0.0F;
  } else if (T11.re == 0.0F) {
    b_re = 0.0F;
    b_im = T11.im / 27.0F;
  } else {
    b_re = T11.re / 27.0F;
    b_im = T11.im / 27.0F;
  }

  c_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  d_T11_re = c_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = c_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  e_T11_re = d_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = d_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (e_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  fc35.re = ((((((fc0.re - T13_re) - T12_re) - b_T11_re) + b_re) + c_T11_re) +
             (b_T12_re * H23 - b_T12_im * 0.0F)) + d_T11_re;
  fc35.im = ((((((fc0.im - T13_im) - T12_im) - T11_im) + b_im) + b_T11_im) +
             (b_T12_re * 0.0F + b_T12_im * H23)) + c_T11_im;
  fc0 = c_mpower(fc35);
  if (fc0.im == 0.0F) {
    b_re = fc0.re / 2.0F;
    b_im = 0.0F;
  } else if (fc0.re == 0.0F) {
    b_re = 0.0F;
    b_im = fc0.im / 2.0F;
  } else {
    b_re = fc0.re / 2.0F;
    b_im = fc0.im / 2.0F;
  }

  a_re = h_a_re * h_a_re;
  a_im = h_a_re * 0.0F + 0.0F * h_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  if (T13_im == 0.0F) {
    T13_re /= 3.0F;
    T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  b_a_re = j_a_re * j_a_re;
  b_a_im = j_a_re * 0.0F + 0.0F * j_a_re;
  if (b_a_im == 0.0F) {
    b_a_re /= 9.0F;
    b_a_im = 0.0F;
  } else if (b_a_re == 0.0F) {
    b_a_re = 0.0F;
    b_a_im /= 9.0F;
  } else {
    b_a_re /= 9.0F;
    b_a_im /= 9.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  if (b_T12_im == 0.0F) {
    b_T12_re /= 3.0F;
    b_T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 3.0F;
  } else {
    b_T12_re /= 3.0F;
    b_T12_im /= 3.0F;
  }

  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  if (b_T13_im == 0.0F) {
    b_T13_re /= 3.0F;
    b_T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 3.0F;
  } else {
    b_T13_re /= 3.0F;
    b_T13_im /= 3.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  if (b_T23_im == 0.0F) {
    b_T23_re /= 3.0F;
    b_T23_im = 0.0F;
  } else if (b_T23_re == 0.0F) {
    b_T23_re = 0.0F;
    b_T23_im /= 3.0F;
  } else {
    b_T23_re /= 3.0F;
    b_T23_im /= 3.0F;
  }

  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 3.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 3.0F;
  } else {
    d_T11_re /= 3.0F;
    c_T11_im /= 3.0F;
  }

  e_T11_re = H11 * H33;
  c_a_re = H11 * 0.0F + 0.0F * H33;
  if (c_a_re == 0.0F) {
    e_T11_re /= 3.0F;
    c_a_re = 0.0F;
  } else if (e_T11_re == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 3.0F;
  } else {
    e_T11_re /= 3.0F;
    c_a_re /= 3.0F;
  }

  c_T22_re = H22 * H33;
  b_T22_im = H22 * 0.0F + 0.0F * H33;
  if (b_T22_im == 0.0F) {
    c_T22_re /= 3.0F;
    b_T22_im = 0.0F;
  } else if (c_T22_re == 0.0F) {
    c_T22_re = 0.0F;
    b_T22_im /= 3.0F;
  } else {
    c_T22_re /= 3.0F;
    b_T22_im /= 3.0F;
  }

  rc_a.re = (((((b_a_re + b_T12_re) + b_T13_re) + b_T23_re) - d_T11_re) -
             e_T11_re) - c_T22_re;
  rc_a.im = (((((b_a_im + b_T12_im) + b_T13_im) + b_T23_im) - c_T11_im) - c_a_re)
    - b_T22_im;
  fc0 = mpower(rc_a);
  sc_a.re = (i_a_re * i_a_re - d_a_im * d_a_im) - fc0.re;
  sc_a.im = (i_a_re * d_a_im + d_a_im * i_a_re) - fc0.im;
  fc0 = b_mpower(sc_a);
  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  c_T13_re = b_T13_re * H22 - b_T13_im * 0.0F;
  b_T13_im = b_T13_re * 0.0F + b_T13_im * H22;
  if (b_T13_im == 0.0F) {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im = 0.0F;
  } else if (c_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 2.0F;
  } else {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im /= 2.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  c_T12_re = b_T12_re * H33 - b_T12_im * 0.0F;
  b_T12_im = b_T12_re * 0.0F + b_T12_im * H33;
  if (b_T12_im == 0.0F) {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im = 0.0F;
  } else if (c_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 2.0F;
  } else {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im /= 2.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  d_T11_re = H11 * b_T23_re - 0.0F * b_T23_im;
  c_T11_im = H11 * b_T23_im + 0.0F * b_T23_re;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re /= 2.0F;
    c_T11_im /= 2.0F;
  }

  vc_T11.re = (H11 + H22) + H33;
  vc_T11.im = 0.0F;
  T11 = mpower(vc_T11);
  if (T11.im == 0.0F) {
    b_a_re = T11.re / 27.0F;
    d_a_re = 0.0F;
  } else if (T11.re == 0.0F) {
    b_a_re = 0.0F;
    d_a_re = T11.im / 27.0F;
  } else {
    b_a_re = T11.re / 27.0F;
    d_a_re = T11.im / 27.0F;
  }

  e_T11_re = (H11 + H22) + H33;
  c_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  c_T22_re = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  b_T22_im = e_T11_re * c_T12_re - 0.0F * c_T22_re;
  c_a_re = e_T11_re * c_T22_re + 0.0F * c_T12_re;
  if (c_a_re == 0.0F) {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re = 0.0F;
  } else if (b_T22_im == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 6.0F;
  } else {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re /= 6.0F;
  }

  c_T12_re = H12 * H13;
  c_T22_re = H12 * 0.0F + 0.0F * H13;
  b_T22_im = H11 * H22;
  d_T11_im = H11 * 0.0F + 0.0F * H22;
  f_T11_re = b_T22_im * H33 - d_T11_im * 0.0F;
  d_T11_im = b_T22_im * 0.0F + d_T11_im * H33;
  if (d_T11_im == 0.0F) {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im = 0.0F;
  } else if (f_T11_re == 0.0F) {
    b_T22_im = 0.0F;
    d_T11_im /= 2.0F;
  } else {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im /= 2.0F;
  }

  fc36.re = ((((((fc0.re - b_T13_re) - b_T12_re) - d_T11_re) + b_a_re) +
              e_T11_re) + (c_T12_re * H23 - c_T22_re * 0.0F)) + b_T22_im;
  fc36.im = ((((((fc0.im - b_T13_im) - b_T12_im) - c_T11_im) + d_a_re) + c_a_re)
             + (c_T12_re * 0.0F + c_T22_re * H23)) + d_T11_im;
  fc0 = c_mpower(fc36);
  f_T11_re = (((((a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  b_T22_im = (((((a_im + T12_im) + T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  d_T11_im = 2.0F * fc0.re;
  c_T22_re = 2.0F * fc0.im;
  if (c_T22_re == 0.0F) {
    if (b_T22_im == 0.0F) {
      a_re = f_T11_re / d_T11_im;
      a_im = 0.0F;
    } else if (f_T11_re == 0.0F) {
      a_re = 0.0F;
      a_im = b_T22_im / d_T11_im;
    } else {
      a_re = f_T11_re / d_T11_im;
      a_im = b_T22_im / d_T11_im;
    }
  } else if (d_T11_im == 0.0F) {
    if (f_T11_re == 0.0F) {
      a_re = b_T22_im / c_T22_re;
      a_im = 0.0F;
    } else if (b_T22_im == 0.0F) {
      a_re = 0.0F;
      a_im = -(f_T11_re / c_T22_re);
    } else {
      a_re = b_T22_im / c_T22_re;
      a_im = -(f_T11_re / c_T22_re);
    }
  } else {
    c_T12_re = (float)fabs(d_T11_im);
    V1_norm = (float)fabs(c_T22_re);
    if (c_T12_re > V1_norm) {
      s = c_T22_re / d_T11_im;
      V1_norm = d_T11_im + s * c_T22_re;
      a_re = (f_T11_re + s * b_T22_im) / V1_norm;
      a_im = (b_T22_im - s * f_T11_re) / V1_norm;
    } else if (V1_norm == c_T12_re) {
      if (d_T11_im > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (c_T22_re > 0.0F) {
        V1_norm = 0.5F;
      } else {
        V1_norm = -0.5F;
      }

      a_re = (f_T11_re * s + b_T22_im * V1_norm) / c_T12_re;
      a_im = (b_T22_im * s - f_T11_re * V1_norm) / c_T12_re;
    } else {
      s = d_T11_im / c_T22_re;
      V1_norm = c_T22_re + s * d_T11_im;
      a_re = (s * f_T11_re + b_T22_im) / V1_norm;
      a_im = (s * b_T22_im - f_T11_re) / V1_norm;
    }
  }

  a_re = ((((T11_re + T22_re) + T33_re) + re) - b_re) - a_re;
  a_im = (im - b_im) - a_im;
  wc_T11.re = (H11 + H22) + H33;
  wc_T11.im = 0.0F;
  fc0 = mpower(wc_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  b_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  b_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  c_a_re = (H11 + H22) + H33;
  d_a_re = (H11 + H22) + H33;
  xc_T11.re = (H11 + H22) + H33;
  xc_T11.im = 0.0F;
  fc0 = mpower(xc_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  e_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  c_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  f_a_re = (H11 + H22) + H33;
  yc_T11.re = (H11 + H22) + H33;
  yc_T11.im = 0.0F;
  fc0 = mpower(yc_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  g_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  d_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  h_a_re = (H11 + H22) + H33;
  i_a_re = (H11 + H22) + H33;
  ad_T11.re = (H11 + H22) + H33;
  ad_T11.im = 0.0F;
  fc0 = mpower(ad_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  j_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  e_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  V23_re = (H11 + H22) + H33;
  V1_norm = a_re * a_re - a_im * a_im;
  a_im = a_re * a_im + a_im * a_re;
  T13_re = H13 * V1_norm - 0.0F * a_im;
  T13_im = H13 * a_im + 0.0F * V1_norm;
  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * H13;
  T11_im = H11 * 0.0F + 0.0F * H13;
  c_T13_re = H13 * H22;
  v21_re = H13 * 0.0F + 0.0F * H22;
  d_T11_im = (((H12 * b_T13_re - 0.0F * b_T13_im) - (H12 * T23_re - 0.0F *
    T23_im)) - (T11_re * H23 - T11_im * 0.0F)) + (c_T13_re * H23 - v21_re * 0.0F);
  c_T22_re = (((H12 * b_T13_im + 0.0F * b_T13_re) - (H12 * T23_im + 0.0F *
    T23_re)) - (T11_re * 0.0F + T11_im * H23)) + (c_T13_re * 0.0F + v21_re * H23);
  if (c_T22_re == 0.0F) {
    if (T13_im == 0.0F) {
      b_T13_re = T13_re / d_T11_im;
      T13_im = 0.0F;
    } else if (T13_re == 0.0F) {
      b_T13_re = 0.0F;
      T13_im /= d_T11_im;
    } else {
      b_T13_re = T13_re / d_T11_im;
      T13_im /= d_T11_im;
    }
  } else if (d_T11_im == 0.0F) {
    if (T13_re == 0.0F) {
      b_T13_re = T13_im / c_T22_re;
      T13_im = 0.0F;
    } else if (T13_im == 0.0F) {
      b_T13_re = 0.0F;
      T13_im = -(T13_re / c_T22_re);
    } else {
      b_T13_re = T13_im / c_T22_re;
      T13_im = -(T13_re / c_T22_re);
    }
  } else {
    c_T12_re = (float)fabs(d_T11_im);
    if (c_T12_re > c_T22_re) {
      s = c_T22_re / d_T11_im;
      V1_norm = d_T11_im + s * c_T22_re;
      b_T13_re = (T13_re + s * T13_im) / V1_norm;
      T13_im = (T13_im - s * T13_re) / V1_norm;
    } else if (c_T22_re == c_T12_re) {
      if (d_T11_im > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (T13_re * s + T13_im * -0.5F > 0.0F) {
        b_T13_re = ((real32_T)rtInf);
      } else if (T13_re * s + T13_im * -0.5F < 0.0F) {
        b_T13_re = ((real32_T)rtMinusInf);
      } else {
        b_T13_re = ((real32_T)rtNaN);
      }

      if (T13_im * s - T13_re * -0.5F > 0.0F) {
        T13_im = ((real32_T)rtInf);
      } else if (T13_im * s - T13_re * -0.5F < 0.0F) {
        T13_im = ((real32_T)rtMinusInf);
      } else {
        T13_im = ((real32_T)rtNaN);
      }
    } else {
      s = d_T11_im / c_T22_re;
      V1_norm = c_T22_re + s * d_T11_im;
      b_T13_re = (s * T13_re + T13_im) / V1_norm;
      T13_im = (s * T13_im - T13_re) / V1_norm;
    }
  }

  fc0 = mpower(T13);
  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * H33;
  T11_im = H11 * 0.0F + 0.0F * H33;
  T12_re = H12 * H33;
  T12_im = H12 * 0.0F + 0.0F * H33;
  T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  b_T11_re = H11 * H13;
  b_T11_im = H11 * 0.0F + 0.0F * H13;
  c_T13_re = H13 * H22;
  v21_re = H13 * 0.0F + 0.0F * H22;
  f_T11_re = ((fc0.re + (H13 * T23_re - 0.0F * T23_im)) - (T11_re * H13 - T11_im
    * 0.0F)) - (T12_re * H23 - T12_im * 0.0F);
  b_T22_im = ((fc0.im + (H13 * T23_im + 0.0F * T23_re)) - (T11_re * 0.0F +
    T11_im * H13)) - (T12_re * 0.0F + T12_im * H23);
  d_T11_im = (((H12 * T13_re - 0.0F * b_T13_im) - (H12 * b_T23_re - 0.0F *
    b_T23_im)) - (b_T11_re * H23 - b_T11_im * 0.0F)) + (c_T13_re * H23 - v21_re *
    0.0F);
  c_T22_re = (((H12 * b_T13_im + 0.0F * T13_re) - (H12 * b_T23_im + 0.0F *
    b_T23_re)) - (b_T11_re * 0.0F + b_T11_im * H23)) + (c_T13_re * 0.0F + v21_re
    * H23);
  if (c_T22_re == 0.0F) {
    if (b_T22_im == 0.0F) {
      re = f_T11_re / d_T11_im;
      im = 0.0F;
    } else if (f_T11_re == 0.0F) {
      re = 0.0F;
      im = b_T22_im / d_T11_im;
    } else {
      re = f_T11_re / d_T11_im;
      im = b_T22_im / d_T11_im;
    }
  } else if (d_T11_im == 0.0F) {
    if (f_T11_re == 0.0F) {
      re = b_T22_im / c_T22_re;
      im = 0.0F;
    } else if (b_T22_im == 0.0F) {
      re = 0.0F;
      im = -(f_T11_re / c_T22_re);
    } else {
      re = b_T22_im / c_T22_re;
      im = -(f_T11_re / c_T22_re);
    }
  } else {
    c_T12_re = (float)fabs(d_T11_im);
    if (c_T12_re > c_T22_re) {
      s = c_T22_re / d_T11_im;
      V1_norm = d_T11_im + s * c_T22_re;
      re = (f_T11_re + s * b_T22_im) / V1_norm;
      im = (b_T22_im - s * f_T11_re) / V1_norm;
    } else if (c_T22_re == c_T12_re) {
      if (d_T11_im > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (f_T11_re * s + b_T22_im * -0.5F > 0.0F) {
        re = ((real32_T)rtInf);
      } else if (f_T11_re * s + b_T22_im * -0.5F < 0.0F) {
        re = ((real32_T)rtMinusInf);
      } else {
        re = ((real32_T)rtNaN);
      }

      if (b_T22_im * s - f_T11_re * -0.5F > 0.0F) {
        im = ((real32_T)rtInf);
      } else if (b_T22_im * s - f_T11_re * -0.5F < 0.0F) {
        im = ((real32_T)rtMinusInf);
      } else {
        im = ((real32_T)rtNaN);
      }
    } else {
      s = d_T11_im / c_T22_re;
      V1_norm = c_T22_re + s * d_T11_im;
      re = (s * f_T11_re + b_T22_im) / V1_norm;
      im = (s * b_T22_im - f_T11_re) / V1_norm;
    }
  }

  T11_re = H11 / 3.0F;
  T22_re = H22 / 3.0F;
  T33_re = H33 / 3.0F;
  a_re = c_a_re * c_a_re;
  a_im = c_a_re * 0.0F + 0.0F * c_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  if (b_T13_im == 0.0F) {
    T13_re /= 3.0F;
    b_T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    b_T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    b_T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  tc_a.re = (((((a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  tc_a.im = (((((a_im + T12_im) + b_T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  fc0 = mpower(tc_a);
  uc_a.re = (b_a_re * b_a_re - b_a_im * b_a_im) - fc0.re;
  uc_a.im = (b_a_re * b_a_im + b_a_im * b_a_re) - fc0.im;
  fc0 = b_mpower(uc_a);
  T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  c_T13_re = T13_re * H22 - b_T13_im * 0.0F;
  b_T13_im = T13_re * 0.0F + b_T13_im * H22;
  if (b_T13_im == 0.0F) {
    T13_re = c_T13_re / 2.0F;
    b_T13_im = 0.0F;
  } else if (c_T13_re == 0.0F) {
    T13_re = 0.0F;
    b_T13_im /= 2.0F;
  } else {
    T13_re = c_T13_re / 2.0F;
    b_T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  b_T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    b_T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    b_T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  bd_T11.re = (H11 + H22) + H33;
  bd_T11.im = 0.0F;
  T11 = mpower(bd_T11);
  if (T11.im == 0.0F) {
    b_re = T11.re / 27.0F;
    b_im = 0.0F;
  } else if (T11.re == 0.0F) {
    b_re = 0.0F;
    b_im = T11.im / 27.0F;
  } else {
    b_re = T11.re / 27.0F;
    b_im = T11.im / 27.0F;
  }

  c_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  d_T11_re = c_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = c_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  e_T11_re = d_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = d_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (e_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  fc37.re = ((((((fc0.re - T13_re) - T12_re) - b_T11_re) + b_re) + c_T11_re) +
             (b_T12_re * H23 - b_T12_im * 0.0F)) + d_T11_re;
  fc37.im = ((((((fc0.im - b_T13_im) - T12_im) - T11_im) + b_im) + b_T11_im) +
             (b_T12_re * 0.0F + b_T12_im * H23)) + c_T11_im;
  fc0 = c_mpower(fc37);
  a_re = d_a_re * d_a_re;
  a_im = d_a_re * 0.0F + 0.0F * d_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  if (b_T13_im == 0.0F) {
    T13_re /= 3.0F;
    b_T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    b_T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    b_T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  b_a_re = f_a_re * f_a_re;
  b_a_im = f_a_re * 0.0F + 0.0F * f_a_re;
  if (b_a_im == 0.0F) {
    b_a_re /= 9.0F;
    b_a_im = 0.0F;
  } else if (b_a_re == 0.0F) {
    b_a_re = 0.0F;
    b_a_im /= 9.0F;
  } else {
    b_a_re /= 9.0F;
    b_a_im /= 9.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  if (b_T12_im == 0.0F) {
    b_T12_re /= 3.0F;
    b_T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 3.0F;
  } else {
    b_T12_re /= 3.0F;
    b_T12_im /= 3.0F;
  }

  c_T13_re = H13 * H13;
  v21_re = H13 * 0.0F + 0.0F * H13;
  if (v21_re == 0.0F) {
    c_T13_re /= 3.0F;
    v21_re = 0.0F;
  } else if (c_T13_re == 0.0F) {
    c_T13_re = 0.0F;
    v21_re /= 3.0F;
  } else {
    c_T13_re /= 3.0F;
    v21_re /= 3.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  if (b_T23_im == 0.0F) {
    b_T23_re /= 3.0F;
    b_T23_im = 0.0F;
  } else if (b_T23_re == 0.0F) {
    b_T23_re = 0.0F;
    b_T23_im /= 3.0F;
  } else {
    b_T23_re /= 3.0F;
    b_T23_im /= 3.0F;
  }

  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 3.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 3.0F;
  } else {
    d_T11_re /= 3.0F;
    c_T11_im /= 3.0F;
  }

  e_T11_re = H11 * H33;
  c_a_re = H11 * 0.0F + 0.0F * H33;
  if (c_a_re == 0.0F) {
    e_T11_re /= 3.0F;
    c_a_re = 0.0F;
  } else if (e_T11_re == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 3.0F;
  } else {
    e_T11_re /= 3.0F;
    c_a_re /= 3.0F;
  }

  c_T22_re = H22 * H33;
  b_T22_im = H22 * 0.0F + 0.0F * H33;
  if (b_T22_im == 0.0F) {
    c_T22_re /= 3.0F;
    b_T22_im = 0.0F;
  } else if (c_T22_re == 0.0F) {
    c_T22_re = 0.0F;
    b_T22_im /= 3.0F;
  } else {
    c_T22_re /= 3.0F;
    b_T22_im /= 3.0F;
  }

  vc_a.re = (((((b_a_re + b_T12_re) + c_T13_re) + b_T23_re) - d_T11_re) -
             e_T11_re) - c_T22_re;
  vc_a.im = (((((b_a_im + b_T12_im) + v21_re) + b_T23_im) - c_T11_im) - c_a_re)
    - b_T22_im;
  T11 = mpower(vc_a);
  wc_a.re = (e_a_re * e_a_re - c_a_im * c_a_im) - T11.re;
  wc_a.im = (e_a_re * c_a_im + c_a_im * e_a_re) - T11.im;
  T11 = b_mpower(wc_a);
  c_T13_re = H13 * H13;
  v21_re = H13 * 0.0F + 0.0F * H13;
  V1_norm = c_T13_re * H22 - v21_re * 0.0F;
  v21_re = c_T13_re * 0.0F + v21_re * H22;
  if (v21_re == 0.0F) {
    c_T13_re = V1_norm / 2.0F;
    v21_re = 0.0F;
  } else if (V1_norm == 0.0F) {
    c_T13_re = 0.0F;
    v21_re /= 2.0F;
  } else {
    c_T13_re = V1_norm / 2.0F;
    v21_re /= 2.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  c_T12_re = b_T12_re * H33 - b_T12_im * 0.0F;
  b_T12_im = b_T12_re * 0.0F + b_T12_im * H33;
  if (b_T12_im == 0.0F) {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im = 0.0F;
  } else if (c_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 2.0F;
  } else {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im /= 2.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  d_T11_re = H11 * b_T23_re - 0.0F * b_T23_im;
  c_T11_im = H11 * b_T23_im + 0.0F * b_T23_re;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re /= 2.0F;
    c_T11_im /= 2.0F;
  }

  cd_T11.re = (H11 + H22) + H33;
  cd_T11.im = 0.0F;
  fc1 = mpower(cd_T11);
  if (fc1.im == 0.0F) {
    b_re = fc1.re / 27.0F;
    b_im = 0.0F;
  } else if (fc1.re == 0.0F) {
    b_re = 0.0F;
    b_im = fc1.im / 27.0F;
  } else {
    b_re = fc1.re / 27.0F;
    b_im = fc1.im / 27.0F;
  }

  e_T11_re = (H11 + H22) + H33;
  c_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  c_T22_re = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  b_T22_im = e_T11_re * c_T12_re - 0.0F * c_T22_re;
  c_a_re = e_T11_re * c_T22_re + 0.0F * c_T12_re;
  if (c_a_re == 0.0F) {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re = 0.0F;
  } else if (b_T22_im == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 6.0F;
  } else {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re /= 6.0F;
  }

  c_T12_re = H12 * H13;
  c_T22_re = H12 * 0.0F + 0.0F * H13;
  b_T22_im = H11 * H22;
  d_T11_im = H11 * 0.0F + 0.0F * H22;
  f_T11_re = b_T22_im * H33 - d_T11_im * 0.0F;
  d_T11_im = b_T22_im * 0.0F + d_T11_im * H33;
  if (d_T11_im == 0.0F) {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im = 0.0F;
  } else if (f_T11_re == 0.0F) {
    b_T22_im = 0.0F;
    d_T11_im /= 2.0F;
  } else {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im /= 2.0F;
  }

  fc38.re = ((((((T11.re - c_T13_re) - b_T12_re) - d_T11_re) + b_re) + e_T11_re)
             + (c_T12_re * H23 - c_T22_re * 0.0F)) + b_T22_im;
  fc38.im = ((((((T11.im - v21_re) - b_T12_im) - c_T11_im) + b_im) + c_a_re) +
             (c_T12_re * 0.0F + c_T22_re * H23)) + d_T11_im;
  T11 = c_mpower(fc38);
  f_T11_re = (((((a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  b_T22_im = (((((a_im + T12_im) + b_T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  if (T11.im == 0.0F) {
    if (b_T22_im == 0.0F) {
      a_re = f_T11_re / T11.re;
      a_im = 0.0F;
    } else if (f_T11_re == 0.0F) {
      a_re = 0.0F;
      a_im = b_T22_im / T11.re;
    } else {
      a_re = f_T11_re / T11.re;
      a_im = b_T22_im / T11.re;
    }
  } else if (T11.re == 0.0F) {
    if (f_T11_re == 0.0F) {
      a_re = b_T22_im / T11.im;
      a_im = 0.0F;
    } else if (b_T22_im == 0.0F) {
      a_re = 0.0F;
      a_im = -(f_T11_re / T11.im);
    } else {
      a_re = b_T22_im / T11.im;
      a_im = -(f_T11_re / T11.im);
    }
  } else {
    c_T12_re = (float)fabs(T11.re);
    V1_norm = (float)fabs(T11.im);
    if (c_T12_re > V1_norm) {
      s = T11.im / T11.re;
      V1_norm = T11.re + s * T11.im;
      a_re = (f_T11_re + s * b_T22_im) / V1_norm;
      a_im = (b_T22_im - s * f_T11_re) / V1_norm;
    } else if (V1_norm == c_T12_re) {
      if (T11.re > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (T11.im > 0.0F) {
        V1_norm = 0.5F;
      } else {
        V1_norm = -0.5F;
      }

      a_re = (f_T11_re * s + b_T22_im * V1_norm) / c_T12_re;
      a_im = (b_T22_im * s - f_T11_re * V1_norm) / c_T12_re;
    } else {
      s = T11.re / T11.im;
      V1_norm = T11.im + s * T11.re;
      a_re = (s * f_T11_re + b_T22_im) / V1_norm;
      a_im = (s * b_T22_im - f_T11_re) / V1_norm;
    }
  }

  b_re = 1.73205078F * (fc0.re - a_re);
  b_im = 1.73205078F * (fc0.im - a_im);
  b_a_re = b_re * 0.0F - b_im;
  b_im = b_re + b_im * 0.0F;
  if (b_im == 0.0F) {
    b_re = b_a_re / 2.0F;
    b_im = 0.0F;
  } else if (b_a_re == 0.0F) {
    b_re = 0.0F;
    b_im /= 2.0F;
  } else {
    b_re = b_a_re / 2.0F;
    b_im /= 2.0F;
  }

  a_re = h_a_re * h_a_re;
  a_im = h_a_re * 0.0F + 0.0F * h_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  if (b_T13_im == 0.0F) {
    T13_re /= 3.0F;
    b_T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    b_T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    b_T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  xc_a.re = (((((a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  xc_a.im = (((((a_im + T12_im) + b_T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  fc0 = mpower(xc_a);
  yc_a.re = (g_a_re * g_a_re - d_a_im * d_a_im) - fc0.re;
  yc_a.im = (g_a_re * d_a_im + d_a_im * g_a_re) - fc0.im;
  fc0 = b_mpower(yc_a);
  T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  c_T13_re = T13_re * H22 - b_T13_im * 0.0F;
  b_T13_im = T13_re * 0.0F + b_T13_im * H22;
  if (b_T13_im == 0.0F) {
    T13_re = c_T13_re / 2.0F;
    b_T13_im = 0.0F;
  } else if (c_T13_re == 0.0F) {
    T13_re = 0.0F;
    b_T13_im /= 2.0F;
  } else {
    T13_re = c_T13_re / 2.0F;
    b_T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  b_T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    b_T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    b_T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  dd_T11.re = (H11 + H22) + H33;
  dd_T11.im = 0.0F;
  T11 = mpower(dd_T11);
  if (T11.im == 0.0F) {
    b_a_re = T11.re / 27.0F;
    d_a_re = 0.0F;
  } else if (T11.re == 0.0F) {
    b_a_re = 0.0F;
    d_a_re = T11.im / 27.0F;
  } else {
    b_a_re = T11.re / 27.0F;
    d_a_re = T11.im / 27.0F;
  }

  c_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  d_T11_re = c_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = c_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  e_T11_re = d_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = d_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (e_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  fc39.re = ((((((fc0.re - T13_re) - T12_re) - b_T11_re) + b_a_re) + c_T11_re) +
             (b_T12_re * H23 - b_T12_im * 0.0F)) + d_T11_re;
  fc39.im = ((((((fc0.im - b_T13_im) - T12_im) - T11_im) + d_a_re) + b_T11_im) +
             (b_T12_re * 0.0F + b_T12_im * H23)) + c_T11_im;
  fc0 = c_mpower(fc39);
  if (fc0.im == 0.0F) {
    b_a_re = fc0.re / 2.0F;
    d_a_re = 0.0F;
  } else if (fc0.re == 0.0F) {
    b_a_re = 0.0F;
    d_a_re = fc0.im / 2.0F;
  } else {
    b_a_re = fc0.re / 2.0F;
    d_a_re = fc0.im / 2.0F;
  }

  a_re = i_a_re * i_a_re;
  a_im = i_a_re * 0.0F + 0.0F * i_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  if (b_T13_im == 0.0F) {
    T13_re /= 3.0F;
    b_T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    b_T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    b_T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  V1_norm = V23_re * V23_re;
  c_T12_re = V23_re * 0.0F + 0.0F * V23_re;
  if (c_T12_re == 0.0F) {
    V23_re = V1_norm / 9.0F;
    c_T12_re = 0.0F;
  } else if (V1_norm == 0.0F) {
    V23_re = 0.0F;
    c_T12_re /= 9.0F;
  } else {
    V23_re = V1_norm / 9.0F;
    c_T12_re /= 9.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  if (b_T12_im == 0.0F) {
    b_T12_re /= 3.0F;
    b_T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 3.0F;
  } else {
    b_T12_re /= 3.0F;
    b_T12_im /= 3.0F;
  }

  c_T13_re = H13 * H13;
  v21_re = H13 * 0.0F + 0.0F * H13;
  if (v21_re == 0.0F) {
    c_T13_re /= 3.0F;
    v21_re = 0.0F;
  } else if (c_T13_re == 0.0F) {
    c_T13_re = 0.0F;
    v21_re /= 3.0F;
  } else {
    c_T13_re /= 3.0F;
    v21_re /= 3.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  if (b_T23_im == 0.0F) {
    b_T23_re /= 3.0F;
    b_T23_im = 0.0F;
  } else if (b_T23_re == 0.0F) {
    b_T23_re = 0.0F;
    b_T23_im /= 3.0F;
  } else {
    b_T23_re /= 3.0F;
    b_T23_im /= 3.0F;
  }

  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 3.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 3.0F;
  } else {
    d_T11_re /= 3.0F;
    c_T11_im /= 3.0F;
  }

  e_T11_re = H11 * H33;
  c_a_re = H11 * 0.0F + 0.0F * H33;
  if (c_a_re == 0.0F) {
    e_T11_re /= 3.0F;
    c_a_re = 0.0F;
  } else if (e_T11_re == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 3.0F;
  } else {
    e_T11_re /= 3.0F;
    c_a_re /= 3.0F;
  }

  c_T22_re = H22 * H33;
  b_T22_im = H22 * 0.0F + 0.0F * H33;
  if (b_T22_im == 0.0F) {
    c_T22_re /= 3.0F;
    b_T22_im = 0.0F;
  } else if (c_T22_re == 0.0F) {
    c_T22_re = 0.0F;
    b_T22_im /= 3.0F;
  } else {
    c_T22_re /= 3.0F;
    b_T22_im /= 3.0F;
  }

  d_V23.re = (((((V23_re + b_T12_re) + c_T13_re) + b_T23_re) - d_T11_re) -
              e_T11_re) - c_T22_re;
  d_V23.im = (((((c_T12_re + b_T12_im) + v21_re) + b_T23_im) - c_T11_im) -
              c_a_re) - b_T22_im;
  fc0 = mpower(d_V23);
  ad_a.re = (j_a_re * j_a_re - e_a_im * e_a_im) - fc0.re;
  ad_a.im = (j_a_re * e_a_im + e_a_im * j_a_re) - fc0.im;
  fc0 = b_mpower(ad_a);
  c_T13_re = H13 * H13;
  v21_re = H13 * 0.0F + 0.0F * H13;
  V1_norm = c_T13_re * H22 - v21_re * 0.0F;
  v21_re = c_T13_re * 0.0F + v21_re * H22;
  if (v21_re == 0.0F) {
    c_T13_re = V1_norm / 2.0F;
    v21_re = 0.0F;
  } else if (V1_norm == 0.0F) {
    c_T13_re = 0.0F;
    v21_re /= 2.0F;
  } else {
    c_T13_re = V1_norm / 2.0F;
    v21_re /= 2.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  c_T12_re = b_T12_re * H33 - b_T12_im * 0.0F;
  b_T12_im = b_T12_re * 0.0F + b_T12_im * H33;
  if (b_T12_im == 0.0F) {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im = 0.0F;
  } else if (c_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 2.0F;
  } else {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im /= 2.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  d_T11_re = H11 * b_T23_re - 0.0F * b_T23_im;
  c_T11_im = H11 * b_T23_im + 0.0F * b_T23_re;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re /= 2.0F;
    c_T11_im /= 2.0F;
  }

  ed_T11.re = (H11 + H22) + H33;
  ed_T11.im = 0.0F;
  T11 = mpower(ed_T11);
  if (T11.im == 0.0F) {
    V1_norm = T11.re / 27.0F;
    s = 0.0F;
  } else if (T11.re == 0.0F) {
    V1_norm = 0.0F;
    s = T11.im / 27.0F;
  } else {
    V1_norm = T11.re / 27.0F;
    s = T11.im / 27.0F;
  }

  e_T11_re = (H11 + H22) + H33;
  c_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  c_T22_re = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  b_T22_im = e_T11_re * c_T12_re - 0.0F * c_T22_re;
  c_a_re = e_T11_re * c_T22_re + 0.0F * c_T12_re;
  if (c_a_re == 0.0F) {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re = 0.0F;
  } else if (b_T22_im == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 6.0F;
  } else {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re /= 6.0F;
  }

  c_T12_re = H12 * H13;
  c_T22_re = H12 * 0.0F + 0.0F * H13;
  b_T22_im = H11 * H22;
  d_T11_im = H11 * 0.0F + 0.0F * H22;
  f_T11_re = b_T22_im * H33 - d_T11_im * 0.0F;
  d_T11_im = b_T22_im * 0.0F + d_T11_im * H33;
  if (d_T11_im == 0.0F) {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im = 0.0F;
  } else if (f_T11_re == 0.0F) {
    b_T22_im = 0.0F;
    d_T11_im /= 2.0F;
  } else {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im /= 2.0F;
  }

  fc40.re = ((((((fc0.re - c_T13_re) - b_T12_re) - d_T11_re) + V1_norm) +
              e_T11_re) + (c_T12_re * H23 - c_T22_re * 0.0F)) + b_T22_im;
  fc40.im = ((((((fc0.im - v21_re) - b_T12_im) - c_T11_im) + s) + c_a_re) +
             (c_T12_re * 0.0F + c_T22_re * H23)) + d_T11_im;
  fc0 = c_mpower(fc40);
  f_T11_re = (((((a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  b_T22_im = (((((a_im + T12_im) + b_T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  d_T11_im = 2.0F * fc0.re;
  c_T22_re = 2.0F * fc0.im;
  if (c_T22_re == 0.0F) {
    if (b_T22_im == 0.0F) {
      a_re = f_T11_re / d_T11_im;
      a_im = 0.0F;
    } else if (f_T11_re == 0.0F) {
      a_re = 0.0F;
      a_im = b_T22_im / d_T11_im;
    } else {
      a_re = f_T11_re / d_T11_im;
      a_im = b_T22_im / d_T11_im;
    }
  } else if (d_T11_im == 0.0F) {
    if (f_T11_re == 0.0F) {
      a_re = b_T22_im / c_T22_re;
      a_im = 0.0F;
    } else if (b_T22_im == 0.0F) {
      a_re = 0.0F;
      a_im = -(f_T11_re / c_T22_re);
    } else {
      a_re = b_T22_im / c_T22_re;
      a_im = -(f_T11_re / c_T22_re);
    }
  } else {
    c_T12_re = (float)fabs(d_T11_im);
    V1_norm = (float)fabs(c_T22_re);
    if (c_T12_re > V1_norm) {
      s = c_T22_re / d_T11_im;
      V1_norm = d_T11_im + s * c_T22_re;
      a_re = (f_T11_re + s * b_T22_im) / V1_norm;
      a_im = (b_T22_im - s * f_T11_re) / V1_norm;
    } else if (V1_norm == c_T12_re) {
      if (d_T11_im > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (c_T22_re > 0.0F) {
        V1_norm = 0.5F;
      } else {
        V1_norm = -0.5F;
      }

      a_re = (f_T11_re * s + b_T22_im * V1_norm) / c_T12_re;
      a_im = (b_T22_im * s - f_T11_re * V1_norm) / c_T12_re;
    } else {
      s = d_T11_im / c_T22_re;
      V1_norm = c_T22_re + s * d_T11_im;
      a_re = (s * f_T11_re + b_T22_im) / V1_norm;
      a_im = (s * b_T22_im - f_T11_re) / V1_norm;
    }
  }

  b_T11_re = (H11 * H13 + H12 * H23) + H13 * H33;
  T11_im = ((H11 * 0.0F + 0.0F * H13) + (H12 * 0.0F + 0.0F * H23)) + (H13 * 0.0F
    + 0.0F * H33);
  T11_re = ((((T11_re + T22_re) + T33_re) + b_re) - b_a_re) - a_re;
  b_T11_im = (b_im - d_a_re) - a_im;
  c_T11_re = b_T11_re * T11_re - T11_im * b_T11_im;
  T11_im = b_T11_re * b_T11_im + T11_im * T11_re;
  T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * H13;
  b_T11_im = H11 * 0.0F + 0.0F * H13;
  c_T13_re = H13 * H22;
  v21_re = H13 * 0.0F + 0.0F * H22;
  d_T11_im = (((H12 * T13_re - 0.0F * b_T13_im) - (H12 * T23_re - 0.0F * T23_im))
              - (T11_re * H23 - b_T11_im * 0.0F)) + (c_T13_re * H23 - v21_re *
    0.0F);
  c_T22_re = (((H12 * b_T13_im + 0.0F * T13_re) - (H12 * T23_im + 0.0F * T23_re))
              - (T11_re * 0.0F + b_T11_im * H23)) + (c_T13_re * 0.0F + v21_re *
    H23);
  if (c_T22_re == 0.0F) {
    if (T11_im == 0.0F) {
      T11_re = c_T11_re / d_T11_im;
      T11_im = 0.0F;
    } else if (c_T11_re == 0.0F) {
      T11_re = 0.0F;
      T11_im /= d_T11_im;
    } else {
      T11_re = c_T11_re / d_T11_im;
      T11_im /= d_T11_im;
    }
  } else if (d_T11_im == 0.0F) {
    if (c_T11_re == 0.0F) {
      T11_re = T11_im / c_T22_re;
      T11_im = 0.0F;
    } else if (T11_im == 0.0F) {
      T11_re = 0.0F;
      T11_im = -(c_T11_re / c_T22_re);
    } else {
      T11_re = T11_im / c_T22_re;
      T11_im = -(c_T11_re / c_T22_re);
    }
  } else {
    c_T12_re = (float)fabs(d_T11_im);
    if (c_T12_re > c_T22_re) {
      s = c_T22_re / d_T11_im;
      V1_norm = d_T11_im + s * c_T22_re;
      T11_re = (c_T11_re + s * T11_im) / V1_norm;
      T11_im = (T11_im - s * c_T11_re) / V1_norm;
    } else if (c_T22_re == c_T12_re) {
      if (d_T11_im > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (c_T11_re * s + T11_im * -0.5F > 0.0F) {
        T11_re = ((real32_T)rtInf);
      } else if (c_T11_re * s + T11_im * -0.5F < 0.0F) {
        T11_re = ((real32_T)rtMinusInf);
      } else {
        T11_re = ((real32_T)rtNaN);
      }

      if (T11_im * s - c_T11_re * -0.5F > 0.0F) {
        T11_im = ((real32_T)rtInf);
      } else if (T11_im * s - c_T11_re * -0.5F < 0.0F) {
        T11_im = ((real32_T)rtMinusInf);
      } else {
        T11_im = ((real32_T)rtNaN);
      }
    } else {
      s = d_T11_im / c_T22_re;
      V1_norm = c_T22_re + s * d_T11_im;
      T11_re = (s * c_T11_re + T11_im) / V1_norm;
      T11_im = (s * T11_im - c_T11_re) / V1_norm;
    }
  }

  V23_re = (b_T13_re - re) - T11_re;
  c_T12_re = (T13_im - im) - T11_im;
  V[0].re = V11_re;
  V[0].im = V11_im;
  V[3].re = V12_re;
  V[3].im = V12_im;
  V[6].re = V13_re;
  V[6].im = V13_im;
  V[1].re = V21_re;
  V[1].im = V21_im;
  V[4].re = V22_re;
  V[4].im = V22_im;
  V[7].re = V23_re;
  V[7].im = c_T12_re;
  for (i0 = 0; i0 < 3; i0++) {
    V[2 + 3 * i0].re = 1.0F;
    V[2 + 3 * i0].im = 0.0F;
  }

  V1_norm = norm(*(creal32_T (*)[3])&V[0]);
  if (V11_im == 0.0F) {
    v11_re = V11_re / V1_norm;
    v11_im = 0.0F;
  } else if (V11_re == 0.0F) {
    v11_re = 0.0F;
    v11_im = V11_im / V1_norm;
  } else {
    v11_re = V11_re / V1_norm;
    v11_im = V11_im / V1_norm;
  }

  if (V21_im == 0.0F) {
    v12_re = V21_re / V1_norm;
    v12_im = 0.0F;
  } else if (V21_re == 0.0F) {
    v12_re = 0.0F;
    v12_im = V21_im / V1_norm;
  } else {
    v12_re = V21_re / V1_norm;
    v12_im = V21_im / V1_norm;
  }

  v13_re = 1.0F / V1_norm;
  V1_norm = norm(*(creal32_T (*)[3])&V[3]);
  if (V12_im == 0.0F) {
    v21_re = V12_re / V1_norm;
    V21_im = 0.0F;
  } else if (V12_re == 0.0F) {
    v21_re = 0.0F;
    V21_im = V12_im / V1_norm;
  } else {
    v21_re = V12_re / V1_norm;
    V21_im = V12_im / V1_norm;
  }

  if (V22_im == 0.0F) {
    v22_re = V22_re / V1_norm;
    v22_im = 0.0F;
  } else if (V22_re == 0.0F) {
    v22_re = 0.0F;
    v22_im = V22_im / V1_norm;
  } else {
    v22_re = V22_re / V1_norm;
    v22_im = V22_im / V1_norm;
  }

  v23_re = 1.0F / V1_norm;
  V1_norm = norm(*(creal32_T (*)[3])&V[6]);
  if (V13_im == 0.0F) {
    V13_re /= V1_norm;
    V13_im = 0.0F;
  } else if (V13_re == 0.0F) {
    V13_re = 0.0F;
    V13_im /= V1_norm;
  } else {
    V13_re /= V1_norm;
    V13_im /= V1_norm;
  }

  if (c_T12_re == 0.0F) {
    V11_re = V23_re / V1_norm;
    V11_im = 0.0F;
  } else if (V23_re == 0.0F) {
    V11_re = 0.0F;
    V11_im = c_T12_re / V1_norm;
  } else {
    V11_re = V23_re / V1_norm;
    V11_im = c_T12_re / V1_norm;
  }

  V21_re = 1.0F / V1_norm;
  fd_T11.re = (H11 + H22) + H33;
  fd_T11.im = 0.0F;
  fc0 = mpower(fd_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
           b_T12_im * 0.0F)) + c_T11_re;
  a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F +
           b_T12_im * H23)) + c_T11_im;
  b_a_re = (H11 + H22) + H33;
  c_a_re = (H11 + H22) + H33;
  gd_T11.re = (H11 + H22) + H33;
  gd_T11.im = 0.0F;
  fc0 = mpower(gd_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  d_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  b_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  e_a_re = (H11 + H22) + H33;
  T11_re = H11 / 3.0F;
  T22_re = H22 / 3.0F;
  T33_re = H33 / 3.0F;
  f_a_re = b_a_re * b_a_re;
  c_a_im = b_a_re * 0.0F + 0.0F * b_a_re;
  if (c_a_im == 0.0F) {
    b_a_re = f_a_re / 9.0F;
    c_a_im = 0.0F;
  } else if (f_a_re == 0.0F) {
    b_a_re = 0.0F;
    c_a_im /= 9.0F;
  } else {
    b_a_re = f_a_re / 9.0F;
    c_a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  if (T13_im == 0.0F) {
    T13_re /= 3.0F;
    T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  bd_a.re = (((((b_a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  bd_a.im = (((((c_a_im + T12_im) + T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  fc0 = mpower(bd_a);
  cd_a.re = (a_re * a_re - a_im * a_im) - fc0.re;
  cd_a.im = (a_re * a_im + a_im * a_re) - fc0.im;
  fc0 = b_mpower(cd_a);
  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  b_T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    b_T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    b_T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  hd_T11.re = (H11 + H22) + H33;
  hd_T11.im = 0.0F;
  T11 = mpower(hd_T11);
  if (T11.im == 0.0F) {
    re = T11.re / 27.0F;
    im = 0.0F;
  } else if (T11.re == 0.0F) {
    re = 0.0F;
    im = T11.im / 27.0F;
  } else {
    re = T11.re / 27.0F;
    im = T11.im / 27.0F;
  }

  c_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  d_T11_re = c_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = c_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  e_T11_re = d_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = d_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (e_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  fc41.re = ((((((fc0.re - T13_re) - T12_re) - b_T11_re) + re) + c_T11_re) +
             (b_T12_re * H23 - b_T12_im * 0.0F)) + d_T11_re;
  fc41.im = ((((((fc0.im - T13_im) - T12_im) - T11_im) + im) + b_T11_im) +
             (b_T12_re * 0.0F + b_T12_im * H23)) + c_T11_im;
  fc0 = c_mpower(fc41);
  a_re = c_a_re * c_a_re;
  a_im = c_a_re * 0.0F + 0.0F * c_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  if (T13_im == 0.0F) {
    T13_re /= 3.0F;
    T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  b_a_re = e_a_re * e_a_re;
  c_a_im = e_a_re * 0.0F + 0.0F * e_a_re;
  if (c_a_im == 0.0F) {
    b_a_re /= 9.0F;
    c_a_im = 0.0F;
  } else if (b_a_re == 0.0F) {
    b_a_re = 0.0F;
    c_a_im /= 9.0F;
  } else {
    b_a_re /= 9.0F;
    c_a_im /= 9.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  if (b_T12_im == 0.0F) {
    b_T12_re /= 3.0F;
    b_T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 3.0F;
  } else {
    b_T12_re /= 3.0F;
    b_T12_im /= 3.0F;
  }

  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  if (b_T13_im == 0.0F) {
    b_T13_re /= 3.0F;
    b_T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 3.0F;
  } else {
    b_T13_re /= 3.0F;
    b_T13_im /= 3.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  if (b_T23_im == 0.0F) {
    b_T23_re /= 3.0F;
    b_T23_im = 0.0F;
  } else if (b_T23_re == 0.0F) {
    b_T23_re = 0.0F;
    b_T23_im /= 3.0F;
  } else {
    b_T23_re /= 3.0F;
    b_T23_im /= 3.0F;
  }

  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 3.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 3.0F;
  } else {
    d_T11_re /= 3.0F;
    c_T11_im /= 3.0F;
  }

  e_T11_re = H11 * H33;
  c_a_re = H11 * 0.0F + 0.0F * H33;
  if (c_a_re == 0.0F) {
    e_T11_re /= 3.0F;
    c_a_re = 0.0F;
  } else if (e_T11_re == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 3.0F;
  } else {
    e_T11_re /= 3.0F;
    c_a_re /= 3.0F;
  }

  c_T22_re = H22 * H33;
  b_T22_im = H22 * 0.0F + 0.0F * H33;
  if (b_T22_im == 0.0F) {
    c_T22_re /= 3.0F;
    b_T22_im = 0.0F;
  } else if (c_T22_re == 0.0F) {
    c_T22_re = 0.0F;
    b_T22_im /= 3.0F;
  } else {
    c_T22_re /= 3.0F;
    b_T22_im /= 3.0F;
  }

  dd_a.re = (((((b_a_re + b_T12_re) + b_T13_re) + b_T23_re) - d_T11_re) -
             e_T11_re) - c_T22_re;
  dd_a.im = (((((c_a_im + b_T12_im) + b_T13_im) + b_T23_im) - c_T11_im) - c_a_re)
    - b_T22_im;
  T11 = mpower(dd_a);
  ed_a.re = (d_a_re * d_a_re - b_a_im * b_a_im) - T11.re;
  ed_a.im = (d_a_re * b_a_im + b_a_im * d_a_re) - T11.im;
  T11 = b_mpower(ed_a);
  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  c_T13_re = b_T13_re * H22 - b_T13_im * 0.0F;
  b_T13_im = b_T13_re * 0.0F + b_T13_im * H22;
  if (b_T13_im == 0.0F) {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im = 0.0F;
  } else if (c_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 2.0F;
  } else {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im /= 2.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  c_T12_re = b_T12_re * H33 - b_T12_im * 0.0F;
  b_T12_im = b_T12_re * 0.0F + b_T12_im * H33;
  if (b_T12_im == 0.0F) {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im = 0.0F;
  } else if (c_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 2.0F;
  } else {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im /= 2.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  d_T11_re = H11 * b_T23_re - 0.0F * b_T23_im;
  c_T11_im = H11 * b_T23_im + 0.0F * b_T23_re;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re /= 2.0F;
    c_T11_im /= 2.0F;
  }

  id_T11.re = (H11 + H22) + H33;
  id_T11.im = 0.0F;
  fc1 = mpower(id_T11);
  if (fc1.im == 0.0F) {
    re = fc1.re / 27.0F;
    im = 0.0F;
  } else if (fc1.re == 0.0F) {
    re = 0.0F;
    im = fc1.im / 27.0F;
  } else {
    re = fc1.re / 27.0F;
    im = fc1.im / 27.0F;
  }

  e_T11_re = (H11 + H22) + H33;
  c_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  c_T22_re = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  b_T22_im = e_T11_re * c_T12_re - 0.0F * c_T22_re;
  c_a_re = e_T11_re * c_T22_re + 0.0F * c_T12_re;
  if (c_a_re == 0.0F) {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re = 0.0F;
  } else if (b_T22_im == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 6.0F;
  } else {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re /= 6.0F;
  }

  c_T12_re = H12 * H13;
  c_T22_re = H12 * 0.0F + 0.0F * H13;
  b_T22_im = H11 * H22;
  d_T11_im = H11 * 0.0F + 0.0F * H22;
  f_T11_re = b_T22_im * H33 - d_T11_im * 0.0F;
  d_T11_im = b_T22_im * 0.0F + d_T11_im * H33;
  if (d_T11_im == 0.0F) {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im = 0.0F;
  } else if (f_T11_re == 0.0F) {
    b_T22_im = 0.0F;
    d_T11_im /= 2.0F;
  } else {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im /= 2.0F;
  }

  fc42.re = ((((((T11.re - b_T13_re) - b_T12_re) - d_T11_re) + re) + e_T11_re) +
             (c_T12_re * H23 - c_T22_re * 0.0F)) + b_T22_im;
  fc42.im = ((((((T11.im - b_T13_im) - b_T12_im) - c_T11_im) + im) + c_a_re) +
             (c_T12_re * 0.0F + c_T22_re * H23)) + d_T11_im;
  T11 = c_mpower(fc42);
  f_T11_re = (((((a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  b_T22_im = (((((a_im + T12_im) + T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  if (T11.im == 0.0F) {
    if (b_T22_im == 0.0F) {
      a_re = f_T11_re / T11.re;
      a_im = 0.0F;
    } else if (f_T11_re == 0.0F) {
      a_re = 0.0F;
      a_im = b_T22_im / T11.re;
    } else {
      a_re = f_T11_re / T11.re;
      a_im = b_T22_im / T11.re;
    }
  } else if (T11.re == 0.0F) {
    if (f_T11_re == 0.0F) {
      a_re = b_T22_im / T11.im;
      a_im = 0.0F;
    } else if (b_T22_im == 0.0F) {
      a_re = 0.0F;
      a_im = -(f_T11_re / T11.im);
    } else {
      a_re = b_T22_im / T11.im;
      a_im = -(f_T11_re / T11.im);
    }
  } else {
    c_T12_re = (float)fabs(T11.re);
    V1_norm = (float)fabs(T11.im);
    if (c_T12_re > V1_norm) {
      s = T11.im / T11.re;
      V1_norm = T11.re + s * T11.im;
      a_re = (f_T11_re + s * b_T22_im) / V1_norm;
      a_im = (b_T22_im - s * f_T11_re) / V1_norm;
    } else if (V1_norm == c_T12_re) {
      if (T11.re > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (T11.im > 0.0F) {
        V1_norm = 0.5F;
      } else {
        V1_norm = -0.5F;
      }

      a_re = (f_T11_re * s + b_T22_im * V1_norm) / c_T12_re;
      a_im = (b_T22_im * s - f_T11_re * V1_norm) / c_T12_re;
    } else {
      s = T11.re / T11.im;
      V1_norm = T11.im + s * T11.re;
      a_re = (s * f_T11_re + b_T22_im) / V1_norm;
      a_im = (s * b_T22_im - f_T11_re) / V1_norm;
    }
  }

  V12_re = (((T11_re + T22_re) + T33_re) + fc0.re) + a_re;
  V12_im = fc0.im + a_im;
  jd_T11.re = (H11 + H22) + H33;
  jd_T11.im = 0.0F;
  fc0 = mpower(jd_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
           b_T12_im * 0.0F)) + c_T11_re;
  a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F +
           b_T12_im * H23)) + c_T11_im;
  b_a_re = (H11 + H22) + H33;
  c_a_re = (H11 + H22) + H33;
  kd_T11.re = (H11 + H22) + H33;
  kd_T11.im = 0.0F;
  fc0 = mpower(kd_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  d_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  b_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  e_a_re = (H11 + H22) + H33;
  ld_T11.re = (H11 + H22) + H33;
  ld_T11.im = 0.0F;
  fc0 = mpower(ld_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  f_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  c_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  g_a_re = (H11 + H22) + H33;
  h_a_re = (H11 + H22) + H33;
  md_T11.re = (H11 + H22) + H33;
  md_T11.im = 0.0F;
  fc0 = mpower(md_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  i_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  d_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  j_a_re = (H11 + H22) + H33;
  T11_re = H11 / 3.0F;
  T22_re = H22 / 3.0F;
  T33_re = H33 / 3.0F;
  V1_norm = b_a_re * b_a_re;
  e_a_im = b_a_re * 0.0F + 0.0F * b_a_re;
  if (e_a_im == 0.0F) {
    b_a_re = V1_norm / 9.0F;
    e_a_im = 0.0F;
  } else if (V1_norm == 0.0F) {
    b_a_re = 0.0F;
    e_a_im /= 9.0F;
  } else {
    b_a_re = V1_norm / 9.0F;
    e_a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  if (T13_im == 0.0F) {
    T13_re /= 3.0F;
    T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  fd_a.re = (((((b_a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  fd_a.im = (((((e_a_im + T12_im) + T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  fc0 = mpower(fd_a);
  gd_a.re = (a_re * a_re - a_im * a_im) - fc0.re;
  gd_a.im = (a_re * a_im + a_im * a_re) - fc0.im;
  fc0 = b_mpower(gd_a);
  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  b_T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    b_T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    b_T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  nd_T11.re = (H11 + H22) + H33;
  nd_T11.im = 0.0F;
  T11 = mpower(nd_T11);
  if (T11.im == 0.0F) {
    re = T11.re / 27.0F;
    im = 0.0F;
  } else if (T11.re == 0.0F) {
    re = 0.0F;
    im = T11.im / 27.0F;
  } else {
    re = T11.re / 27.0F;
    im = T11.im / 27.0F;
  }

  c_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  d_T11_re = c_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = c_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  e_T11_re = d_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = d_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (e_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  fc43.re = ((((((fc0.re - T13_re) - T12_re) - b_T11_re) + re) + c_T11_re) +
             (b_T12_re * H23 - b_T12_im * 0.0F)) + d_T11_re;
  fc43.im = ((((((fc0.im - T13_im) - T12_im) - T11_im) + im) + b_T11_im) +
             (b_T12_re * 0.0F + b_T12_im * H23)) + c_T11_im;
  fc0 = c_mpower(fc43);
  a_re = c_a_re * c_a_re;
  a_im = c_a_re * 0.0F + 0.0F * c_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  if (T13_im == 0.0F) {
    T13_re /= 3.0F;
    T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  b_a_re = e_a_re * e_a_re;
  e_a_im = e_a_re * 0.0F + 0.0F * e_a_re;
  if (e_a_im == 0.0F) {
    b_a_re /= 9.0F;
    e_a_im = 0.0F;
  } else if (b_a_re == 0.0F) {
    b_a_re = 0.0F;
    e_a_im /= 9.0F;
  } else {
    b_a_re /= 9.0F;
    e_a_im /= 9.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  if (b_T12_im == 0.0F) {
    b_T12_re /= 3.0F;
    b_T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 3.0F;
  } else {
    b_T12_re /= 3.0F;
    b_T12_im /= 3.0F;
  }

  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  if (b_T13_im == 0.0F) {
    b_T13_re /= 3.0F;
    b_T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 3.0F;
  } else {
    b_T13_re /= 3.0F;
    b_T13_im /= 3.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  if (b_T23_im == 0.0F) {
    b_T23_re /= 3.0F;
    b_T23_im = 0.0F;
  } else if (b_T23_re == 0.0F) {
    b_T23_re = 0.0F;
    b_T23_im /= 3.0F;
  } else {
    b_T23_re /= 3.0F;
    b_T23_im /= 3.0F;
  }

  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 3.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 3.0F;
  } else {
    d_T11_re /= 3.0F;
    c_T11_im /= 3.0F;
  }

  e_T11_re = H11 * H33;
  c_a_re = H11 * 0.0F + 0.0F * H33;
  if (c_a_re == 0.0F) {
    e_T11_re /= 3.0F;
    c_a_re = 0.0F;
  } else if (e_T11_re == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 3.0F;
  } else {
    e_T11_re /= 3.0F;
    c_a_re /= 3.0F;
  }

  c_T22_re = H22 * H33;
  b_T22_im = H22 * 0.0F + 0.0F * H33;
  if (b_T22_im == 0.0F) {
    c_T22_re /= 3.0F;
    b_T22_im = 0.0F;
  } else if (c_T22_re == 0.0F) {
    c_T22_re = 0.0F;
    b_T22_im /= 3.0F;
  } else {
    c_T22_re /= 3.0F;
    b_T22_im /= 3.0F;
  }

  hd_a.re = (((((b_a_re + b_T12_re) + b_T13_re) + b_T23_re) - d_T11_re) -
             e_T11_re) - c_T22_re;
  hd_a.im = (((((e_a_im + b_T12_im) + b_T13_im) + b_T23_im) - c_T11_im) - c_a_re)
    - b_T22_im;
  T11 = mpower(hd_a);
  id_a.re = (d_a_re * d_a_re - b_a_im * b_a_im) - T11.re;
  id_a.im = (d_a_re * b_a_im + b_a_im * d_a_re) - T11.im;
  T11 = b_mpower(id_a);
  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  c_T13_re = b_T13_re * H22 - b_T13_im * 0.0F;
  b_T13_im = b_T13_re * 0.0F + b_T13_im * H22;
  if (b_T13_im == 0.0F) {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im = 0.0F;
  } else if (c_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 2.0F;
  } else {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im /= 2.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  c_T12_re = b_T12_re * H33 - b_T12_im * 0.0F;
  b_T12_im = b_T12_re * 0.0F + b_T12_im * H33;
  if (b_T12_im == 0.0F) {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im = 0.0F;
  } else if (c_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 2.0F;
  } else {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im /= 2.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  d_T11_re = H11 * b_T23_re - 0.0F * b_T23_im;
  c_T11_im = H11 * b_T23_im + 0.0F * b_T23_re;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re /= 2.0F;
    c_T11_im /= 2.0F;
  }

  od_T11.re = (H11 + H22) + H33;
  od_T11.im = 0.0F;
  fc1 = mpower(od_T11);
  if (fc1.im == 0.0F) {
    re = fc1.re / 27.0F;
    im = 0.0F;
  } else if (fc1.re == 0.0F) {
    re = 0.0F;
    im = fc1.im / 27.0F;
  } else {
    re = fc1.re / 27.0F;
    im = fc1.im / 27.0F;
  }

  e_T11_re = (H11 + H22) + H33;
  c_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  c_T22_re = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  b_T22_im = e_T11_re * c_T12_re - 0.0F * c_T22_re;
  c_a_re = e_T11_re * c_T22_re + 0.0F * c_T12_re;
  if (c_a_re == 0.0F) {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re = 0.0F;
  } else if (b_T22_im == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 6.0F;
  } else {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re /= 6.0F;
  }

  c_T12_re = H12 * H13;
  c_T22_re = H12 * 0.0F + 0.0F * H13;
  b_T22_im = H11 * H22;
  d_T11_im = H11 * 0.0F + 0.0F * H22;
  f_T11_re = b_T22_im * H33 - d_T11_im * 0.0F;
  d_T11_im = b_T22_im * 0.0F + d_T11_im * H33;
  if (d_T11_im == 0.0F) {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im = 0.0F;
  } else if (f_T11_re == 0.0F) {
    b_T22_im = 0.0F;
    d_T11_im /= 2.0F;
  } else {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im /= 2.0F;
  }

  fc44.re = ((((((T11.re - b_T13_re) - b_T12_re) - d_T11_re) + re) + e_T11_re) +
             (c_T12_re * H23 - c_T22_re * 0.0F)) + b_T22_im;
  fc44.im = ((((((T11.im - b_T13_im) - b_T12_im) - c_T11_im) + im) + c_a_re) +
             (c_T12_re * 0.0F + c_T22_re * H23)) + d_T11_im;
  T11 = c_mpower(fc44);
  f_T11_re = (((((a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  b_T22_im = (((((a_im + T12_im) + T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  if (T11.im == 0.0F) {
    if (b_T22_im == 0.0F) {
      a_re = f_T11_re / T11.re;
      a_im = 0.0F;
    } else if (f_T11_re == 0.0F) {
      a_re = 0.0F;
      a_im = b_T22_im / T11.re;
    } else {
      a_re = f_T11_re / T11.re;
      a_im = b_T22_im / T11.re;
    }
  } else if (T11.re == 0.0F) {
    if (f_T11_re == 0.0F) {
      a_re = b_T22_im / T11.im;
      a_im = 0.0F;
    } else if (b_T22_im == 0.0F) {
      a_re = 0.0F;
      a_im = -(f_T11_re / T11.im);
    } else {
      a_re = b_T22_im / T11.im;
      a_im = -(f_T11_re / T11.im);
    }
  } else {
    c_T12_re = (float)fabs(T11.re);
    V1_norm = (float)fabs(T11.im);
    if (c_T12_re > V1_norm) {
      s = T11.im / T11.re;
      V1_norm = T11.re + s * T11.im;
      a_re = (f_T11_re + s * b_T22_im) / V1_norm;
      a_im = (b_T22_im - s * f_T11_re) / V1_norm;
    } else if (V1_norm == c_T12_re) {
      if (T11.re > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (T11.im > 0.0F) {
        V1_norm = 0.5F;
      } else {
        V1_norm = -0.5F;
      }

      a_re = (f_T11_re * s + b_T22_im * V1_norm) / c_T12_re;
      a_im = (b_T22_im * s - f_T11_re * V1_norm) / c_T12_re;
    } else {
      s = T11.re / T11.im;
      V1_norm = T11.im + s * T11.re;
      a_re = (s * f_T11_re + b_T22_im) / V1_norm;
      a_im = (s * b_T22_im - f_T11_re) / V1_norm;
    }
  }

  re = 1.73205078F * (fc0.re - a_re);
  im = 1.73205078F * (fc0.im - a_im);
  b_re = re * 0.0F - im;
  im = re + im * 0.0F;
  if (im == 0.0F) {
    re = b_re / 2.0F;
    im = 0.0F;
  } else if (b_re == 0.0F) {
    re = 0.0F;
    im /= 2.0F;
  } else {
    re = b_re / 2.0F;
    im /= 2.0F;
  }

  a_re = g_a_re * g_a_re;
  a_im = g_a_re * 0.0F + 0.0F * g_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  if (T13_im == 0.0F) {
    T13_re /= 3.0F;
    T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  jd_a.re = (((((a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  jd_a.im = (((((a_im + T12_im) + T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  fc0 = mpower(jd_a);
  kd_a.re = (f_a_re * f_a_re - c_a_im * c_a_im) - fc0.re;
  kd_a.im = (f_a_re * c_a_im + c_a_im * f_a_re) - fc0.im;
  fc0 = b_mpower(kd_a);
  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  b_T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    b_T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    b_T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  pd_T11.re = (H11 + H22) + H33;
  pd_T11.im = 0.0F;
  T11 = mpower(pd_T11);
  if (T11.im == 0.0F) {
    b_re = T11.re / 27.0F;
    b_im = 0.0F;
  } else if (T11.re == 0.0F) {
    b_re = 0.0F;
    b_im = T11.im / 27.0F;
  } else {
    b_re = T11.re / 27.0F;
    b_im = T11.im / 27.0F;
  }

  c_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  d_T11_re = c_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = c_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  e_T11_re = d_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = d_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (e_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  fc45.re = ((((((fc0.re - T13_re) - T12_re) - b_T11_re) + b_re) + c_T11_re) +
             (b_T12_re * H23 - b_T12_im * 0.0F)) + d_T11_re;
  fc45.im = ((((((fc0.im - T13_im) - T12_im) - T11_im) + b_im) + b_T11_im) +
             (b_T12_re * 0.0F + b_T12_im * H23)) + c_T11_im;
  fc0 = c_mpower(fc45);
  if (fc0.im == 0.0F) {
    b_re = fc0.re / 2.0F;
    b_im = 0.0F;
  } else if (fc0.re == 0.0F) {
    b_re = 0.0F;
    b_im = fc0.im / 2.0F;
  } else {
    b_re = fc0.re / 2.0F;
    b_im = fc0.im / 2.0F;
  }

  a_re = h_a_re * h_a_re;
  a_im = h_a_re * 0.0F + 0.0F * h_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  if (T13_im == 0.0F) {
    T13_re /= 3.0F;
    T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  b_a_re = j_a_re * j_a_re;
  b_a_im = j_a_re * 0.0F + 0.0F * j_a_re;
  if (b_a_im == 0.0F) {
    b_a_re /= 9.0F;
    b_a_im = 0.0F;
  } else if (b_a_re == 0.0F) {
    b_a_re = 0.0F;
    b_a_im /= 9.0F;
  } else {
    b_a_re /= 9.0F;
    b_a_im /= 9.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  if (b_T12_im == 0.0F) {
    b_T12_re /= 3.0F;
    b_T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 3.0F;
  } else {
    b_T12_re /= 3.0F;
    b_T12_im /= 3.0F;
  }

  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  if (b_T13_im == 0.0F) {
    b_T13_re /= 3.0F;
    b_T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 3.0F;
  } else {
    b_T13_re /= 3.0F;
    b_T13_im /= 3.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  if (b_T23_im == 0.0F) {
    b_T23_re /= 3.0F;
    b_T23_im = 0.0F;
  } else if (b_T23_re == 0.0F) {
    b_T23_re = 0.0F;
    b_T23_im /= 3.0F;
  } else {
    b_T23_re /= 3.0F;
    b_T23_im /= 3.0F;
  }

  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 3.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 3.0F;
  } else {
    d_T11_re /= 3.0F;
    c_T11_im /= 3.0F;
  }

  e_T11_re = H11 * H33;
  c_a_re = H11 * 0.0F + 0.0F * H33;
  if (c_a_re == 0.0F) {
    e_T11_re /= 3.0F;
    c_a_re = 0.0F;
  } else if (e_T11_re == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 3.0F;
  } else {
    e_T11_re /= 3.0F;
    c_a_re /= 3.0F;
  }

  c_T22_re = H22 * H33;
  b_T22_im = H22 * 0.0F + 0.0F * H33;
  if (b_T22_im == 0.0F) {
    c_T22_re /= 3.0F;
    b_T22_im = 0.0F;
  } else if (c_T22_re == 0.0F) {
    c_T22_re = 0.0F;
    b_T22_im /= 3.0F;
  } else {
    c_T22_re /= 3.0F;
    b_T22_im /= 3.0F;
  }

  ld_a.re = (((((b_a_re + b_T12_re) + b_T13_re) + b_T23_re) - d_T11_re) -
             e_T11_re) - c_T22_re;
  ld_a.im = (((((b_a_im + b_T12_im) + b_T13_im) + b_T23_im) - c_T11_im) - c_a_re)
    - b_T22_im;
  fc0 = mpower(ld_a);
  md_a.re = (i_a_re * i_a_re - d_a_im * d_a_im) - fc0.re;
  md_a.im = (i_a_re * d_a_im + d_a_im * i_a_re) - fc0.im;
  fc0 = b_mpower(md_a);
  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  c_T13_re = b_T13_re * H22 - b_T13_im * 0.0F;
  b_T13_im = b_T13_re * 0.0F + b_T13_im * H22;
  if (b_T13_im == 0.0F) {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im = 0.0F;
  } else if (c_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 2.0F;
  } else {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im /= 2.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  c_T12_re = b_T12_re * H33 - b_T12_im * 0.0F;
  b_T12_im = b_T12_re * 0.0F + b_T12_im * H33;
  if (b_T12_im == 0.0F) {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im = 0.0F;
  } else if (c_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 2.0F;
  } else {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im /= 2.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  d_T11_re = H11 * b_T23_re - 0.0F * b_T23_im;
  c_T11_im = H11 * b_T23_im + 0.0F * b_T23_re;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re /= 2.0F;
    c_T11_im /= 2.0F;
  }

  qd_T11.re = (H11 + H22) + H33;
  qd_T11.im = 0.0F;
  T11 = mpower(qd_T11);
  if (T11.im == 0.0F) {
    b_a_re = T11.re / 27.0F;
    d_a_re = 0.0F;
  } else if (T11.re == 0.0F) {
    b_a_re = 0.0F;
    d_a_re = T11.im / 27.0F;
  } else {
    b_a_re = T11.re / 27.0F;
    d_a_re = T11.im / 27.0F;
  }

  e_T11_re = (H11 + H22) + H33;
  c_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  c_T22_re = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  b_T22_im = e_T11_re * c_T12_re - 0.0F * c_T22_re;
  c_a_re = e_T11_re * c_T22_re + 0.0F * c_T12_re;
  if (c_a_re == 0.0F) {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re = 0.0F;
  } else if (b_T22_im == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 6.0F;
  } else {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re /= 6.0F;
  }

  c_T12_re = H12 * H13;
  c_T22_re = H12 * 0.0F + 0.0F * H13;
  b_T22_im = H11 * H22;
  d_T11_im = H11 * 0.0F + 0.0F * H22;
  f_T11_re = b_T22_im * H33 - d_T11_im * 0.0F;
  d_T11_im = b_T22_im * 0.0F + d_T11_im * H33;
  if (d_T11_im == 0.0F) {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im = 0.0F;
  } else if (f_T11_re == 0.0F) {
    b_T22_im = 0.0F;
    d_T11_im /= 2.0F;
  } else {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im /= 2.0F;
  }

  fc46.re = ((((((fc0.re - b_T13_re) - b_T12_re) - d_T11_re) + b_a_re) +
              e_T11_re) + (c_T12_re * H23 - c_T22_re * 0.0F)) + b_T22_im;
  fc46.im = ((((((fc0.im - b_T13_im) - b_T12_im) - c_T11_im) + d_a_re) + c_a_re)
             + (c_T12_re * 0.0F + c_T22_re * H23)) + d_T11_im;
  fc0 = c_mpower(fc46);
  f_T11_re = (((((a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  b_T22_im = (((((a_im + T12_im) + T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  d_T11_im = 2.0F * fc0.re;
  c_T22_re = 2.0F * fc0.im;
  if (c_T22_re == 0.0F) {
    if (b_T22_im == 0.0F) {
      a_re = f_T11_re / d_T11_im;
      a_im = 0.0F;
    } else if (f_T11_re == 0.0F) {
      a_re = 0.0F;
      a_im = b_T22_im / d_T11_im;
    } else {
      a_re = f_T11_re / d_T11_im;
      a_im = b_T22_im / d_T11_im;
    }
  } else if (d_T11_im == 0.0F) {
    if (f_T11_re == 0.0F) {
      a_re = b_T22_im / c_T22_re;
      a_im = 0.0F;
    } else if (b_T22_im == 0.0F) {
      a_re = 0.0F;
      a_im = -(f_T11_re / c_T22_re);
    } else {
      a_re = b_T22_im / c_T22_re;
      a_im = -(f_T11_re / c_T22_re);
    }
  } else {
    c_T12_re = (float)fabs(d_T11_im);
    V1_norm = (float)fabs(c_T22_re);
    if (c_T12_re > V1_norm) {
      s = c_T22_re / d_T11_im;
      V1_norm = d_T11_im + s * c_T22_re;
      a_re = (f_T11_re + s * b_T22_im) / V1_norm;
      a_im = (b_T22_im - s * f_T11_re) / V1_norm;
    } else if (V1_norm == c_T12_re) {
      if (d_T11_im > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (c_T22_re > 0.0F) {
        V1_norm = 0.5F;
      } else {
        V1_norm = -0.5F;
      }

      a_re = (f_T11_re * s + b_T22_im * V1_norm) / c_T12_re;
      a_im = (b_T22_im * s - f_T11_re * V1_norm) / c_T12_re;
    } else {
      s = d_T11_im / c_T22_re;
      V1_norm = c_T22_re + s * d_T11_im;
      a_re = (s * f_T11_re + b_T22_im) / V1_norm;
      a_im = (s * b_T22_im - f_T11_re) / V1_norm;
    }
  }

  V22_re = ((((T11_re + T22_re) + T33_re) - re) - b_re) - a_re;
  V22_im = ((0.0F - im) - b_im) - a_im;
  rd_T11.re = (H11 + H22) + H33;
  rd_T11.im = 0.0F;
  fc0 = mpower(rd_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
           b_T12_im * 0.0F)) + c_T11_re;
  a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F +
           b_T12_im * H23)) + c_T11_im;
  b_a_re = (H11 + H22) + H33;
  c_a_re = (H11 + H22) + H33;
  sd_T11.re = (H11 + H22) + H33;
  sd_T11.im = 0.0F;
  fc0 = mpower(sd_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  d_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  b_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  e_a_re = (H11 + H22) + H33;
  td_T11.re = (H11 + H22) + H33;
  td_T11.im = 0.0F;
  fc0 = mpower(td_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  f_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  c_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  g_a_re = (H11 + H22) + H33;
  h_a_re = (H11 + H22) + H33;
  ud_T11.re = (H11 + H22) + H33;
  ud_T11.im = 0.0F;
  fc0 = mpower(ud_T11);
  if (fc0.im == 0.0F) {
    re = fc0.re / 27.0F;
    im = 0.0F;
  } else if (fc0.re == 0.0F) {
    re = 0.0F;
    im = fc0.im / 27.0F;
  } else {
    re = fc0.re / 27.0F;
    im = fc0.im / 27.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (T11_re == 0.0F) {
    T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  b_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  c_T11_re = b_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = b_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    b_T11_re = c_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  c_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  d_T11_re = c_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = c_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    c_T11_re = d_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  i_a_re = (((((re - T13_re) - T12_re) - T11_re) + b_T11_re) + (b_T12_re * H23 -
             b_T12_im * 0.0F)) + c_T11_re;
  d_a_im = (((((im - T13_im) - T12_im) - T11_im) + b_T11_im) + (b_T12_re * 0.0F
             + b_T12_im * H23)) + c_T11_im;
  j_a_re = (H11 + H22) + H33;
  T11_re = H11 / 3.0F;
  T22_re = H22 / 3.0F;
  T33_re = H33 / 3.0F;
  V1_norm = b_a_re * b_a_re;
  e_a_im = b_a_re * 0.0F + 0.0F * b_a_re;
  if (e_a_im == 0.0F) {
    b_a_re = V1_norm / 9.0F;
    e_a_im = 0.0F;
  } else if (V1_norm == 0.0F) {
    b_a_re = 0.0F;
    e_a_im /= 9.0F;
  } else {
    b_a_re = V1_norm / 9.0F;
    e_a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  if (T13_im == 0.0F) {
    T13_re /= 3.0F;
    T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  nd_a.re = (((((b_a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  nd_a.im = (((((e_a_im + T12_im) + T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  fc0 = mpower(nd_a);
  od_a.re = (a_re * a_re - a_im * a_im) - fc0.re;
  od_a.im = (a_re * a_im + a_im * a_re) - fc0.im;
  fc0 = b_mpower(od_a);
  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  b_T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    b_T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    b_T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  vd_T11.re = (H11 + H22) + H33;
  vd_T11.im = 0.0F;
  T11 = mpower(vd_T11);
  if (T11.im == 0.0F) {
    re = T11.re / 27.0F;
    im = 0.0F;
  } else if (T11.re == 0.0F) {
    re = 0.0F;
    im = T11.im / 27.0F;
  } else {
    re = T11.re / 27.0F;
    im = T11.im / 27.0F;
  }

  c_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  d_T11_re = c_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = c_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  e_T11_re = d_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = d_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (e_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  fc47.re = ((((((fc0.re - T13_re) - T12_re) - b_T11_re) + re) + c_T11_re) +
             (b_T12_re * H23 - b_T12_im * 0.0F)) + d_T11_re;
  fc47.im = ((((((fc0.im - T13_im) - T12_im) - T11_im) + im) + b_T11_im) +
             (b_T12_re * 0.0F + b_T12_im * H23)) + c_T11_im;
  fc0 = c_mpower(fc47);
  a_re = c_a_re * c_a_re;
  a_im = c_a_re * 0.0F + 0.0F * c_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  if (T13_im == 0.0F) {
    T13_re /= 3.0F;
    T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  b_a_re = e_a_re * e_a_re;
  e_a_im = e_a_re * 0.0F + 0.0F * e_a_re;
  if (e_a_im == 0.0F) {
    b_a_re /= 9.0F;
    e_a_im = 0.0F;
  } else if (b_a_re == 0.0F) {
    b_a_re = 0.0F;
    e_a_im /= 9.0F;
  } else {
    b_a_re /= 9.0F;
    e_a_im /= 9.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  if (b_T12_im == 0.0F) {
    b_T12_re /= 3.0F;
    b_T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 3.0F;
  } else {
    b_T12_re /= 3.0F;
    b_T12_im /= 3.0F;
  }

  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  if (b_T13_im == 0.0F) {
    b_T13_re /= 3.0F;
    b_T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 3.0F;
  } else {
    b_T13_re /= 3.0F;
    b_T13_im /= 3.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  if (b_T23_im == 0.0F) {
    b_T23_re /= 3.0F;
    b_T23_im = 0.0F;
  } else if (b_T23_re == 0.0F) {
    b_T23_re = 0.0F;
    b_T23_im /= 3.0F;
  } else {
    b_T23_re /= 3.0F;
    b_T23_im /= 3.0F;
  }

  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 3.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 3.0F;
  } else {
    d_T11_re /= 3.0F;
    c_T11_im /= 3.0F;
  }

  e_T11_re = H11 * H33;
  c_a_re = H11 * 0.0F + 0.0F * H33;
  if (c_a_re == 0.0F) {
    e_T11_re /= 3.0F;
    c_a_re = 0.0F;
  } else if (e_T11_re == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 3.0F;
  } else {
    e_T11_re /= 3.0F;
    c_a_re /= 3.0F;
  }

  c_T22_re = H22 * H33;
  b_T22_im = H22 * 0.0F + 0.0F * H33;
  if (b_T22_im == 0.0F) {
    c_T22_re /= 3.0F;
    b_T22_im = 0.0F;
  } else if (c_T22_re == 0.0F) {
    c_T22_re = 0.0F;
    b_T22_im /= 3.0F;
  } else {
    c_T22_re /= 3.0F;
    b_T22_im /= 3.0F;
  }

  pd_a.re = (((((b_a_re + b_T12_re) + b_T13_re) + b_T23_re) - d_T11_re) -
             e_T11_re) - c_T22_re;
  pd_a.im = (((((e_a_im + b_T12_im) + b_T13_im) + b_T23_im) - c_T11_im) - c_a_re)
    - b_T22_im;
  T11 = mpower(pd_a);
  qd_a.re = (d_a_re * d_a_re - b_a_im * b_a_im) - T11.re;
  qd_a.im = (d_a_re * b_a_im + b_a_im * d_a_re) - T11.im;
  T11 = b_mpower(qd_a);
  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  c_T13_re = b_T13_re * H22 - b_T13_im * 0.0F;
  b_T13_im = b_T13_re * 0.0F + b_T13_im * H22;
  if (b_T13_im == 0.0F) {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im = 0.0F;
  } else if (c_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 2.0F;
  } else {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im /= 2.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  c_T12_re = b_T12_re * H33 - b_T12_im * 0.0F;
  b_T12_im = b_T12_re * 0.0F + b_T12_im * H33;
  if (b_T12_im == 0.0F) {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im = 0.0F;
  } else if (c_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 2.0F;
  } else {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im /= 2.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  d_T11_re = H11 * b_T23_re - 0.0F * b_T23_im;
  c_T11_im = H11 * b_T23_im + 0.0F * b_T23_re;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re /= 2.0F;
    c_T11_im /= 2.0F;
  }

  wd_T11.re = (H11 + H22) + H33;
  wd_T11.im = 0.0F;
  fc1 = mpower(wd_T11);
  if (fc1.im == 0.0F) {
    re = fc1.re / 27.0F;
    im = 0.0F;
  } else if (fc1.re == 0.0F) {
    re = 0.0F;
    im = fc1.im / 27.0F;
  } else {
    re = fc1.re / 27.0F;
    im = fc1.im / 27.0F;
  }

  e_T11_re = (H11 + H22) + H33;
  c_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  c_T22_re = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  b_T22_im = e_T11_re * c_T12_re - 0.0F * c_T22_re;
  c_a_re = e_T11_re * c_T22_re + 0.0F * c_T12_re;
  if (c_a_re == 0.0F) {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re = 0.0F;
  } else if (b_T22_im == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 6.0F;
  } else {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re /= 6.0F;
  }

  c_T12_re = H12 * H13;
  c_T22_re = H12 * 0.0F + 0.0F * H13;
  b_T22_im = H11 * H22;
  d_T11_im = H11 * 0.0F + 0.0F * H22;
  f_T11_re = b_T22_im * H33 - d_T11_im * 0.0F;
  d_T11_im = b_T22_im * 0.0F + d_T11_im * H33;
  if (d_T11_im == 0.0F) {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im = 0.0F;
  } else if (f_T11_re == 0.0F) {
    b_T22_im = 0.0F;
    d_T11_im /= 2.0F;
  } else {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im /= 2.0F;
  }

  fc48.re = ((((((T11.re - b_T13_re) - b_T12_re) - d_T11_re) + re) + e_T11_re) +
             (c_T12_re * H23 - c_T22_re * 0.0F)) + b_T22_im;
  fc48.im = ((((((T11.im - b_T13_im) - b_T12_im) - c_T11_im) + im) + c_a_re) +
             (c_T12_re * 0.0F + c_T22_re * H23)) + d_T11_im;
  T11 = c_mpower(fc48);
  f_T11_re = (((((a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  b_T22_im = (((((a_im + T12_im) + T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  if (T11.im == 0.0F) {
    if (b_T22_im == 0.0F) {
      a_re = f_T11_re / T11.re;
      a_im = 0.0F;
    } else if (f_T11_re == 0.0F) {
      a_re = 0.0F;
      a_im = b_T22_im / T11.re;
    } else {
      a_re = f_T11_re / T11.re;
      a_im = b_T22_im / T11.re;
    }
  } else if (T11.re == 0.0F) {
    if (f_T11_re == 0.0F) {
      a_re = b_T22_im / T11.im;
      a_im = 0.0F;
    } else if (b_T22_im == 0.0F) {
      a_re = 0.0F;
      a_im = -(f_T11_re / T11.im);
    } else {
      a_re = b_T22_im / T11.im;
      a_im = -(f_T11_re / T11.im);
    }
  } else {
    c_T12_re = (float)fabs(T11.re);
    V1_norm = (float)fabs(T11.im);
    if (c_T12_re > V1_norm) {
      s = T11.im / T11.re;
      V1_norm = T11.re + s * T11.im;
      a_re = (f_T11_re + s * b_T22_im) / V1_norm;
      a_im = (b_T22_im - s * f_T11_re) / V1_norm;
    } else if (V1_norm == c_T12_re) {
      if (T11.re > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (T11.im > 0.0F) {
        V1_norm = 0.5F;
      } else {
        V1_norm = -0.5F;
      }

      a_re = (f_T11_re * s + b_T22_im * V1_norm) / c_T12_re;
      a_im = (b_T22_im * s - f_T11_re * V1_norm) / c_T12_re;
    } else {
      s = T11.re / T11.im;
      V1_norm = T11.im + s * T11.re;
      a_re = (s * f_T11_re + b_T22_im) / V1_norm;
      a_im = (s * b_T22_im - f_T11_re) / V1_norm;
    }
  }

  re = 1.73205078F * (fc0.re - a_re);
  im = 1.73205078F * (fc0.im - a_im);
  b_re = re * 0.0F - im;
  im = re + im * 0.0F;
  if (im == 0.0F) {
    re = b_re / 2.0F;
    im = 0.0F;
  } else if (b_re == 0.0F) {
    re = 0.0F;
    im /= 2.0F;
  } else {
    re = b_re / 2.0F;
    im /= 2.0F;
  }

  a_re = g_a_re * g_a_re;
  a_im = g_a_re * 0.0F + 0.0F * g_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  if (T13_im == 0.0F) {
    T13_re /= 3.0F;
    T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  rd_a.re = (((((a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  rd_a.im = (((((a_im + T12_im) + T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  fc0 = mpower(rd_a);
  sd_a.re = (f_a_re * f_a_re - c_a_im * c_a_im) - fc0.re;
  sd_a.im = (f_a_re * c_a_im + c_a_im * f_a_re) - fc0.im;
  fc0 = b_mpower(sd_a);
  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  b_T13_re = T13_re * H22 - T13_im * 0.0F;
  T13_im = T13_re * 0.0F + T13_im * H22;
  if (T13_im == 0.0F) {
    T13_re = b_T13_re / 2.0F;
    T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 2.0F;
  } else {
    T13_re = b_T13_re / 2.0F;
    T13_im /= 2.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  b_T12_re = T12_re * H33 - T12_im * 0.0F;
  T12_im = T12_re * 0.0F + T12_im * H33;
  if (T12_im == 0.0F) {
    T12_re = b_T12_re / 2.0F;
    T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 2.0F;
  } else {
    T12_re = b_T12_re / 2.0F;
    T12_im /= 2.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  b_T11_re = H11 * T23_re - 0.0F * T23_im;
  T11_im = H11 * T23_im + 0.0F * T23_re;
  if (T11_im == 0.0F) {
    b_T11_re /= 2.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 2.0F;
  } else {
    b_T11_re /= 2.0F;
    T11_im /= 2.0F;
  }

  xd_T11.re = (H11 + H22) + H33;
  xd_T11.im = 0.0F;
  T11 = mpower(xd_T11);
  if (T11.im == 0.0F) {
    b_re = T11.re / 27.0F;
    b_im = 0.0F;
  } else if (T11.re == 0.0F) {
    b_re = 0.0F;
    b_im = T11.im / 27.0F;
  } else {
    b_re = T11.re / 27.0F;
    b_im = T11.im / 27.0F;
  }

  c_T11_re = (H11 + H22) + H33;
  b_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  b_T12_im = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  d_T11_re = c_T11_re * b_T12_re - 0.0F * b_T12_im;
  b_T11_im = c_T11_re * b_T12_im + 0.0F * b_T12_re;
  if (b_T11_im == 0.0F) {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 6.0F;
  } else {
    c_T11_re = d_T11_re / 6.0F;
    b_T11_im /= 6.0F;
  }

  b_T12_re = H12 * H13;
  b_T12_im = H12 * 0.0F + 0.0F * H13;
  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  e_T11_re = d_T11_re * H33 - c_T11_im * 0.0F;
  c_T11_im = d_T11_re * 0.0F + c_T11_im * H33;
  if (c_T11_im == 0.0F) {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im = 0.0F;
  } else if (e_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re = e_T11_re / 2.0F;
    c_T11_im /= 2.0F;
  }

  fc49.re = ((((((fc0.re - T13_re) - T12_re) - b_T11_re) + b_re) + c_T11_re) +
             (b_T12_re * H23 - b_T12_im * 0.0F)) + d_T11_re;
  fc49.im = ((((((fc0.im - T13_im) - T12_im) - T11_im) + b_im) + b_T11_im) +
             (b_T12_re * 0.0F + b_T12_im * H23)) + c_T11_im;
  fc0 = c_mpower(fc49);
  if (fc0.im == 0.0F) {
    b_re = fc0.re / 2.0F;
    b_im = 0.0F;
  } else if (fc0.re == 0.0F) {
    b_re = 0.0F;
    b_im = fc0.im / 2.0F;
  } else {
    b_re = fc0.re / 2.0F;
    b_im = fc0.im / 2.0F;
  }

  a_re = h_a_re * h_a_re;
  a_im = h_a_re * 0.0F + 0.0F * h_a_re;
  if (a_im == 0.0F) {
    a_re /= 9.0F;
    a_im = 0.0F;
  } else if (a_re == 0.0F) {
    a_re = 0.0F;
    a_im /= 9.0F;
  } else {
    a_re /= 9.0F;
    a_im /= 9.0F;
  }

  T12_re = H12 * H12;
  T12_im = H12 * 0.0F + 0.0F * H12;
  if (T12_im == 0.0F) {
    T12_re /= 3.0F;
    T12_im = 0.0F;
  } else if (T12_re == 0.0F) {
    T12_re = 0.0F;
    T12_im /= 3.0F;
  } else {
    T12_re /= 3.0F;
    T12_im /= 3.0F;
  }

  T13_re = H13 * H13;
  T13_im = H13 * 0.0F + 0.0F * H13;
  if (T13_im == 0.0F) {
    T13_re /= 3.0F;
    T13_im = 0.0F;
  } else if (T13_re == 0.0F) {
    T13_re = 0.0F;
    T13_im /= 3.0F;
  } else {
    T13_re /= 3.0F;
    T13_im /= 3.0F;
  }

  T23_re = H23 * H23;
  T23_im = H23 * 0.0F + 0.0F * H23;
  if (T23_im == 0.0F) {
    T23_re /= 3.0F;
    T23_im = 0.0F;
  } else if (T23_re == 0.0F) {
    T23_re = 0.0F;
    T23_im /= 3.0F;
  } else {
    T23_re /= 3.0F;
    T23_im /= 3.0F;
  }

  b_T11_re = H11 * H22;
  T11_im = H11 * 0.0F + 0.0F * H22;
  if (T11_im == 0.0F) {
    b_T11_re /= 3.0F;
    T11_im = 0.0F;
  } else if (b_T11_re == 0.0F) {
    b_T11_re = 0.0F;
    T11_im /= 3.0F;
  } else {
    b_T11_re /= 3.0F;
    T11_im /= 3.0F;
  }

  c_T11_re = H11 * H33;
  b_T11_im = H11 * 0.0F + 0.0F * H33;
  if (b_T11_im == 0.0F) {
    c_T11_re /= 3.0F;
    b_T11_im = 0.0F;
  } else if (c_T11_re == 0.0F) {
    c_T11_re = 0.0F;
    b_T11_im /= 3.0F;
  } else {
    c_T11_re /= 3.0F;
    b_T11_im /= 3.0F;
  }

  b_T22_re = H22 * H33;
  T22_im = H22 * 0.0F + 0.0F * H33;
  if (T22_im == 0.0F) {
    b_T22_re /= 3.0F;
    T22_im = 0.0F;
  } else if (b_T22_re == 0.0F) {
    b_T22_re = 0.0F;
    T22_im /= 3.0F;
  } else {
    b_T22_re /= 3.0F;
    T22_im /= 3.0F;
  }

  b_a_re = j_a_re * j_a_re;
  b_a_im = j_a_re * 0.0F + 0.0F * j_a_re;
  if (b_a_im == 0.0F) {
    b_a_re /= 9.0F;
    b_a_im = 0.0F;
  } else if (b_a_re == 0.0F) {
    b_a_re = 0.0F;
    b_a_im /= 9.0F;
  } else {
    b_a_re /= 9.0F;
    b_a_im /= 9.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  if (b_T12_im == 0.0F) {
    b_T12_re /= 3.0F;
    b_T12_im = 0.0F;
  } else if (b_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 3.0F;
  } else {
    b_T12_re /= 3.0F;
    b_T12_im /= 3.0F;
  }

  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  if (b_T13_im == 0.0F) {
    b_T13_re /= 3.0F;
    b_T13_im = 0.0F;
  } else if (b_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 3.0F;
  } else {
    b_T13_re /= 3.0F;
    b_T13_im /= 3.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  if (b_T23_im == 0.0F) {
    b_T23_re /= 3.0F;
    b_T23_im = 0.0F;
  } else if (b_T23_re == 0.0F) {
    b_T23_re = 0.0F;
    b_T23_im /= 3.0F;
  } else {
    b_T23_re /= 3.0F;
    b_T23_im /= 3.0F;
  }

  d_T11_re = H11 * H22;
  c_T11_im = H11 * 0.0F + 0.0F * H22;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 3.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 3.0F;
  } else {
    d_T11_re /= 3.0F;
    c_T11_im /= 3.0F;
  }

  e_T11_re = H11 * H33;
  c_a_re = H11 * 0.0F + 0.0F * H33;
  if (c_a_re == 0.0F) {
    e_T11_re /= 3.0F;
    c_a_re = 0.0F;
  } else if (e_T11_re == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 3.0F;
  } else {
    e_T11_re /= 3.0F;
    c_a_re /= 3.0F;
  }

  c_T22_re = H22 * H33;
  b_T22_im = H22 * 0.0F + 0.0F * H33;
  if (b_T22_im == 0.0F) {
    c_T22_re /= 3.0F;
    b_T22_im = 0.0F;
  } else if (c_T22_re == 0.0F) {
    c_T22_re = 0.0F;
    b_T22_im /= 3.0F;
  } else {
    c_T22_re /= 3.0F;
    b_T22_im /= 3.0F;
  }

  td_a.re = (((((b_a_re + b_T12_re) + b_T13_re) + b_T23_re) - d_T11_re) -
             e_T11_re) - c_T22_re;
  td_a.im = (((((b_a_im + b_T12_im) + b_T13_im) + b_T23_im) - c_T11_im) - c_a_re)
    - b_T22_im;
  fc0 = mpower(td_a);
  ud_a.re = (i_a_re * i_a_re - d_a_im * d_a_im) - fc0.re;
  ud_a.im = (i_a_re * d_a_im + d_a_im * i_a_re) - fc0.im;
  fc0 = b_mpower(ud_a);
  b_T13_re = H13 * H13;
  b_T13_im = H13 * 0.0F + 0.0F * H13;
  c_T13_re = b_T13_re * H22 - b_T13_im * 0.0F;
  b_T13_im = b_T13_re * 0.0F + b_T13_im * H22;
  if (b_T13_im == 0.0F) {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im = 0.0F;
  } else if (c_T13_re == 0.0F) {
    b_T13_re = 0.0F;
    b_T13_im /= 2.0F;
  } else {
    b_T13_re = c_T13_re / 2.0F;
    b_T13_im /= 2.0F;
  }

  b_T12_re = H12 * H12;
  b_T12_im = H12 * 0.0F + 0.0F * H12;
  c_T12_re = b_T12_re * H33 - b_T12_im * 0.0F;
  b_T12_im = b_T12_re * 0.0F + b_T12_im * H33;
  if (b_T12_im == 0.0F) {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im = 0.0F;
  } else if (c_T12_re == 0.0F) {
    b_T12_re = 0.0F;
    b_T12_im /= 2.0F;
  } else {
    b_T12_re = c_T12_re / 2.0F;
    b_T12_im /= 2.0F;
  }

  b_T23_re = H23 * H23;
  b_T23_im = H23 * 0.0F + 0.0F * H23;
  d_T11_re = H11 * b_T23_re - 0.0F * b_T23_im;
  c_T11_im = H11 * b_T23_im + 0.0F * b_T23_re;
  if (c_T11_im == 0.0F) {
    d_T11_re /= 2.0F;
    c_T11_im = 0.0F;
  } else if (d_T11_re == 0.0F) {
    d_T11_re = 0.0F;
    c_T11_im /= 2.0F;
  } else {
    d_T11_re /= 2.0F;
    c_T11_im /= 2.0F;
  }

  yd_T11.re = (H11 + H22) + H33;
  yd_T11.im = 0.0F;
  T11 = mpower(yd_T11);
  if (T11.im == 0.0F) {
    b_a_re = T11.re / 27.0F;
    d_a_re = 0.0F;
  } else if (T11.re == 0.0F) {
    b_a_re = 0.0F;
    d_a_re = T11.im / 27.0F;
  } else {
    b_a_re = T11.re / 27.0F;
    d_a_re = T11.im / 27.0F;
  }

  e_T11_re = (H11 + H22) + H33;
  c_T12_re = ((((H12 * H12 + H13 * H13) + H23 * H23) - H11 * H22) - H11 * H33) -
    H22 * H33;
  c_T22_re = (((((H12 * 0.0F + 0.0F * H12) + (H13 * 0.0F + 0.0F * H13)) + (H23 *
    0.0F + 0.0F * H23)) - (H11 * 0.0F + 0.0F * H22)) - (H11 * 0.0F + 0.0F * H33))
    - (H22 * 0.0F + 0.0F * H33);
  b_T22_im = e_T11_re * c_T12_re - 0.0F * c_T22_re;
  c_a_re = e_T11_re * c_T22_re + 0.0F * c_T12_re;
  if (c_a_re == 0.0F) {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re = 0.0F;
  } else if (b_T22_im == 0.0F) {
    e_T11_re = 0.0F;
    c_a_re /= 6.0F;
  } else {
    e_T11_re = b_T22_im / 6.0F;
    c_a_re /= 6.0F;
  }

  c_T12_re = H12 * H13;
  c_T22_re = H12 * 0.0F + 0.0F * H13;
  b_T22_im = H11 * H22;
  d_T11_im = H11 * 0.0F + 0.0F * H22;
  f_T11_re = b_T22_im * H33 - d_T11_im * 0.0F;
  d_T11_im = b_T22_im * 0.0F + d_T11_im * H33;
  if (d_T11_im == 0.0F) {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im = 0.0F;
  } else if (f_T11_re == 0.0F) {
    b_T22_im = 0.0F;
    d_T11_im /= 2.0F;
  } else {
    b_T22_im = f_T11_re / 2.0F;
    d_T11_im /= 2.0F;
  }

  fc50.re = ((((((fc0.re - b_T13_re) - b_T12_re) - d_T11_re) + b_a_re) +
              e_T11_re) + (c_T12_re * H23 - c_T22_re * 0.0F)) + b_T22_im;
  fc50.im = ((((((fc0.im - b_T13_im) - b_T12_im) - c_T11_im) + d_a_re) + c_a_re)
             + (c_T12_re * 0.0F + c_T22_re * H23)) + d_T11_im;
  fc0 = c_mpower(fc50);
  f_T11_re = (((((a_re + T12_re) + T13_re) + T23_re) - b_T11_re) - c_T11_re) -
    b_T22_re;
  b_T22_im = (((((a_im + T12_im) + T13_im) + T23_im) - T11_im) - b_T11_im) -
    T22_im;
  d_T11_im = 2.0F * fc0.re;
  c_T22_re = 2.0F * fc0.im;
  if (c_T22_re == 0.0F) {
    if (b_T22_im == 0.0F) {
      a_re = f_T11_re / d_T11_im;
      a_im = 0.0F;
    } else if (f_T11_re == 0.0F) {
      a_re = 0.0F;
      a_im = b_T22_im / d_T11_im;
    } else {
      a_re = f_T11_re / d_T11_im;
      a_im = b_T22_im / d_T11_im;
    }
  } else if (d_T11_im == 0.0F) {
    if (f_T11_re == 0.0F) {
      a_re = b_T22_im / c_T22_re;
      a_im = 0.0F;
    } else if (b_T22_im == 0.0F) {
      a_re = 0.0F;
      a_im = -(f_T11_re / c_T22_re);
    } else {
      a_re = b_T22_im / c_T22_re;
      a_im = -(f_T11_re / c_T22_re);
    }
  } else {
    c_T12_re = (float)fabs(d_T11_im);
    V1_norm = (float)fabs(c_T22_re);
    if (c_T12_re > V1_norm) {
      s = c_T22_re / d_T11_im;
      V1_norm = d_T11_im + s * c_T22_re;
      a_re = (f_T11_re + s * b_T22_im) / V1_norm;
      a_im = (b_T22_im - s * f_T11_re) / V1_norm;
    } else if (V1_norm == c_T12_re) {
      if (d_T11_im > 0.0F) {
        s = 0.5F;
      } else {
        s = -0.5F;
      }

      if (c_T22_re > 0.0F) {
        V1_norm = 0.5F;
      } else {
        V1_norm = -0.5F;
      }

      a_re = (f_T11_re * s + b_T22_im * V1_norm) / c_T12_re;
      a_im = (b_T22_im * s - f_T11_re * V1_norm) / c_T12_re;
    } else {
      s = d_T11_im / c_T22_re;
      V1_norm = c_T22_re + s * d_T11_im;
      a_re = (s * f_T11_re + b_T22_im) / V1_norm;
      a_im = (s * b_T22_im - f_T11_re) / V1_norm;
    }
  }

  V23_re = ((((T11_re + T22_re) + T33_re) + re) - b_re) - a_re;
  c_T12_re = (im - b_im) - a_im;

  /*  Returning Eigen-values sorted ASCENDING! |El1|<=|El2|<=|El3| */
  if ((rt_hypotf_snf(V12_re, V12_im) <= rt_hypotf_snf(V22_re, V22_im)) &&
      (rt_hypotf_snf(V22_re, V22_im) <= rt_hypotf_snf(V23_re, c_T12_re))) {
    /*  Returning Eigen-values */
    V1_norm = V12_re;
    b_sign(&V1_norm);
    *El1 = V1_norm * rt_hypotf_snf(V12_re, V12_im);
    V1_norm = V22_re;
    b_sign(&V1_norm);
    *El2 = V1_norm * rt_hypotf_snf(V22_re, V22_im);
    V1_norm = V23_re;
    b_sign(&V1_norm);
    *El3 = V1_norm * rt_hypotf_snf(V23_re, c_T12_re);

    /*  Returning Eigen-vectors */
    V1_norm = v11_re;
    b_sign(&V1_norm);
    *Ev11 = V1_norm * rt_hypotf_snf(v11_re, v11_im);
    V1_norm = v12_re;
    b_sign(&V1_norm);
    *Ev12 = V1_norm * rt_hypotf_snf(v12_re, v12_im);
    V1_norm = v13_re;
    b_sign(&V1_norm);
    *Ev13 = V1_norm * rt_hypotf_snf(v13_re, 0.0F);
    V1_norm = v21_re;
    b_sign(&V1_norm);
    *Ev21 = V1_norm * rt_hypotf_snf(v21_re, V21_im);
    V1_norm = v22_re;
    b_sign(&V1_norm);
    *Ev22 = V1_norm * rt_hypotf_snf(v22_re, v22_im);
    V1_norm = v23_re;
    b_sign(&V1_norm);
    *Ev23 = V1_norm * rt_hypotf_snf(v23_re, 0.0F);
    V1_norm = V13_re;
    b_sign(&V1_norm);
    *Ev31 = V1_norm * rt_hypotf_snf(V13_re, V13_im);
    V1_norm = V11_re;
    b_sign(&V1_norm);
    *Ev32 = V1_norm * rt_hypotf_snf(V11_re, V11_im);
    V1_norm = V21_re;
    b_sign(&V1_norm);
    *Ev33 = V1_norm * rt_hypotf_snf(V21_re, 0.0F);
  } else if ((rt_hypotf_snf(V12_re, V12_im) <= rt_hypotf_snf(V23_re, c_T12_re)) &&
             (rt_hypotf_snf(V23_re, c_T12_re) <= rt_hypotf_snf(V22_re, V22_im)))
  {
    /*  Returning Eigen-values */
    V1_norm = V12_re;
    b_sign(&V1_norm);
    *El1 = V1_norm * rt_hypotf_snf(V12_re, V12_im);
    V1_norm = V23_re;
    b_sign(&V1_norm);
    *El2 = V1_norm * rt_hypotf_snf(V23_re, c_T12_re);
    V1_norm = V22_re;
    b_sign(&V1_norm);
    *El3 = V1_norm * rt_hypotf_snf(V22_re, V22_im);

    /*  Returning Eigen-vectors */
    V1_norm = v11_re;
    b_sign(&V1_norm);
    *Ev11 = V1_norm * rt_hypotf_snf(v11_re, v11_im);
    V1_norm = v12_re;
    b_sign(&V1_norm);
    *Ev12 = V1_norm * rt_hypotf_snf(v12_re, v12_im);
    V1_norm = v13_re;
    b_sign(&V1_norm);
    *Ev13 = V1_norm * rt_hypotf_snf(v13_re, 0.0F);
    V1_norm = V13_re;
    b_sign(&V1_norm);
    *Ev21 = V1_norm * rt_hypotf_snf(V13_re, V13_im);
    V1_norm = V11_re;
    b_sign(&V1_norm);
    *Ev22 = V1_norm * rt_hypotf_snf(V11_re, V11_im);
    V1_norm = V21_re;
    b_sign(&V1_norm);
    *Ev23 = V1_norm * rt_hypotf_snf(V21_re, 0.0F);
    V1_norm = v21_re;
    b_sign(&V1_norm);
    *Ev31 = V1_norm * rt_hypotf_snf(v21_re, V21_im);
    V1_norm = v22_re;
    b_sign(&V1_norm);
    *Ev32 = V1_norm * rt_hypotf_snf(v22_re, v22_im);
    V1_norm = v23_re;
    b_sign(&V1_norm);
    *Ev33 = V1_norm * rt_hypotf_snf(v23_re, 0.0F);
  } else if ((rt_hypotf_snf(V22_re, V22_im) <= rt_hypotf_snf(V12_re, V12_im)) &&
             (rt_hypotf_snf(V12_re, V12_im) <= rt_hypotf_snf(V23_re, c_T12_re)))
  {
    /*  Returning Eigen-values */
    V1_norm = V22_re;
    b_sign(&V1_norm);
    *El1 = V1_norm * rt_hypotf_snf(V22_re, V22_im);
    V1_norm = V12_re;
    b_sign(&V1_norm);
    *El2 = V1_norm * rt_hypotf_snf(V12_re, V12_im);
    V1_norm = V23_re;
    b_sign(&V1_norm);
    *El3 = V1_norm * rt_hypotf_snf(V23_re, c_T12_re);

    /*  Returning Eigen-vectors */
    V1_norm = v21_re;
    b_sign(&V1_norm);
    *Ev11 = V1_norm * rt_hypotf_snf(v21_re, V21_im);
    V1_norm = v22_re;
    b_sign(&V1_norm);
    *Ev12 = V1_norm * rt_hypotf_snf(v22_re, v22_im);
    V1_norm = v23_re;
    b_sign(&V1_norm);
    *Ev13 = V1_norm * rt_hypotf_snf(v23_re, 0.0F);
    V1_norm = v11_re;
    b_sign(&V1_norm);
    *Ev21 = V1_norm * rt_hypotf_snf(v11_re, v11_im);
    V1_norm = v12_re;
    b_sign(&V1_norm);
    *Ev22 = V1_norm * rt_hypotf_snf(v12_re, v12_im);
    V1_norm = v13_re;
    b_sign(&V1_norm);
    *Ev23 = V1_norm * rt_hypotf_snf(v13_re, 0.0F);
    V1_norm = V13_re;
    b_sign(&V1_norm);
    *Ev31 = V1_norm * rt_hypotf_snf(V13_re, V13_im);
    V1_norm = V11_re;
    b_sign(&V1_norm);
    *Ev32 = V1_norm * rt_hypotf_snf(V11_re, V11_im);
    V1_norm = V21_re;
    b_sign(&V1_norm);
    *Ev33 = V1_norm * rt_hypotf_snf(V21_re, 0.0F);
  } else if ((rt_hypotf_snf(V22_re, V22_im) <= rt_hypotf_snf(V23_re, c_T12_re)) &&
             (rt_hypotf_snf(V23_re, c_T12_re) <= rt_hypotf_snf(V12_re, V12_im)))
  {
    /*  Returning Eigen-values */
    V1_norm = V22_re;
    b_sign(&V1_norm);
    *El1 = V1_norm * rt_hypotf_snf(V22_re, V22_im);
    V1_norm = V23_re;
    b_sign(&V1_norm);
    *El2 = V1_norm * rt_hypotf_snf(V23_re, c_T12_re);
    V1_norm = V12_re;
    b_sign(&V1_norm);
    *El3 = V1_norm * rt_hypotf_snf(V12_re, V12_im);

    /*  Returning Eigen-vectors */
    V1_norm = v21_re;
    b_sign(&V1_norm);
    *Ev11 = V1_norm * rt_hypotf_snf(v21_re, V21_im);
    V1_norm = v22_re;
    b_sign(&V1_norm);
    *Ev12 = V1_norm * rt_hypotf_snf(v22_re, v22_im);
    V1_norm = v23_re;
    b_sign(&V1_norm);
    *Ev13 = V1_norm * rt_hypotf_snf(v23_re, 0.0F);
    V1_norm = V13_re;
    b_sign(&V1_norm);
    *Ev21 = V1_norm * rt_hypotf_snf(V13_re, V13_im);
    V1_norm = V11_re;
    b_sign(&V1_norm);
    *Ev22 = V1_norm * rt_hypotf_snf(V11_re, V11_im);
    V1_norm = V21_re;
    b_sign(&V1_norm);
    *Ev23 = V1_norm * rt_hypotf_snf(V21_re, 0.0F);
    V1_norm = v11_re;
    b_sign(&V1_norm);
    *Ev31 = V1_norm * rt_hypotf_snf(v11_re, v11_im);
    V1_norm = v12_re;
    b_sign(&V1_norm);
    *Ev32 = V1_norm * rt_hypotf_snf(v12_re, v12_im);
    V1_norm = v13_re;
    b_sign(&V1_norm);
    *Ev33 = V1_norm * rt_hypotf_snf(v13_re, 0.0F);
  } else if ((rt_hypotf_snf(V23_re, c_T12_re) <= rt_hypotf_snf(V22_re, V22_im)) &&
             (rt_hypotf_snf(V22_re, V22_im) <= rt_hypotf_snf(V12_re, V12_im))) {
    /*  Returning Eigen-values */
    V1_norm = V23_re;
    b_sign(&V1_norm);
    *El1 = V1_norm * rt_hypotf_snf(V23_re, c_T12_re);
    V1_norm = V22_re;
    b_sign(&V1_norm);
    *El2 = V1_norm * rt_hypotf_snf(V22_re, V22_im);
    V1_norm = V12_re;
    b_sign(&V1_norm);
    *El3 = V1_norm * rt_hypotf_snf(V12_re, V12_im);

    /*  Returning Eigen-vectors */
    V1_norm = V13_re;
    b_sign(&V1_norm);
    *Ev11 = V1_norm * rt_hypotf_snf(V13_re, V13_im);
    V1_norm = V11_re;
    b_sign(&V1_norm);
    *Ev12 = V1_norm * rt_hypotf_snf(V11_re, V11_im);
    V1_norm = V21_re;
    b_sign(&V1_norm);
    *Ev13 = V1_norm * rt_hypotf_snf(V21_re, 0.0F);
    if (v21_re < 0.0F) {
      b_v21_re = -1.0F;
    } else if (v21_re > 0.0F) {
      b_v21_re = 1.0F;
    } else if (v21_re == 0.0F) {
      b_v21_re = 0.0F;
    } else {
      b_v21_re = v21_re;
    }

    *Ev21 = b_v21_re * rt_hypotf_snf(v21_re, V21_im);
    if (v22_re < 0.0F) {
      b_v22_re = -1.0F;
    } else if (v22_re > 0.0F) {
      b_v22_re = 1.0F;
    } else if (v22_re == 0.0F) {
      b_v22_re = 0.0F;
    } else {
      b_v22_re = v22_re;
    }

    *Ev22 = b_v22_re * rt_hypotf_snf(v22_re, v22_im);
    if (v23_re < 0.0F) {
      b_v23_re = -1.0F;
    } else if (v23_re > 0.0F) {
      b_v23_re = 1.0F;
    } else if (v23_re == 0.0F) {
      b_v23_re = 0.0F;
    } else {
      b_v23_re = v23_re;
    }

    *Ev23 = b_v23_re * rt_hypotf_snf(v23_re, 0.0F);
    if (v11_re < 0.0F) {
      c_v11_re = -1.0F;
    } else if (v11_re > 0.0F) {
      c_v11_re = 1.0F;
    } else if (v11_re == 0.0F) {
      c_v11_re = 0.0F;
    } else {
      c_v11_re = v11_re;
    }

    *Ev31 = c_v11_re * rt_hypotf_snf(v11_re, v11_im);
    if (v12_re < 0.0F) {
      c_v12_re = -1.0F;
    } else if (v12_re > 0.0F) {
      c_v12_re = 1.0F;
    } else if (v12_re == 0.0F) {
      c_v12_re = 0.0F;
    } else {
      c_v12_re = v12_re;
    }

    *Ev32 = c_v12_re * rt_hypotf_snf(v12_re, v12_im);
    if (v13_re < 0.0F) {
      c_v13_re = -1.0F;
    } else if (v13_re > 0.0F) {
      c_v13_re = 1.0F;
    } else if (v13_re == 0.0F) {
      c_v13_re = 0.0F;
    } else {
      c_v13_re = v13_re;
    }

    *Ev33 = c_v13_re * rt_hypotf_snf(v13_re, 0.0F);
  } else {
    /*  abs(l3) <= abs(l1) && abs(l1)<= abs(l2) */
    /*  Returning Eigen-values */
    if (V23_re < 0.0F) {
      b_V23_re = -1.0F;
    } else if (V23_re > 0.0F) {
      b_V23_re = 1.0F;
    } else if (V23_re == 0.0F) {
      b_V23_re = 0.0F;
    } else {
      b_V23_re = V23_re;
    }

    *El1 = b_V23_re * rt_hypotf_snf(V23_re, c_T12_re);
    if (V12_re < 0.0F) {
      b_V12_re = -1.0F;
    } else if (V12_re > 0.0F) {
      b_V12_re = 1.0F;
    } else if (V12_re == 0.0F) {
      b_V12_re = 0.0F;
    } else {
      b_V12_re = V12_re;
    }

    *El2 = b_V12_re * rt_hypotf_snf(V12_re, V12_im);
    if (V22_re < 0.0F) {
      b_V22_re = -1.0F;
    } else if (V22_re > 0.0F) {
      b_V22_re = 1.0F;
    } else if (V22_re == 0.0F) {
      b_V22_re = 0.0F;
    } else {
      b_V22_re = V22_re;
    }

    *El3 = b_V22_re * rt_hypotf_snf(V22_re, V22_im);

    /*  Returning Eigen-vectors */
    if (V13_re < 0.0F) {
      b_V13_re = -1.0F;
    } else if (V13_re > 0.0F) {
      b_V13_re = 1.0F;
    } else if (V13_re == 0.0F) {
      b_V13_re = 0.0F;
    } else {
      b_V13_re = V13_re;
    }

    *Ev11 = b_V13_re * rt_hypotf_snf(V13_re, V13_im);
    if (V11_re < 0.0F) {
      b_V11_re = -1.0F;
    } else if (V11_re > 0.0F) {
      b_V11_re = 1.0F;
    } else if (V11_re == 0.0F) {
      b_V11_re = 0.0F;
    } else {
      b_V11_re = V11_re;
    }

    *Ev12 = b_V11_re * rt_hypotf_snf(V11_re, V11_im);
    if (V21_re < 0.0F) {
      b_V21_re = -1.0F;
    } else if (V21_re > 0.0F) {
      b_V21_re = 1.0F;
    } else if (V21_re == 0.0F) {
      b_V21_re = 0.0F;
    } else {
      b_V21_re = V21_re;
    }

    *Ev13 = b_V21_re * rt_hypotf_snf(V21_re, 0.0F);
    if (v11_re < 0.0F) {
      b_v11_re = -1.0F;
    } else if (v11_re > 0.0F) {
      b_v11_re = 1.0F;
    } else if (v11_re == 0.0F) {
      b_v11_re = 0.0F;
    } else {
      b_v11_re = v11_re;
    }

    *Ev21 = b_v11_re * rt_hypotf_snf(v11_re, v11_im);
    if (v12_re < 0.0F) {
      b_v12_re = -1.0F;
    } else if (v12_re > 0.0F) {
      b_v12_re = 1.0F;
    } else if (v12_re == 0.0F) {
      b_v12_re = 0.0F;
    } else {
      b_v12_re = v12_re;
    }

    *Ev22 = b_v12_re * rt_hypotf_snf(v12_re, v12_im);
    if (v13_re < 0.0F) {
      b_v13_re = -1.0F;
    } else if (v13_re > 0.0F) {
      b_v13_re = 1.0F;
    } else if (v13_re == 0.0F) {
      b_v13_re = 0.0F;
    } else {
      b_v13_re = v13_re;
    }

    *Ev23 = b_v13_re * rt_hypotf_snf(v13_re, 0.0F);
    if (v21_re < 0.0F) {
      c_v21_re = -1.0F;
    } else if (v21_re > 0.0F) {
      c_v21_re = 1.0F;
    } else if (v21_re == 0.0F) {
      c_v21_re = 0.0F;
    } else {
      c_v21_re = v21_re;
    }

    *Ev31 = c_v21_re * rt_hypotf_snf(v21_re, V21_im);
    if (v22_re < 0.0F) {
      c_v22_re = -1.0F;
    } else if (v22_re > 0.0F) {
      c_v22_re = 1.0F;
    } else if (v22_re == 0.0F) {
      c_v22_re = 0.0F;
    } else {
      c_v22_re = v22_re;
    }

    *Ev32 = c_v22_re * rt_hypotf_snf(v22_re, v22_im);
    if (v23_re < 0.0F) {
      c_v23_re = -1.0F;
    } else if (v23_re > 0.0F) {
      c_v23_re = 1.0F;
    } else if (v23_re == 0.0F) {
      c_v23_re = 0.0F;
    } else {
      c_v23_re = v23_re;
    }

    *Ev33 = c_v23_re * rt_hypotf_snf(v23_re, 0.0F);
  }

  sortElvs( El1, El2, El3,
            Ev11,Ev12,Ev13,
            Ev21,Ev22,Ev23,
            Ev31,Ev32,Ev33);

}

