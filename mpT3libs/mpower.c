/* Include Files */
#include "rt_nonfinite.h"
#include "mpower.h"
#include "T3_eig.h"
#include "T3_eig_rtwutil.h"

/* Function Declarations */
static float rt_atan2f_snf(float u0, float u1);
static float rt_powf_snf(float u0, float u1);

/* Function Definitions */
static float rt_atan2f_snf(float u0, float u1)
{
  float y;
  int b_u0;
  int b_u1;
  if (rtIsNaNF(u0) || rtIsNaNF(u1)) {
    y = ((real32_T)rtNaN);
  } else if (rtIsInfF(u0) && rtIsInfF(u1)) {
    if (u0 > 0.0F) {
      b_u0 = 1;
    } else {
      b_u0 = -1;
    }

    if (u1 > 0.0F) {
      b_u1 = 1;
    } else {
      b_u1 = -1;
    }

    y = (float)atan2((float)b_u0, (float)b_u1);
  } else if (u1 == 0.0F) {
    if (u0 > 0.0F) {
      y = RT_PIF / 2.0F;
    } else if (u0 < 0.0F) {
      y = -(RT_PIF / 2.0F);
    } else {
      y = 0.0F;
    }
  } else {
    y = (float)atan2(u0, u1);
  }

  return y;
}

static float rt_powf_snf(float u0, float u1)
{
  float y;
  float f0;
  float f1;
  if (rtIsNaNF(u0) || rtIsNaNF(u1)) {
    y = ((real32_T)rtNaN);
  } else {
    f0 = (float)fabs(u0);
    f1 = (float)fabs(u1);
    if (rtIsInfF(u1)) {
      if (f0 == 1.0F) {
        y = ((real32_T)rtNaN);
      } else if (f0 > 1.0F) {
        if (u1 > 0.0F) {
          y = ((real32_T)rtInf);
        } else {
          y = 0.0F;
        }
      } else if (u1 > 0.0F) {
        y = 0.0F;
      } else {
        y = ((real32_T)rtInf);
      }
    } else if (f1 == 0.0F) {
      y = 1.0F;
    } else if (f1 == 1.0F) {
      if (u1 > 0.0F) {
        y = u0;
      } else {
        y = 1.0F / u0;
      }
    } else if (u1 == 2.0F) {
      y = u0 * u0;
    } else if ((u1 == 0.5F) && (u0 >= 0.0F)) {
      y = (float)sqrt(u0);
    } else if ((u0 < 0.0F) && (u1 > (float)floor(u1))) {
      y = ((real32_T)rtNaN);
    } else {
      y = (float)pow(u0, u1);
    }
  }

  return y;
}

creal32_T b_mpower(const creal32_T a)
{
  creal32_T c;
  float absxi;
  float absxr;
  if (a.im == 0.0F) {
    if (a.re < 0.0F) {
      absxi = 0.0F;
      absxr = (float)sqrt((float)fabs(a.re));
    } else {
      absxi = (float)sqrt(a.re);
      absxr = 0.0F;
    }
  } else if (a.re == 0.0F) {
    if (a.im < 0.0F) {
      absxi = (float)sqrt(-a.im / 2.0F);
      absxr = -absxi;
    } else {
      absxi = (float)sqrt(a.im / 2.0F);
      absxr = absxi;
    }
  } else if (rtIsNaNF(a.re) || rtIsNaNF(a.im)) {
    absxi = ((real32_T)rtNaN);
    absxr = ((real32_T)rtNaN);
  } else if (rtIsInfF(a.im)) {
    absxi = ((real32_T)rtInf);
    absxr = a.im;
  } else if (rtIsInfF(a.re)) {
    if (a.re < 0.0F) {
      absxi = 0.0F;
      absxr = ((real32_T)rtInf);
    } else {
      absxi = ((real32_T)rtInf);
      absxr = 0.0F;
    }
  } else {
    absxr = (float)fabs(a.re);
    absxi = (float)fabs(a.im);
    if ((absxr > 8.50705867E+37F) || (absxi > 8.50705867E+37F)) {
      absxr *= 0.5F;
      absxi *= 0.5F;
      absxi = rt_hypotf_snf(absxr, absxi);
      if (absxi > absxr) {
        absxi = (float)sqrt(absxi) * (float)sqrt(1.0F + absxr / absxi);
      } else {
        absxi = (float)sqrt(absxi) * 1.41421354F;
      }
    } else {
      absxi = (float)sqrt((rt_hypotf_snf(absxr, absxi) + absxr) * 0.5F);
    }

    if (a.re > 0.0F) {
      absxr = 0.5F * (a.im / absxi);
    } else {
      if (a.im < 0.0F) {
        absxr = -absxi;
      } else {
        absxr = absxi;
      }

      absxi = 0.5F * (a.im / absxr);
    }
  }

  c.re = absxi;
  c.im = absxr;
  return c;
}

creal32_T c_mpower(const creal32_T a)
{
  creal32_T c;
  float r;
  float x_im;
  if ((a.im == 0.0F) && (a.re >= 0.0F)) {
    c.re = rt_powf_snf(a.re, 0.333333343F);
    c.im = 0.0F;
  } else {
    r = a.re;
    x_im = a.im;
    if ((a.im == 0.0F) && rtIsNaNF(a.re)) {
    } else if (((float)fabs(a.re) > 1.70141173E+38F) || ((float)fabs(a.im) >
                1.70141173E+38F)) {
      r = (float)log(rt_hypotf_snf(a.re / 2.0F, a.im / 2.0F)) + 0.693147182F;
      x_im = rt_atan2f_snf(a.im, a.re);
    } else {
      r = (float)log(rt_hypotf_snf(a.re, a.im));
      x_im = rt_atan2f_snf(a.im, a.re);
    }

    r *= 0.333333343F;
    x_im *= 0.333333343F;
    if (x_im == 0.0F) {
      c.re = (float)exp(r);
      c.im = 0.0F;
    } else if (rtIsInfF(x_im) && rtIsInfF(r) && (r < 0.0F)) {
      c.re = 0.0F;
      c.im = 0.0F;
    } else {
      r = (float)exp(r / 2.0F);
      c.re = r * (r * (float)cos(x_im));
      c.im = r * (r * (float)sin(x_im));
    }
  }

  return c;
}

/*
 * Arguments    : const creal32_T a
 * Return Type  : creal32_T
 */
creal32_T mpower(const creal32_T a)
{
  creal32_T c;
  float r;
  float x_im;
  if ((a.im == 0.0F) && (a.re >= 0.0F)) {
    c.re = rt_powf_snf(a.re, 3.0F);
    c.im = 0.0F;
  } else if (a.re == 0.0F) {
    c.re = 0.0F;
    c.im = -rt_powf_snf(a.im, 3.0F);
  } else {
    r = a.re;
    x_im = a.im;
    if ((a.im == 0.0F) && rtIsNaNF(a.re)) {
    } else if (((float)fabs(a.re) > 1.70141173E+38F) || ((float)fabs(a.im) >
                1.70141173E+38F)) {
      r = (float)log(rt_hypotf_snf(a.re / 2.0F, a.im / 2.0F)) + 0.693147182F;
      x_im = rt_atan2f_snf(a.im, a.re);
    } else {
      r = (float)log(rt_hypotf_snf(a.re, a.im));
      x_im = rt_atan2f_snf(a.im, a.re);
    }

    r *= 3.0F;
    x_im *= 3.0F;
    if (x_im == 0.0F) {
      c.re = (float)exp(r);
      c.im = 0.0F;
    } else if (rtIsInfF(x_im) && rtIsInfF(r) && (r < 0.0F)) {
      c.re = 0.0F;
      c.im = 0.0F;
    } else {
      r = (float)exp(r / 2.0F);
      c.re = r * (r * (float)cos(x_im));
      c.im = r * (r * (float)sin(x_im));
    }
  }

  return c;
}
