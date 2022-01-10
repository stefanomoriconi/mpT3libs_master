#ifndef MPOWER_H
#define MPOWER_H

/* Include Files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include "rt_defines.h"
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "T3_eig_types.h"

/* Function Declarations */
extern creal32_T b_mpower(const creal32_T a);
extern creal32_T c_mpower(const creal32_T a);
extern creal32_T mpower(const creal32_T a);

#endif