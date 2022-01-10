#ifndef T3_EIG_H
#define T3_EIG_H

/* Include Files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include "rt_defines.h"
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "T3_eig_types.h"

/* Function Declarations */
extern void T3_eig(float *H11_in, float *H12_in, float *H13_in,
				   float *H22_in, float *H23_in, float *H33_in,
				   float *El1, float *El2, float *El3,
				   float *Ev11, float *Ev12, float *Ev13,
				   float *Ev21, float *Ev22, float *Ev23,
                   float *Ev31, float *Ev32, float *Ev33);

#endif
