#ifndef T3LIC_TO_T3_H
#define T3LIC_TO_T3_H

/* Include Files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include "rt_defines.h"
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "T3LIC_to_T3_types.h"

/* Function Declarations */
extern void T3LIC_to_T3(float *T11_in, float *T12_in, float *T13_in,
						float *T22_in, float *T23_in, float *T33_in,
						float *El1, float *El2, float *El3,
						float *Ev11, float *Ev12, float *Ev13,
						float *Ev21, float *Ev22, float *Ev23,
						float *Ev31, float *Ev32, float *Ev33);

#endif
