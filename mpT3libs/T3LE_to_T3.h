#ifndef T3LE_TO_T3_H
#define T3LE_TO_T3_H

/* Include Files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include "rt_defines.h"
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "T3LE_to_T3_types.h"

/* Function Declarations */
extern void T3LE_to_T3(float *T11LE_in, float *T12LE_in, float *T13LE_in,
					   float *T22LE_in, float *T23LE_in, float *T33LE_in,
					   float *El1, float *El2, float *El3,
					   float *Ev11, float *Ev12, float *Ev13,
					   float *Ev21, float *Ev22, float *Ev23,
					   float *Ev31, float *Ev32, float *Ev33);

#endif
