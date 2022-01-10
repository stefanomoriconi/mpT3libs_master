#ifndef T3_TO_T3LE_H
#define T3_TO_T3LE_H

/* Include Files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "T3_to_T3LE_types.h"

/* Function Declarations */
extern void T3_to_T3LE(float *El1_in, float *El2_in, float *El3_in,
					   float *Ev11_in, float *Ev12_in, float *Ev13_in,
					   float *Ev21_in, float *Ev22_in, float *Ev23_in,
					   float *Ev31_in, float *Ev32_in, float *Ev33_in,
					   float *T11LE, float *T12LE, float *T13LE,
					   float *T22LE, float *T23LE, float *T33LE);

#endif
