#ifndef T3_TO_T3LIC_H
#define T3_TO_T3LIC_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "T3_to_T3LIC_types.h"

/* Function Declarations */
extern void T3_to_T3LIC(float *El1_in, float *El2_in, float *El3_in,
						float *Ev11_in, float *Ev12_in, float *Ev13_in,
						float *Ev21_in, float *Ev22_in, float *Ev23_in,
						float *Ev31_in, float *Ev32_in, float *Ev33_in,
						float *T11, float *T12, float *T13,
						float *T22, float *T23, float *T33);

#endif
