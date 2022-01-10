#include <omp.h>
#include "T3_to_T3LIC.h"
#include "mpT3_to_T3LIC.h"

void mpT3_to_T3LIC(float *El1_in, float *El2_in, float *El3_in,
                   float *Ev11_in, float *Ev12_in, float *Ev13_in,
                   float *Ev21_in, float *Ev22_in, float *Ev23_in,
                   float *Ev31_in, float *Ev32_in, float *Ev33_in,
                   float *T11, float *T12, float *T13,
                   float *T22, float *T23, float *T33,
                   int *idx, int *numel)
{
	int i;
	#pragma omp parallel for
    for (i=0; i<numel[0]; ++i)
      {
        T3_to_T3LIC( &El1_in[idx[i]] , &El2_in[idx[i]] , &El3_in[idx[i]] ,
                     &Ev11_in[idx[i]] , &Ev12_in[idx[i]] , &Ev13_in[idx[i]] ,
                     &Ev21_in[idx[i]] , &Ev22_in[idx[i]] , &Ev23_in[idx[i]] ,
                     &Ev31_in[idx[i]] , &Ev32_in[idx[i]] , &Ev33_in[idx[i]] ,
                     &T11[idx[i]] , &T12[idx[i]] , &T13[idx[i]],
                     &T22[idx[i]] , &T23[idx[i]] , &T33[idx[i]] );
      }
}