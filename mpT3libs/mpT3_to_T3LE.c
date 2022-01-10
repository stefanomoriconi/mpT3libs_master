#include <omp.h>
#include "T3_to_T3LE.h"
#include "mpT3_to_T3LE.h"

void mpT3_to_T3LE(float *El1_in, float *El2_in, float *El3_in,
                  float *Ev11_in, float *Ev12_in, float *Ev13_in,
                  float *Ev21_in, float *Ev22_in, float *Ev23_in,
                  float *Ev31_in, float *Ev32_in, float *Ev33_in,
                  float *T11LE, float *T12LE, float *T13LE,
                  float *T22LE, float *T23LE, float *T33LE,
                  int *idx, int *numel)
{
	int i;
	#pragma omp parallel for
    for (i=0; i<numel[0]; ++i)
      {
        T3_to_T3LE( &El1_in[idx[i]] , &El2_in[idx[i]] , &El3_in[idx[i]] ,
                    &Ev11_in[idx[i]] , &Ev12_in[idx[i]] , &Ev13_in[idx[i]] ,
                    &Ev21_in[idx[i]] , &Ev22_in[idx[i]] , &Ev23_in[idx[i]] ,
                    &Ev31_in[idx[i]] , &Ev32_in[idx[i]] , &Ev33_in[idx[i]] ,
                    &T11LE[idx[i]] , &T12LE[idx[i]] , &T13LE[idx[i]],
                    &T22LE[idx[i]] , &T23LE[idx[i]] , &T33LE[idx[i]] );
      }
}