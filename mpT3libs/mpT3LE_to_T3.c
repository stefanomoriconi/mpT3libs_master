#include <omp.h>
#include <stdbool.h>
#include "rt_nonfinite.h"
#include "T3_eig.h"
#include "regOrthog.h"
#include "sortElvs.h"
#include "mpT3LE_to_T3.h"

void mpT3LE_to_T3(float *T11LE_in, float *T12LE_in, float *T13LE_in,
                  float *T22LE_in, float *T23LE_in, float *T33LE_in,
                  float *El1, float *El2, float *El3,
                  float *Ev11, float *Ev12, float *Ev13,
                  float *Ev21, float *Ev22, float *Ev23,
                  float *Ev31, float *Ev32, float *Ev33,
                  bool *mskValid, int *idx, int *numel)
{
	int i;
	#pragma omp parallel for
    for (i=0; i<numel[0]; ++i)
      {
        T3_eig( &T11LE_in[idx[i]] , &T12LE_in[idx[i]] , &T13LE_in[idx[i]] ,
        		    &T22LE_in[idx[i]] , &T23LE_in[idx[i]] , &T33LE_in[idx[i]] ,
        		    &El1[idx[i]]    , &El2[idx[i]]    , &El3[idx[i]]    ,
        		    &Ev11[idx[i]]   , &Ev12[idx[i]]   , &Ev13[idx[i]]   ,
        		    &Ev21[idx[i]]   , &Ev22[idx[i]]   , &Ev23[idx[i]]   ,
        		    &Ev31[idx[i]]   , &Ev32[idx[i]]   , &Ev33[idx[i]]   );

        regOrthog( &Ev11[idx[i]]   , &Ev12[idx[i]]   , &Ev13[idx[i]]   ,
                   &Ev21[idx[i]]   , &Ev22[idx[i]]   , &Ev23[idx[i]]   ,
                   &Ev31[idx[i]]   , &Ev32[idx[i]]   , &Ev33[idx[i]]   );

        El1[idx[i]] = (float)exp(El1[idx[i]]);
        El2[idx[i]] = (float)exp(El2[idx[i]]);
        El3[idx[i]] = (float)exp(El3[idx[i]]);

        sortElvs( &El1[idx[i]]    , &El2[idx[i]]    , &El3[idx[i]]    ,
                  &Ev11[idx[i]]   , &Ev12[idx[i]]   , &Ev13[idx[i]]   ,
                  &Ev21[idx[i]]   , &Ev22[idx[i]]   , &Ev23[idx[i]]   ,
                  &Ev31[idx[i]]   , &Ev32[idx[i]]   , &Ev33[idx[i]]   );

        /*  Check for Numerical Degenerate Case - Re-Set as Isotropic */
  		  if ((!((!rtIsInfF(El1[idx[i]])) && (!rtIsNaNF(El1[idx[i]])))) || (!((!rtIsInfF(El2[idx[i]])) &&
        		  (!rtIsNaNF(El2[idx[i]])))) || (!((!rtIsInfF(El3[idx[i]])) && (!rtIsNaNF(El3[idx[i]])))) ||
      		  (!((!rtIsInfF(Ev11[idx[i]])) && (!rtIsNaNF(Ev11[idx[i]])))) || (!((!rtIsInfF(Ev12[idx[i]])) &&
        		  (!rtIsNaNF(Ev12[idx[i]])))) || (!((!rtIsInfF(Ev13[idx[i]])) && (!rtIsNaNF(Ev13[idx[i]])))) ||
      		  (!((!rtIsInfF(Ev21[idx[i]])) && (!rtIsNaNF(Ev21[idx[i]])))) || (!((!rtIsInfF(Ev22[idx[i]])) &&
        		  (!rtIsNaNF(Ev22[idx[i]])))) || (!((!rtIsInfF(Ev23[idx[i]])) && (!rtIsNaNF(Ev23[idx[i]])))) ||
      		  (!((!rtIsInfF(Ev31[idx[i]])) && (!rtIsNaNF(Ev31[idx[i]])))) || (!((!rtIsInfF(Ev32[idx[i]])) &&
        		  (!rtIsNaNF(Ev32[idx[i]])))) || (!((!rtIsInfF(Ev33[idx[i]])) && (!rtIsNaNF(Ev33[idx[i]]))))) {
    		  El1[idx[i]] = 1.0F;
    		  El2[idx[i]] = 1.0F;
    		  El3[idx[i]] = 1.0F;
    		  Ev11[idx[i]] = 1.0F;
    		  Ev12[idx[i]] = 0.0F;
    		  Ev13[idx[i]] = 0.0F;
    		  Ev21[idx[i]] = 0.0F;
    		  Ev22[idx[i]] = 1.0F;
    		  Ev23[idx[i]] = 0.0F;
    		  Ev31[idx[i]] = 0.0F;
    		  Ev32[idx[i]] = 0.0F;
    		  Ev33[idx[i]] = 1.0F;
          mskValid[idx[i]] = false;
  		  }

      }
}