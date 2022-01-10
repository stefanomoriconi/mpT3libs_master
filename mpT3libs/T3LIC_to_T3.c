/* Include Files */
#include "rt_nonfinite.h"
#include "T3LIC_to_T3.h"
#include "T3_eig.h"

/* Function Definitions */

void T3LIC_to_T3(float *T11_in, float *T12_in, float *T13_in,
                 float *T22_in, float *T23_in, float *T33_in,
                 float *El1, float *El2, float *El3,
                 float *Ev11, float*Ev12, float *Ev13,
                 float *Ev21, float *Ev22, float *Ev23,
                 float *Ev31, float *Ev32, float *Ev33)
{
  T3_eig(T11_in, T12_in, T13_in, T22_in, T23_in, T33_in, El1, El2, El3, Ev11, Ev12, Ev13, Ev21,
         Ev22, Ev23, Ev31, Ev32, Ev33);

  /*  Check for Numerical Degenerate Case - Re-Set as Isotropic */
  if ((!((!rtIsInfF(*El1)) && (!rtIsNaNF(*El1)))) || (!((!rtIsInfF(*El2)) &&
        (!rtIsNaNF(*El2)))) || (!((!rtIsInfF(*El3)) && (!rtIsNaNF(*El3)))) ||
      (!((!rtIsInfF(*Ev11)) && (!rtIsNaNF(*Ev11)))) || (!((!rtIsInfF(*Ev12)) &&
        (!rtIsNaNF(*Ev12)))) || (!((!rtIsInfF(*Ev13)) && (!rtIsNaNF(*Ev13)))) ||
      (!((!rtIsInfF(*Ev21)) && (!rtIsNaNF(*Ev21)))) || (!((!rtIsInfF(*Ev22)) &&
        (!rtIsNaNF(*Ev22)))) || (!((!rtIsInfF(*Ev23)) && (!rtIsNaNF(*Ev23)))) ||
      (!((!rtIsInfF(*Ev31)) && (!rtIsNaNF(*Ev31)))) || (!((!rtIsInfF(*Ev32)) &&
        (!rtIsNaNF(*Ev32)))) || (!((!rtIsInfF(*Ev33)) && (!rtIsNaNF(*Ev33))))) {
    *El1 = 1.0F;
    *El2 = 1.0F;
    *El3 = 1.0F;
    *Ev11 = 1.0F;
    *Ev12 = 0.0F;
    *Ev13 = 0.0F;
    *Ev21 = 0.0F;
    *Ev22 = 1.0F;
    *Ev23 = 0.0F;
    *Ev31 = 0.0F;
    *Ev32 = 0.0F;
    *Ev33 = 1.0F;
  }
}

