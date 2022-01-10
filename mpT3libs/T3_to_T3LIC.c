/* Include Files */
#include "rt_nonfinite.h"
#include "T3_to_T3LIC.h"

/* Function Definitions */

void T3_to_T3LIC(float *El1_in, float *El2_in, float *El3_in,
                 float *Ev11_in, float *Ev12_in, float *Ev13_in,
                 float *Ev21_in, float *Ev22_in, float *Ev23_in,
                 float *Ev31_in, float *Ev32_in, float *Ev33_in,
                 float *T11, float *T12, float *T13,
                 float *T22, float *T23, float *T33)
{
  float El1, El2, El3, Ev11, Ev12, Ev13, Ev21, Ev22, Ev23, Ev31, Ev32, Ev33;
  El1 = *El1_in;
  El2 = *El2_in;
  El3 = *El3_in;

  Ev11 = *Ev11_in;
  Ev12 = *Ev12_in;
  Ev13 = *Ev13_in;

  Ev21 = *Ev21_in;
  Ev22 = *Ev22_in;
  Ev23 = *Ev23_in;

  Ev31 = *Ev31_in;
  Ev32 = *Ev32_in;
  Ev33 = *Ev33_in;

  *T11 = (El1 * (Ev11 * Ev11) + El2 * (Ev21 * Ev21)) + El3 * (Ev31 * Ev31);
  *T12 = (El1 * Ev11 * Ev12 + El2 * Ev21 * Ev22) + El3 * Ev31 * Ev32;
  *T13 = (El1 * Ev11 * Ev13 + El2 * Ev21 * Ev23) + El3 * Ev31 * Ev33;
  *T22 = (El1 * (Ev12 * Ev12) + El2 * (Ev22 * Ev22)) + El3 * (Ev32 * Ev32);
  *T23 = (El1 * Ev12 * Ev13 + El2 * Ev22 * Ev23) + El3 * Ev32 * Ev33;
  *T33 = (El1 * (Ev13 * Ev13) + El2 * (Ev23 * Ev23)) + El3 * (Ev33 * Ev33);
}

