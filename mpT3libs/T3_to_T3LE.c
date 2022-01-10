/* Include Files */
#include "rt_nonfinite.h"
#include "T3_to_T3LE.h"

/* Function Definitions */

void T3_to_T3LE(float *El1_in, float *El2_in, float *El3_in,
                float *Ev11_in, float *Ev12_in, float *Ev13_in,
                float *Ev21_in, float *Ev22_in, float *Ev23_in,
                float *Ev31_in, float *Ev32_in, float *Ev33_in,
                float *T11LE, float *T12LE, float *T13LE,
                float *T22LE, float *T23LE, float *T33LE)
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

  *T11LE = (Ev11 * (float)log(El1) * Ev11 + Ev21 * (float)log(El2) * Ev21) +
    Ev31 * (float)log(El3) * Ev31;

  *T12LE = (Ev12 * (float)log(El1) * Ev11 + Ev22 * (float)log(El2) * Ev21) +
    Ev32 * (float)log(El3) * Ev31;

  *T13LE = (Ev13 * (float)log(El1) * Ev11 + Ev23 * (float)log(El2) * Ev21) +
    Ev33 * (float)log(El3) * Ev31;

  *T22LE = (Ev12 * (float)log(El1) * Ev12 + Ev22 * (float)log(El2) * Ev22) +
    Ev32 * (float)log(El3) * Ev32;

  *T23LE = (Ev13 * (float)log(El1) * Ev12 + Ev23 * (float)log(El2) * Ev22) +
    Ev33 * (float)log(El3) * Ev32;

  *T33LE = (Ev13 * (float)log(El1) * Ev13 + Ev23 * (float)log(El2) * Ev23) +
    Ev33 * (float)log(El3) * Ev33;

}
