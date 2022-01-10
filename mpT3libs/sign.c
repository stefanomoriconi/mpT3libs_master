/* Include Files */
#include "rt_nonfinite.h"
#include "T3_eig.h"
#include "sign.h"

/* Function Definitions */

void b_sign(float *x)
{
  if (*x < 0.0F) {
    *x = -1.0F;
  } else if (*x > 0.0F) {
    *x = 1.0F;
  } else {
    if (*x == 0.0F) {
      *x = 0.0F;
    }
  }
}
