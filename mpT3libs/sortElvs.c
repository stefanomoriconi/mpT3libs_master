//#include <stdio.h>
//#include <stdlib.h>
#include <math.h>
#include "sortElvs.h"

void sortElvs(float *El1, float *El2, float *El3, 
			  float *Ev11, float *Ev12, float *Ev13,
              float *Ev21, float *Ev22, float *Ev23,
              float *Ev31, float *Ev32, float *Ev33)
{
	float El1_temp,El2_temp,El3_temp;
	float Ev11_temp, Ev12_temp, Ev13_temp, Ev21_temp, Ev22_temp, Ev23_temp, Ev31_temp, Ev32_temp, Ev33_temp;

	El1_temp = El1[0];
	El2_temp = El2[0];
	El3_temp = El3[0];

	Ev11_temp = Ev11[0];
	Ev12_temp = Ev12[0];
	Ev13_temp = Ev13[0];
	Ev21_temp = Ev21[0];
	Ev22_temp = Ev22[0];
	Ev23_temp = Ev23[0];
	Ev31_temp = Ev31[0];
	Ev32_temp = Ev32[0];
	Ev33_temp = Ev33[0];

	if (fabs(El1_temp) != 1.0 && fabs(El2_temp) != 1.0 && fabs(El3_temp) != 1.0 ){

		if ( (fabs(El1_temp) <= fabs(El2_temp)) && (fabs(El2_temp) <= fabs(El3_temp)) ){
			//printf("El1_temp = %f < El2_temp = %f < El3_temp = %f \n", El1_temp,El2_temp,El3_temp);
			El1[0] = El1_temp;
			El2[0] = El2_temp;
			El3[0] = El3_temp;
			Ev11[0] = Ev11_temp;
			Ev12[0] = Ev12_temp;
			Ev13[0] = Ev13_temp;
			Ev21[0] = Ev21_temp;
			Ev22[0] = Ev22_temp;
			Ev23[0] = Ev23_temp;
			Ev31[0] = Ev31_temp;
			Ev32[0] = Ev32_temp;
			Ev33[0] = Ev33_temp;
			//printf("El1 = %f , El2 = %f , El3 = %f \n", El1[0],El2[0],El3[0]);
		}
		else if ( (fabs(El1_temp) <= fabs(El3_temp)) && (fabs(El3_temp) <= fabs(El2_temp)) ){
			//printf("El1_temp = %f < El3_temp = %f < El2_temp = %f \n", El1_temp,El3_temp,El2_temp);
			El1[0] = El1_temp;
			El2[0] = El3_temp;
			El3[0] = El2_temp;
			Ev11[0] = Ev11_temp;
			Ev12[0] = Ev12_temp;
			Ev13[0] = Ev13_temp;
			Ev21[0] = Ev31_temp;
			Ev22[0] = Ev32_temp;
			Ev23[0] = Ev33_temp;
			Ev31[0] = Ev21_temp;
			Ev32[0] = Ev22_temp;
			Ev33[0] = Ev23_temp;
			//printf("El1 = %f , El2 = %f , El3 = %f \n", El1[0],El2[0],El3[0]);
		}
		else if ( (fabs(El2_temp) <= fabs(El1_temp)) && (fabs(El1_temp) <= fabs(El3_temp)) ){
			//printf("El2_temp = %f < El1_temp = %f < El3_temp = %f \n", El2_temp,El1_temp,El3_temp);
			El1[0] = El2_temp;
			El2[0] = El1_temp;
			El3[0] = El3_temp;
			Ev11[0] = Ev21_temp;
			Ev12[0] = Ev22_temp;
			Ev13[0] = Ev23_temp;
			Ev21[0] = Ev11_temp;
			Ev22[0] = Ev12_temp;
			Ev23[0] = Ev13_temp;
			Ev31[0] = Ev31_temp;
			Ev32[0] = Ev32_temp;
			Ev33[0] = Ev33_temp;
			//printf("El1 = %f , El2 = %f , El3 = %f \n", El1[0],El2[0],El3[0]);
		}
		else if ( (fabs(El2_temp) <= fabs(El3_temp)) && (fabs(El3_temp) <= fabs(El1_temp)) ){
			//printf("El2_temp = %f < El3_temp = %f < El1_temp = %f \n", El2_temp,El3_temp,El1_temp);
			El1[0] = El2_temp;
			El2[0] = El3_temp;
			El3[0] = El1_temp;
			Ev11[0] = Ev21_temp;
			Ev12[0] = Ev22_temp;
			Ev13[0] = Ev23_temp;
			Ev21[0] = Ev31_temp;
			Ev22[0] = Ev32_temp;
			Ev23[0] = Ev33_temp;
			Ev31[0] = Ev11_temp;
			Ev32[0] = Ev12_temp;
			Ev33[0] = Ev13_temp;
			//printf("El1 = %f , El2 = %f , El3 = %f \n", El1[0],El2[0],El3[0]);
		}
		else if ( (fabs(El3_temp) <= fabs(El1_temp)) && (fabs(El1_temp) <= fabs(El2_temp)) ){
			//printf("El3_temp = %f < El1_temp = %f < El2_temp = %f \n", El3_temp,El1_temp,El2_temp);
			El1[0] = El3_temp;
			El2[0] = El1_temp;
			El3[0] = El2_temp;
			Ev11[0] = Ev31_temp;
			Ev12[0] = Ev32_temp;
			Ev13[0] = Ev33_temp;
			Ev21[0] = Ev11_temp;
			Ev22[0] = Ev12_temp;
			Ev23[0] = Ev13_temp;
			Ev31[0] = Ev21_temp;
			Ev32[0] = Ev22_temp;
			Ev33[0] = Ev23_temp;
			//printf("El1 = %f , El2 = %f , El3 = %f \n", El1[0],El2[0],El3[0]);
		}
		else if ( (fabs(El3_temp) <= fabs(El2_temp)) && (fabs(El2_temp) <= fabs(El1_temp)) ){
			//printf("El3_temp = %f < El2_temp = %f < El1_temp = %f \n", El3_temp,El2_temp,El1_temp);
			El1[0] = El3_temp;
			El2[0] = El2_temp;
			El3[0] = El1_temp;
			Ev11[0] = Ev31_temp;
			Ev12[0] = Ev32_temp;
			Ev13[0] = Ev33_temp;
			Ev21[0] = Ev21_temp;
			Ev22[0] = Ev22_temp;
			Ev23[0] = Ev23_temp;
			Ev31[0] = Ev11_temp;
			Ev32[0] = Ev12_temp;
			Ev33[0] = Ev13_temp;
			//printf("El1 = %f , El2 = %f , El3 = %f \n", El1[0],El2[0],El3[0]);
		}
		else
		{
			//printf("Unrecognised Case: \n");
			//printf("El1_temp = %f , El2_temp = %f , El3_temp = %f \n", El1_temp,El2_temp,El3_temp);
			//printf("|El1_temp| = %f , |El2_temp| = %f , |El3_temp| = %f \n", (float)fabs(El1_temp),(float)fabs(El2_temp),(float)fabs(El3_temp));
			//printf("El1 = %f , El2 = %f , El3 = %f \n", El1[0],El2[0],El3[0]);
		}

	}

}