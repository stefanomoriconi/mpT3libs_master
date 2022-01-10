#include <math.h>
#include "regOrthog.h"

void regOrthog(float *Ev11, float *Ev12, float *Ev13,
               float *Ev21, float *Ev22, float *Ev23,
               float *Ev31, float *Ev32, float *Ev33)
{

	float Ev11_out,Ev12_out,Ev13_out,Ev21_out,Ev22_out,Ev23_out,Ev31_out,Ev32_out,Ev33_out;
	float Ev1_norm,Ev2_norm,Ev3_norm,dot_Ev1Ev2,dot_Ev3Ev3_out;

	//Ev1
	Ev1_norm = (float)sqrt(Ev11[0]*Ev11[0] + Ev12[0]*Ev12[0] + Ev13[0]*Ev13[0]);
	Ev11_out = Ev11[0]/Ev1_norm;
	Ev12_out = Ev12[0]/Ev1_norm;
	Ev13_out = Ev13[0]/Ev1_norm;

	//Ev2
	dot_Ev1Ev2 = Ev11_out*Ev21[0] + Ev12_out*Ev22[0] + Ev13_out*Ev23[0];
	Ev21_out = Ev21[0] - Ev11_out*dot_Ev1Ev2;
	Ev22_out = Ev22[0] - Ev12_out*dot_Ev1Ev2;
	Ev23_out = Ev23[0] - Ev13_out*dot_Ev1Ev2;

	Ev2_norm = (float)sqrt(Ev21_out*Ev21_out + Ev22_out*Ev22_out + Ev23_out*Ev23_out);
	Ev21_out = Ev21_out/Ev2_norm;
	Ev22_out = Ev22_out/Ev2_norm;
	Ev23_out = Ev23_out/Ev2_norm;

	//Ev3
	Ev31_out = Ev12_out*Ev23_out - Ev13_out*Ev22_out;
	Ev32_out = Ev13_out*Ev21_out - Ev11_out*Ev23_out;
	Ev33_out = Ev11_out*Ev22_out - Ev12_out*Ev21_out;

	dot_Ev3Ev3_out = Ev31[0]*Ev31_out + Ev32[0]*Ev32_out + Ev33[0]*Ev33_out;
	Ev31_out = Ev31_out*dot_Ev3Ev3_out;
	Ev32_out = Ev32_out*dot_Ev3Ev3_out;
	Ev33_out = Ev33_out*dot_Ev3Ev3_out;

	Ev3_norm = (float)sqrt(Ev31_out*Ev31_out + Ev32_out*Ev32_out + Ev33_out*Ev33_out);
	Ev31_out = Ev31_out/Ev3_norm;
	Ev32_out = Ev32_out/Ev3_norm;
	Ev33_out = Ev33_out/Ev3_norm;

	//Assignment
	Ev11[0] = Ev11_out;
	Ev12[0] = Ev12_out;
	Ev13[0] = Ev13_out;
	Ev21[0] = Ev21_out;
	Ev22[0] = Ev22_out;
	Ev23[0] = Ev23_out;
	Ev31[0] = Ev31_out;
	Ev32[0] = Ev32_out;
	Ev33[0] = Ev33_out;

}