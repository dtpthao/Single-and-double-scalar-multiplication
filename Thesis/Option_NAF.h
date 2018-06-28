#ifndef _OPTION_NAF_H
#define _OPTION_NAF_H

#include "Ellipse.h"

/*
Computing the NAF of a positive integer
*/
DWORD GenNAF(big k, char *rNAF);

//Single scalar multiplication
void ScalarMul_NAF(big k, epoint* P, epoint* R);

//Double scalar multiplication
void PreMul(pepoint P, pepoint Q, pepoint *listpoint);
void ShamirMul_NAF(big a, pepoint P, big b, pepoint Q, pepoint R);

#endif

