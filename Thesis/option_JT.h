#ifndef _OPTION_JT_H
#define _OPTION_JT_H

#include <iostream>
#include "Ellipse.h"
#include "Duration.h"
using namespace std;

/*
Computing the JT{+/-1; +/-3} of an positive integer
*/
DWORD GenJT(big k, char* rJT);
int checkGenJT(char* rJT, DWORD len, big k);

//JT option of single scalar multiplication
void ScalarMul_JT(big k, pepoint P, pepoint R);

//JT option of double scalar multiplication
void ShamirMul_JT(big a, pepoint P, big b, pepoint Q, pepoint R);

/*
*/
void Test_SclMul_JT(csprng &Rng, pepoint P);
void Test_ShrMul_JT(csprng &Rng, pepoint P);

#endif
