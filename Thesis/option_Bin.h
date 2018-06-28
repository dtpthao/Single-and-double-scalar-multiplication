#ifndef _OPTION_BIN_H
#define _OPTION_BIN_H

extern "C" {
#include "miracl.h"
}

typedef epoint* pepoint;

//ScalarMul
void ScalarMul_Bin(big k, pepoint P, pepoint R);

//ShamirMul
void PreMul_Bin(pepoint P, pepoint Q, pepoint *plist);

void ShamirMul_Bin(big a, pepoint P, big b, pepoint Q, pepoint R);

#endif
