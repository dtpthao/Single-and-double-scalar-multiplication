#ifndef _OPTION_AKDAC_H
#define _OPTION_AKDAC_H

#include <iostream>
#include "Ellipse.h"

void R1(big d, big e, pepoint Ru, pepoint Rv);

void R2(big d, pepoint Ru);

void ShamirMul_AKDAC(big a, pepoint P, big b, pepoint Q, pepoint R);

void TestShrMulAKDAC(csprng &Rng, pepoint P);

#endif
