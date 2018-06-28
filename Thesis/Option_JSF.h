#ifndef _OPTION_JSF_H
#define _OPTION_JSF_H

#include "Ellipse.h"

void SubGenJSF(big Var1, big Var2, char *JSF, DWORD &lenJSF, bool &d, big RS);

DWORD GenJSF(big R, big S, char *JSFr, char *JSFs);

void ShamirMul_JSF(big a, pepoint P, big b, pepoint Q, pepoint R);

#endif
