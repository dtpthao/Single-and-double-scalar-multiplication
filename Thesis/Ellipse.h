#ifndef _ELLIPSE_H
#define _ELLIPSE_H
//#include <iostream>
//FILE _iob[] = { *stdin, *stdout, *stderr };
//extern "C" FILE * __cdecl __iob_func(void) { return _iob; }

#include <Windows.h>
extern "C" {
#include"miracl.h"
}

#define EC163 0
#define EC233 1
#define EC283 2
#define EC409 3
#define EC571 4

#define MAX_M 571
#define LEN_HEX_F2m (((MAX_M+3)/4)+1)	// 143+1
#define LEN_DWORD_F2m (((MAX_M+31)/32)*3+1)	// 18∙3+1=55

typedef struct {	//	E: y2 + x∙y = (x3 + a∙x2 + b) mod f (t)  (#E — порядок ЭК)
	DWORD m, k3, k2, k1,	// f (t) = 
		a,	// коэффициент ЭК (a  GF(2m), a = 1)
		h;	// кофактор (#E = n∙h) (h = 2)
	CHAR b[LEN_HEX_F2m],	// коэффициент ЭК b  GF(2m))
		n[LEN_HEX_F2m],	// порядок точки G (#E = n∙h) (простое число)
		Gx[LEN_HEX_F2m], Gy[LEN_HEX_F2m];	// базовая точка G = (xG , yG)
} EC_CONSTANTS_F2m_POLY;//EC_163, EC_233, EC_283, EC_409, EC_571;



typedef epoint* pepoint;

void InitStrongRNG(csprng *Rng);

int GenEC(EC_CONSTANTS_F2m_POLY EC, big a, big b, pepoint G, big x, big y, big n);

void cotnumEp(pepoint P);

void GetConstainsEC(EC_CONSTANTS_F2m_POLY &EC, unsigned int m);

#endif
