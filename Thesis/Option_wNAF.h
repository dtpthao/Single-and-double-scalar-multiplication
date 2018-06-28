#ifndef _OPTION_WNAF_H
#define _OPTION_WNAF_H

#include <Windows.h>
#include "Ellipse.h"

#define wNAF 4
#define Pow2wNAF (1<<wNAF)
#define HalfPow2wNAF (Pow2wNAF>>1)
#define MAXeSDR (Pow2wNAF-1)

/*
Computing the wNAF of a positive integer 
*/
DWORD GenwNAF(big K, char *SDR);

/*
Computing the wNAF of a positive integer.
This option is much better than above for point multiplication when w = 6, 7.
But in the report we used w = 4 as the best value, so it wasn't mentioned.
*/
DWORD GenwNAF(big K, char *SDR, char *eSDRout, DWORD &lene);

/*
Those are reported options for single and double scalar multiplication.
They come with GenwNAF(big K, char *SDR).
*/
void PreMul_wNAF(pepoint P, pepoint *listP);
void ScalarMul_wNAF(big k, pepoint P, pepoint R);
void ShamirMul_wNAF(big a, pepoint P, big b, pepoint Q, pepoint R);

/*
Тhose are unreported options for single and double scalar multiplication.
They come with GenwNAF(big K, char *SDR, char *eSDRout, DWORD &lene).
*/
void PreMul_wNAF(pepoint P, pepoint *listP, char *eSDR, DWORD lene);
void ScalarMul_wNAF1(big k, pepoint P, pepoint R);
void ShamirMul_wNAF1(big a, pepoint P, big b, pepoint Q, pepoint R);


/*******************************************************************
Those functions and definitions below are (~)completely copied above
for testing different values of w (3 - 7). Ignore them, please!!!
********************************************************************/

////////////////////////////////////////////
#define w3NAF 3
#define Pow2w3NAF (1<<w3NAF)
#define HalfPow2w3NAF (Pow2w3NAF>>1)
#define MAXe3SDR (Pow2w3NAF-1)

DWORD Genw3NAF(big K, char *SDR);
void ScalarMul_w3NAF(big k, pepoint P, pepoint R);
void ShamirMul_w3NAF(big a, pepoint P, big b, pepoint Q, pepoint R);

////////////////////////////////////////////
#define w4NAF 4
#define Pow2w4NAF (1<<w4NAF)
#define HalfPow2w4NAF (Pow2w4NAF>>1)
#define MAXe4SDR (Pow2w4NAF-1)

DWORD Genw4NAF(big K, char *SDR);
void ScalarMul_w4NAF(big k, pepoint P, pepoint R);
void ShamirMul_w4NAF(big a, pepoint P, big b, pepoint Q, pepoint R);

////////////////////////////////////////////
#define w5NAF 5
#define Pow2w5NAF (1<<w5NAF)
#define HalfPow2w5NAF (Pow2w5NAF>>1)
#define MAXe5SDR (Pow2w5NAF-1)

DWORD Genw5NAF(big K, char *SDR);
void ScalarMul_w5NAF(big k, pepoint P, pepoint R);
void ShamirMul_w5NAF(big a, pepoint P, big b, pepoint Q, pepoint R);

////////////////////////////////////////////
#define w6NAF 6
#define Pow2w6NAF (1<<w6NAF)
#define HalfPow2w6NAF (Pow2w6NAF>>1)
#define MAXe6SDR (Pow2w6NAF-1)

DWORD Genw6NAF(big K, char *SDR);
void ScalarMul_w6NAF(big k, pepoint P, pepoint R);
void ShamirMul_w6NAF(big a, pepoint P, big b, pepoint Q, pepoint R);

////////////////////////////////////////////
#define w7NAF 7
#define Pow2w7NAF (1<<w7NAF)
#define HalfPow2w7NAF (Pow2w7NAF>>1)
#define MAXe7SDR (Pow2w7NAF-1)

DWORD Genw7NAF(big K, char *SDR);
void ScalarMul_w7NAF(big k, pepoint P, pepoint R);
void ShamirMul_w7NAF(big a, pepoint P, big b, pepoint Q, pepoint R);
#endif

