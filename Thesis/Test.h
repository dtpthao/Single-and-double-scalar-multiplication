#ifndef _TEST_H
#define _TEST_H

#include <iostream>
#include <math.h>
#include <iomanip>
using namespace std;

#include "Ellipse.h"
#include "Duration.h"

#include "ShamirMul.h"
#include "option_Bin.h"
#include "option_NAF.h"
#include "option_wNAF.h"
#include "option_JSF.h"
#include "option_JT.h"
#include "option_BNBC.h"
#include "option_AKDAC.h"

#define TESTALL 1000
#define TESTSCL 1100
#define TESTSHR 1200
#define TESTwNAF 1300
#define TESTGF	1400

#define SclBin	0
#define ShrBin	2
#define ShrJSF	4
#define SclNAF	5
#define ShrNAF	6
#define SclwNAF	7
#define ShrwNAF	8
#define SclJT	9
#define ShrJT	10
#define ShrBNBC 12
#define ShrAKDAC 14
#define LIB		16

/*
Test all single scalar multiplication options 
(normal scalar and scalar via Shamir method)
*/
void TestAllScalar(big, pepoint, csprng &, big,
	double tmin[LIB + 1], double per[LIB + 1], unsigned int correct[LIB]);

/**
Test only normal single scalar multiplication options 
(without via Shamir method options)
(this test already occurs in TestAllScalar)
**/
void TestScalar(big, pepoint, csprng &, big,
	double tmin[LIB + 1], double per[LIB + 1], unsigned int correct[LIB]);

/*
Test all double scalar multiplication options
*/
void TestShamir(big, pepoint, csprng &, big,
	double tmin[LIB + 1], double per[LIB + 1], unsigned int correct[LIB]);

/*
Test values of w (from 3 to 7) for the wNAF single scalar
*/
void TestScalarwNAF(big k, pepoint P, csprng &Rng, big n,
	double tmin[7], double per[7], unsigned int correct[7]);

/*
Test values of w (from 3 to 7) for the wNAF double scalar
The result of this test is the same as TestScalarwNAF
*/
void TestShrwNAF(big k, pepoint P, csprng &Rng, big n,
	double tmin[7], double per[7], unsigned int correct[7]);

/*
Get running time of calculating number representations
*/
void TestGenForm(big k, csprng &Rng, big n, double tmin[LIB + 1], DWORD tick[LIB + 1]);

#endif
