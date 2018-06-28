#ifndef _COUNT_TIME_H
#define _COUNT_TIME_H

#include <Windows.h>
#include <time.h>
extern "C" {
#include "miracl.h"
}

typedef epoint* pepoint;

typedef struct {
	LARGE_INTEGER start;
	LARGE_INTEGER stop;
} stopWatch;

void startTimer(stopWatch *timer);
void stopTimer(stopWatch *timer);
double LIToSecs(LARGE_INTEGER *L);
double getElapsedTime(stopWatch *timer);
DWORD getQuadPart(stopWatch *timer);


void SclDuration(big k, pepoint P, pepoint R,
	void(*func) (big, pepoint, pepoint), double &min);

void ShrDuration(big k, pepoint P, pepoint R, 
	void(*func) (big, pepoint, big, pepoint, pepoint), double &min);

void ShrDuration(big a, pepoint P, big b, pepoint Q, pepoint R, 
	void(*func) (big, pepoint, big, pepoint, pepoint), double &min);

void GenFormSclDuration(big k, DWORD (*func) (big, char*), double &min, DWORD&);

void GenFormShrDuration(big a, big b, DWORD(*func) (big, char*), double &min, DWORD&);

#endif
