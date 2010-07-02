
#ifndef	_WAVELET_TRANSFORM_H_
#define _WAVELET_TRANSFORM_H_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

void WaveletTransform2D(double *pSrc, int row, int col, double *pLoFilter, double *pHiFilter, int filterLen, double *pLow, double *pVer, double *pHor, double *pDiag);



#endif

