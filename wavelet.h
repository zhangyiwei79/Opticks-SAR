
#ifndef	_WAVELET_TRANSFORM_H_
#define _WAVELET_TRANSFORM_H_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct _WaveletNode WaveletNode;
struct _WaveletNode
{

  struct _WaveletNode *sibling; //Point to the peer next
    
	double *pVer;  //Vertical wavelet coefficiet
	double *pHor;  //Horizontal wavelet coefficiet
	double *pDiag; //Diagonal wavelet coefficiet
	double *pLow;  //Approximate wavelet coefficiet

	int coeffRow;  // Number of rows
	int coeffCol;  //Number of columns

};

void WaveletTransform2D(double *pSrc, int row, int col, double *pLoFilter, double *pHiFilter, int filterLen, double *pLow, double *pVer, double *pHor, double *pDiag);
void InverseWaveletTransform2D(double *pLow, double *pHor, double *pVer, double *pDiag,  int row, int col, int orig_row, int orig_col, double *pLoFilter, double *pHiFilter, int filterLen, double *pResult);

void ShiftInvariantWaveletTransform(double *pSrc, int row, int col, double *pLoFilter, double *pHiFilter, int filterLen, int nLayer, WaveletNode *pNodeList);
void ShiftInvariantInverseWaveletTransform(int row, int col, double *pLoFilter, double *pHiFilter, int filterLen, int nLayer, WaveletNode *pNodeList);

void WaveletDenoise(WaveletNode *pNodeList, double *pBuffer, int nLayer);
void ReleaseList(WaveletNode *pNodeList, int nLayer);

#endif

