

#include "wavelet.h"

//Perform column convlolution
void ColConvolution2D(double *pSrc, double *pFilter, int filterLength, int row, int col, double *pRes, bool bAddZero)
{
	double temp, a, b;
	int i,j,k;
	double *pData;

	if ((pSrc == NULL) || (pRes == NULL))
	{
		return;
	}

	pData = (double *)malloc(sizeof(double)*(row+2*(filterLength-1)));

	for (k=0; k<col; k++)
	{
		for (i=0; i<(filterLength-1); i++)
		{		
			if (!bAddZero)
			{
			    *(pData + i) = *(pSrc + (filterLength-2)*col -i*col);
			}
			else
			{
				*(pData + i) = 0;
			}
		}
		pData = pData + filterLength - 1;

		for (i=0; i<row; i++)
		{
			*(pData + i) = *(pSrc+i*col);
		}
		pData = pData + row;

		for (i=0; i<(filterLength-1); i++)
		{		
			if (!bAddZero)
			{
			    *(pData + i) = *(pSrc+col*(row-1)-i*col);
			}
			else
			{
				*(pData + i) = 0;
			}
		}

		pData = pData - row -filterLength + 1;

	    for (i=0; i<(row+filterLength-1); i++)
	    {
		    temp = 0;

		    for (j=0; j<filterLength; j++)
		    {
			    a = *(pFilter + j);
			    b = *(pData + i + j);
			    temp = temp + a*b;
		    }
		
		    *(pRes + k + i*col) = temp;
		}

		pSrc++;

	}

	free(pData);
}

//Perform row convlolution
void RowConvolution2D(double *pSrc, double *pFilter, int filterLength, int row, int col, double *pRes, bool bAddZero)
{
	double temp, a, b;
	double *pData;
	int i,j,k;

	if ((pSrc == NULL) || (pRes == NULL))
	{
		return;
	}

	pData = (double *)malloc(sizeof(double)*(col+2*(filterLength-1)));

	for (k=0; k<row; k++)
	{
		//Fill the outer part
		for (i=0; i<(filterLength-1); i++)
		{
			if (!bAddZero)
			{
			    *(pData + i) = *(pSrc + filterLength - 2 - i);
			}
			else
			{
				*(pData + i) = 0;
			}
		}
		pData = pData + filterLength - 1;

		for (i=0; i<col; i++)
		{
			*(pData + i) = *(pSrc+i);
		}
		pData = pData + col;

		//Fill the outer part
		for (i=0; i<(filterLength-1); i++)
		{		
			if (!bAddZero)
			{
			    *(pData + i) = *(pSrc+col-1-i);
			}
			else
			{
				*(pData + i) = 0;
			}
		}

        pData = pData - col -filterLength + 1;

	    for (i=0; i<(col+filterLength-1); i++)
	    {
		    temp = 0;

		    for (j=0; j<filterLength; j++)
		    {
			    a = *(pFilter + j);
				b = *(pData + i + j);

			    temp = temp + a*b;
		    }
		
		    *pRes = temp;
			pRes++;
		}

		pSrc = pSrc + col;
	}

	free(pData);
}

void WaveletTransform2D(double *pSrc, int row, int col, double *pLoFilter, double *pHiFilter, int filterLen, double *pLow, double *pVer, double *pHor, double *pDiag)
{
	int nIndex = 0;
	int coeff_row = row + filterLen -1;
	int coeff_col = col + filterLen - 1;
	
	double *pTemp1 = (double *)malloc(sizeof(double)*row*coeff_col);
	double *pTemp2 = (double *)malloc(sizeof(double)*coeff_row*coeff_col);
	
	RowConvolution2D(pSrc, pLoFilter, filterLen, row, col, pTemp1, false);
	
	//Get the Low frquency part
	ColConvolution2D(pTemp1, pLoFilter, filterLen, row, col+filterLen-1, pTemp2, false);
	//Downsample the coefficient
	nIndex = 0;
	for (int i=0; i<coeff_row; i++)
	{
		for(int j=0; j<coeff_col; j++)
		{
			if ((i%2 == 1) && (j%2 == 1))
			{
		        pLow[nIndex] = pTemp2[i*coeff_col+j];
				nIndex++;
			}
		}
	}
	
	//Get Horizontal coefficients
	ColConvolution2D(pTemp1, pHiFilter, filterLen, row, col+filterLen-1, pTemp2, false);
	nIndex = 0;
	for (int i=0; i<coeff_row; i++)
	{
		for(int j=0; j<coeff_col; j++)
		{
			if ((i%2 == 1) && (j%2 == 1))
			{
		        pHor[nIndex] = pTemp2[i*coeff_col+j];
				nIndex++;
			}
		}
	}

	RowConvolution2D(pSrc, pHiFilter, filterLen, row, col, pTemp1, false);
	//Get vertical coefficients
	ColConvolution2D(pTemp2, pLoFilter, filterLen, row, col+filterLen-1, pTemp2, false);
	nIndex = 0;
	for (int i=0; i<coeff_row; i++)
	{
		for(int j=0; j<coeff_col; j++)
		{
			if ((i%2 == 1) && (j%2 == 1))
			{
		        pVer[nIndex] = pTemp2[i*coeff_col+j];
				nIndex++;
			}
		}
	}
	
	//Get the diagonal coefficients
	ColConvolution2D(pTemp1, pHiFilter, filterLen, row, col+filterLen-1, pTemp2, false);
	nIndex = 0;
	for (int i=0; i<coeff_row; i++)
	{
		for(int j=0; j<coeff_col; j++)
		{
			if ((i%2 == 1) && (j%2 == 1))
			{
		        pDiag[nIndex] = pTemp2[i*coeff_col+j];
				nIndex++;
			}
		}
	}
	
	free(pTemp1);
	free(pTemp2);
	
	return;
}
