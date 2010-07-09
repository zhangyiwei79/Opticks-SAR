

#include "wavelet.h"

//Daubechies 2 wavelet basis
double pLoFilter[4] = {0.4830, 0.8365, 0.2241, -0.1294};
double pHiFilter[4] = {-0.1294,-0.2241, 0.8365, -0.4830};
double pRecHiFilter[4] = {-0.4830, 0.8365, -0.2241, -0.1294};
double pRecLoFilter[4] = {-0.1294,0.2241, 0.8365, 0.4830};

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
void RowConvolution2D(double *pSrc, double *pFilter, int filterLength, int row, int col, double *pRes, bool bAddZero, bool bAddResult)
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
		
			if (bAddResult)
			{
				*pRes = *pRes + temp;
			}
			else
			{
		        *pRes = temp;
			}
			pRes++;
		}

		pSrc = pSrc + col;
	}

	free(pData);
}




void IDWT_Partial(double *pSrc, int row, int col, int orig_row, int orig_col, double *pFilter1, double *pFilter2, int filterLen, double *pResult)
{
	int count = 0;

	double *pTemp1 = (double *)malloc(sizeof(double)*(2*row-1)*col);
	double *pTemp2 = (double *)malloc(sizeof(double)*(2*row-1+filterLen-1)*col);
	double *pTemp3 = (double *)malloc(sizeof(double)*(2*row-1+filterLen-1)*(2*col-1));
	

	//Up Sample rows
	for (int i=0; i<2*row-1; i++)
	{
		for (int j=0; j<col; j++)
		{
			if (i%2 == 1)
			{
				*(pTemp1+i*col+j) = 0;
			}
			else
			{
				*(pTemp1+i*col+j) = *(pSrc+i/2*col+j);
			}
		}
	}

	//ShowResult(pTemp1, 2*row-1, col);
	ColConvolution2D(pTemp1, pFilter1, filterLen, 2*row-1, col, pTemp2, true);
	//ShowResult(pTemp2, 2*row-1+filterLen-1,col);

	//Up Sample cols
	count = 2*row-1+filterLen-1;
	for (int i=0; i<count; i++)
	{
		for (int j=0; j<(2*col-1); j++)
		{
			if (j%2 == 1)
			{
				*(pTemp3+i*(col*2-1)+j) = 0;
			}
			else
			{
				*(pTemp3+i*(col*2-1)+j) = *(pTemp2 + i*col + j/2);
			}
		}
	}

	RowConvolution2D(pTemp3, pFilter2, filterLen, count, (2*col-1), pResult, true, true);

	free(pTemp1);
	free(pTemp2);
	free(pTemp3);

}

//Perform discrete 2 dimensinal inverse wavelet transform
void InverseWaveletTransform2D(double *pLow, double *pHor, double *pVer, double *pDiag,  int row, int col, int orig_row, int orig_col, double *pLoFilter, double *pHiFilter, int filterLen, double *pResult)
{
	int row_shift = 0;
	int col_shift = 0;
	int index = 0;

	double *pTemp = (double *)malloc(sizeof(double)*(2*row-1+filterLen-1)*(2*col-1+filterLen-1));
	memset(pTemp, 0, sizeof(double)*(2*row-1+filterLen-1)*(2*col-1+filterLen-1));

	IDWT_Partial(pLow,  row, col, orig_row, orig_col, pLoFilter, pLoFilter, filterLen, pTemp);
	IDWT_Partial(pHor,  row, col, orig_row, orig_col, pHiFilter, pLoFilter, filterLen, pTemp);
	IDWT_Partial(pVer,  row, col, orig_row, orig_col, pLoFilter, pHiFilter, filterLen, pTemp);
	IDWT_Partial(pDiag, row, col, orig_row, orig_col, pHiFilter, pHiFilter, filterLen, pTemp);

	ShowResultInt(pTemp, 2*row-1+filterLen-1, 2*col-1+filterLen-1);

	if ((orig_row + filterLen - 1)%2 == 0)
	{
		row_shift = filterLen - 2;
	}
	else
	{
		if (row*2 < (orig_row + filterLen - 1))
		{
			row_shift = filterLen - 2;
		}
		else
		{
			row_shift = filterLen - 1;
		}
	}

	if ((orig_col + filterLen - 1)%2 == 0)
	{
		col_shift = filterLen - 2;
	}
	else
	{
		if (col*2 < (orig_col + filterLen - 1))
		{
			col_shift = filterLen - 2;
		}
		else
		{
			col_shift = filterLen - 1;
		}
	}
	col = 2*col-1+filterLen-1;
	index = col*row_shift + col_shift;

	//Keep the center part as result
	for (int i=0; i<orig_row; i++)
	{
		for (int j=0; j<orig_col; j++)
		{
			*(pResult+i*orig_col+j) = *(pResult+i*orig_col+j) + *(pTemp+index);
			index++;
		}
		index = index + (col-orig_col);
	}

	free(pTemp);
}

//Perform discrete 2 dimensinal wavelet transform
void WaveletTransform2D(double *pSrc, int row, int col, double *pLoFilter, double *pHiFilter, int filterLen, double *pLow, double *pVer, double *pHor, double *pDiag)
{
	int nIndex = 0;
	int coeff_row = row + filterLen -1;
	int coeff_col = col + filterLen - 1;
	
	double *pTemp1 = (double *)malloc(sizeof(double)*row*coeff_col);
	double *pTemp2 = (double *)malloc(sizeof(double)*coeff_row*coeff_col);
	
	RowConvolution2D(pSrc, pLoFilter, filterLen, row, col, pTemp1, false, false);
	
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

	RowConvolution2D(pSrc, pHiFilter, filterLen, row, col, pTemp1, false, false);
	//Get vertical coefficients
	ColConvolution2D(pTemp1, pLoFilter, filterLen, row, col+filterLen-1, pTemp2, false);
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
