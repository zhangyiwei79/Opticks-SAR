#include "ImageSegmentation.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.1415926

#define X(ix,iy) (ix)*iNy+ (iy)
#define Xe(ix,iy) (ix)*iNye+ (iy)

#define SQR(x) (x)*(x)

#ifndef HAVE_RINT 
#define rint(A) floor((A)+(((A) < 0)? -0.5 : 0.5)) 
#endif



/****************************************/
/* This function computes the texture feature */
void vComputeTextureFeature(float *pfIm0, int iNx, int iNy, float *pfTextFeat, int ip, int iNbBin)
{
    int    ix, iy, iX, iNxe, iNye, ip2, ixe, iye, ik;
    int    iyp2, ixp2, iy2, ix2, iXe, iNyxe;
    float  *pfI, *pfIx, *pfIy, *pfIxe, *pfIye, fG11, fG22, fG12;
    float  fMaxI, fMinI, fNbBin;
    
    
    
    /* Allocate memory */
    pfI = (float *) calloc( (unsigned)(iNx*iNy), sizeof(float) );
 
    /* Normalize Im0 */
    fMaxI = 0.0;
    fMinI = 1e10;
    for (ix=0; ix< iNx; ix++)
        for (iy=0; iy< iNy; iy++)
    {
        if (pfIm0[X(ix,iy)]>fMaxI) fMaxI= pfIm0[X(ix,iy)];
        if (pfIm0[X(ix,iy)]<fMinI) fMinI= pfIm0[X(ix,iy)];
        }
    for (ix=0; ix< iNx; ix++)
        for (iy=0; iy< iNy; iy++)
            pfI[X(ix,iy)] = (pfIm0[X(ix,iy)]-fMinI)/ (fMaxI-fMinI);
    
    
    /* Allocate memory */
    pfIx = (float *) calloc( (unsigned)(iNx*iNy), sizeof(float) );  
    pfIy = (float *) calloc( (unsigned)(iNx*iNy), sizeof(float) );

    
    /* Compute partial_x I=Ix, partial_y I=Iy */
    for (ix=0; ix< iNx-1; ix++)
        for (iy=0; iy< iNy-1; iy++)
    {
        iX = X(ix,iy);
        pfIx[iX] = pfI[X(ix+1,iy)] - pfI[iX];
        pfIy[iX] = pfI[X(ix,iy+1)] - pfI[iX];
        }
    ix = iNx-1;
    for (iy=0; iy< iNy-1; iy++)
    {
        iX = X(ix,iy);
        pfIx[iX] = pfI[iX] - pfI[X(ix-1,iy)];
        pfIy[iX] = pfI[X(ix,iy+1)] - pfI[iX];
    }
    iy = iNy-1;
    for (ix=0; ix< iNx-1; ix++)
    {
        iX = X(ix,iy);
        pfIx[iX] = pfI[X(ix+1,iy)] - pfI[iX];
        pfIy[iX] = pfI[iX] - pfI[X(ix,iy-1)];
    }
    ix = iNx-1;
    iy = iNy-1;
    iX = X(ix,iy);
    pfIx[iX] = pfI[iX] - pfI[X(ix-1,iy)];
    pfIy[iX] = pfI[iX] - pfI[X(ix,iy-1)];
    
    
    
    /* Compute extended Ix, Iy with mirror conditions */
    ip2 = (ip-1)/2;
    iNxe = iNx + 2* ip2;
    iNye = iNy + 2* ip2;
    iNyxe = iNxe*iNye;
    
    
    /* Allocate memory */
    pfIxe = (float *) calloc( (unsigned)(iNxe*iNye), sizeof(float) );
    pfIye = (float *) calloc( (unsigned)(iNxe*iNye), sizeof(float) );
    
    for (ixe=0,ik=0;ixe<iNxe;ixe++)
    {
        if (ixe<ip2) ix=ip2-ixe;
        else if (ixe>iNx+ip2-1) ix=2*iNx+ip2-ixe-2;
        else ix=ixe-ip2;
        
        for (iye=0;iye<iNye;iye++,ik++)
        {
            if (iye<ip2)  iy=ip2-iye;
            else if (iye>iNy+ip2-1) iy=2*iNy+ip2-iye-2;
            else iy=iye-ip2;
            
            pfIxe[ik]=pfIx[X(ix,iy)];
            pfIye[ik]=pfIy[X(ix,iy)];
        }
    }
    
    
    /* Compute Texture Feature */
    for (iye=ip2;iye<iNye-ip2;iye++)
        for (ixe=ip2;ixe<iNxe-ip2;ixe++)
    {
        fG11 = 0.0;
        fG22 = 0.0;
        fG12 = 0.0;
        for (iyp2=-ip2;iyp2<ip2;iyp2++)
            for (ixp2=-ip2;ixp2<ip2;ixp2++)
        {
            iy2 = iye+iyp2;
            ix2 = ixe+ixp2;
            iXe = Xe(ix2,iy2);
            fG11 += SQR(pfIxe[iXe]);
            fG22 += SQR(pfIye[iXe]);
            fG12 += pfIxe[iXe]*pfIye[iXe];
            }
        
        iy=iye-ip2;
        ix=ixe-ip2;
        pfTextFeat[X(ix,iy)] = 1.0/( (1.0+fG11)*(1.0+fG22) - SQR(fG12) );
        }
    
    
    /* Normalize Texture Feature */
    fMaxI = 0.0;
    fMinI = 1e10;
    for (ix=0; ix< iNx; ix++)
        for (iy=0; iy< iNy; iy++)
    {
        if (pfTextFeat[X(ix,iy)]>fMaxI) fMaxI= pfTextFeat[X(ix,iy)];
        if (pfTextFeat[X(ix,iy)]<fMinI) fMinI= pfTextFeat[X(ix,iy)];
        }
    
    /* Range Texture Feature is in [1 NbBin] */
    fNbBin = (float)iNbBin;
    for (ix=0; ix< iNx; ix++)
        for (iy=0; iy< iNy; iy++)
            pfTextFeat[X(ix,iy)] = 1.0 + (pfTextFeat[X(ix,iy)]-fMinI)/ (fMaxI-fMinI)* fNbBin;
    
    
    
    /* Free memory */
    free( (float *) pfI );
    free( (float *) pfIx );
    free( (float *) pfIy );
    free( (float *) pfIxe );
    free( (float *) pfIye );
    
}
/****************************************/

/****************************************/
/* This function computes a Gaussian function */
void vComputeGaussian(float *pfGaussian, float fStdParzenWindow, int iNGaussian)
{
    int    ic, i;
    float  f1, f2;
    
    ic = (int) ((iNGaussian-1)/2);
    
    f1 = 2.0* SQR(fStdParzenWindow);
    f2 = fStdParzenWindow* sqrt(2.0* PI);
    
    for (i=0; i< iNGaussian; i++)
        pfGaussian[i] = (float)exp( -SQR(i-ic)/ f1 ) / f2;
    
}
/****************************************/

/****************************************/
/* This function estimates the intensity distribution from an histogram */
void vComputeParzenWindFunction_p(float *pfp, float *pfHisto, float fAera, int iNbBin, float *pfGaussian, int iNGaussian)
{
  int    iNg, iNg2, ic, i, i2, i_shift;
  float  fSumGauss;

  iNg = iNGaussian;
  iNg2 = (int) ((iNGaussian-1)/2);
  ic = (int) ((iNGaussian-1)/2);
  
  
  /* Pre-processing to avoid division by zero */
  for (i=0; i< iNbBin; i++)
      if( pfHisto[i] < 1 )
          pfHisto[i] = 1;
  
  
  /* Center of p */
  for (i=iNg2; i< iNbBin- iNg2; i++)
  {
      for (i2=0; i2< iNg; i2++)
          pfp[i] += pfHisto[i+(i2-ic)]* pfGaussian[i2];
      pfp[i] /= fAera;
  }
  
  
  /* Borders of p */
  for (i=0; i< iNg2; i++)
  {
      fSumGauss = 0.0;
      for (i2=0; i2< iNg; i2++)
      {
          i_shift = i+(i2-ic);
          if ( i_shift>= 0  &&  i_shift< iNbBin )
              pfp[i] += pfHisto[i_shift]* pfGaussian[i2];
          fSumGauss += pfGaussian[i2];
      }
      pfp[i] /= fAera;
      pfp[i] /= fSumGauss;
  }
  
  for (i=iNbBin- iNg2; i< iNbBin; i++)
  {
      fSumGauss = 0.0;
      for (i2=0; i2< iNg; i2++)
      {
          i_shift = i+(i2-ic);
          if ( i_shift>= 0  &&  i_shift< iNbBin )
              pfp[i] += pfHisto[i_shift]* pfGaussian[i2];
          fSumGauss += pfGaussian[i2];
      }
      pfp[i] /= fAera;
      pfp[i] /= fSumGauss;
  }

  
}
/****************************************/

/****************************************/
void vComputeParzenWindFunction_E(float *pfpIn, float *pfLogpIn, float *pfpOut, float *pfLogpOut, float *pfE1, float *pfE2, int iNbBin, float *pfGaussian, int iNGaussian)
{
    int    iNg, iNg2, ic, i, i2, i_shift;
    float  fSumGauss;
    
    
    iNg = iNGaussian;
    iNg2 = (int) ((iNGaussian-1)/2);
    ic = (int) ((iNGaussian-1)/2);
    
    
    /* Init */
    for (i=0; i< iNbBin; i++) { pfE1[i] = 0.0; pfE2[i] = 0.0;  }
    
    
    /* Center of E1 and E2 */
    for (i=iNg2; i< iNbBin- iNg2; i++){
        for (i2=0; i2< iNg; i2++)
        {
            pfE1[i] += ( pfLogpIn[i+(i2-ic)] - pfLogpOut[i+(i2-ic)] + 1.0 - pfpOut[i+(i2-ic)]/pfpIn[i+(i2-ic)])* pfGaussian[i2];
            pfE2[i] += ( pfLogpIn[i+(i2-ic)] - pfLogpOut[i+(i2-ic)] - 1.0 + pfpIn[i+(i2-ic)]/pfpOut[i+(i2-ic)])* pfGaussian[i2];
        }
    }
    
    
    /* Boders of E1 and E2 */
    for (i=0; i< iNg2; i++)
    {
        fSumGauss = 0.0;
        for (i2=0; i2< iNg; i2++)
        {
            i_shift = i+(i2-ic);
            if ( i_shift>= 0  &&  i_shift< iNbBin )
            {
                pfE1[i] += ( pfLogpIn[i_shift] - pfLogpOut[i_shift] + 1.0 - pfpOut[i_shift]/pfpIn[i_shift])* pfGaussian[i2];
                pfE2[i] += ( pfLogpIn[i_shift] - pfLogpOut[i_shift] - 1.0 + pfpIn[i_shift]/pfpOut[i_shift])* pfGaussian[i2];
                fSumGauss += pfGaussian[i2];
            }
        }
        pfE1[i] /= fSumGauss;
        pfE2[i] /= fSumGauss;
    }
    
    for (i= iNbBin- iNg2; i< iNbBin; i++)
    {
        fSumGauss = 0.0;
        for (i2=0; i2< iNg; i2++)
        {
            i_shift = i+(i2-ic);
            if ( i_shift>= 0  &&  i_shift< iNbBin )
            {
                pfE1[i] += ( pfLogpIn[i_shift] - pfLogpOut[i_shift] + 1.0 - pfpOut[i_shift]/pfpIn[i_shift])* pfGaussian[i2];
                pfE2[i] += ( pfLogpIn[i_shift] - pfLogpOut[i_shift] - 1.0 + pfpIn[i_shift]/pfpOut[i_shift])* pfGaussian[i2];
                fSumGauss += pfGaussian[i2];
            }
        }
        pfE1[i] /= fSumGauss;
        pfE2[i] /= fSumGauss;
    }
    
    
}
/****************************************/

/****************************************/
void vComputeParzenWindFunction_F(float *pfpIn, float *pfLogpIn, float *pfpOut, float *pfLogpOut, float *pfF1, float *pfF2, int iNbBin, float *pfGaussian, int iNGaussian)
{
    int    iNg, iNg2, ic, i;
    float  fF1, fF2;
    
    iNg = iNGaussian;
    iNg2 = (int) ((iNGaussian-1)/2);
    ic = (int) ((iNGaussian-1)/2);
    
    fF1 = 0.0;
    fF2 = 0.0;
    
    for (i=0; i< iNbBin; i++)
    {
        fF1 += (pfpIn[i]*pfLogpIn[i] - pfpIn[i]*pfLogpOut[i] + pfpIn[i] - pfpOut[i]);
        fF2 += (pfpOut[i]*pfLogpIn[i] - pfpOut[i]*pfLogpOut[i] - pfpOut[i] + pfpIn[i]);
    }
    
    *pfF1 = fF1;
    *pfF2 = fF2;
    
}
/****************************************/


/****************************************/
void vComputeKL (float *pfIm0, float *pfu, int iNx, int iNy, float *pfHr, int iNbItersUpdateHr,
                 float *pfpIn, float *pfpOut, float *pfHistoIn, float *pfHistoOut, float *pfLogpIn, 
                 float *pfLogpOut, int iNI, int iNbBin, float *pfGaussian, int iNGaussian, 
                 float *pfE1, float *pfE2)
{
    float  fNormalizationIn, fNormalizationOut;
    float  fEps, fLogEps, fF1, fF2, fEvolIn, fEvolOut;
    int    ix, iy, i, iI;
    
    
    /* Update the region function Hr every "iNbItersUpdateHr" iterations */
    if ( (iNI==0) || (iNI%iNbItersUpdateHr==0) )
    {
        
        /* Init */
        for (i=0; i< iNbBin; i++) { pfHistoIn[i] = 0.0; pfHistoOut[i] = 0.0; }
        fNormalizationIn = 0.0;
        fNormalizationOut = 0.0;
        
        
        /* Compute histograms inside and outside the boundary of {u>0.5} */
        for ( ix=0; ix< iNx; ix++)
            for (iy=0; iy< iNy; iy++)
        {
            iI = (int)(rint(pfIm0[X(ix,iy)]));
            if (iI<=0) iI=0;
            if (pfu[X(ix,iy)] >= 0.5 )
                pfHistoIn[iI] += 1.0;
            else
                pfHistoOut[iI] += 1.0;

            fNormalizationIn += pfu[X(ix,iy)];
            fNormalizationOut += 1.0-pfu[X(ix,iy)];
            }
        
        
        /* Estimate the intensity distributions inside and outside from histograms, using the Parzen estimation method */
        vComputeParzenWindFunction_p(pfpIn,pfHistoIn,fNormalizationIn,iNbBin,pfGaussian,iNGaussian);
        vComputeParzenWindFunction_p(pfpOut,pfHistoOut,fNormalizationOut,iNbBin,pfGaussian,iNGaussian);
        
        
        /* Compute the log of distributions */
        fEps = 1e-32;
        fLogEps = log10(fEps);
        for (i=0; i< iNbBin; i++)
        {
            if ( pfpIn[i] > fEps ) pfLogpIn[i] = log10(pfpIn[i]);
            else  pfLogpIn[i] = fLogEps;
            if ( pfpOut[i] > fEps ) pfLogpOut[i] = log10(pfpOut[i]);
            else  pfLogpOut[i] = fLogEps;
        }
        
        
        vComputeParzenWindFunction_E(pfpIn,pfLogpIn,pfpOut,pfLogpOut,pfE1,pfE2,iNbBin,pfGaussian,iNGaussian);
        vComputeParzenWindFunction_F(pfpIn,pfLogpIn,pfpOut,pfLogpOut,&fF1,&fF2,iNbBin,pfGaussian,iNGaussian);
        
        
        for (ix=0; ix< iNx; ix++)
            for (iy=0; iy< iNy; iy++)
        {
            iI =(int)rint(pfIm0[X(ix,iy)]);
            if (iI<=0) iI=0;
            if (fNormalizationIn>1e-4) fEvolIn = (pfE1[iI] - fF1)/ fNormalizationIn; else fEvolIn = 0.0;
            if (fNormalizationOut>1e-4) fEvolOut =(pfE2[iI] - fF2)/fNormalizationOut; else fEvolOut = 0.0;
            pfHr[X(ix,iy)] = -(fEvolIn + fEvolOut);
            }
        
        
    }
    
    
}

/**********************************************/
/************** END SUB FUNCTIONS *************/
/**********************************************/




/**********************************************/
/************** MAIN FUNCTION *****************/
/**********************************************/

/****************************************/
void ImageSegmentation(float *pOrigIm, float *pfu, float *pfVecParameters)
{
    float   *pfdx, *pfdy, *pfuOld;
    float   *pfbx, *pfby, *pfIm0 = pOrigIm;
    float   fLambda, fMu, fct1, fct2, fctST, fG, fDxu, fDyu, fs, fTemp, fSumDiff;
    float   fct1b, fct2b, fct1c, fct2c, fInvMu, fInvMu2, f1, f2;
    float   *pfHr, fMaxImRef, *pfpIn, *pfpOut, *pfHistoIn, *pfHistoOut, *pfLogpIn, *pfLogpOut;
    float   *pfE1, *pfE2, *pfGaussian, *pfuOld2, *pfTextFeat;
    float   fStdParzenWindow, fRangeHr, fMinHr, fMaxHr, fNyx;
    float   fDiffNew, fDiffOld, fDiffNew2, fSumU, fError;
    float   fSumUold, fDiffFirst2, fStopThres;
    int     iNy, iNx, ix, iy;
    int     iNI = 0, iX, iGS, nIterations, iNbBin, iNGaussian;
    int     iNbItersUpdateHr, iMeanGS, iCptGS, ip;
  
    /* Size */
    iNy = (int) pfVecParameters[0];
    iNx = (int) pfVecParameters[1];   
    
    /* Choice of region segmentation model */
    nIterations = (int) pfVecParameters[2]; /* Maximum number of loops */
    iNbBin = (int) pfVecParameters[5]; /* Number of bins for histograms */
    fStdParzenWindow = pfVecParameters[6]; /* standard deviation in Parzen estimation */
    ip = (int)pfVecParameters[7]; /* patch size */
    
    /* Memory allocation */
    pfdx = (float *) calloc( (unsigned)(iNy*iNx), sizeof(float) );
    pfdy = (float *) calloc( (unsigned)(iNy*iNx), sizeof(float) ); 
    pfbx = (float *) calloc( (unsigned)(iNy*iNx), sizeof(float) ); 
    pfby = (float *) calloc( (unsigned)(iNy*iNx), sizeof(float) ); 
    
    pfuOld = (float *) calloc( (unsigned)(iNy*iNx), sizeof(float) );  
    pfHr = (float *) calloc( (unsigned)(iNy*iNx), sizeof(float) ); 
    pfuOld2 = (float *) calloc( (unsigned)(iNy*iNx), sizeof(float) );  
    

    iNGaussian = (int)rint(6.0* fStdParzenWindow);
    iNGaussian = (iNGaussian%2 == 1) ? iNGaussian : iNGaussian+1;
    pfGaussian = (float *) calloc( (unsigned)(iNGaussian), sizeof(float) );
 
    pfpIn = (float *) calloc( (unsigned)(iNbBin*2), sizeof(float) );
    pfpOut = (float *) calloc( (unsigned)(iNbBin*2), sizeof(float) );
    
    pfHistoIn = (float *) calloc( (unsigned)(iNbBin*2), sizeof(float) );
    pfHistoOut = (float *) calloc( (unsigned)(iNbBin*2), sizeof(float) );

    pfLogpIn = (float *) calloc( (unsigned)(iNbBin*2), sizeof(float) );
    pfLogpOut = (float *) calloc( (unsigned)(iNbBin*2), sizeof(float) );

    pfE1 = (float *) calloc( (unsigned)(iNbBin*2), sizeof(float) );
    pfE2 = (float *) calloc( (unsigned)(iNbBin*2), sizeof(float) );
    pfTextFeat = (float *) calloc( (unsigned)(iNy*iNx), sizeof(float) );

    /* Compute Gaussian kernel for the KL model */
    vComputeGaussian(pfGaussian,fStdParzenWindow,iNGaussian);
        
    /* Compute Texture Feature */
    vComputeTextureFeature((float *)pfIm0,iNx,iNy,pfTextFeat,ip,iNbBin);
    pfIm0 = pfTextFeat;

    /* Compute init function u (any non-zero function is fine) */
    fMaxImRef = 0.0;
    for (ix=0; ix< iNx; ix++)
        for (iy=0; iy< iNy; iy++)
            if (abs(pfIm0[X(ix,iy)])>fMaxImRef) fMaxImRef= abs(pfIm0[X(ix,iy)]);
    for (ix=0; ix< iNx; ix++)
        for (iy=0; iy< iNy; iy++)
            pfu[X(ix,iy)] = pfIm0[X(ix,iy)]/fMaxImRef;
      
    /* Parameters for the segmentation code */
    fLambda = pfVecParameters[3];
    fMu = pfVecParameters[4];
    

    /* Normalize lambda value for the computation of Gauss-Seidel iterations */
    /* Estimate of g_r */
	iNbItersUpdateHr = 1;

    vComputeKL(pfIm0,pfu,iNx,iNy,pfHr,iNbItersUpdateHr,pfpIn,pfpOut,pfHistoIn,pfHistoOut,
               pfLogpIn,pfLogpOut,0,iNbBin,pfGaussian,iNGaussian,pfE1,pfE2);

    /* Estimate of the range of g_r */
    fMinHr = 1e10;
    fMaxHr = -1e10;
    for (ix=1; ix< iNx-1; ix++)
        for (iy=1; iy< iNy-1; iy++)
    {
        if ( pfHr[X(ix,iy)]>fMaxHr ) fMaxHr = pfHr[X(ix,iy)];
        if ( pfHr[X(ix,iy)]<fMinHr ) fMinHr = pfHr[X(ix,iy)];
        }
    fRangeHr = fMaxHr-fMinHr;
        
    /* Normalize lambda with respect to the range of hr */
    fLambda /= fRangeHr;     
    
    
    /* Constants */
    iNbItersUpdateHr = 1;  /* number of iterations to update the region function g_r */
    
    fInvMu = 1./ fMu;
    fInvMu2 = SQR(fInvMu);
    
    fct1 = 1./4.;
    fct2 = fLambda/(4.0*fMu);
    fct1b = 1./3.;
    fct2b = fLambda/(3.0*fMu);
    fct1c = 1./2.;
    fct2c = fLambda/(2.0*fMu);

    fNyx = (float)(iNy*iNx);
    
    fStopThres = 1e-10;       
    
    
    // Iterative minimization scheme 
    fDiffOld = 1e10; fDiffNew = 1e11;
    iMeanGS = 0; iCptGS = 0;
    iNI=0; /* number of iterations (outer iterations) */
    while ( abs(fDiffNew-fDiffOld)>fStopThres && iNI<nIterations ) 
    {
        
        /* Store u^old for outer iterations */
        for (ix=0; ix< iNx; ix++)
            for (iy=0; iy< iNy; iy++)
                pfuOld[X(ix,iy)] = pfu[X(ix,iy)];
        
        
        /* Update region function hr */
        vComputeKL(pfIm0,pfu,iNx,iNy,pfHr,iNbItersUpdateHr,pfpIn,pfpOut,pfHistoIn,pfHistoOut,
                   pfLogpIn,pfLogpOut,iNI,iNbBin,pfGaussian,iNGaussian,pfE1,pfE2);
      
        /* Compute u^{k+1} with Gauss-Seidel */
        /* Solve u^k+1 = arg min int lambda h_r u + mu/2 |d - grad u - b^k|^2 */
        /* Euler-Lagrange is  mu Laplacian u = lambda hr + mu div (b^k-d^k), u in [0,1] */
        iGS=0; /* number of iterations for Gauss-Seidel (inner iterations) */
        fError = 1e10;
        while ( fError>1e-2 && iGS<50 )  
        {
            
            /* Store u^old for inner iterations (Gauss-Seidel) */
            for (ix=0; ix< iNx; ix++)
                for (iy=0; iy< iNy; iy++)
                    pfuOld2[X(ix,iy)] = pfu[X(ix,iy)];
            
            /* Center */
            for (ix=1; ix< iNx-1; ix++)
                for (iy=1; iy< iNy-1; iy++)
            {
                iX = X(ix,iy);
                fG = pfu[X(ix+1,iy)] + pfu[X(ix-1,iy)] + pfu[X(ix,iy+1)] + pfu[X(ix,iy-1)];
                fG += pfdx[X(ix-1,iy)] - pfdx[iX] - pfbx[X(ix-1,iy)] + pfbx[iX];
                fG += pfdy[X(ix,iy-1)] - pfdy[iX] - pfby[X(ix,iy-1)] + pfby[iX];
                fG *= fct1;
                fG -= fct2* pfHr[iX];
                if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                pfu[iX] = fG;
                }
            
            /* Borders */
            ix=0;
            for (iy=1; iy< iNy-1; iy++)
            {
                iX = X(ix,iy);
                fG = pfu[X(ix+1,iy)] + pfu[X(ix,iy+1)] + pfu[X(ix,iy-1)];
                fG += - pfdx[iX] + pfbx[iX];
                fG += pfdy[X(ix,iy-1)] - pfdy[iX] - pfby[X(ix,iy-1)] + pfby[iX];
                fG *= fct1b;
                fG -= fct2b* pfHr[iX];
                if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                pfu[iX] = fG;
            }
            
            ix=iNx-1;
            for (iy=1; iy< iNy-1; iy++)
            {
                iX = X(ix,iy);
                fG = pfu[X(ix-1,iy)] + pfu[X(ix,iy+1)] + pfu[X(ix,iy-1)];
                fG += pfdx[X(ix-1,iy)] - pfdx[iX] - pfbx[X(ix-1,iy)] + pfbx[iX];
                fG += pfdy[X(ix,iy-1)] - pfdy[iX] - pfby[X(ix,iy-1)] + pfby[iX];
                fG *= fct1b;
                fG -= fct2b* pfHr[iX];
                if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                pfu[iX] = fG;
            }
            
            iy=0;
            for (ix=1; ix< iNx-1; ix++)
            {
                iX = X(ix,iy);
                fG = pfu[X(ix+1,iy)] + pfu[X(ix-1,iy)] + pfu[X(ix,iy+1)];
                fG += pfdx[X(ix-1,iy)] - pfdx[iX] - pfbx[X(ix-1,iy)] + pfbx[iX];
                fG += - pfdy[iX] + pfby[iX];
                fG *= fct1b;
                fG -= fct2b* pfHr[iX];
                if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                pfu[iX] = fG;
            }
            
            iy=iNy-1;
            for (ix=1; ix< iNx-1; ix++)
            {
                iX = X(ix,iy);
                fG = pfu[X(ix+1,iy)] + pfu[X(ix-1,iy)] + pfu[X(ix,iy-1)];
                fG += pfdx[X(ix-1,iy)] - pfdx[iX] - pfbx[X(ix-1,iy)] + pfbx[iX];
                fG += pfdy[X(ix,iy-1)] - pfdy[iX] - pfby[X(ix,iy-1)] + pfby[iX];
                fG *= fct1b;
                fG -= fct2b* pfHr[iX];
                if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                pfu[iX] = fG;
            }
            
            ix=0; iy=0;
            iX = X(ix,iy);
            fG = pfu[X(ix+1,iy)] + pfu[X(ix,iy+1)];
            fG += - pfdx[iX] + pfbx[iX];
            fG += - pfdy[iX] + pfby[iX];
            fG *= fct1c;
            fG -= fct2c* pfHr[iX];
            if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
            pfu[iX] = fG;
            
            ix=iNx-1; iy=0;
            iX = X(ix,iy);
            fG = pfu[X(ix-1,iy)] + pfu[X(ix,iy+1)];
            fG += pfdx[X(ix-1,iy)] - pfdx[iX] - pfbx[X(ix-1,iy)] + pfbx[iX];
            fG += - pfdy[iX] + pfby[iX];
            fG *= fct1c;
            fG -= fct2c* pfHr[iX];
            if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
            pfu[iX] = fG;
            
            ix=0; iy=iNy-1;
            iX = X(ix,iy);
            fG = pfu[X(ix-1,iy)] + pfu[X(ix,iy+1)];
            fG += pfdx[X(ix-1,iy)] - pfdx[iX] - pfbx[X(ix-1,iy)] + pfbx[iX];
            fG += - pfdy[iX] + pfby[iX];
            fG *= fct1c;
            fG -= fct2c* pfHr[iX];
            if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
            pfu[iX] = fG;
            
            ix=iNx-1; iy=iNy-1;
            iX = X(ix,iy);
            fG = pfu[X(ix-1,iy)] + pfu[X(ix,iy-1)];
            fG += pfdx[X(ix-1,iy)] - pfdx[iX] - pfbx[X(ix-1,iy)] + pfbx[iX];
            fG += pfdy[X(ix,iy-1)] - pfdy[iX] - pfby[X(ix,iy-1)] + pfby[iX];
            fG *= fct1c;
            fG -= fct2c* pfHr[iX];
            if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
            pfu[iX] = fG;
            /* end Borders */
            
            
            /* Compute diff ( u - uold ) */
            fSumDiff = 0.0;
            fSumU = 0.0;
            fSumUold = 0.0;
            for (ix=0; ix< iNx; ix++)
                for (iy=0; iy< iNy; iy++)
                    fSumDiff += SQR(pfu[X(ix,iy)]-pfuOld2[X(ix,iy)]);
            fDiffNew2 = fSumDiff/ fNyx;
            if ( iGS==0 ) 
            {
                fDiffFirst2 = fDiffNew2;
                fDiffNew2 = 1e10;
                fError = 1e10;
            }
            else
                fError = 1.0 - abs(fDiffNew2-fDiffFirst2)/fDiffFirst2;
            iGS++;
            
        }
        
        
        iMeanGS += iGS;
        iCptGS++;
        

        /* Compute d^{k+1} (Soft-Thresholding) and b^{k+1} (Bregman function) */
        /* d^k+1 = arg min int g_b |d| + mu/2 |d - grad u - b^k|^2 */
        /* d^k+1 = (grad u^k+1 + b^k)/ |grad u^k+1 + b^k| max(|grad u^k+1 + b^k|-1/mu,0) */
        /* b^k+1 = b^k + grad u^k+1 - d^k+1 */
        /* Center */
        for (ix=0; ix< iNx-1; ix++)
            for (iy=0; iy< iNy-1; iy++)
        {
            iX = X(ix,iy);
            /* d */
            fDxu = pfu[X(ix+1,iy)] - pfu[iX];
            fDyu = pfu[X(ix,iy+1)] - pfu[iX];
            f1 = fDxu+pfbx[iX];
            f2 = fDyu+pfby[iX];
            fs = SQR(f1)+SQR(f2);
            fctST = fInvMu2;
            if ( fs<fctST ) { pfdx[iX]=0.0; pfdy[iX]=0.0; }
            else {
                fs = sqrt(fs);
                fctST = sqrt(fctST);
                fTemp = fs-fctST; fTemp /= fs;
                pfdx[iX] = fTemp* f1;
                pfdy[iX] = fTemp* f2; }
            /* b */
            pfbx[iX] += fDxu - pfdx[iX];
            pfby[iX] += fDyu - pfdy[iX];
            }
        
        /* Borders */
        ix=iNx-1;
        for (iy=1; iy< iNy-1; iy++)
        {
            iX = X(ix,iy);
            fDxu = 0.0;
            fDyu = pfu[X(ix,iy+1)] - pfu[iX];
            f1 = fDxu+pfbx[iX]; f2 = fDyu+pfby[iX];
            fs = SQR(f1)+SQR(f2); fctST = fInvMu2;
            if ( fs<fctST ) { pfdx[iX]=0.0; pfdy[iX]=0.0; }
            else {
                fs = sqrt(fs); fctST = sqrt(fctST);
                fTemp = fs-fctST; fTemp /= fs;
                pfdx[iX] = fTemp* f1;
                pfdy[iX] = fTemp* f2; }
            pfbx[iX] += fDxu - pfdx[iX];
            pfby[iX] += fDyu - pfdy[iX];
        }
        
        iy=iNy-1;
        for (ix=1; ix< iNx-1; ix++)
        {
            iX = X(ix,iy);
            fDxu = pfu[X(ix+1,iy)] - pfu[iX];
            fDyu = 0.0;
            f1 = fDxu+pfbx[iX]; f2 = fDyu+pfby[iX];
            fs = SQR(f1)+SQR(f2); fctST = fInvMu2;
            if ( fs<fctST ) { pfdx[iX]=0.0; pfdy[iX]=0.0; }
            else {
                fs = sqrt(fs); fctST = sqrt(fctST);
                fTemp = fs-fctST; fTemp /= fs;
                pfdx[iX] = fTemp* f1;
                pfdy[iX] = fTemp* f2; }
            pfbx[iX] += fDxu - pfdx[iX];
            pfby[iX] += fDyu - pfdy[iX];
        }
        
        ix=iNx-1; iy=0;
        iX = X(ix,iy);
        fDxu = 0.0;
        fDyu = pfu[X(ix,iy+1)] - pfu[iX];
        f1 = fDxu+pfbx[iX]; f2 = fDyu+pfby[iX];
        fs = SQR(f1)+SQR(f2); fctST = fInvMu2;
        if ( fs<fctST ) { pfdx[iX]=0.0; pfdy[iX]=0.0; }
        else {
            fs = sqrt(fs); fctST = sqrt(fctST);
            fTemp = fs-fctST; fTemp /= fs;
            pfdx[iX] = fTemp* f1;
            pfdy[iX] = fTemp* f2; }
        pfbx[iX] += fDxu - pfdx[iX];
        pfby[iX] += fDyu - pfdy[iX];
        
        ix=0; iy=iNy-1;
        iX = X(ix,iy);
        fDxu = pfu[X(ix+1,iy)] - pfu[iX];
        fDyu = 0.0;
        f1 = fDxu+pfbx[iX]; f2 = fDyu+pfby[iX];
        fs = SQR(f1)+SQR(f2); fctST = fInvMu2;
        if ( fs<fctST ) { pfdx[iX]=0.0; pfdy[iX]=0.0; }
        else {
            fs = sqrt(fs); fctST = sqrt(fctST);
            fTemp = fs-fctST; fTemp /= fs;
            pfdx[iX] = fTemp* f1;
            pfdy[iX] = fTemp* f2; }
        pfbx[iX] += fDxu - pfdx[iX];
        pfby[iX] += fDyu - pfdy[iX];
        
        ix=iNx-1; iy=iNy-1;
        iX = X(ix,iy);
        fDxu = 0.0;
        fDyu = 0.0;
        f1 = fDxu+pfbx[iX]; f2 = fDyu+pfby[iX];
        fs = SQR(f1)+SQR(f2); fctST = fInvMu2;
        if ( fs<fctST ) { pfdx[iX]=0.0; pfdy[iX]=0.0; }
        else {
            fs = sqrt(fs); fctST = sqrt(fctST);
            fTemp = fs-fctST; fTemp /= fs;
            pfdx[iX] = fTemp* f1;
            pfdy[iX] = fTemp* f2; }
        pfbx[iX] += fDxu - pfdx[iX];
        pfby[iX] += fDyu - pfdy[iX];
        /* Borders */
        
        
        
        fSumDiff = 0.0;
        fSumU = 0.0;
        fSumUold = 0.0;
        for (ix=0; ix< iNx; ix++)
            for (iy=0; iy< iNy; iy++)
        {
            fSumDiff += SQR(pfu[X(ix,iy)]-pfuOld[X(ix,iy)]);
            fSumU += SQR(pfu[X(ix,iy)]);
            fSumUold += SQR(pfuOld[X(ix,iy)]);
            }
        fDiffOld = fDiffNew;
        fDiffNew = fSumDiff/ (fSumU*fSumUold);
        iNI++;

    }

       
    /* Free memory */
    free( (float *) pfdx );
    free( (float *) pfdy );
    free( (float *) pfbx );
    free( (float *) pfby );
    free( (float *) pfHr );
    free( (float *) pfuOld );
    free( (float *) pfuOld2 );

    free( (float *) pfpIn );
    free( (float *) pfpOut );
    free( (float *) pfHistoIn );
    free( (float *) pfHistoOut );
    free( (float *) pfLogpIn );
    free( (float *) pfLogpOut );
    free( (float *) pfE1 );
    free( (float *) pfE2 );
    free( (float *) pfGaussian );
    free( (float *) pfTextFeat );
    
}
/****************************************/


/**********************************************/
/************** END MAIN FUNCTION *************/
/**********************************************/