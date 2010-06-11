/*
 * The information in this file is
 * Copyright(c) 2007 Ball Aerospace & Technologies Corporation
 * and is subject to the terms and conditions of the
 * GNU Lesser General Public License Version 2.1
 * The license text is available from   
 * http://www.gnu.org/licenses/lgpl.html
 */

/*
The following is aid to remove the speckle noise using refined Lee filter
based on local statistcs
*/
#include "DataAccessor.h"
#include "DataAccessorImpl.h"
#include "DataRequest.h"
#include "DesktopServices.h"
#include "MessageLogResource.h"
#include "ObjectResource.h"
#include "PlugInArgList.h"
#include "PlugInManagerServices.h"
#include "PlugInRegistration.h"
#include "Progress.h"
#include "RasterDataDescriptor.h"
#include "RasterElement.h"
#include "RasterUtilities.h"
#include "SpatialDataView.h"
#include "SpatialDataWindow.h"
#include "switchOnEncoding.h"
#include "SpeckleRemove.h"
#include <limits>

REGISTER_PLUGIN_BASIC(OpticksTutorial, SpeckleRemove);

namespace
{

   #define OVERLAP_WINDOW_SIZE 3        //Window length
   #define GRADIENT_DIFF_THRESHOLD 100  //threshold for edge detection
   #define NUM_OF_VAR 5                 

	//Calculate the local mean and variance when horizantal/vertical edge is present
   double CaculatePartialMeanHV(int nrowStart, int rowEnd, int colStart, int colEnd, double pMatrix[][OVERLAP_WINDOW_SIZE*2+1], double *var)
   {
	   int i,j, nCount = 0;;
	   double meanVal = 0;
	   double varVal = 0;

	   for (i=nrowStart; i<rowEnd; i++)
	   {
		   for (j=colStart; j<colEnd; j++)
		   {
			   meanVal = meanVal + pMatrix[i][j];
			   nCount++;
		   }
	   }

	   meanVal = meanVal/nCount;

	   for (i=nrowStart; i<rowEnd; i++)
	   {
		   for (j=colStart; j<colEnd; j++)
		   {
			   varVal = varVal + (pMatrix[i][j]-meanVal)*(pMatrix[i][j]-meanVal);
		   }
	   }

	   varVal = varVal/nCount;
	   *var = varVal;

	   return meanVal;
   }

   //Calculate the local mean for diagonal edge is present
   double CaculatePartialMeanDiag(int colStart, int colEnd, int nIndex, double pMatrix[][OVERLAP_WINDOW_SIZE*2+1], double *var)
   {
	   int i,j, col1, col2, nCount = 0;;
	   double meanVal = 0;
	   double varVal = 0;

	   col1 = colStart;
	   col2 = colEnd;

	   for (i=0; i<OVERLAP_WINDOW_SIZE*2+1; i++)
	   {
		   for (j=col1; j<col2; j++)
		   {
			   meanVal = meanVal + pMatrix[i][j];
			   nCount++;
		   }

		   if (1 == nIndex)  //case 1
		   {
			   col1++;
		   }
		   else if (3 == nIndex) //case 3
		   {
			   col2--;
		   }
		   else if (5 == nIndex) //case 5
		   {
			   col2++;
		   }
		   else                 //case 7
		   {
			   col1--;
		   }
	   }

	   meanVal =  meanVal/nCount;

	   col1 = colStart;
	   col2 = colEnd;

	   for (i=0; i<OVERLAP_WINDOW_SIZE*2+1; i++)
	   {
		   for (j=col1; j<col2; j++)
		   {
			  varVal = varVal + (pMatrix[i][j]-meanVal)*(pMatrix[i][j]-meanVal);
		   }

		   if (1 == nIndex)  //case 1
		   {
			   col1++;
		   }
		   else if (3 == nIndex) //case 3
		   {
			   col2--;
		   }
		   else if (5 == nIndex) //case 5
		   {
			   col2++;
		   }
		   else                 //case 7
		   {
			   col1--;
		   }
	   }

	   varVal = varVal/nCount;
	   *var = varVal;

	   return meanVal;
   }

   //Calculate the local mean within the sliding window (the neighborhood of current pixel)
   double CaculatePartialMean(int nIndex, double pMatrix[][OVERLAP_WINDOW_SIZE*2+1], double *var)
   {
	   int i,j, m, n;
	   double meanVal = 0;

	   if (nIndex == 0)
	   {
		   i = 0;
		   m = OVERLAP_WINDOW_SIZE*2 + 1;
		   j = 4;
		   n = OVERLAP_WINDOW_SIZE*2 + 1;

		   meanVal = CaculatePartialMeanHV(i, m, j, n, pMatrix, var);
	   }
	   else if (nIndex == 2)
	   {
		   i = 0;
		   m = OVERLAP_WINDOW_SIZE + 1;
		   j = 0;
		   n = OVERLAP_WINDOW_SIZE*2 + 1;

		   meanVal = CaculatePartialMeanHV(i, m, j, n, pMatrix, var);
	   }
	   else if (nIndex == 4)
	   {
		   i = 0;
		   m = OVERLAP_WINDOW_SIZE*2 + 1;
		   j = 0;
		   n = OVERLAP_WINDOW_SIZE + 1;

		   meanVal = CaculatePartialMeanHV(i, m, j, n, pMatrix, var);
	   }
	   else if (nIndex == 6)
	   {
		   i = OVERLAP_WINDOW_SIZE;
		   m = OVERLAP_WINDOW_SIZE*2 + 1;
		   j = 0;
		   n = OVERLAP_WINDOW_SIZE*2 + 1;

		   meanVal = CaculatePartialMeanHV(i, m, j, n, pMatrix, var);
	   }
	   else if (nIndex == 1)
	   {
		   m = 0;
		   n = OVERLAP_WINDOW_SIZE*2 + 1;

		   meanVal = CaculatePartialMeanDiag(m, n, nIndex, pMatrix, var);
	   }
	   else if (nIndex == 3)
	   {
		   m = 0;
		   n = OVERLAP_WINDOW_SIZE*2 + 1;

		   meanVal = CaculatePartialMeanDiag(m, n, nIndex, pMatrix, var);
	   }
	   else if (nIndex == 5)
	   {
		   m = 0;
		   n = 1;

		   meanVal = CaculatePartialMeanDiag(m, n, nIndex,  pMatrix, var);
	   }
	   else 
	   {
		   m = OVERLAP_WINDOW_SIZE*2;
		   n = OVERLAP_WINDOW_SIZE*2 + 1;

		   meanVal = CaculatePartialMeanDiag(m, n, nIndex, pMatrix, var);
	   }

	   return meanVal;
   }

   //Estimate the local noise variance
   double EstimateNoise(double windowVal[][OVERLAP_WINDOW_SIZE*2+1])
   {
       double noiseVector[OVERLAP_WINDOW_SIZE*OVERLAP_WINDOW_SIZE];

	   double meanVal;
	   double var;
	   double noiseVar = 0;

       int count = 0;
       int currentRow = 0;
       int currentCol = 0;

	   int i,j,m,n,index;

       for (i=0;i<OVERLAP_WINDOW_SIZE;i++)
	   {
           currentCol = 0;
           for (j=0;j<OVERLAP_WINDOW_SIZE;j++)
		   {

               meanVal = 0;
               for (m=0; m<OVERLAP_WINDOW_SIZE; m++)
			   {
                   for (n=0; n<OVERLAP_WINDOW_SIZE; n++)
				   {
                       meanVal = meanVal +  windowVal[currentRow+m][currentCol+n];
				   }
			   }     
               meanVal = meanVal/(OVERLAP_WINDOW_SIZE*OVERLAP_WINDOW_SIZE);
      
               var = 0;
			   for (m=0; m<OVERLAP_WINDOW_SIZE; m++)
			   {
                   for (n=0; n<OVERLAP_WINDOW_SIZE; n++)
				   {
                       var = var + (windowVal[currentRow+m][currentCol+n]-meanVal)*(windowVal[currentRow+m][currentCol+n]-meanVal);
				   }
			   }
               
               var = var/(OVERLAP_WINDOW_SIZE*OVERLAP_WINDOW_SIZE);
      
			   //Sort the var in asending
			   var = var/(meanVal*meanVal);
			   index = count;
			   for (m=0; m<count; m++)
			   {
				   if (noiseVector[m] > var)
				   {
					   index = m;
					   break;
				   }
			   }
			   if (index != count)
			   {
				   for (m = count; m>index; m--)
				   {
					   noiseVector[m] = noiseVector[m-1];
				   }
			   }
               noiseVector[index] = var;
               currentCol = currentCol + 2;
               count = count+1;
		   }

		   currentRow = currentRow+2;
	   }

       for (i=0; i<NUM_OF_VAR; i++)
	   {
           noiseVar = noiseVector[i] + noiseVar;
	   }
	   noiseVar = noiseVar/NUM_OF_VAR;

	  return noiseVar;
   }


   template<typename T>
   void speckleNoiseRemove(T* pData, DataAccessor pSrcAcc, int row, int col, int rowSize, int colSize)
   {
	  int i, j, m, n;
	  int distX, distY;
	  double meanVal = 0.0;
	  double coeff = 0.0;
	  double gradientVector[4];
	  double pixelVal = 0.0;
	  double noiseVar = 64;
	  double localVar = 0.0;

	  T tempVal;
	  double windowVal[OVERLAP_WINDOW_SIZE*2+1][OVERLAP_WINDOW_SIZE*2+1];
	  double meanWindow[OVERLAP_WINDOW_SIZE][OVERLAP_WINDOW_SIZE];

	  if ((col-OVERLAP_WINDOW_SIZE < 0) || (col+OVERLAP_WINDOW_SIZE > colSize - 1))
	  {
		  pSrcAcc->toPixel(i, j);
          VERIFYNRV(pSrcAcc.isValid());
          *pData = *reinterpret_cast<T*>(pSrcAcc->getColumn()); 
		  return;
	  }

	  if ((row-OVERLAP_WINDOW_SIZE < 0) || (row+OVERLAP_WINDOW_SIZE > rowSize - 1))
	  {
		  pSrcAcc->toPixel(i, j);
          VERIFYNRV(pSrcAcc.isValid());
          *pData = *reinterpret_cast<T*>(pSrcAcc->getColumn()); 
		  return;
	  }

	  //Get the pixels in the sliding window
	  m = 0;
	  for (i=row - OVERLAP_WINDOW_SIZE; i<= row + OVERLAP_WINDOW_SIZE; i++)
	  {
		  n = 0;
		  for (j=col - OVERLAP_WINDOW_SIZE; j<= col + OVERLAP_WINDOW_SIZE; j++)
		  {
			 pSrcAcc->toPixel(i, j);
             VERIFYNRV(pSrcAcc.isValid());
             tempVal = *reinterpret_cast<T*>(pSrcAcc->getColumn()); 
			 windowVal[m][n] = tempVal;
			 n++;
		  }
		  m++;
	  }
	  pixelVal = windowVal[OVERLAP_WINDOW_SIZE][OVERLAP_WINDOW_SIZE];

	  //Calculate the mean in each 3 by 3 window
	  distX = 0;
	  for (i = 0; i< OVERLAP_WINDOW_SIZE; i++)
	  {
		  distY = 0;

		  for (j=0; j<OVERLAP_WINDOW_SIZE; j++)
		  {
			  meanVal = 0.0;

			  for (m = 0; m < OVERLAP_WINDOW_SIZE; m++)
			  {
				  for (n = 0 ; n < OVERLAP_WINDOW_SIZE; n++)
				  {
					  meanVal = meanVal + windowVal[distX + m][distY + n];
				  }
			  }
			  meanWindow[i][j] = meanVal/(OVERLAP_WINDOW_SIZE * OVERLAP_WINDOW_SIZE);

			  distY += 2;
		  }
		  distX += 2;
	  }

	  //Calculate the gradient at 4 directions
	  gradientVector[0] = abs((meanWindow[0][0] - meanWindow[2][0]) + 2*(meanWindow[0][1] - meanWindow[2][1]) + (meanWindow[0][2] - meanWindow[2][2])); // horizantal
	  gradientVector[1] = abs((meanWindow[0][0] - meanWindow[0][2]) + 2*(meanWindow[1][0] - meanWindow[1][2]) + (meanWindow[2][0] - meanWindow[2][2])); // vertical
	  gradientVector[2] = abs((meanWindow[0][1] - meanWindow[1][2]) + 2*(meanWindow[0][0] - meanWindow[2][2]) + (meanWindow[1][0] - meanWindow[2][1])); // digonal 45
	  gradientVector[3] = abs((meanWindow[0][1] - meanWindow[1][0]) + 2*(meanWindow[0][2] - meanWindow[2][0]) + (meanWindow[1][2] - meanWindow[2][1])); // digonal 135

	  //Get the maxium and minimum gradient
	  double maxGradient = -1;
	  double minGradient = 65535;
	  int    nIndex = 0;
	  for (i=0; i<4; i++)
	  {
		  if (gradientVector[i] > maxGradient)
		  {
			  maxGradient = gradientVector[i];
			  nIndex = i;
		  }

		  if (gradientVector[i] < minGradient)
		  {
			  minGradient = gradientVector[i];
		  }
	  }

   //If there is an edge, we need to determine to which sub-area the pixel belongs
   // and caluclate the local mean and variance of the sub-area
	  if (abs(maxGradient) > GRADIENT_DIFF_THRESHOLD)
	  {
		  if (0 == nIndex)
		  {
			  if (abs(meanWindow[0][1] - meanWindow[1][1]) < abs(meanWindow[2][1] - meanWindow[1][1]))
			  {
			      meanVal = CaculatePartialMean(2, windowVal, &localVar);
			  }
			  else
			  {
			      meanVal = CaculatePartialMean(6, windowVal, &localVar);
			  }
		  }
		  else if (1 == nIndex)
		  {
			  if (abs(meanWindow[1][1] - meanWindow[1][2]) < abs(meanWindow[1][0] - meanWindow[1][1]))
			  {
			      meanVal = CaculatePartialMean(0, windowVal, &localVar);
			  }
			  else
			  {
			      meanVal = CaculatePartialMean(4, windowVal, &localVar);
			  }
		  }
		  if (2 == nIndex)
		  {
			  if (abs(meanWindow[0][0] - meanWindow[1][1]) < abs(meanWindow[2][2] - meanWindow[1][1]))
			  {
			      meanVal = CaculatePartialMean(3, windowVal, &localVar);
			  }
			  else
			  {
			      meanVal = CaculatePartialMean(7, windowVal, &localVar);
			  }
		  }
		  else
		  {
			  if (abs(meanWindow[0][2] - meanWindow[1][1]) < abs(meanWindow[2][0] - meanWindow[1][1]))
			  {
			      meanVal = CaculatePartialMean(1, windowVal, &localVar);
			  }
			  else
			  {
			      meanVal = CaculatePartialMean(5, windowVal, &localVar);
			  }
		  }

	  }
	  else
	  {
		  meanVal = CaculatePartialMeanHV(0, OVERLAP_WINDOW_SIZE*2+1, 0, OVERLAP_WINDOW_SIZE*2+1, windowVal, &localVar);

	  }

	  //Use MMSEE filter to denoise
	  if (localVar > 0)
	  {
	      double L = 1/EstimateNoise(windowVal); //Estimate noise
          double  deltax = (L*localVar-meanVal*meanVal)/(L+1);
          if (deltax >= 0)
	      {
              pixelVal = deltax*(pixelVal - meanVal)/(deltax + (meanVal*meanVal + deltax)/L);
	      }
	      else
	      {
              pixelVal = meanVal;
	      }
	  }
	  else
	  {
		   pixelVal = meanVal;
	  }

      *pData = static_cast<T>(pixelVal);
   }
};

SpeckleRemove::SpeckleRemove()
{
   setDescriptorId("{73AC829B-C975-48a1-9A97-C697CFA709E6}");
   setName("Speckle Remove");
   setDescription("Remove speckle noise for SAR");
   setCreator("Yiwei Zhang");
   setVersion("Sample");
   setCopyright("Copyright (C) 2008, Ball Aerospace & Technologies Corp.");
   setProductionStatus(false);
   setType("Sample");
   setSubtype("Denoise");
   setMenuLocation("[Tutorial]/Speckle Remove");
   setAbortSupported(true);
}

SpeckleRemove::~SpeckleRemove()
{
}

bool SpeckleRemove::getInputSpecification(PlugInArgList*& pInArgList)
{
   VERIFY(pInArgList = Service<PlugInManagerServices>()->getPlugInArgList());
   pInArgList->addArg<Progress>(Executable::ProgressArg(), NULL, "Progress reporter");
   pInArgList->addArg<RasterElement>(Executable::DataElementArg(), "Perform speckle remove on this data element");
   return true;
}

bool SpeckleRemove::getOutputSpecification(PlugInArgList*& pOutArgList)
{
   VERIFY(pOutArgList = Service<PlugInManagerServices>()->getPlugInArgList());
   pOutArgList->addArg<RasterElement>("Result", NULL);
   return true;
}

bool SpeckleRemove::execute(PlugInArgList* pInArgList, PlugInArgList* pOutArgList)
{
   StepResource pStep("Speckle remove", "app", "94F01F0D-EDA9-4f2c-AB00-1BA1448E8783");
   if (pInArgList == NULL || pOutArgList == NULL)
   {
      return false;
   }
   Progress* pProgress = pInArgList->getPlugInArgValue<Progress>(Executable::ProgressArg());
   RasterElement* pCube = pInArgList->getPlugInArgValue<RasterElement>(Executable::DataElementArg());
   if (pCube == NULL)
   {
      std::string msg = "A raster cube must be specified.";
      pStep->finalize(Message::Failure, msg);
      if (pProgress != NULL) 
      {
         pProgress->updateProgress(msg, 0, ERRORS);
      }
      return false;
   }
   RasterDataDescriptor* pDesc = static_cast<RasterDataDescriptor*>(pCube->getDataDescriptor());
   VERIFY(pDesc != NULL);
   if (pDesc->getDataType() == INT4SCOMPLEX || pDesc->getDataType() == FLT8COMPLEX)
   {
      std::string msg = "Speckle remove cannot be performed on complex types.";
      pStep->finalize(Message::Failure, msg);
      if (pProgress != NULL) 
      {
         pProgress->updateProgress(msg, 0, ERRORS);
      }
      return false;
   }

   FactoryResource<DataRequest> pRequest;
   pRequest->setInterleaveFormat(BSQ);
   DataAccessor pSrcAcc = pCube->getDataAccessor(pRequest.release());

   ModelResource<RasterElement> pResultCube(RasterUtilities::createRasterElement(pCube->getName() +
      "_Speckle_Remove_Result", pDesc->getRowCount(), pDesc->getColumnCount(), pDesc->getDataType()));
   if (pResultCube.get() == NULL)
   {
      std::string msg = "A raster cube could not be created.";
      pStep->finalize(Message::Failure, msg);
      if (pProgress != NULL) 
      {
         pProgress->updateProgress(msg, 0, ERRORS);
      }
      return false;
   }
   FactoryResource<DataRequest> pResultRequest;
   pResultRequest->setWritable(true);
   DataAccessor pDestAcc = pResultCube->getDataAccessor(pResultRequest.release());

   for (unsigned int row = 0; row < pDesc->getRowCount(); ++row)
   {
      if (pProgress != NULL)
      {
         pProgress->updateProgress("Remove result", row * 100 / pDesc->getRowCount(), NORMAL);
      }
      if (isAborted())
      {
         std::string msg = getName() + " has been aborted.";
         pStep->finalize(Message::Abort, msg);
         if (pProgress != NULL)
         {
            pProgress->updateProgress(msg, 0, ABORT);
         }
         return false;
      }
      if (!pDestAcc.isValid())
      {
         std::string msg = "Unable to access the cube data.";
         pStep->finalize(Message::Failure, msg);
         if (pProgress != NULL) 
         {
            pProgress->updateProgress(msg, 0, ERRORS);
         }
         return false;
      }
      for (unsigned int col = 0; col < pDesc->getColumnCount(); ++col)
      {
         switchOnEncoding(pDesc->getDataType(), speckleNoiseRemove, pDestAcc->getColumn(), pSrcAcc, row, col,
            pDesc->getRowCount(), pDesc->getColumnCount());
         pDestAcc->nextColumn();
      }

      pDestAcc->nextRow();
   }

   if (!isBatch())
   {
      Service<DesktopServices> pDesktop;

      SpatialDataWindow* pWindow = static_cast<SpatialDataWindow*>(pDesktop->createWindow(pResultCube->getName(),
         SPATIAL_DATA_WINDOW));

      SpatialDataView* pView = (pWindow == NULL) ? NULL : pWindow->getSpatialDataView();
      if (pView == NULL)
      {
         std::string msg = "Unable to create view.";
         pStep->finalize(Message::Failure, msg);
         if (pProgress != NULL) 
         {
            pProgress->updateProgress(msg, 0, ERRORS);
         }
         return false;
      }

      pView->setPrimaryRasterElement(pResultCube.get());
      pView->createLayer(RASTER, pResultCube.get());
   }

   if (pProgress != NULL)
   {
      pProgress->updateProgress("Speckle remove is compete.", 100, NORMAL);
   }

   pOutArgList->setPlugInArgValue("Speckle Remove Result", pResultCube.release());

   pStep->finalize();
   return true;
}
