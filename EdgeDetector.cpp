/*
 * The information in this file is
 * Copyright(c) 2007 Ball Aerospace & Technologies Corporation
 * and is subject to the terms and conditions of the
 * GNU Lesser General Public License Version 2.1
 * The license text is available from   
 * http://www.gnu.org/licenses/lgpl.html
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
#include "EdgeDetector.h"
#include "EdgeRatioThresholdDlg.h"
#include <limits>
#include <algorithm>

REGISTER_PLUGIN_BASIC(OpticksTutorial, EdgeDetector);

namespace
{

    #define SMALL_WINDOW_SIZE 3
    #define MEDIAN_WINDOW_SIZE 5
    #define LARGE_WINDOW_SIZE 7
    
    #define SMALL_WINDOW_THRESHOLD 0.65
	#define MEDIAN_WINDOW_THRESHOLD 0.5
	#define LARGE_WINDOW_THRESHOLD 0.4

	//Calculate the block means for horizantal/vertical edge case
   void CaculateHVMean(int tag, int SUB_WINDOW_SIZE, double **subWindow, double &lineMean, double &blockMean1, double &blockMean2)
   {
	   int currentRow = 0;
       int currentCol = 0;
	   int i,j;

	   lineMean = 0;
       blockMean1 = 0;
       blockMean2 = 0;

	   //tag is 0 corresponds to horizontal case, tag is 1 corresponds to vertical case
       if (tag == 0)
	   {
           currentRow = (SUB_WINDOW_SIZE-1)/2;
           currentCol = 0;
	   }
       else
	   {
           currentRow = 0;
           currentCol = (SUB_WINDOW_SIZE-1)/2;
	   }

	   //First, we calculate the mean value of the middle vertical/horizontal line
       for (i=1; i<=SUB_WINDOW_SIZE; i++)
	   {
           lineMean = lineMean + *(*(subWindow+currentRow)+currentCol);
 
           if (tag == 0)
               currentCol = currentCol + 1;
           else
               currentRow = currentRow + 1;
	   }
       lineMean = lineMean/SUB_WINDOW_SIZE;

	   //Then we calculate the mean of the blocks besides the line
       currentRow = 0;
       currentCol = 0;
       if (tag == 0)
	   {
		   for (i=1; i<=(SUB_WINDOW_SIZE-1)/2; i++)
		   {
               currentCol = 0;

			   for (j=1; j<=SUB_WINDOW_SIZE; j++)
			   {
                   blockMean1 = blockMean1 + *(*(subWindow+currentRow)+currentCol);
                   blockMean2 = blockMean2 + *(*(subWindow+currentRow+(SUB_WINDOW_SIZE+1)/2)+currentCol);
                   currentCol = currentCol+1;
			   }
               currentRow = currentRow + 1;
		   }
	   }
       else
	   {
           for (i=1; i<=SUB_WINDOW_SIZE; i++)
		   {
               currentCol = 0;

               for (j=1; j<=(SUB_WINDOW_SIZE-1)/2;j++)
			   {
                   blockMean1 = blockMean1 + *(*(subWindow+currentRow)+currentCol);
                   blockMean2 = blockMean2 + *(*(subWindow+currentRow)+currentCol+(SUB_WINDOW_SIZE+1)/2);
                   currentCol = currentCol+1;
			   }
               currentRow = currentRow + 1;
		   }
	   }

       blockMean1 = blockMean1/(SUB_WINDOW_SIZE*(SUB_WINDOW_SIZE-1)/2);
       blockMean2 = blockMean2/(SUB_WINDOW_SIZE*(SUB_WINDOW_SIZE-1)/2);
   }

   //Calculate the block means for diagonal edge case
   void CaculateDiagMean(int tag, int SUB_WINDOW_SIZE, double **subWindow, double &lineMean, double &blockMean1, double &blockMean2)
   {
	   int colStart = 0;
       int colEnd = 0;
	   int i,j;

	   lineMean = 0;
       blockMean1 = 0;
       blockMean2 = 0;


	   //tag is 1 corresponds to 45 degree case, tag is 0 corresponds to 135 degree case
       if (tag == 1)
	   {
           j = SUB_WINDOW_SIZE-1;
	   }
       else
	   {
           j = 0;
	   }

	   //First, we calculate the mean value of the diagonal line
	   for (i=0; i<SUB_WINDOW_SIZE; i++)
	   {
           lineMean = lineMean + *(*(subWindow+i)+j);
 
           if (tag == 1)
               j = j - 1;
           else      
               j = j + 1;
	   }
       lineMean = lineMean/SUB_WINDOW_SIZE;

	   //Then we calculate the mean of the blocks besides the line
       if (tag == 1)
	   {
           colStart = 0;
           colEnd = SUB_WINDOW_SIZE-2;
    
           for (i=0; i<=SUB_WINDOW_SIZE-2;i++)
		   {
               for (j=colStart; j<=colEnd; j++)
			   {
                   blockMean1 = blockMean1 + *(*(subWindow+i)+j);
			   }
               colEnd = colEnd - 1;
		   }
    
           colStart = SUB_WINDOW_SIZE-1;
           colEnd = SUB_WINDOW_SIZE-1;
           for (i=1; i<=SUB_WINDOW_SIZE-1; i++)
		   {
               for (j=colStart; j<=colEnd; j++)
			   {
                   blockMean2 = blockMean2 + *(*(subWindow+i)+j);
			   }
               colStart = colStart-1;
		   }
	   }
       else
	   {
           colStart = 0;
           colEnd = 0;
    
           for (i=1;i<=SUB_WINDOW_SIZE-1;i++)
		   {
               for (j=colStart; j<=colEnd;j++)
			   {
                   blockMean1 = blockMean1 + *(*(subWindow+i)+j);
			   }
               colEnd = colEnd + 1;
		   }
    
           colStart = 1;
           colEnd = SUB_WINDOW_SIZE-1;
           for (i=0;i<=SUB_WINDOW_SIZE-2;i++)
		   {
               for (j=colStart; j<=colEnd;j++)
			   {
                   blockMean2 = blockMean2 + *(*(subWindow+i)+j);
			   }
               colStart = colStart+1;
		   }
	   }

       blockMean1 = blockMean1/(SUB_WINDOW_SIZE*(SUB_WINDOW_SIZE-1)/2);
       blockMean2 = blockMean2/(SUB_WINDOW_SIZE*(SUB_WINDOW_SIZE-1)/2);
   }

   //Calculate the edge ratio of current pixel's neighbor in all directions and get the minimum ratio
   //If the ratio is less than threshold, then there is an edge crossing current pixel
   bool SubWindowEdgeDetect(double **subWindow, int SUB_WINDOW_SIZE, double threshold)
   {
	   //Vertical direction
       double leftBlockMean = 0;
       double rightBlockMean = 0;
	   double verticalLineMean = 0;

	   //horizontal direction
       double upperBlockMean = 0;
       double lowerBlockMean = 0;     
       double horizantalLineMean = 0;

	   //digonal 45 degree direction
       double leftUpperDiagBlockMean = 0;
       double rightLowerDiagBlockMean = 0;
	   double diag45LineMean = 0;

       //digonal 135 degree direction
       double leftLowerDiagBlockMean = 0;
       double rightUpperDiagBlockMean = 0;     
       double diag135LineMean = 0;


       double EdgeRatio[4]; //Edge ratio at 4 directions
	   double P1, P2, P3;
	   int index = 0;
	   int i;

       //vertical
	   CaculateHVMean(1, SUB_WINDOW_SIZE, subWindow, verticalLineMean, leftBlockMean, rightBlockMean);
       P1 = std::min(leftBlockMean/rightBlockMean, rightBlockMean/leftBlockMean);
       P2 = std::min(leftBlockMean/verticalLineMean, verticalLineMean/leftBlockMean);
       P3 = std::min(rightBlockMean/verticalLineMean, verticalLineMean/rightBlockMean);
       EdgeRatio[index] = std::min(P1, P2);
       EdgeRatio[index] = std::min(EdgeRatio[index], P3);

       //horizantal
       index = index + 1;
	   CaculateHVMean(0, SUB_WINDOW_SIZE, subWindow, horizantalLineMean, upperBlockMean, lowerBlockMean);
       P1 = std::min(upperBlockMean/lowerBlockMean, lowerBlockMean/upperBlockMean);
       P2 = std::min(upperBlockMean/horizantalLineMean, horizantalLineMean/upperBlockMean);
       P3 = std::min(lowerBlockMean/horizantalLineMean, horizantalLineMean/lowerBlockMean);
       EdgeRatio[index] = std::min(P1, P2);
       EdgeRatio[index] = std::min(EdgeRatio[index], P3);

       //diagonal 45 degree
       index = index + 1;
	   CaculateDiagMean(1, SUB_WINDOW_SIZE, subWindow, diag45LineMean, leftUpperDiagBlockMean, rightLowerDiagBlockMean);
       P1 = std::min(leftUpperDiagBlockMean/rightLowerDiagBlockMean, rightLowerDiagBlockMean/leftUpperDiagBlockMean);
       P2 = std::min(leftUpperDiagBlockMean/diag45LineMean, diag45LineMean/leftUpperDiagBlockMean);
       P3 = std::min(rightLowerDiagBlockMean/diag45LineMean, diag45LineMean/rightLowerDiagBlockMean);
       EdgeRatio[index] = std::min(P1, P2);
       EdgeRatio[index] = std::min(EdgeRatio[index], P3);

       //diagonal 135 degree
       index = index + 1;
	   CaculateDiagMean(0, SUB_WINDOW_SIZE, subWindow, diag135LineMean, leftLowerDiagBlockMean, rightUpperDiagBlockMean);
       P1 = std::min(leftLowerDiagBlockMean/rightUpperDiagBlockMean, rightUpperDiagBlockMean/leftLowerDiagBlockMean);
       P2 = std::min(leftLowerDiagBlockMean/diag135LineMean, diag135LineMean/leftLowerDiagBlockMean);
       P3 = std::min(rightUpperDiagBlockMean/diag135LineMean, diag135LineMean/rightUpperDiagBlockMean);
       EdgeRatio[index] = std::min(P1, P2);
       EdgeRatio[index] = std::min(EdgeRatio[index], P3);


       //Get the minimum ratio
       double minRatio = EdgeRatio[0];
       index = 0;
       for (i=1; i<4; i++)
	   {
           if (EdgeRatio[i] < minRatio)
		   {
		       minRatio = EdgeRatio[i];
		       index = i;
		   }
	   }

	   //Compare with threshold to determine if there is edge
       if (EdgeRatio[index] < threshold)
	   {
           return true;
	   }
       else
	   {
           return false;
	   }

   }


   //Estimate the local noise variance
   void RatioEdgeDetect(DataAccessor pSrcAcc, int row, int col, int rowSize, int colSize, EncodingType type, int SUB_WINDOW_SIZE, double threshold, bool &bHasEdge)
   {
	   int HALF_WINDOW_SIZE = (SUB_WINDOW_SIZE-1)/2;
	   int i,j,m,n;
	   double **subWindow;

       if ((col-HALF_WINDOW_SIZE < 0) || (col+HALF_WINDOW_SIZE > colSize - 1))
	   {
		   bHasEdge = false;
		   return;
	   }

	   if ((row-HALF_WINDOW_SIZE < 0) || (row+HALF_WINDOW_SIZE > rowSize - 1))
	   {
		   bHasEdge = false;
		   return;
	   }

	   subWindow = (double **)malloc(sizeof(double *)*SUB_WINDOW_SIZE);
	   for (i=0; i<SUB_WINDOW_SIZE;i++)
	   {
		   *(subWindow+i) = (double *)malloc(sizeof(double)*SUB_WINDOW_SIZE);
	   }

	   //Get the pixels in the window
	   m = 0;
	   for (i=row - HALF_WINDOW_SIZE; i<= row + HALF_WINDOW_SIZE; i++)
	   {
	   	   n = 0;
		   for (j=col - HALF_WINDOW_SIZE; j<= col + HALF_WINDOW_SIZE; j++)
		   {
		       pSrcAcc->toPixel(i, j);
               VERIFYNRV(pSrcAcc.isValid());
			   *(*(subWindow+m)+n) = Service<ModelServices>()->getDataValue(type, pSrcAcc->getColumn(), COMPLEX_MAGNITUDE, 0)+0.000001;
			   n++;
		   }
		   m++;
	   }

	   bHasEdge =  SubWindowEdgeDetect(subWindow, SUB_WINDOW_SIZE, threshold);

	   for (i=0; i<SUB_WINDOW_SIZE;i++)
	   {
		   free(*(subWindow+i));
	   }
	   free(subWindow);

	   return;
   }


   template<typename T>
   void EdgeDetectSAR(T* pData, DataAccessor pSrcAcc, int row, int col, int rowSize, int colSize, EncodingType type, double t1, double t2, double t3)
   {
	   bool bHasEdge = false;
	   unsigned char pixelVal = 255;


	   //Test edge for different window size and threshold
	   RatioEdgeDetect(pSrcAcc, row,  col, rowSize, colSize, type, SMALL_WINDOW_SIZE, t1, bHasEdge);
       if (bHasEdge)  //For small window, noise may be detected as edge, so need to check with larger window sizes
	   {
		   //pixelVal = 0;
	  
           RatioEdgeDetect(pSrcAcc, row,  col, rowSize, colSize, type, MEDIAN_WINDOW_SIZE, t2, bHasEdge);

           if (bHasEdge)
		   {
			   pixelVal = 0;

		   }
		   else
		   {
			   RatioEdgeDetect(pSrcAcc, row,  col, rowSize, colSize, type, LARGE_WINDOW_SIZE, t3, bHasEdge);
			   if (bHasEdge)
		       {
			       pixelVal = 0;
		       }
		   }
	   }

	   *pData = static_cast<T>(pixelVal);

   }
};

EdgeDetector::EdgeDetector()
{
   setDescriptorId("{C72E631E-3B4A-4b40-93B0-D2AA01E23315}");
   setName("Edge Detector");
   setDescription("Edge Detector for SAR");
   setCreator("Yiwei Zhang");
   setVersion("Sample");
   setCopyright("Copyright (C) 2008, Ball Aerospace & Technologies Corp.");
   setProductionStatus(false);
   setType("Sample");
   setSubtype("SAR Edge");
   setMenuLocation("[SAR]/SAR Edge Detect");
   setAbortSupported(true);
}

EdgeDetector::~EdgeDetector()
{
}

bool EdgeDetector::getInputSpecification(PlugInArgList*& pInArgList)
{
   VERIFY(pInArgList = Service<PlugInManagerServices>()->getPlugInArgList());
   pInArgList->addArg<Progress>(Executable::ProgressArg(), NULL, "Progress reporter");
   pInArgList->addArg<RasterElement>(Executable::DataElementArg(), "Perform edge detect on this data element");
   return true;
}

bool EdgeDetector::getOutputSpecification(PlugInArgList*& pOutArgList)
{
   VERIFY(pOutArgList = Service<PlugInManagerServices>()->getPlugInArgList());
   pOutArgList->addArg<RasterElement>("Result", NULL);
   return true;
}

bool EdgeDetector::execute(PlugInArgList* pInArgList, PlugInArgList* pOutArgList)
{
   StepResource pStep("Edge Detector", "app", "37C57772-DD49-4532-8BC6-9CFB8587D0C9");
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
   EncodingType ResultType = INT1UBYTE;

   FactoryResource<DataRequest> pRequest;
   pRequest->setInterleaveFormat(BSQ);
   DataAccessor pSrcAcc = pCube->getDataAccessor(pRequest.release());

   ModelResource<RasterElement> pResultCube(RasterUtilities::createRasterElement(pCube->getName() +
      "_Edge_Detect_Result", pDesc->getRowCount(), pDesc->getColumnCount(), ResultType));
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

   Service<DesktopServices> pDesktop;
   EdgeRatioThresholdDlg dlg(pDesktop->getMainWidget(), SMALL_WINDOW_THRESHOLD, MEDIAN_WINDOW_THRESHOLD, LARGE_WINDOW_THRESHOLD);
   int stat = dlg.exec();
   if (stat == QDialog::Accepted)
   {
      for (unsigned int row = 0; row < pDesc->getRowCount(); ++row)
      {
         if (pProgress != NULL)
         {
            pProgress->updateProgress("Edge detect ", row * 100 / pDesc->getRowCount(), NORMAL);
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
            switchOnEncoding(ResultType, EdgeDetectSAR, pDestAcc->getColumn(), pSrcAcc, row, col,
                             pDesc->getRowCount(), pDesc->getColumnCount(), pDesc->getDataType(), 
			                 dlg.getSmallThreshold(), dlg.getMedianThreshold(), dlg.getLargeThreshold());
            pDestAcc->nextColumn();
         }

         pDestAcc->nextRow();
      }

      if (!isBatch())
      {
         //Service<DesktopServices> pDesktop;

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
         pProgress->updateProgress("Edge detect compete.", 100, NORMAL);
      }

      pOutArgList->setPlugInArgValue("Edge detect result", pResultCube.release());

      pStep->finalize();
   }
   return true;
}
