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
#include "WaveletSpeckleRemove.h"
#include "Wavelet.h"
#include <limits>

REGISTER_PLUGIN_BASIC(OpticksTutorial, WaveletSpeckleRemove);

namespace
{
    #define LAYERS 4
	WaveletNode NodeList[LAYERS+1];

	int filterLen = 4;
	double pLoFilter[4] = {0.4830, 0.8365, 0.2241, -0.1294};
	double pHiFilter[4] = {-0.1294,-0.2241, 0.8365, -0.4830};
	double pRecHiFilter[4] = {-0.4830, 0.8365, -0.2241, -0.1294};
	double pRecLoFilter[4] = {-0.1294,0.2241, 0.8365, 0.4830};

   void ProcessData(DataAccessor pSrcAcc, double *pBuffer, unsigned int row, unsigned int col, unsigned int rowBlocks, unsigned int colBlocks, EncodingType type)
   {
	  unsigned int nCount = 0;
	  double pixelVal;

	  for (unsigned int i=0; i<rowBlocks; i++)
	  {
		  for (unsigned int j=0; j<colBlocks; j++)
		  {
			  pSrcAcc->toPixel(row+i, col+j);
              VERIFYNRV(pSrcAcc.isValid());
		      pixelVal = Service<ModelServices>()->getDataValue(type, pSrcAcc->getColumn(), COMPLEX_MAGNITUDE, 0);

			  *(pBuffer + nCount) = log(1+pixelVal);
			  nCount++;
		  }
	  }

	  ShiftInvariantWaveletTransform(pBuffer, rowBlocks, colBlocks, pLoFilter, pHiFilter, filterLen, LAYERS, NodeList);

	  WaveletDenoise(NodeList, pBuffer, LAYERS);

	  ShiftInvariantInverseWaveletTransform(rowBlocks, colBlocks, pRecLoFilter, pRecHiFilter, filterLen, LAYERS, NodeList);

	  for (nCount=0; nCount<rowBlocks*colBlocks; nCount++)
	  {
		  pixelVal = *(pBuffer + nCount);
		  *(pBuffer + nCount) = exp(pixelVal) - 1;
	  }

	  ReleaseList(NodeList, LAYERS);
   }

   template<typename T>
   void speckleNoiseRemove(T* pData, double *pSrc)
   {
	   double pixelVal;
	   if (pSrc != NULL)
	   {
		   pixelVal = *pSrc;
	       *pData = static_cast<T>(pixelVal);
	   }
   }
};

WaveletSpeckleRemove::WaveletSpeckleRemove()
{
   setDescriptorId("{0183FA2C-18E9-4516-8931-9EFBCD3616B9}");
   setName("Wavelet Speckle Remove");
   setDescription("Remove speckle noise for SAR");
   setCreator("Yiwei Zhang");
   setVersion("Sample");
   setCopyright("Copyright (C) 2008, Ball Aerospace & Technologies Corp.");
   setProductionStatus(false);
   setType("Sample");
   setSubtype("Denoise");
   setMenuLocation("[SAR]/Wavelet Speckle Remove");
   setAbortSupported(true);

   rowBlocks = 128;
   colBlocks = 128;
   pBuffer = (double *)malloc(sizeof(double)*(10+rowBlocks)*(10+colBlocks));  
}

WaveletSpeckleRemove::~WaveletSpeckleRemove()
{
	free(pBuffer);
}

bool WaveletSpeckleRemove::getInputSpecification(PlugInArgList*& pInArgList)
{
   VERIFY(pInArgList = Service<PlugInManagerServices>()->getPlugInArgList());
   pInArgList->addArg<Progress>(Executable::ProgressArg(), NULL, "Progress reporter");
   pInArgList->addArg<RasterElement>(Executable::DataElementArg(), "Perform speckle remove on this data element");
   return true;
}

bool WaveletSpeckleRemove::getOutputSpecification(PlugInArgList*& pOutArgList)
{
   VERIFY(pOutArgList = Service<PlugInManagerServices>()->getPlugInArgList());
   pOutArgList->addArg<RasterElement>("Result", NULL);
   return true;
}

bool WaveletSpeckleRemove::execute(PlugInArgList* pInArgList, PlugInArgList* pOutArgList)
{
   StepResource pStep("Wavelet Speckle remove", "app", "9B815808-7E87-4449-8549-4B6AB45F816B");
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
   EncodingType ResultType = pDesc->getDataType();
   if (pDesc->getDataType() == INT4SCOMPLEX)
   {
      ResultType = INT4SBYTES;
   }
   else if (pDesc->getDataType() == FLT8COMPLEX)
   {
      ResultType = FLT8BYTES;
   }

   FactoryResource<DataRequest> pRequest;
   pRequest->setInterleaveFormat(BSQ);
   DataAccessor pSrcAcc = pCube->getDataAccessor(pRequest.release());

   ModelResource<RasterElement> pResultCube(RasterUtilities::createRasterElement(pCube->getName() +
      "_Speckle_Remove_Result", pDesc->getRowCount(), pDesc->getColumnCount(), ResultType));
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

   unsigned int rowLoops;
   unsigned int colLoops;
   unsigned int rowIndex = 0;
   unsigned int colIndex = 0;
   
   if (0 == pDesc->getRowCount()%rowBlocks)
   {
	   rowLoops = pDesc->getRowCount()/rowBlocks;
   }
   else
   {
	   rowLoops = pDesc->getRowCount()/rowBlocks + 1;
   }

   if (0 == pDesc->getColumnCount()%colBlocks)
   {
	   colLoops = pDesc->getColumnCount()/colBlocks;
   }
   else
   {
	   colLoops = pDesc->getColumnCount()/colBlocks + 1;
   }

   for (unsigned int i = 0; i < rowLoops; i++)
   {
	   if ( rowIndex + rowBlocks > pDesc->getRowCount())
	   {
		   rowIndex = pDesc->getRowCount() - rowBlocks;
	   }

	   colIndex = 0;

	   for (unsigned int j = 0; j < colLoops; j++)
	   {
		   if ( colIndex + colBlocks > pDesc->getColumnCount())
	       {
		       colIndex = pDesc->getColumnCount() - colBlocks;
	       }

		   if (pProgress != NULL)
           {
               pProgress->updateProgress("Remove result", (i*colLoops+j) / (rowLoops*colLoops), NORMAL);
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
      
           //Process the data in current block
		   ProcessData(pSrcAcc, pBuffer, rowIndex, colIndex, rowBlocks, colBlocks, pDesc->getDataType());

		   //Output the value 
           for (unsigned int m = 0; m < rowBlocks; m++)
		   {
			   for (unsigned int n = 0; n < colBlocks; n++)
			   {
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

				   pDestAcc->toPixel(rowIndex+m, colIndex+n);
				   
				   switchOnEncoding(ResultType, speckleNoiseRemove, pDestAcc->getColumn(), (pBuffer+m*colBlocks+n));
			   }
		   }
		   colIndex += colBlocks;
	   }
	   rowIndex += rowBlocks;
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
