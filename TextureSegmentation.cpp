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
#include "TextureSegmentation.h"
#include "ImageSegmentation.h"
#include <limits>

REGISTER_PLUGIN_BASIC(OpticksTutorial, TextureSegmentation);

namespace
{
   #define EDGE_THRESHOLD_VALUE 0.8
   #define MAXIMUM_LOOPS 30
   #define LAMDA_VALUE  40
   #define MU_VALUE     1000
   #define NUMBER_OF_BINS  256
   #define PARZEN_WINDOW_SIZE 2 
   #define PATCH_SIZE 7 

   void MakeSegmentation(DataAccessor pSrcAcc, float *pBuffer, float *pfu, unsigned int rows, unsigned int cols, EncodingType type)
   {
	  double pixelVal;
	  double maxVal = 0;
	  float pfVecParameters[8];

	  unsigned int nCount = 0;

	  pfVecParameters[0] = rows;
      pfVecParameters[1] = cols;  
	  pfVecParameters[2] = MAXIMUM_LOOPS; /* Maximum number of loops */    
	  pfVecParameters[3] = LAMDA_VALUE; //Lamda value
      pfVecParameters[4] = MU_VALUE;    // Mu value
	  pfVecParameters[5] = NUMBER_OF_BINS; /* Number of bins for histograms */
      pfVecParameters[6] = PARZEN_WINDOW_SIZE; /* standard deviation in Parzen estimation */
      pfVecParameters[7] = PATCH_SIZE; /* patch size */

	  for (unsigned int j=0; j<cols; j++)
	  {
		  for (unsigned int i=0; i<rows; i++)
		  {
			  pSrcAcc->toPixel(i, j);
              VERIFYNRV(pSrcAcc.isValid());
		      pixelVal = Service<ModelServices>()->getDataValue(type, pSrcAcc->getColumn(), COMPLEX_MAGNITUDE, 0);

			  *(pBuffer + nCount) = pixelVal;
			  nCount++;

			  if (maxVal < pixelVal)
				  maxVal = pixelVal;
		  }
	  }

	  for (unsigned int i=0; i<nCount; i++)
	  {
		  *(pBuffer + i) = 1 + (*(pBuffer + i)/maxVal)*NUMBER_OF_BINS;
	  }

	  ImageSegmentation(pBuffer, pfu, pfVecParameters);
   }

   template<typename T>
   void restoreSegmentationValue(T* pData, float *pSrc)
   {
	   unsigned char pixelVal = 0;
	   float dVal = *pSrc;

	   if (pSrc != NULL)
	   {
		   if (*pSrc > EDGE_THRESHOLD_VALUE)
		   {
			   pixelVal = 1;
		   }
	   }

	   *pData = static_cast<T>(pixelVal);
   }
};

TextureSegmentation::TextureSegmentation()
{
   setDescriptorId("{21860EDB-3761-4e7f-B4DF-169369576749}");
   setName("SAR Image Segmentation");
   setDescription("Segmentation for SAR");
   setCreator("Yiwei Zhang");
   setVersion("Sample");
   setCopyright("Copyright (C) 2008, Ball Aerospace & Technologies Corp.");
   setProductionStatus(false);
   setType("Sample");
   setSubtype("Segmentation");
   setMenuLocation("[SAR]/Segmentation");
   setAbortSupported(true);

   pBuffer = NULL;
}

TextureSegmentation::~TextureSegmentation()
{
	if (pBuffer != NULL)
		free(pBuffer);
}

bool TextureSegmentation::getInputSpecification(PlugInArgList*& pInArgList)
{
   VERIFY(pInArgList = Service<PlugInManagerServices>()->getPlugInArgList());
   pInArgList->addArg<Progress>(Executable::ProgressArg(), NULL, "Progress reporter");
   pInArgList->addArg<RasterElement>(Executable::DataElementArg(), "Perform segmentation on this data element");
   return true;
}

bool TextureSegmentation::getOutputSpecification(PlugInArgList*& pOutArgList)
{
   VERIFY(pOutArgList = Service<PlugInManagerServices>()->getPlugInArgList());
   pOutArgList->addArg<RasterElement>("Result", NULL);
   return true;
}

bool TextureSegmentation::execute(PlugInArgList* pInArgList, PlugInArgList* pOutArgList)
{
   StepResource pStep("SAR image segmentation", "app", "CC430C1A-9D8C-4bb5-9254-FCF7EECAFA3C");
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
      "_Segmentation_Result", pDesc->getRowCount(), pDesc->getColumnCount(), ResultType));
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

   if (NULL != pBuffer)
   {
	   free(pBuffer);
   }
   pBuffer = (float *)malloc(sizeof(float)*pDesc->getRowCount()*pDesc->getColumnCount());
  
   MakeSegmentation(pSrcAcc, pBuffer, pBuffer, pDesc->getRowCount(), pDesc->getColumnCount(), pDesc->getDataType());

   //Output the value 
   unsigned int nCount = 0;
   for (unsigned int j = 0; j < pDesc->getColumnCount(); j++)
   {
       for (unsigned int i = 0; i < pDesc->getRowCount(); i++)		   
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
			   
		   pDestAcc->toPixel(i, j);		   
		   switchOnEncoding(ResultType, restoreSegmentationValue, pDestAcc->getColumn(), (pBuffer+nCount));
		   nCount++;
	   }
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
      pProgress->updateProgress("Image segmentation is compete.", 100, NORMAL);
   }

   pOutArgList->setPlugInArgValue("Image segmentation result", pResultCube.release());

   pStep->finalize();
   return true;
}
