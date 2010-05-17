/*
 * The information in this file is
 * Copyright(c) 2007 Ball Aerospace & Technologies Corporation
 * and is subject to the terms and conditions of the
 * GNU Lesser General Public License Version 2.1
 * The license text is available from   
 * http://www.gnu.org/licenses/lgpl.html
 */

#include "AppConfig.h"
#include "AppVerify.h"
#include "DataAccessor.h"
#include "DataAccessorImpl.h"
#include "DataRequest.h"
#include "MessageLogResource.h"
#include "ObjectResource.h"
#include "PlugInArgList.h"
#include "PlugInManagerServices.h"
#include "PlugInRegistration.h"
#include "PlugInResource.h"
#include "Progress.h"
#include "RasterDataDescriptor.h"
#include "RasterElement.h"
#include "StringUtilities.h"
#include "switchOnEncoding.h"
#include "MyFirstTest.h"
#include <limits>

REGISTER_PLUGIN_BASIC(OpticksTutorial, MyFirstTest);

namespace
{
   template<typename T>
   void updateStatistics(T* pData, double& total)
   {
      total += *pData;
   }
};

MyFirstTest::MyFirstTest()
{
   setDescriptorId("{C5B6AF80-1E0B-4751-A5DC-B77837797A75}");
   setName("My First Test");
   setDescription("Calculate mean.");
   setCreator("Yiwei Zhang");
   setVersion("Sample");
   setCopyright("Copyright (C) 2008, Ball Aerospace & Technologies Corp.");
   setProductionStatus(false);
   setType("Sample");
   setSubtype("Statistics");
   setMenuLocation("[Tutorial]/My First Test");
   setAbortSupported(true);
}

MyFirstTest::~MyFirstTest()
{
}

bool MyFirstTest::getInputSpecification(PlugInArgList* &pInArgList)
{
   VERIFY(pInArgList = Service<PlugInManagerServices>()->getPlugInArgList());
   pInArgList->addArg<Progress>(Executable::ProgressArg(), NULL, "Progress reporter");
   pInArgList->addArg<RasterElement>(Executable::DataElementArg(), "Calculate mean for this raster element");
   return true;
}

bool MyFirstTest::getOutputSpecification(PlugInArgList*& pOutArgList)
{
   VERIFY(pOutArgList = Service<PlugInManagerServices>()->getPlugInArgList());
   pOutArgList->addArg<double>("Mean", "The average value");
   return true;
}

bool MyFirstTest::execute(PlugInArgList* pInArgList, PlugInArgList* pOutArgList)
{
   StepResource pStep("My First Test", "app", "0FD3C564-041D-4f8f-BBF8-96A7A165AB61");
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

   FactoryResource<DataRequest> pRequest;
   pRequest->setInterleaveFormat(BSQ);
   DataAccessor pAcc = pCube->getDataAccessor(pRequest.release());

   double total = 0.0;
   for (unsigned int row = 0; row < pDesc->getRowCount(); ++row)
   {
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
      if (!pAcc.isValid())
      {
         std::string msg = "Unable to access the cube data.";
         pStep->finalize(Message::Failure, msg);
         if (pProgress != NULL)
         {
            pProgress->updateProgress(msg, 0, ERRORS);
         }

         return false;
      }

      if (pProgress != NULL)
      {
         pProgress->updateProgress("Calculating statistics", row * 100 / pDesc->getRowCount(), NORMAL);
      }

      for (unsigned int col = 0; col < pDesc->getColumnCount(); ++col)
      {
         switchOnEncoding(pDesc->getDataType(), updateStatistics, pAcc->getColumn(), total);
         pAcc->nextColumn();
      }
      pAcc->nextRow();
   }
   unsigned int count = pDesc->getColumnCount() * pDesc->getRowCount();
   double mean = total / count;

   if (pProgress != NULL)
   {
      std::string msg =  "Average: " + StringUtilities::toDisplayString(mean);
      pProgress->updateProgress(msg, 100, NORMAL);
   }

   pStep->addProperty("Mean", mean);
   
   pOutArgList->setPlugInArgValue("Mean", &mean);

   pStep->finalize();
   return true;
}