#ifndef ZYW_Edge_Detector_H
#define ZYW_Edge_Detector_H

#include "ExecutableShell.h"

class EdgeDetector : public ExecutableShell
{
public:
   EdgeDetector();
   virtual ~EdgeDetector();

   virtual bool getInputSpecification(PlugInArgList*& pInArgList);
   virtual bool getOutputSpecification(PlugInArgList*& pOutArgList);
   virtual bool execute(PlugInArgList* pInArgList, PlugInArgList* pOutArgList);
};

#endif
