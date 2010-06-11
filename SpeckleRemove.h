#ifndef ZYW_SPECKLE_REMOVE_H
#define ZYW_SPECKLE_REMOVE_H

#include "ExecutableShell.h"

class SpeckleRemove : public ExecutableShell
{
public:
   SpeckleRemove();
   virtual ~SpeckleRemove();

   virtual bool getInputSpecification(PlugInArgList*& pInArgList);
   virtual bool getOutputSpecification(PlugInArgList*& pOutArgList);
   virtual bool execute(PlugInArgList* pInArgList, PlugInArgList* pOutArgList);
};

#endif
