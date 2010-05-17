#ifndef ZYW_TEST_H
#define ZYW_TEST_H

#include "ExecutableShell.h"

class MyFirstTest : public ExecutableShell
{
public:
   MyFirstTest();
   virtual ~MyFirstTest();

   virtual bool getInputSpecification(PlugInArgList*& pInArgList);
   virtual bool getOutputSpecification(PlugInArgList*& pOutArgList);
   virtual bool execute(PlugInArgList* pInArgList, PlugInArgList* pOutArgList);
};

#endif
