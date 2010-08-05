#ifndef ZYW_TEXTURE_SEGMENTATION_H
#define ZYW_TEXTURE_SEGMENTATION_H

#include "ExecutableShell.h"

class TextureSegmentation: public ExecutableShell
{
public:
   TextureSegmentation();
   virtual ~TextureSegmentation();

   virtual bool getInputSpecification(PlugInArgList*& pInArgList);
   virtual bool getOutputSpecification(PlugInArgList*& pOutArgList);
   virtual bool execute(PlugInArgList* pInArgList, PlugInArgList* pOutArgList);

   float *pBuffer;
};

#endif
