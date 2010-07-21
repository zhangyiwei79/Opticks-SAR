#ifndef ZYW_WAVELET_SPECKLE_REMOVE_H
#define ZYW_WAVELET_SPECKLE_REMOVE_H

#include "ExecutableShell.h"

class WaveletSpeckleRemove : public ExecutableShell
{
public:
   WaveletSpeckleRemove();
   virtual ~WaveletSpeckleRemove();

   virtual bool getInputSpecification(PlugInArgList*& pInArgList);
   virtual bool getOutputSpecification(PlugInArgList*& pOutArgList);
   virtual bool execute(PlugInArgList* pInArgList, PlugInArgList* pOutArgList);

   double *pBuffer;
   unsigned int rowBlocks;
   unsigned int colBlocks;
};

#endif
