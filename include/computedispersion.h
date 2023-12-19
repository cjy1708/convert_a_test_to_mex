#ifndef __computedispersion_h
#define __computedispersion_h
#include "fiberbundle.h"


int
computeDispersion(fiberbundle& bundle,
                  double scale,
                  unsigned int numberOfSamplingDirections,
                  const std::string& outputFilename,
                  unsigned int tractSubSampling = 1,
                  unsigned int fiberPointSubSampling = 1);

#endif // __computedispersion_h

