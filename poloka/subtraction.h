// This may look like C code, but it is really -*- C++ -*-
#ifndef SUBTRACTION__H
#define SUBTRACTION__H

#include <poloka/reducedimage.h>
class KernelFitter;

void MakeSubtractionWeight(const ReducedImage& Ref, ReducedImage& Im, const KernelFitter& Kfit);
void MakeSubtraction(const ReducedImage& Ref, ReducedImage& Im, const KernelFitter& Kfit);
void MakeSubtraction(const ReducedImage& Ref, ReducedImage& Im);

struct Subtractor {
  bool overWrite, noSwap, doSub, doDetect;
  ReducedImageRef Ref;
  string subName;
  Subtractor()
    : overWrite(false), noSwap(false), doSub(true), doDetect(false) {}
  bool operator () (ReducedImageRef Im) const;
};
  
#endif // SUBTRACTION__H
