#include <poloka/fitsimage.h>
#include <poloka/reducedimage.h>

#include <poloka/imagepair.h>


ImagePair::ImagePair(const ReducedImageRef &Best, const ReducedImageRef &Worst) : best(Best), worst(Worst)
{
  bestSeeing = best->Seeing();
  worstSeeing = worst->Seeing();
  bestGain = best->Gain();
  worstGain = worst->Gain();
  commonFrame = best->UsablePart()*worst->UsablePart();
  for (unsigned k=0; k<4; ++k) images[k] = NULL;
}


const Image& ImagePair::BestImage()
{
  if (images[0] == NULL) images[0] = new FitsImage(best->FitsName());
  if (images[0] == NULL) { 
    FitsImage im(best->FitsName());
    im -= best->BackLevel();    
    images[0] = new Image(im);
  }  
  return *images[0];
}

const Image& ImagePair::WorstImage()
{
  if (images[1] == NULL) images[1] = new FitsImage(worst->FitsName());
  if (images[1] == NULL) {
    FitsImage im(worst->FitsName());
    im -= worst->BackLevel();    
    images[1] = new Image(im);
  }  
  return *images[1];
}

static Image* weight_with_satur(const ReducedImageRef &Ri)
{
  Image *res = new FitsImage(Ri->FitsWeightName());
  if (Ri->HasSatur()) {
    FitsImage satur(Ri->FitsSaturName());
    *res *= (1.-satur);
  }
  return res;
}

const Image& ImagePair::BestWeight()
{
  if (images[2] == NULL) images[2] = weight_with_satur(best);
  return *images[2];
}

const Image& ImagePair::WorstWeight()
{
  if (images[3] == NULL) images[3] = weight_with_satur(worst);
  return *images[3];
}

ImagePair::~ImagePair()
{
  //"delete NULL" is safe .
  for (unsigned k=0; k<4; ++k) delete images[k];
}
