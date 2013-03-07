#include <limits>

#include <poloka/subtraction.h>
#include <poloka/kernelfitter.h>
#include <poloka/polokaexception.h>
#include <poloka/fitsimage.h>
#include <poloka/polokaconf.h>
#include <poloka/datacards.h>
#include <poloka/imageback.h>
#include <poloka/apersestar.h>
#include <poloka/detection.h>

static double sqr(const double& x) { return x*x; }

static void addStarVariance(const SEStar *Star, const double& Gain, Image &Variance) {
  double mxx, myy, mxy;

  const AperSEStar* ap = dynamic_cast<const AperSEStar*>(Star);
  if (ap) {
    mxx = ap->gmxx;
    myy = ap->gmyy;
    mxy = ap->gmxy;
  } else {
    mxx = Star->Mxx();
    myy = Star->Myy();
    mxy = Star->Mxy();
  }
  
  // account for object noise using a gaussian "model"
  double det = mxx * myy - sqr(mxy);
  double wxx = myy/det;
  double wyy = mxx/det;
  double wxy = mxy/det;
  double factor = Star->flux/Gain * (2*M_PI*sqrt(det));
  double hwidth = 3*Star->A() + 2;
  int xstart = max(0, int(floor(Star->x - hwidth)));
  int ystart = max(0, int(floor(Star->y - hwidth)));
  int xend = min(int(ceil(Star->x + hwidth)), Variance.Nx());
  int yend = min(int(ceil(Star->y + hwidth)), Variance.Ny());

  for (int j=ystart; j<yend; ++j) {
    double y = j - Star->y;
    for (int i=xstart; i<xend; ++i) {
      double x = i - Star->x;
      Variance(i,j) += factor * exp(-0.5*(wxx*x*x + wyy*y*y + 2*wxy*x*y));
    }
  }
}

static Image* loadVariance(const ReducedImage& Im, const int VarianceType) {

  // initialize with current weight
  FitsImage *var = new FitsImage(Im.FitsWeightName());

  // add a small constant to weights so that variances of 
  // zero weight pixels remain finite and go to variances
  const Pixel *pvend = var->end();
  Pixel eps = numeric_limits<Pixel>::min();

  // mask saturated pixels
  if (Im.HasSatur()) {
    FitsImage satur(Im.FitsSaturName());
    *var *= (1 - satur);
  }

  for (Pixel *pv = var->begin(); pv < pvend; ++pv)
    *pv = 1. / (*pv + eps);

  // Variance types:
  //   0: no extra variance var(sub) = sky(best) * k^2 + sky(worst) default
  //   1: var(sub) = k^2 * (sky(best) + best) + sky(worst) + worst
  //      includes Poisson noise from objects as "data": bias photometry but useful when many bright stars
  //   2: var(sub) = k^2 * (sky(best) + model(best)) + sky(worst) + model(worst)
  //      includes poisson noise from objects as Gaussian model. Might be less biased than 
  //	  previous, but usually very wrong for crowded fields and sextractor based catalogs

  Pixel factfmax = 10.;
  if (VarianceType == 1) {
    cout << " loadVariance: adding pixel intensity as Poisson noise for " << Im.Name() << endl;
    Pixel fact = 1. / Im.Gain();
    Pixel back = Im.BackLevel();
    FitsImage im(Im.FitsName());
    for (Pixel *pim=im.begin(), *pv=var->begin(); pv < pvend; ++pim, ++pv)
      *pv += (*pim - back) * fact;
  } else if (VarianceType == 2) {
    cout << " loadVariance: adding Gaussian star models as Poisson noise for " << Im.Name() << endl;
    Pixel gain = Im.Gain();
    Pixel minfmax = factfmax * sqr(Im.SigmaBack());
    SEStarList imList;
    if (Im.HasAperCatalog()) 
      imList.read(Im.AperCatalogName());
    else 
      imList.read(Im.CatalogName());
    for (SEStarCIterator it=imList.begin(); it != imList.end(); ++it)
      if ((*it)->Fluxmax() > minfmax) 
	addStarVariance(*it, gain, *var);
  }

  return var;
}

void varToWeight(Image& Var, const Frame& Borders) {

  Pixel threshold = numeric_limits<Pixel>::max()*0.9;
  int nx = Var.Nx();
  int ny = Var.Ny();

  for (int j=0; j<ny; ++j)
    for (int i=0; i<nx; ++i) {
      Pixel *pw = &Var(i,j);
      if (j<Borders.yMin || 
	  j>Borders.yMax || 
	  i<Borders.xMin || 
	  i>Borders.xMax || 
	  *pw> threshold)
	*pw = 0;
      else
	*pw = 1. / *pw;      
    }
}

// mask = satur * common_frame * mask_bright_stars * cosmics
// weight = mask / ( kern^2*var(ref) + var(im) )
void MakeSubtractionWeight(const ReducedImage& Ref, ReducedImage& Im, const KernelFitter& Kfit) {

  string fileName = Im.FitsSubWeightName();
  if (FileExists(fileName)) return;

  cout << " MakeSubtractionWeight: making weight image " << fileName << endl;

  int subVarianceType = 0;
  int subMaskDilate = 0;
  Pixel subMaskSatur = 1;    
  string datacards = DefaultDatacards("sub.conf");
  if (FileExists(datacards)) {
    DataCards cards(datacards);
    if (cards.HasKey("SUB_VARIANCE_TYPE"))
      subVarianceType = cards.IParam("SUB_VARIANCE_TYPE");
    if (cards.HasKey("SUB_MASK_DILATE"))
      subMaskDilate = cards.IParam("SUB_MASK_DILATE");
    if (cards.HasKey("SUB_MASK_SATUR"))
      subMaskSatur = cards.DParam("SUB_MASK_SATUR");
  }

  Image *var = loadVariance(Im, subVarianceType);
  {
    Image *refVar = loadVariance(Ref, subVarianceType);
    Image refVarConv(refVar->Nx(), refVar->Ny());
    Kfit.VarianceConvolve(*refVar, refVarConv);
    *var += refVarConv;
    delete refVar;
  }
  

  FitsHeader head(Im.FitsName());
  FitsImage subWeight(fileName, head, *var);
  delete var;
  subWeight.PreserveZeros();
  // set to 0 the "side bands" where the convolution did not go
  // size of the dead band come from variance convolution
  Frame borders(Kfit.CommonFrame());
  Kernel kern;
  Kfit.KernAllocateAndCompute(kern, subWeight.Nx()/2., subWeight.Ny()/2);
  borders.CutMargin(kern.Nx()/2 -1, kern.Ny()/2 -1);   
  varToWeight(subWeight, borders);
#if 0
  FitsImage subWeightFits(FitsWeightName(), RW);
  if (!subWeightFits.IsValid())
    {
      cerr << " ImageSubtraction: missing weight to count cosmics\n";
      return false;
    }
  subFits *= subWeightFits;
  Pixel mean, sigma;
  subFits.SkyLevel(&mean, &sigma);
  FitsImage subCosmicFits(FitsCosmicName(), (FitsHeader&)subFits);;
  subFits.Cosmics(sigma, mean, Seeing(), subCosmicFits);
  subCosmicFits.AddOrModKey("BITPIX",8);
  // update weights
  subWeightFits *= (1.- subCosmicFits);
  subWeightFits.AddOrModKey("COSMPIXS", true, "This weight accounts for (identified) cosmics");
  return true;
  // check weight validity
  float wmin, wmax;
  subWeight.MinMaxValue(&wmin,&wmax);
  if (fabs(wmin - wmax) < 1e-10) 
    cerr << " MakeSubtractionWeight: " << Im.Name() << " weight seems to be zero everywhere\n";
#endif
}

void MakeSubtraction(const ReducedImage& Ref, ReducedImage& Im) {

  KernelFitter kernFit(&Ref, &Im, true);
  if (!kernFit.ReadKernel()) {
    kernFit.DoTheFit();
    kernFit.WriteKernel();
  }
  MakeSubtraction(Ref, Im, kernFit);
}

void MakeSubtraction(const ReducedImage& Ref, ReducedImage& Im, const KernelFitter& KernFit) {

  string fileName = Im.FitsSubName();
  if (FileExists(fileName)) return;

  if (!Ref.HasImage())
    throw(PolokaException("MakeSubtraction: missing FITS image for reference " + Ref.Name()));

  if (!Im.MakeFits())
    throw(PolokaException("MakeSubtraction: could not produce " + Im.Name()));

  cout << " MakeSubtraction: " << Im.Name()  << " - Kernel*" << Ref.Name() << endl;

  // produce subtraction
  FitsImage im(Im.FitsName());
  {
    FitsImage ref(Ref.FitsName());
    Image refConv(ref.Nx(), ref.Ny());
    KernFit.ImageConvolve(ref, refConv);
    im -= refConv;
  }
  
  MakeSubtractionWeight(Ref, Im, KernFit);
  int backMesh = 64;

  {
    FitsImage subWeight(Im.FitsSubWeightName());

    // subtract background
    string datacards = DefaultDatacards("sub.conf");
    if (FileExists(datacards)) {
      DataCards cards(datacards);
      if (cards.HasKey("SUB_BACK_MESHSTEP"))
	backMesh = cards.IParam("SUB_BACK_MESHSTEP");
    }
    {
      ImageBack b(im, backMesh, &subWeight);
      Image *back = b.BackgroundImage();
      im -= *back;
      // back is replaced
      FitsImage backFits(Im.FitsBackName(), im, *back);
      backFits.AddOrModKey("BITPIX", 32);
      delete back;
    }
  }
  
  // update header
  FitsImage sub(fileName, (const FitsHeader&)im, im);
  KernFit.CommonFrame().WriteInHeader(sub);
  double photomRatio = KernFit.PhotomRatio();
  sub.ModKey("BITPIX", 16);
  sub.AddOrModKey("KERNREF" , Ref.Name(), "name of the seeing reference image");
  sub.AddOrModKey("KERNCHI2", KernFit.Chi2(), "chi2/dof of the kernel fit");
  sub.AddOrModKey("PHORATIO", photomRatio, "flux = PHORATIO * flux(KERNREF)");
  sub.AddOrModKey("BACK_SUB", true, "subtracted weighed background");
  sub.AddOrModKey("BACKMESH", backMesh, "mesh size used for back computation ");
  sub.AddOrModKey("BACKLEV" , 0, "");
  sub.AddOrModKey("SKYSIGMA", sqrt(sqr(Ref.SigmaBack()*photomRatio) + sqr(Im.SigmaBack())), "");
  sub.AddOrModKey("TOADRDON", (sqrt(sqr(photomRatio*Ref.ReadoutNoise()) + sqr(Im.ReadoutNoise()))),"");
  sub.AddOrModKey("SKYLEV"  , (photomRatio*Ref.OriginalSkyLevel() + Im.OriginalSkyLevel()),
		  "sum of sky levels of the subtraction terms");
}

bool Subtractor::operator () (ReducedImageRef Im) const {

  // read or redo kernel fit
  bool status = true;
  try {
    KernelFitter kernFit(Ref, Im, noSwap);
    if (overWrite || !kernFit.ReadKernel())
      if (kernFit.DoTheFit())
	kernFit.WriteKernel(overWrite);
      else {
	cerr << " Subtractor: kernel fit between "
	     << Ref->Name() << " and "
	     << Im->Name() << " failed\n";
	return false;
      }
    if (!doSub) return true;
    
    if (overWrite) {
      string fileName = Im->FitsSubName();
      if (FileExists(fileName)) remove(fileName.c_str());
      fileName = Im->FitsSubWeightName();
      if (FileExists(fileName)) remove(fileName.c_str());
      fileName = Im->ImageCatalogName(Subtraction);
      if (FileExists(fileName)) remove(fileName.c_str());
    }

    MakeSubtraction(*Ref, *Im, kernFit);
    if (doDetect) {
      DetectionList detections;
      double sigma = Im->Seeing();
      status = (ImageDetect(Im->FitsSubName(), Im->FitsSubWeightName(), sigma, detections, Ref) &&
	detections.write(Im->ImageCatalogName(Subtraction)) == 1);
    }
  } catch (PolokaException p) {
    p.PrintMessage(cerr);
    status = false;
  }
  return status;
}

