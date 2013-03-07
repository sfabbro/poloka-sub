#include <algorithm>

#include <poloka/polokaexception.h>
#include <poloka/subtraction.h>

static void usage(const char *progname) {
  cerr << "Usage: " << progname << " [OPTION]... DBIMAGE DBIMAGE...\n"
       << "PSF match and subtract first DBIMAGE (reference) to each successive DBIMAGE\n\n"
       << "   -d        : perform candidate detection on each subtracted image\n"
       << "   -f        : force overwrite\n"
       << "   -n        : do not subtract, do PSF match only\n"
       << "   -o DBIMAGE: output subtraction images in a different DBIMAGE\n"
       << "   -r        : reference image will be forced to be one to convolve (default is best seeing one)\n\n";
  exit(EXIT_FAILURE);
}

int main(int nargs, char **args) {

  if (nargs < 2) usage(args[0]);

  ReducedImageList imList;
  Subtractor imSubtract;
  bool status = true;

  for (int i=1; i<nargs; ++i) {
    char *arg = args[i];
    if (arg[0] != '-') {
      ReducedImageRef im = new ReducedImage(arg);
      if (!im || !im->IsValid()) { 
	cerr << args[0] << ": " << arg << " is not a valid dbimage\n";
	status = false;
	continue;
      }
      im->Execute(DoFits|DoCatalog);
      imList.push_back(im);
      continue;
    }
    switch (arg[1]) {
    case 'd': imSubtract.doDetect = true; break;
    case 'r': imSubtract.noSwap = true; break;
    case 'n': imSubtract.doSub = false; break;
    case 'f': imSubtract.overWrite = true; break;
    case 'o': imSubtract.subName = args[++i]; break;
    default : usage(args[0]);
    }
  }

  if (imList.size() < 2) {
    cerr << args[0] << ": needs an image to subtract reference from\n";
    return EXIT_FAILURE;
  }

  imSubtract.Ref = imList.front();
  imList.pop_front();

  for_each(imList.begin(), imList.end(), imSubtract);

  return status? EXIT_SUCCESS : EXIT_FAILURE;
}
