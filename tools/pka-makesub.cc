#include<iostream>

#include <poloka/polokaexception.h>
#include <poloka/sub.h>
 
static void usage(const char *progname) {
  cerr << "Usage: " << progname << " FILE\n"
       << "Resample, stack, fit convolution kernels, and produce subtractions\n\n"
       << "    -o :  overwrite\n\n";
  exit(EXIT_FAILURE);
}

int main(int argc, char **argv) {
  
  if (argc < 2) usage(argv[0]);

  string subfile;
  bool overwrite = false;

  for (int i=1; i< argc; i++)  {
    char *arg = argv[i];
    if (arg[0] != '-')  subfile = arg;
    switch (arg[1]) {
    case 'o': overwrite = true; break;
    default: usage(argv[0]);
    }
  }

  // if sub/matcheddet.list exists the subtraction is normally well finished
  // You don't want to rerun the whole process
  // except when explicitly required with overwrite
  bool status = true;
  try {
    if (!FileExists("sub/matcheddet.list") || overwrite) {
      Sub sub(subfile, overwrite);
      status = sub.DoIt();
    } else
      cout << argv[0] << ": subtractions from " << subfile << " already done\n";
  } catch(PolokaException p) {
    p.PrintMessage(cerr);
    status = false;
  }
  
  return status? EXIT_SUCCESS : EXIT_FAILURE;
}
