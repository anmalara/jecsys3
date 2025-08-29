// Purpose: Run Global Fit
// Authors: Andrea Malara,  Mikko Voutilainen

#include "GlobalFit.hpp"
using namespace std;

int main(int argc, char **argv) {

  if (argc != 2)
    throw invalid_argument("Only 1 parameter needed: json file path.");

  unique_ptr<GlobalFit> GF(new GlobalFit(argv[1]));
  GF->Run();

  return 0;
}
