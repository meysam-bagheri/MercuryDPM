#include "muICal3DParameters.h"


int main(int argc, char **argv) {
  if (argc == 1)
    muICal3DParameters::presets()->write(stdout);
  else {
    auto paramsFileName = argv[1];
    auto params = muICal3DParameters::read(paramsFileName);
    params->write(stdout);
  }
}
