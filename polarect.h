#ifndef POLARECT_H
#define POLARECT_H

#include "libOrsa/libNumerics/matrix.h"
#include "libOrsa/match.hpp"
#include <utility>

bool orientedEpipoles(const std::vector<Match>& m,
                      const libNumerics::matrix<double>& F,
                      libNumerics::vector<double>& eL,
                      libNumerics::vector<double>& eR);

void polarect(const libNumerics::matrix<double>& F,
              const libNumerics::vector<double>& eL,
              const libNumerics::vector<double>& eR,
              int& wL, int& hL,int& wR, int& hR,
              std::pair<double,double>*& pullbackL,
              std::pair<double,double>*& pullbackR);

#endif
