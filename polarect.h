#ifndef POLARECT_H
#define POLARECT_H

#include "libOrsa/libNumerics/matrix.h"
#include "libOrsa/match.hpp"
#include <utility>

bool orientedEpipoles(const std::vector<Match>& m,
                      const libNumerics::matrix<double>& F,
                      libNumerics::vector<double>& eL,
                      libNumerics::vector<double>& eR);

/// Transform image in rectangle [0,w]x[0,h] to polar.
class Polarizer {
     /// Center of polar transform (3-vector, homogeneous coordinates)
    libNumerics::vector<double> c;
    double r, R; ///< Min/max radii
    double t, T; ///< Min/max angles
    bool rotate; ///< Image needs rotation

    bool inf() const { return c(2)==0; } ///< Center at infinity?
    double dist_rect(int w, int h) const;
    std::pair<int,int> region(int w, int h) const;
    std::pair<double,double> transfer_angle(const Polarizer& P,
        double theta, const libNumerics::matrix<double>* F) const;
    void pullback(int w, int h, std::pair<double,double>*& pb,
        const Polarizer& pol, const libNumerics::matrix<double>* F) const;
public:
    Polarizer(const libNumerics::vector<double>& center, int w, int h);
    int width() const;
    int height() const;
    std::pair<double,double> polar(double x, double y) const;
    void restrict_angles(const Polarizer& polR,
                         const libNumerics::matrix<double>& F);

    /// Generate pullback map to transform image to polar.
    void pullback_map(int w, int h, std::pair<double,double>*& p) const
    { pullback(w,h,p, *this, 0); }
    void pullback_map(int w, int h, std::pair<double,double>*& p,
        const Polarizer& pol, const libNumerics::matrix<double>& F) const
    { pullback(w,h,p, pol, &F); }
};

class Polarectifyer {
    Polarizer polL;
    Polarizer polR;
    libNumerics::matrix<double> F;
    std::pair<double,double> *pbMapL, *pbMapR;
public:
    Polarectifyer(const libNumerics::matrix<double>& F,
                  const libNumerics::vector<double>& eL,
                  const libNumerics::vector<double>& eR,
                  int wL, int hL, int wR, int hR);
    ~Polarectifyer();
    int widthL() const { return polL.width(); }
    int widthR() const { return polR.width(); }
    int height() const { return polL.height(); }
    const std::pair<double,double>* pullback_map(bool left);
};

#endif
