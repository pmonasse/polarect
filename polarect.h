/**
 * SPDX-License-Identifier: GPL-2.0+
 * @file polarect.h
 * @brief Image pair polar rectification.
 *
 * Copyright (c) 2020 Pascal Monasse
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

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
    friend class Polarectifyer;
     /// Center of polar transform (3-vector, homogeneous coordinates)
    libNumerics::vector<double> c;
    double r, R; ///< Min/max radii
    double t, T; ///< Min/max angles
    bool rotate; ///< Image needs rotation

    bool inf() const { return c(2)==0; } ///< Center at infinity?
    double dist_rect(int w, int h) const;
    std::pair<int,int> region(int w, int h) const;
    std::pair<double,double> transfer_theta(double theta,
        const Polarizer& P, const libNumerics::matrix<double>* F) const;
    std::pair<double,double>* pullback(int w, int h,
        const Polarizer& pol, const libNumerics::matrix<double>* F);
public:
    Polarizer(const libNumerics::vector<double>& center, int w, int h);
    int width() const;
    int height() const;
    std::pair<double,double> polar(double x, double y) const;
    double sample_rho(double rho) const;
    double sample_theta(double theta) const;
    void restrict_angles(const Polarizer& polR,
                         const libNumerics::matrix<double>& F);

    /// Generate pullback map to transform image to polar.
    std::pair<double,double>* pullback_map(int w, int h)
    { return pullback(w,h,*this, 0); }
    std::pair<double,double>* pullback_map(int w, int h,
        const Polarizer& pol, const libNumerics::matrix<double>& F)
    { return pullback(w,h,pol, &F); }
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
    const std::pair<double,double>* pullback_map(bool left) const;
    void pushforward(std::pair<double,double>& xy, bool left) const;
};

#endif
