/**
 * SPDX-License-Identifier: GPL-2.0+
 * @file polarect_test.cpp
 * @brief Unit test for polarect.
 *
 * Copyright (c) 2021 Pascal Monasse
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

// Unit test for polar rectification:
// - start from dummy square images.
// - generate horizontal matches, so that we start from rectified situation
// - Apply homographies to left and right images so that epipoles move to
//   positions at mid-height of image
// - Generate random rotations arount image center to each image
// - Rectify and measure the average y-error of matches after rectification

#include <cstdlib>
#include <ctime>
#include <iostream>
#include "libOrsa/homography_model.hpp"
#include "polarect.h"

const int DIM=1000; ///< Pixel dimensions of dummy image
const int NMATCHES=100; ///< Number of matches

/// Build homography mapping infinity point (1,0,0) to \a e.
/// It is assumed that e has abscissa less than DIM. Corners at x=DIM are
/// fixed and corners at x=0 map to points such that lines y=0 and y=DIM are
/// mapped to lines containing e.
libNumerics::matrix<double> build_h(const libNumerics::vector<double>& e) {
    assert(e(0)<DIM);
    std::vector<Match> m;
    m.push_back(Match(DIM,0,DIM,0));
    m.push_back(Match(DIM,DIM,DIM,DIM));
    m.push_back(Match(0,0,(DIM+e(0))/2,(0+e(1))/2));
    if(m.back().x2<0) {
        m.back().x2=0;
        m.back().y2=DIM*e(1)/(DIM-e(0));
    }
    m.push_back(Match(0,DIM,(DIM+e(0))/2,(DIM+e(1))/2));
    if(m.back().x2<0) {
        m.back().x2=0;
        m.back().y2=DIM*(e(1)-e(0))/(DIM-e(0));
    }
    orsa::HomographyModel H(m);
    std::vector<orsa::ModelEstimator::Mat> v;
    float ind[4]={0,1,2,3};
    H.Fit(std::vector<int>(ind,ind+4), &v);
    if(v.size()!=1) {
        std::cerr << "Error in building homography" << std::endl;
        std::exit(1);
    }
    return v.back();
}

/// Random rotation matrix around image center
libNumerics::matrix<double> randRot() {
    libNumerics::matrix<double> T = libNumerics::matrix<double>::eye(3);
    T(0,2)=T(1,2)=DIM/2;
    libNumerics::matrix<double> R(3,3);
    R.fill(0);
    double theta=rand()/(double)RAND_MAX*2*M_PI;
    R(0,0)=R(1,1)=cos(theta);
    R(0,1)=-(R(1,0)=sin(theta));
    R(2,2)=1;
    return T*R*T.inv();
}

/// Map by H1 the left point of each match of m and by H2 the right point.
void transform(std::vector<Match>& m,
               const libNumerics::matrix<double>& H1,
               const libNumerics::matrix<double>& H2) {
    for(std::vector<Match>::iterator it=m.begin(); it!=m.end(); ++it) {
        libNumerics::vector<double> p(3);
        p(0)=it->x1; p(1)=it->y1; p(2)=1;
        p = H1*p;
        it->x1=p(0)/p(2); it->y1=p(1)/p(2);
        p(0)=it->x2; p(1)=it->y2; p(2)=1;
        p = H2*p;
        it->x2=p(0)/p(2); it->y2=p(1)/p(2);
    }
}

/// Rectify matches and return mean y-error after rectification.
double test_rect(const std::vector<Match>& m,
                 const libNumerics::matrix<double>& F) {
    libNumerics::vector<double> eL(3), eR(3);
    if(! orientedEpipoles(m, F, eL, eR)) {
        std::cerr << "Error in finding oriented epipoles" << std::endl;
        std::exit(1);
    }
    Polarectifyer pol(F, eL, eR, DIM, DIM, DIM, DIM);

    double error=0;
    for(std::vector<Match>::const_iterator it=m.begin(); it!=m.end(); ++it) {
        std::pair<double,double> p1(it->x1, it->y1);
        std::pair<double,double> p2(it->x2, it->y2);
        pol.pushforward(p1,true);
        pol.pushforward(p2,false);
        error += std::abs(p1.second-p2.second);
    }
    return error/m.size();
}

int main() {
    srand((unsigned int)time(0));

    // Generate random matches in rectified images
    std::vector<Match> matchings;
    for(int i=0; i<NMATCHES; i++) {
        float x1= rand()/(float)RAND_MAX*DIM;
        float x2= rand()/(float)RAND_MAX*DIM;
        float y = rand()/(float)RAND_MAX*DIM;
        matchings.push_back(Match(x1,y,x2,y));
    }

    // Basic rectified fundamental matrix
    libNumerics::matrix<double> F(3,3);
    F.fill(0); F(1,2)=1; F(2,1)=-1;

    // Generate homographies on left and right images and test rectification
    for(int i=0; i<4; i++) {
        libNumerics::vector<double> e1(2); // Left epipole (will be rotated)
        e1(0)=DIM/2-pow(100,i)*DIM/4; e1(1)=DIM/2;
        for(int j=0; j<4; j++) {
            libNumerics::vector<double> e2(2); // Right epipole
            e2(0)=DIM/2-pow(100,j)*DIM/4; e2(1)=DIM/2;
            libNumerics::matrix<double> H1=randRot()*build_h(e1);
            libNumerics::matrix<double> H2=randRot()*build_h(e2);
            std::vector<Match> m(matchings);
            transform(m,H1,H2);
            double error = test_rect(m, H1.inv().t()*F*H2.inv());
            std::cout << "Mean error: " << error << std::endl;
            if(error > 1.5) {
                std::cerr << "ERROR TOO LARGE:" << std::endl;
                std::cerr << "H1=" << H1 << std::endl;
                std::cerr << "H2=" << H2 << std::endl;
                return 1;
            }
        }
    }

    return 0;
}
