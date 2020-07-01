/**
 * @file main.cpp
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

#include <cstdlib>
#include <ctime>

#include <iostream>

#include "libImageIO/image_io.hpp"
#include "libOrsa/eval_model.hpp"
extern "C" {
#include "libSplinter/splinter.h"
}

#include "cmdLine.h"
#include "siftMatch.hpp"
#include "polarect.h"

/// Number of random samples in ORSA
static const int ITER_ORSA=10000;

inline unsigned char clip_col(double c) {
    return (c<0.)? 0: (c>=255.)? 255: (unsigned char)c;
}

inline RGBColor clip_col(double c[3]) {
    return RGBColor(clip_col(c[0]), clip_col(c[1]), clip_col(c[2]));
}

Image<RGBColor>* sample(const Image<RGBColor>& im, int w, int h,
                        const std::pair<double,double>* pullback) {
    size_t s = im.Width()*im.Height();
    double* in = new double[s*3];
    double* pix = in;
    for(size_t y=0; y<im.Height(); y++)
        for(size_t x=0; x<im.Width(); x++, pix++) {
            RGBColor c = im(y,x);
            pix[0*s] = c.r;
            pix[1*s] = c.g;
            pix[2*s] = c.b;
        }

    splinter_plan_t plan = splinter_plan(in, im.Width(), im.Height(), 3,
                                         3, BOUNDARY_PERIODIC, 1.0e-1, 1);
    delete [] in;
    Image<RGBColor>* out = new Image<RGBColor>(w,h,WHITE);
    RGBColor *o=out->data();
    const std::pair<double,double>* p = pullback;
    for(int y=0; y<h; y++)
        for(int x=0; x<w; x++, p++,o++) {
            double x0 = p->first;
            double y0 = p->second;
            if(im.Contains(y0,x0)) {
                double c[3];
                splinter(c, x0, y0, plan);
                *o = clip_col(c);
            }
        }
    splinter_destroy_plan(plan);
    return out;
}

int main(int argc, char **argv) {
    double precision=0;
    float fSiftRatio=0.6f;

    CmdLine cmd;
    cmd.add( make_option('s',fSiftRatio, "sift")
             .doc("SIFT distance ratio of descriptors") );
    cmd.add( make_switch('r', "read")
             .doc("Read file of matches allMatches.txt, do not use SIFT"));
    try {
        cmd.process(argc, argv);
    } catch(const std::string& s) {
        std::cerr << s << std::endl;
        return 1;
    }
    if(argc!=5 && argc!=7) {
        std::cerr << "Usage: " << argv[0] << " [options] imgInA imgInB "
                  << "allMatches.txt orsaMatches.txt [imgOutA imgOutB]\n"
                  << "- imgInA, imgInB: two input images (JPG or PNG format)\n"
                  << "- allMatches.txt: output (input if option -r) text file"
                  << "of format \"x1 y1 x2 y2\"\n"
                  << "- orsaMatches.txt: output, but only with inliers.\n"
                  << "- imgOutA, imgOutB (optional): output rectified image\n"
                  << "\tOptions:\n" << cmd;
        return 1;
    }

    // Init random seed
    srand((unsigned int)time(0));

    // Read images
    Image<RGBColor> image1, image2;
    if(! (libs::ReadImage(argv[1],&image1) && libs::ReadImage(argv[2],&image2)))
        return 1;
    Image<unsigned char> image1Gray, image2Gray;
    libs::convertImage(image1, &image1Gray);
    libs::convertImage(image2, &image2Gray);

    // Find matches with SIFT or read correspondence file
    std::vector<Match> matchings;
    if(cmd.used('r')) {
        if(Match::loadMatch(argv[3], matchings))
            std::cout << "Read " <<matchings.size()<< " matches" <<std::endl;
        else {
            std::cerr << "Failed reading matches from " << argv[3] <<std::endl;
            return 1;
        }
    } else
        SIFT(image1Gray, image2Gray, matchings, fSiftRatio);
    rm_duplicates(matchings); // Remove duplicates (frequent with SIFT)

    // Save match files
    if(! cmd.used('r') && ! Match::saveMatch(argv[3], matchings)) {
        std::cerr << "Failed saving matches into " <<argv[3] <<std::endl;
        return 1;
    }

    // Estimation of fundamental matrix with ORSA
    libNumerics::matrix<double> F(3,3);
    std::vector<int> vec_inliers;
    int w1 = image1Gray.Width(), h1 = image1Gray.Height();
    int w2 = image2Gray.Width(), h2 = image2Gray.Height();
    bool ok = orsa::orsa_fundamental(matchings,w1,h1,w2,h2,precision,ITER_ORSA,
                                     F, vec_inliers);
    if(ok)
        std::cout << "F=" << F <<std::endl;

    // Save inliers
    std::vector<Match> good_match;
    std::vector<int>::const_iterator it = vec_inliers.begin();
    for(; it != vec_inliers.end(); it++)
        good_match.push_back(matchings[*it]);
    if(! Match::saveMatch(argv[4], good_match)) {
        std::cerr << "Failed saving matches into " <<argv[4] <<std::endl;
        return 1;
    }

    // XXXXXXXXXXXX DEBUG
    //    double coeff[3*3] = { -3.87567e-06, 0.000195535, -0.0359479,  -0.000200956, -1.64802e-06, 0.0596173,  0.0380931, -0.0340711, -4.77871};
    //    F.read(coeff);
    // XXXXXXXXXXXX DEBUG
    libNumerics::vector<double> eL(3), eR(3);
    orientedEpipoles(good_match, F, eL, eR);
    std::pair<double,double> *pullbackL, *pullbackR;
    polarect(F, eL, eR, w1, h1, w2, h2, pullbackL, pullbackR);

    if(argc>6) { // Output images
        const char* fileL = argv[5];
        const char* fileR = argv[6];
        Image<RGBColor>* out1 = sample(image1, w1, h1, pullbackL);
        Image<RGBColor>* out2 = sample(image2, w2, h2, pullbackR);
        libs::WriteImage(fileL, *out1);
        libs::WriteImage(fileR, *out2);
        delete out1;
        delete out2;
    }
    delete [] pullbackL;
    delete [] pullbackR;

    if(! ok) {
        std::cerr << "Failed to estimate the fundamental matrix" << std::endl;
        return 1;
    }

    return 0;
}
