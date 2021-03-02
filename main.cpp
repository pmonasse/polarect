/**
 * SPDX-License-Identifier: GPL-2.0+
 * @file main.cpp
 * @brief Image pair polar rectification.
 *
 * Copyright (c) 2020-2021 Pascal Monasse
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
#include "libOrsa/fundamental_model.hpp"

#include "cmdLine.h"
#include "siftMatch.hpp"
#include "polarect.h"
#include "draw.h"

/// Number of random samples in ORSA
static const int ITER_ORSA=10000;

int main(int argc, char **argv) {
    double precision=0;
    float fSiftRatio=0.6f;
    std::string fileMatchesIn, fileMatchesOut;
    bool findInliers=false;
    bool annotate=false;

    CmdLine cmd;
    cmd.add( make_option('s',fSiftRatio, "sift")
             .doc("SIFT distance ratio of descriptors") );
    cmd.add( make_option('r',fileMatchesIn, "read")
             .doc("Read file of matches, do not use SIFT"));
    cmd.add( make_option('w',fileMatchesOut, "write")
             .doc("Write file of inlier matches from SIFT algorithm"));
    cmd.add( make_option('i',findInliers, "inliers")
             .doc("Find inlier matches from file of matches (with option -r)"));
    cmd.add( make_option('a',annotate, "annotate")
             .doc("Draw annotations (guidelines and SIFT matches)"));
    try {
        cmd.process(argc, argv);
    } catch(const std::string& s) {
        std::cerr << s << std::endl;
        return 1;
    }
    if(argc!=5) {
        std::cerr << "Usage: " << argv[0] << " [options] imgInA imgInB "
                  << "imgOutA imgOutB\n"
                  << "- imgInA, imgInB: two input images (JPG or PNG format)\n"
                  << "- imgOutA, imgOutB: output rectified images\n"
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
    const int w1 = image1Gray.Width(), h1 = image1Gray.Height();
    const int w2 = image2Gray.Width(), h2 = image2Gray.Height();

    // Find matches with SIFT or read correspondence file
    std::vector<Match> matchings;
    if(fileMatchesIn.empty())
        SIFT(image1Gray, image2Gray, matchings, fSiftRatio);
    else if(Match::loadMatch(fileMatchesIn.c_str(), matchings))
        std::cout << "Read " <<matchings.size()<< " matches" <<std::endl;
    else {
        std::cerr<<"Failed reading matches from "<<fileMatchesIn<<std::endl;
        return 1;
    }

    // Estimation of fundamental matrix...
    libNumerics::matrix<double> F(3,3);
    std::vector<int> vec_inliers;
    bool ok=true;
    if(fileMatchesIn.empty() || findInliers) { // ...with ORSA
        rm_duplicates(matchings); // Remove duplicates (frequent with SIFT)
        if((ok=orsa::orsa_fundamental(matchings,w1,h1,w2,h2,precision,ITER_ORSA,
                                      F, vec_inliers))) {
            std::vector<Match> good_match;
            std::vector<int>::const_iterator it = vec_inliers.begin();
            for(; it != vec_inliers.end(); it++)
                good_match.push_back(matchings[*it]);
            matchings = good_match;
        }
    } else {                                   // ...with all matches
        orsa::FundamentalModel model(matchings);
        std::vector<orsa::ModelEstimator::Model> models;
        std::vector<int> indices;
        int i=0;
        for(std::vector<Match>::const_iterator it=matchings.begin();
            it!=matchings.end(); ++it)
            indices.push_back(i++);
        model.Fit(indices, &models);
        if((ok = (models.size()==1)))
            F = models.front();
    }
    if(! ok) {
        std::cerr << "Failed to estimate the fundamental matrix" << std::endl;
        return 1;
    }
    std::cout << "F=" << F <<std::endl;

    // Save inliers
    if(! fileMatchesOut.empty() &&
       ! Match::saveMatch(fileMatchesOut.c_str(), matchings)) {
        std::cerr << "Failed saving matches into " <<fileMatchesOut <<std::endl;
        return 1;
    }

    libNumerics::vector<double> eL(3), eR(3);
    orientedEpipoles(matchings, F, eL, eR); // Find epipole in each image
    Polarectifyer pol(F, eL, eR, w1, h1, w2, h2); // Compute rectification maps

    const char* fileL = argv[3];
    const char* fileR = argv[4];
    Image<RGBColor>* out1 = sample(image1, pol.widthL(), pol.height(),
                                   pol.pullback_map(true));
    Image<RGBColor>* out2 = sample(image2, pol.widthR(), pol.height(),
                                   pol.pullback_map(false));
    if(annotate) {
        draw_guidelines(*out1, CYAN);
        draw_sift(*out1, matchings, true, pol, RED);
        draw_guidelines(*out2, CYAN);
        draw_sift(*out2, matchings, false, pol, RED);
    }
    libs::WriteImage(fileL, *out1);
    libs::WriteImage(fileR, *out2);
    delete out1;
    delete out2;

    return 0;
}
