/**
 * SPDX-License-Identifier: GPL-2.0+
 * @file draw.cpp
 * @brief Draw annotations on image
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

#include "draw.h"
#include "polarect.h"

extern "C" {
#include "libSplinter/splinter.h"
}

inline unsigned char clip_col(double c) {
    return (c<0.)? 0: (c>=255.)? 255: (unsigned char)c;
}

inline RGBColor clip_col(double c[3]) {
    return RGBColor(clip_col(c[0]), clip_col(c[1]), clip_col(c[2]));
}

/// Transform image based on \a pullback.
/// Pixels of the target image, of size wxh, are taken from \a im according to
/// the \a pullback map.
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
                                         3, BOUNDARY_HSYMMETRIC, 1.0e-1, 1);
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

/// Draw a dashed horizontal line in image.
void draw_horizontal_dashed_line(Image<RGBColor>& im, int y, RGBColor c, 
                                 int length, int gap) {
    assert(0<=y && (size_t)y<im.Height());
    size_t w = im.Width();
    RGBColor* p = im.data() + (size_t)y*w;
    for(size_t i=0; i<w; i+=gap)
        for(int j=0; j<length && i<w; j++, i++)
            p[i] = c;
}

/// Draw a cross in image.
void draw_cross(Image<RGBColor>& im, int x,int y, RGBColor c, int halfLength) {
    size_t w=im.Width(), h=im.Height();
    if(0<=y && (size_t)y<h) { // horizontal
        RGBColor* p = im.data()+y*w;
        for(int i=-halfLength; i<=halfLength; i++)
            if(0<=x+i && size_t(x+i)<w)
                p[x+i] = c;
    }
    if(0<=x && (size_t)x<w) { // vertical
        RGBColor* p = im.data()+x;
        for(int i=-halfLength; i<=halfLength; i++)
            if(0<=y+i && size_t(y+i)<h)
                p[size_t(y+i)*w] = c;
    }
}

/// Draw \a num horizontal dashed lines in image \a im.
void draw_guidelines(Image<RGBColor>& im, RGBColor c, int num,
                     int length, int gap){
    assert(num>=0);
    double delta = im.Height()/double(num+1);
    if(delta==0) delta=1;
    for(double i=delta; i<(double)im.Height(); i+=delta)
        draw_horizontal_dashed_line(im, (int)i, c, length, gap);
}

/// Draw SIFT locations as crosses of color \a c in image \a im.
void draw_sift(Image<RGBColor>& im, const std::vector<Match>& m, bool bLeft,
               const Polarectifyer& pol, RGBColor c, int halfLength) {
    for(std::vector<Match>::const_iterator it=m.begin(); it!=m.end(); ++it) {
        std::pair<double,double> xy = std::make_pair(it->x1,it->y1);
        if(! bLeft) xy = std::make_pair(it->x2,it->y2);
        pol.pushforward(xy, bLeft);
        int ix = static_cast<int>(std::floor(xy.first+0.5));
        int iy = static_cast<int>(std::floor(xy.second+0.5));
        draw_cross(im, ix, iy, c, halfLength);
    }
}
