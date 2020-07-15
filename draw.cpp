#include "draw.h"
#include "polarect.h"

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
