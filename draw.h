#ifndef DRAW_H
#define DRAW_H

#include "libImageIO/image.hpp"
#include "libImageIO/pixelTypes.hpp"
#include <vector>

class Polarectifyer;
class Match;

void draw_horizontal_dashed_line(Image<RGBColor>& im, int y, RGBColor c, 
                                 int length=7, int gap=3);
void draw_cross(Image<RGBColor>& im, int x,int y, RGBColor c, int halfLength=3);

void draw_guidelines(Image<RGBColor>& im, RGBColor c, int num=15,
                     int length=7, int gap=3);
void draw_sift(Image<RGBColor>& im, const std::vector<Match>& m, bool bLeft,
               const Polarectifyer& pol, RGBColor c, int halfLength=3);

#endif
