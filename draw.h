/**
 * SPDX-License-Identifier: GPL-2.0+
 * @file draw.h
 * @brief Draw annotations on image
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
