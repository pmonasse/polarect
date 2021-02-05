Polarect - Polar rectification of a stereo image pair

Pascal Monasse pascal.monasse@enpc.fr IMAGINE/LIGM, Universite Paris Est

The program relies on the following fine packages (IPOL publications):
  - Fundamental matrix (https://doi.org/10.5201/ipol.2016.147)
  - B-spline interpolation (https://doi.org/10.5201/ipol.2018.221)
  - SIFT Anatomy (https://doi.org/10.5201/ipol.2014.82)

# Future releases and updates:
https://github.com/pmonasse/polarect.git

# Description:
Perform epipolar rectification of a stereo pair based on polar transform around
each epipole. The advantage over regular homographies is that it can handle
an epipole located inside the image.

# Licensing: See LICENSE.txt file

# Build
The build procedure is handled by CMake (https://cmake.org/). Under Linux, the
following commands are enough:
$ mkdir build
$ cd build
$ cmake -DCMAKE_BUILD_TYPE=Release /path/to/source
$ make

# Usage
./polarect [options] imgInA imgInB imgOutA imgOutB
- imgInA, imgInB: two input images (JPG or PNG format)
- imgOutA, imgOutB: output rectified images
        Options:
-s, --sift=ARG SIFT distance ratio of descriptors (0.6)
-r, --read=ARG Read file of matches, do not use SIFT
-w, --write=ARG Write file of inlier matches from SIFT algorithm
-i, --inliers Find inlier matches from file of matches (with option -r)
-a, --annotate Draw annotations (guidelines and SIFT matches)
Options -r, -w and -i are not relevant for the demo, only used to avoid SIFT
extraction. Option -a draws guidelines and SIFT points to inspect the result.

# Example
Try with the images of folder data:
$ ./polarect -a ../data/corridor0.jpg ../data/corridor1.jpg out0.jpg out1.jpg
Compare visually the output files out0.jpg and out1.jpg with
../data/corridor0_rect.jpg and ../data/corridor1_rect.jpg.

# Reviewed files in IPOL:
  - polarect.{h,cpp}
  - main.cpp
