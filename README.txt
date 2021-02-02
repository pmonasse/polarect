Polarect - Polar rectification of a stereo image pair

Pascal Monasse pascal.monasse@enpc.fr IMAGINE/LIGM, Universite Paris Est

The program relies on the following fine packages (IPOL publications):
  - Fundamental matrix (https://doi.org/10.5201/ipol.2016.147)
  - B-spline interpolation (https://doi.org/10.5201/ipol.2018.221)
  - SIFT Anatomy (https://doi.org/10.5201/ipol.2014.82)

Future releases and updates:
https://github.com/pmonasse/polarect.git

Description:
Perform epipolar rectification of a stereo pair based on polar transform around
each epipole. The advantage over regular homographies is that it can handle
an epipole located inside the image.

Licensing: See LICENSE.txt file

Reviewed files in IPOL:
  - polarect.{h,cpp}
  - main.cpp
