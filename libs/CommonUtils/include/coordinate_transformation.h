#ifndef COORDTRANSFORM_H
#define COORDTRANSFORM_H

#include <Matrix/matrix_inv_class.h>
#include<cmath>


MatrixInv<float> Geodetic2Ecef(float lat, float lon, float height);
MatrixInv<float> Geodetic2Ned(float lat, float lon, float height, float lat_ref, float lon_ref, float height_ref);

#endif