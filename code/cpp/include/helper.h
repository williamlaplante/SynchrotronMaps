#ifndef HELPER_H
#define HELPER_H

#include <healpix_cxx/healpix_map.h>

enum angular_units {DEGREES, RADIANS, ARCMIN};

class Angle {

    public:
    Angle();
    Angle(double, angular_units);
    double value;
    angular_units unit;
    void deg_to_rad();
    void rad_to_deg();
    void arcmin_to_deg();
    void deg_to_arcmin();

};

arr<double> pw_dot(const arr<double> & arr1, const arr<double> & arr2);
double mean(const arr<double> & array);
double stdev(const arr<double> & array);
double mean(const arr<int> & array);
double stdev(const arr<int> & array);
#endif