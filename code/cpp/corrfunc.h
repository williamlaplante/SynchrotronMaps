#ifndef CORRFUNC_H
#define CORRFUNC_H
#include <tuple>
#include "helper.h"
#include <healpix_cxx/healpix_map.h>
/*
Two-Point Correlation Function
--------------------------------
Returns a tuple (result, error) of the correlation function. Error is spatial.

Parameters : 
    Healpix_Map<double> & map1 : first healpix map.
    Healpix_Map<double> & map2 : second healpix map.
    const double R : angular distance input of correlation function.
    const double dr : thickness of annuli. 


*/
std::tuple<double,double> compute_corr(Healpix_Map<double> & map1, Healpix_Map<double> & map2, Angle & R, Angle & dr);

#endif