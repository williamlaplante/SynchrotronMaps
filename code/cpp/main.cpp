#include "corrfunc.h"
#include "helper.h"
#include <healpix_cxx/healpix_map_fitsio.h>
#include <string>
//to compile, run g++ *.cpp -lhealpix_cxx -o myprogram

/*
To create a map with Nside resolution and ordering scheme: 
----------------------
int nside = 128;
Heapix_Ordering_Scheme scheme = RING;
const nside_dummy dummy;
Healpix_Map<double> map(nside, scheme, dummy);
-------------------
Now one can use map.fill(value) or just iterate over map and assign values.
*/

void print_tuple(std::tuple<double,double> out){
        std::cout << "Result is : " << std::get<0>(out) << " +/- " << std::get<1>(out) << std::endl;
}
int main() {
    Healpix_Ordering_Scheme scheme = RING;
    const nside_dummy dummy;
    Healpix_Map<double> map(128, scheme, dummy);
    map.fill(10);
    
    Angle R(1, DEGREES);
    Angle dr(0.1, DEGREES);
    
    print_tuple(compute_corr(map,map,R,dr));
    
    return 0;
}