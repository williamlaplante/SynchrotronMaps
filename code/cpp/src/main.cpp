#include "../include/corrfunc.h"
#include "../include/helper.h"
#include <healpix_cxx/healpix_map_fitsio.h>
#include <string>
#include <ctime>
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

    time_t start, end;
    Healpix_Ordering_Scheme scheme = RING;
    const nside_dummy dummy;
    Healpix_Map<double> map(128, scheme, dummy);
    map.fill(1);

    Angle dr(0.5, DEGREES);

    time(&start);
    for (double i=1; i<=5; i+=1){
        Angle R(i, DEGREES);
        print_tuple(compute_corr(map,map,R,dr));
    }
    time(&end);
    std::cout << "program ended in " << (double)(end-start) << " second(s)." << std::endl;
    return 0;
}
