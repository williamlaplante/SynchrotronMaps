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

    time_t t_start, t_end;
    time(&t_start);

    //Setting up parameters
    //-----------------------------------------
    int nside = 2048;
    Healpix_Ordering_Scheme scheme = RING;
    const nside_dummy dummy;
    //-------------------------------------------


    //Reading in maps
    //--------------------------------------------------
    
    Healpix_Map<double> ref_map(nside, scheme, dummy);
    Healpix_Map<double> dust_map(nside, scheme, dummy);

    //----------------------------------------------------


    //Get resolution in degrees
    //--------------------------------------------------
    Angle resol(dust_map.max_pixrad(), RADIANS);
    resol.rad_to_deg();
    std::cout << "resolution of " << resol.value << std::endl;
    //----------------------------------------------------------


    //Setting up correlation function's input : x and dr vectors.
    //--------------------------------------------------------------
    double start = std::log(resol.value);
    double end = std::log(2.6); //Degrees
    int num_points = std::log2(nside);
    double step = (end-start)/(num_points-1);

    std::vector<Angle> x_vector;
    std::vector<Angle> dr_vector;

    for (int i=0; i<num_points; i++){
        x_vector.push_back(Angle(std::exp(start + i*step), DEGREES));
    }

    for (int i=0; i<num_points-1; i++){
        dr_vector.push_back(Angle(x_vector[i+1].value - x_vector[i].value, DEGREES));
    }
    dr_vector.push_back(Angle(x_vector[num_points-1].value - x_vector[num_points-2].value, DEGREES));
//-------------------------------------------------------------------------------



    time(&t_end);
    std::cout << "program ended in " << (double)(t_end-t_start) << " second(s)." << std::endl;
    return 0;
}
