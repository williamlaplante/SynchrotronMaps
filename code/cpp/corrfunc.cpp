#include <healpix_cxx/pointing.h>
#include "corrfunc.h"
#include "helper.h"
#include <cassert>

std::tuple<double,double> compute_corr(Healpix_Map<double> & map1, Healpix_Map<double> & map2, Angle & R, Angle & dr) 
{
    
    double result,error; //Initialize the return values
    Healpix_Ordering_Scheme scheme = RING;
    assert(map1.Nside() == map2.Nside() && "mismatch in map resolution.");//Check if matching resolution
    assert((map1.Order()!=scheme || map2.Order()!=scheme) && "Ordering must be RING."); // Check if correct scheme
    assert((R.unit==DEGREES && dr.unit==DEGREES) && "units of R and dr must be DEGREES.");
    
    R.deg_to_rad();
    dr.deg_to_rad();

    const int Nside = map1.Nside();
    const int Npix = map1.Npix();
    Angle resol;
    resol.value = map1.max_pixrad(); 
    resol.unit = RADIANS;
    
    assert(dr.value < resol.value && "Input thickness smaller than resolution");

    if (R.value==0) {
        arr<double> res_arr = pw_dot(map1.Map(), map2.Map());
        result = mean(res_arr);
        error = stdev(res_arr)/std::sqrt(Npix);
        return std::make_tuple(result, error);
    }
    assert(resol.value < R.value && "Radius smaller than resolution.");


    arr<int> ring_len(Npix);
    arr<double> result_arr(Npix);
    pointing ptg;
    rangeset<int> pix_range;

    //Iterate over pixels
    for (int i=0; i < Npix; i++)
    {
        ptg = map1.pix2ang(i); //Get theta,phi of ith pixel
        pix_range = map2.query_disc(ptg, R.value+dr.value).op_xor(map2.query_disc(ptg,R.value)); //Get the rangeset of pixels on ith annulus
        
        double sum=0;
        for (auto pix : pix_range.toVector()){sum+=map2[pix];} //sum over values at pixels in the ith annulus

        result_arr[i] = sum * map1[i]; // scale by value of map1 at ith pixel
        ring_len[i] = pix_range.size(); // store rangeset size for further usage. (we summed instead of taking average before, must be fixed afterwards)
    }

    double len = mean(ring_len); // get the average rangeset size.
    result = mean(result_arr)/len; // get the average result, then scale by 1/size
    error = stdev(result_arr)/(std::sqrt(Npix) * len); //Get error

    return std::make_tuple(result, error);
    
}
