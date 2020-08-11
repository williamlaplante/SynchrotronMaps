#include <healpix_cxx/pointing.h>
#include "../include/corrfunc.h"
#include "../include/helper.h"
#include <cassert>
#include <cmath>

bool isUnseen(double pixel);

std::tuple<double,double> compute_corr(Healpix_Map<double> & map1, Healpix_Map<double> & map2, Angle R, Angle dr) 
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
    
    assert(resol.value < dr.value && "Input thickness smaller than resolution");

    if (R.value==0) {
        arr<double> res_arr = pw_dot(map1.Map(), map2.Map());
        result = mean(res_arr);
        error = stdev(res_arr)/std::sqrt(Npix);
        return std::make_tuple(result, error);
    }
    assert(resol.value < R.value && "Radius smaller than resolution.");

    int N;
    double sum;
    int numUnseen;
    arr<int> result_len(Npix);
    arr<double> result_arr(Npix);
    pointing ptg;
    rangeset<int> pix_range;
    
    //Iterate over pixels
    for (int i=0; i < Npix; i++)
    {
        if (isUnseen(map1[i])){
            continue;
        }
        ptg = map1.pix2ang(i); //Get theta,phi of ith pixel
        pix_range = map2.query_disc(ptg, R.value+dr.value).op_xor(map2.query_disc(ptg,R.value)); //Get the rangeset of pixels on ith annulus
        N = pix_range.toVector().size();
        assert(N>0 && "queried annulus has null size.");
        sum=0;
        numUnseen=0;
        for (int pix : pix_range.toVector()){ //sum over values at non-unseen pixels in the ith annulus
            if (isUnseen(map2[i])) {
                numUnseen++;
                }
            else {
                sum += map2[pix];
            }
            } 
        result_arr[i] = map1[i] * sum; // store ith correlation value
        result_len[i] = N - numUnseen;
    }

    double len = mean(result_len);
    result = mean(result_arr)/len; // get the average correlation value
    error = stdev(result_arr)/((std::sqrt(Npix))*len); //Get error in averaging

    return std::make_tuple(result, error);
    
}

const double UNSEEN = -1.6375 * std::pow(10,30);

bool isUnseen(double pixel){
    return std::abs(pixel - UNSEEN) < 10.0;
}