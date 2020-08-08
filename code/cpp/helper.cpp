#include "helper.h"
#include <cassert>
#ifndef PI
#define PI 3.1415926535897932
#endif

Angle::Angle(double val, angular_units u){
    value = val;
    unit = u;
}

Angle::Angle(){
    value = 0;
    unit = DEGREES;
}

void Angle::rad_to_deg() {
    value = (value*180)/PI;
    unit = DEGREES;
}

void Angle::deg_to_rad() {
    value = (value*PI)/180;
    unit = RADIANS;
}

void Angle::deg_to_arcmin() {
    value = value/60;
    unit = ARCMIN;
}

void Angle::arcmin_to_deg() {
    value = value*60;
    unit = DEGREES;
}

arr<double> pw_dot(const arr<double> & arr1, const arr<double> & arr2)
{
    assert(arr1.size()!=arr2.size());
    int size = arr1.size();
    arr<double> result(size);

    for (int i=0; i<arr1.size(); i++){
        result[i]=arr1[i]*arr2[i];
    }

    return result;
}

double mean(const arr<int> & array)
{
    int size = array.size();
    assert(size!=0);
    double sum=0;
    
    for (int i=0; i<size; i++){
        sum+=array[i];
    }
    return sum/size;
}
double mean(const arr<double> & array)
{
    int size = array.size();
    assert(size!=0);
    double sum=0;
    for (int i=0; i<size; i++){
        sum+=array[i];
    }
    return sum/size;
}

double stdev(const arr<int> & array){

    double avg = mean(array);
    double sum = 0;
    int N = array.size();
    if (N<=1){
        return 0;
    }
    for (int i =0; i<N; i++){
        sum+= std::pow(((double)array[i] - avg),2);
    }
    return std::sqrt(sum/(N-1));
}

double stdev(const arr<double> & array){

    double avg = mean(array);
    double sum = 0;
    int N = array.size();
    if (N<=1){
        return 0;
    }
    for (int i =0; i<N; i++){
        sum+= std::pow((array[i] - avg),2);
    }
    return std::sqrt(sum/(N-1));
}

