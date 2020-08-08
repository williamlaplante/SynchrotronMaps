#include "helper.h"
#ifndef PI
#define PI 3.1415926535897932
#endif

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
    if (arr1.size()!=arr2.size()){
        throw "Size mismatch. Cannot execute pointwise multiplication.";
    }
    int size = arr1.size();
    arr<double> result(size);

    for (int i=0; i<arr1.size(); i++){
        result[i]=arr1[i]*arr2[i];
    }

    return result;
}

double sum(const arr<double> & array){
    double sum=0;
    for (int i=0; i<array.size(); i++){
        sum+=array[i];
    }
    return sum;
}

int sum(const arr<int> & array){
    int sum=0;
    for (int i=0; i<array.size(); i++){
        sum+=array[i];
    }
    return sum;
}

double mean(const arr<int> & array)
{
    int size = array.size();
    if (size==0){throw "Array is of null size. Cannot take the mean of array.";}
    int res = sum(array);
    return res/size;
}
double mean(const arr<double> & array)
{
    int size = array.size();
    if (size==0){throw "Array is of null size. Cannot take the mean of array.";}
    int res = sum(array);
    return res/size;
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

