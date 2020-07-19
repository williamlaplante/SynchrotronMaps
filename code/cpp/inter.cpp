#include "inter.h"


int * intersection(int *p1, int *p2) {

    static int r[10];
    for (int i=0; i<10; i++){
        r[i] = i;
    }

    return r;
}

