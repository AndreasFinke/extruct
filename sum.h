#pragma once 

#include "xsum.h"

class Sum {
public:
    Sum() {
        xsum_large_init(&state);
    }

    void add(double * vec, int size) {
        xsum_large_addv(&state,vec,size);
    }

    operator double() {
        return xsum_large_round(&state);
    }

    /* for performance reasons this returns 0! */
    double operator+=(double t) {
        xsum_large_addv(&state,&t,1);
        return 0;
    }

private:

    xsum_large_accumulator state;

};
