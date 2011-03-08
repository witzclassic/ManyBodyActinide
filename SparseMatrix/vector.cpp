/*
 *  vector.cpp
 *  SparseMatrix
 *
 *  Created by Brad Marston on 4/1/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#import <libc.h>

#import "vector.h"


const double reciprocolLarge = 1.0/((double)(pow(2, 31) - 1));
double ran(void) // random number between 0 and 1
{
    return reciprocolLarge * random();
}

void random(vector *a)
{
    for (long i = 0; i < a->d; i++) {
        a->component[i] = ran() - 0.5;
    }
    
    normalize(a);
}

void multiplyAdd(const vector *a, const vector *b, const complex c, vector *d)
{
    d->eigenvalue = a->eigenvalue;
    for (long i = 0; i < a->d; i++) {
        d->component[i] = a->component[i] + c * b->component[i];
    }
}
