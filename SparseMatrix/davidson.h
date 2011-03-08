/*
 *  davidson.h
 *  SparseMatrix
 *
 *  Created by Brad Marston on 4/14/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#import <iostream>
#import <Accelerate/accelerate.h>

#import "vector.h"
#import "sparseMatrix.h"

static const long DIM = 10; // size of variational subspace

static const double DISCREPANCY = 1.0e-8; /* maximum error */
static const double EPSILON = 1.0e-4; /* cutoff in preconditioner */

bool davidson(const sparseMatrix *h, vector v[DIM], vector Av[DIM], vector *q, const long k, const vector eigenvector[]);
void orthonormalize(vector *x, const long k, const vector v[]);
