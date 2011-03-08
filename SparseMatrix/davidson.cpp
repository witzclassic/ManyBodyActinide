		/* DAVIDSON.C: low lying eigenvectors and values */
	      /* of complex, sparse, Hermitian Hamiltonian */
	 	   /* using the Davidson-Liu algorithm */

/*
	see:  Ernest R. Davidson, J. Comp. Phys. 17, 87 -- 94 (1975) 	
	      and Comp. in Phys. 7, 519 -- 522 (1993);
	      Ronald B. Morgan, J. Comp. Phys. 101, 287 -- 291 (1992);
	      and Albert Booten and Henk van der Vorst, Comp. in Phys. 10,
	      239 -- 242 and 331 -- 334 (1996).
*/	
#import <complex>
typedef std::complex<double> complex;

#import "davidson.h"

bool davidson(const sparseMatrix *h, vector v[DIM], vector Av[DIM], vector *q, const long k, const vector eigenvector[])  
{	
    int m = 1; /* initial dimension of subspace */
    
    /* orthonormalize starting vector w.r.t. eigenvectors already found */
    orthonormalize(&v[0], k, eigenvector);
    
    /* (A) compute matrix elements A in reduced basis and diagonalize */
    sparseMultiply(h, &v[0], &Av[0]);
        
    /* initial eigenvalue and eigenvector */
    double lambda[DIM];
    memset(lambda, 0, DIM*sizeof(double));
    complex alpha[DIM*DIM];
    memset(alpha, 0, DIM*DIM*sizeof(complex));
    
    complex vAv[DIM][DIM]; 
    vAv[0][0] = v[0] * Av[0];
    lambda[0] = real(vAv[0][0]);
    alpha[0] = 1.0;
    
    std::cout << "\n\nsize \t\teigenvalue \t\terror\n";
    std::cout.flush();
    
    /* begin iterating: initial values of control parameters */
    double qnorm = 1.0; /* large initial value */
    
    /* iteration loop */
    while ((qnorm > DISCREPANCY) && (m < DIM)) {
        
        /* (B) form q */
        (*q).zero(); 
        for (long i = 0; i < m; i++) {
            multiplyAdd(q, &Av[i], alpha[i], q);
            multiplyAdd(q, &v[i], -alpha[i] * lambda[0], q);
        }
        
        /* (C) calculate and report error */
        qnorm = sqrt(real((*q) * (*q)));
        std::cout << m << "\t\t" << lambda[0] << "\t\t" << qnorm << "\n\n";
        std::cout.flush();

        /* (D) form v[m] using preconditioner */
        preconditioner(h, lambda[0], q, &v[m]);
        
        /* (E1) orthogonalize v[m] w.r.t. existing v's */
        orthonormalize(&v[m], m, v);
        
        /* (E2) orthonormalize v[m] w.r.t. eigenvectors already found */
        orthonormalize(&v[m], k, eigenvector);
            
        /* (F) add new row and column to matrix A */  
        sparseMultiply(h, &v[m], &Av[m]);
        for (long i = 0; i < m + 1; i++) {
            vAv[i][m] = v[i] * Av[m];
        }
        
        m = m + 1; /* increase size of subspace */
        
        /* (G) diagonalize A and return to step (B) */
        complex a[DIM*DIM];
        for (long i = 0; i < m; i++) { 
            for (long j = i; j < m; j++) {
                a[i+(j+1)*j/2] = vAv[i][j];
            }
        }
        
        char jobz, uplow;
        jobz = 'v';
        uplow = 'u';
        complex work[2*DIM];
        double rwork[3*DIM];
        int info;        
        int status = zhpev_(&jobz, &uplow, (__CLPK_integer *)(&m), (__CLPK_doublecomplex *)a, lambda, 
                            (__CLPK_doublecomplex *)alpha, (__CLPK_integer *)(&m), (__CLPK_doublecomplex *)work, rwork, (__CLPK_integer *)(&info));

        if (info != 0) {
            std::cerr << "zhpev failed.  status = " << status << "\n";
            exit(EXIT_FAILURE);
        }        
    }
    
    /* update v[0] */
    (*q).zero();
    for (long i = 0; i < m; i++) {
        multiplyAdd(q, &v[i], alpha[i], q);
    }
    normalize(q);
    
    v[0] = *q;
    eigenvalue(&v[0], lambda[0]);
        
    /* finished? */
    if (m == DIM) return false;
    else return true;
}

void orthonormalize(vector *x, const long k, const vector v[])
{
    complex c;
    for (long i = 0; i < k; i++) {
        c = v[i] * (*x);
        multiplyAdd(x, &v[i], -c, x);
    }
    normalize(x);  
}

void sparseMultiply(const sparseMatrix *h, const vector *v, vector *Av)
{
    (*Av).zero();
    
    if ((v->d != h->rows) || (Av->d != h->columns) || (v->d != Av->d)) {
        std::cerr << "Mismatch in rows or columns in matrix-vector multiply\n";
        exit(EXIT_FAILURE);
    }
    
    int r, c;
    double hr;
    for (long i = 0; i < h->realElements; i++) { // off-diagonal elements
        
        r = h->realOffdiagonalMatrixElement[i].rowIndex;
        c = h->realOffdiagonalMatrixElement[i].columnIndex;
        
        hr = h->realOffdiagonalMatrixElement[i].value;
        
        Av->component[c] += hr * v->component[r];
        Av->component[r] += hr * v->component[c];
    }
    
    complex hu, hl;
    for (long i = 0; i < h->complexElements; i++) { // off-diagonal elements
        
        r = h->complexOffdiagonalMatrixElement[i].rowIndex;
        c = h->complexOffdiagonalMatrixElement[i].columnIndex;
        
        hu = complex(h->complexOffdiagonalMatrixElement[i].value.r, h->complexOffdiagonalMatrixElement[i].value.i);
        hl = conj(hu);
        
        Av->component[c] += hu * v->component[r];
        Av->component[r] += hl * v->component[c];
    }
    
    for (long i = 0; i < h->rows; i++) { // diagonal elements
        Av->component[i] += double(h->diagonalMatrixElement[i]) * v->component[i];
    }
    Av->eigenvalue = v->eigenvalue;
}

void preconditioner(const sparseMatrix *h, double lambda, const vector *q, vector *v)
{
    double precond;
    for (long i = 0; i < q->d; i++) {
        precond = fabs(h->diagonalMatrixElement[i] - lambda) + EPSILON;
        v->component[i] = q->component[i] / precond;
    }
    v->eigenvalue = lambda;
}   

