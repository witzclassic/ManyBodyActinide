/*
 *  vector.h
 *  SparseMatrix
 *
 *  Created by Brad Marston on 4/1/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#import <iostream>
#import <complex>
typedef std::complex<double> complex;

class sparseMatrix;

class vector {
    
    friend void eigenvalue(vector *a, double eigenvalue);
    friend double eigenvalue(vector *a);
    friend void normalize(vector *a);
    friend void random(vector *a);
    friend void multiplyAdd(const vector *a, const vector *b, const complex c, vector *d);
    friend complex operator*(const vector &a, const vector &b);
    friend std::ostream &operator<<(std::ostream &os, vector *v);
    friend std::istream &operator>>(std::istream &is, vector *v);

    friend void sparseMultiply(const sparseMatrix *h, const vector *v, vector *Av);
    friend void preconditioner(const sparseMatrix *h, double lambda, const vector *q, vector *v);

public:
    static long d;    
    vector(); // the default constructor
    ~vector();
    void zero();
    vector &operator=(const vector &v);
    
private:
    double eigenvalue;
    complex *component;

};

inline vector::vector()
{
    eigenvalue = 0.0;
    component = new complex[d];
    memset(component, 0, d*sizeof(complex));
}

inline vector::~vector()
{
    if (component) delete[] component;
}   

inline void vector::zero()
{
    eigenvalue = 0.0;
    memset(component, 0, d*sizeof(complex));
}

inline vector &vector::operator=(const vector &v)
{
    eigenvalue = (&v)->eigenvalue;
    memcpy(component, (&v)->component, d*sizeof(complex));
    return *this;
}

inline void eigenvalue(vector *a, double eigenvalue)
{
    a->eigenvalue = eigenvalue;
}

inline double eigenvalue(vector *a)
{
    return a->eigenvalue;
}

inline void normalize(vector *a)
{
    complex inner = *a * *a;
    double norm = pow(real(inner), -0.5);
    long i;
    for (i = 0; i < a->d; i++) {
        a->component[i] = norm * a->component[i];
    }
}

inline complex operator*(const vector &a, const vector &b)
{
    long i;
    complex inner;
    inner = 0.0;
    for (i = 0; i < a.d; i++) {
        inner += conj(a.component[i]) * b.component[i];
    }
    
    return inner;
}

inline std::ostream &operator<<(std::ostream &os, vector *v)
{
    os << v->d << "\t\t" << v->eigenvalue << "\n";
    for (long i = 0; i < v->d; i++) {
        os << real(v->component[i]) << "\t" << imag(v->component[i]) << "\n";
    }
    
    return os;
}

inline std::istream &operator>>(std::istream &is, vector *v)
{
    while (is.peek() == '%') {
        is.ignore(1024, '\n');
    }    
    
    is >> v->d >> v->eigenvalue;

    double real, imag;
    for (long i = 0; i < v->d; i++) {
        if (is >> real >> imag) {
            v->component[i] = complex(real, imag);
        } else {
            std::cerr << "Mismatch in number of vector components";
            exit(EXIT_FAILURE);
        }
    }
        
    return is;
}

