#!/usr/bin/env python

# CG.py - Clebsch Gordan Coefficients
#
# Imports: math
# Requires: Adapted from matlab code ClebschGordan.m by David Terr, 1994
#
# Adapted by: Steve Horowitz, 2009, Brown University Physics (Prof. J.B. Marston)
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# David Terr owns the intellectual property right for this file and reserves the
# right to distribute it under a license other than LGPL

#
# <======================================================================>
# - factorial
# - binomialCoeff: binomial coefficient
# - Winger3j
# - ClebschGordan
# <======================================================================>
# <======================================================================>


import math
# 
def factorial (n,limit=-1):

#    print n, limit
    if n < 0:
        return -1
    if n == 0 or n == 1 or limit == 0:
        return 1
    return n*factorial(n-1,limit-1)

# See wikipedia page for binomial coefficient definition
def binomialCoeff(n,k):  

    return factorial(n,k)/factorial(k)


# Wigner3j.m by David Terr, Raytheon, 6-17-04

# Compute the Wigner 3j symbol using the Racah formula [1]. 

def Wigner3j(j1,j2,j3,m1,m2,m3):

# error checking
    if ( 2*j1 != math.floor(2*j1) or 2*j2 != math.floor(2*j2) or 2*j3 != math.floor(2*j3) 
         or 2*m1 != math.floor(2*m1) or 2*m2 != math.floor(2*m2) or 2*m3 != math.floor(2*m3) ):
#    error('All arguments must be integers or half-integers.');
        return 0 

    if ( j1 - m1 != math.floor ( j1 - m1 ) ):
#    error('2*j1 and 2*m1 must have the same parity');
        return 0

    if ( j2 - m2 != math.floor ( j2 - m2 ) ):
    # error('2*j2 and 2*m2 must have the same parity')
        return 0


    if ( j3 - m3 != math.floor ( j3 - m3 ) ) :
    # error('2*j3 and 2*m3 must have the same parity')
        return 0


    if j3 > j1 + j2 or j3 < math.fabs(j1 - j2) :
    # error('j3 is out of bounds.')
        return 0 


    if math.fabs(m1) > j1:
    # error('m1 is out of bounds.')
        return 0

    if math.fabs(m2) > j2:
        # error('m2 is out of bounds.')
        return 0

    if math.fabs(m3) > j3 :
    # error('m3 is out of bounds.')
        return 0


    t1 = j2 - m1 - j3
    t2 = j1 + m2 - j3
    t3 = j1 + j2 - j3
    t4 = j1 - m1
    t5 = j2 + m2

    tmin = int(max( 0, max( t1, t2 ) ))
    tmax = int(min( t3, min( t4, t5 ) ))

    wigner = 0

    for t in range(tmin,tmax+1):
        wigner = (wigner + math.pow(-1,t) / ( factorial(t)  
                                      * factorial(t-t1)  
                                      * factorial(t-t2) 
                                      * factorial(t3-t) 
                                      * factorial(t4-t) 
                                      * factorial(t5-t) ))   


    return (wigner * math.pow(-1,j1-j2-m3) 
            * math.sqrt(factorial(j1+j2-j3)  
                        * factorial(j1-j2+j3) 
                        * factorial(-j1+j2+j3) 
                        / factorial(j1+j2+j3+1) 
                        * factorial(j1+m1) 
                        * factorial(j1-m1) 
                        * factorial(j2+m2) 
                        * factorial(j2-m2) 
                        * factorial(j3+m3) 
                        * factorial(j3-m3) ) )


# Reference: Wigner 3j-Symbol entry of Eric Weinstein's Mathworld: http://mathworld.wolfram.com/Wigner3j-Symbol.html


def ClebschGordan(j1, j2, j, m1, m2, m):

    cg = 0

    # warning checking
    if ( ( 2*j1 != math.floor(2*j1) 
           or 2*j2 != math.floor(2*j2) 
           or 2*j != math.floor(2*j) 
           or 2*m1 != math.floor(2*m1) 
           or 2*m2 != math.floor(2*m2) 
           or 2*m != math.floor(2*m) )):
#        warning('All arguments must be integers or half-integers.')
        return 0


    if m1 + m2 != m:
#    warning('m1 + m2 must equal m.')
        return 0


    if ( j1 - m1 != math.floor ( j1 - m1 ) ):
#    warning('2*j1 and 2*m1 must have the same parity')
        return 0

    if ( j2 - m2 != math.floor ( j2 - m2 ) ):
#    warning('2*j2 and 2*m2 must have the same parity')
        return 0


    if ( j - m != math.floor ( j - m ) ):
#    warning('2*j and 2*m must have the same parity')
        return 0

    if j > j1 + j2 or j < math.fabs(j1 - j2):
#    warning('j is out of bounds.')
        return 0 

    if math.fabs(m1) > j1:
#        warning('m1 is out of bounds.')
        return 0

    if math.fabs(m2) > j2:
#    warning('m2 is out of bounds.')
        return 0 

    if math.fabs(m) > j:
#  warning('m(#d) > j(#d) is out of bounds.', m, j)
        return 0 

    return math.pow(-1,j1-j2+m) * math.sqrt(2*j + 1) * Wigner3j(j1,j2,j,m1,m2,-m)



