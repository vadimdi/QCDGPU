/******************************************************************************
 * @file     algebra_su2.cpp
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.6
 *
 * @brief    [QCDGPU]
 *           Defines required matrix operations for the SU(2) group
 *
 * @section  LICENSE
 *
 * Copyright (c) 2013-2016 Vadim Demchik, Natalia Kolomoyets
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 *    Redistributions of source code must retain the above copyright notice,
 *      this list of conditions and the following disclaimer.
 *
 *    Redistributions in binary form must reproduce the above copyright notice,
 *      this list of conditions and the following disclaimer in the documentation
 *      and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 *****************************************************************************/

#include "algebra_su2.h"

void        lattice_unity(su_2 *a){
        (*a).u1.re = 1.0;  (*a).u1.im = 0.0;
        (*a).u2.re = 0.0;  (*a).u2.im = 0.0;
        (*a).v1.re = 0.0;  (*a).v1.im = 0.0;
        (*a).v2.re = 1.0;  (*a).v2.im = 0.0;
}

void        lattice_zero(su_2 *a){
        (*a).u1.re = 0.0;  (*a).u1.im = 0.0;
        (*a).u2.re = 0.0;  (*a).u2.im = 0.0;
        (*a).v1.re = 0.0;  (*a).v1.im = 0.0;
        (*a).v2.re = 0.0;  (*a).v2.im = 0.0;
}

void lattice_matrixGID(su_2 *a, int gid, int dir, int sites){
    (*a).u1.re = sin(((dir + 1) / 1000.0 + 10.0 / sites) * gid) / 2.0;
    (*a).u2.im = cos(((dir + 1) / 1000.0 - 10.0 / sites) * gid) / 3.0;
    (*a).u2.re = sin(((dir + 1) / 1000.0 - 18.0 / sites) * gid) / 4.0;
    (*a).u1.im = sqrt(1 - (*a).u1.re * (*a).u1.re - (*a).u2.im * (*a).u2.im - (*a).u2.re * (*a).u2.re);
    
    (*a).v1.re = - (*a).u2.re;
    (*a).v1.im = (*a).u2.im;
    (*a).v2.re = (*a).u1.re;
    (*a).v2.im = - (*a).u1.im;
}

su_2 operator + (su_2 a, su_2 b){
    su_2 c;
    
    c.u1.re = a.u1.re + b.u1.re;    c.u1.im = a.u1.im + b.u1.im;
    c.u2.re = a.u2.re + b.u2.re;    c.u2.im = a.u2.im + b.u2.im;
    
    c.v1.re = a.v1.re + b.v1.re;    c.v1.im = a.v1.im + b.v1.im;
    c.v2.re = a.v2.re + b.v2.re;    c.v2.im = a.v2.im + b.v2.im;
    
    return c;
}

su_2 operator - (su_2 a, su_2 b){
    su_2 c;
    
    c.u1.re = a.u1.re - b.u1.re;    c.u1.im = a.u1.im - b.u1.im;
    c.u2.re = a.u2.re - b.u2.re;    c.u2.im = a.u2.im - b.u2.im;
    
    c.v1.re = a.v1.re - b.v1.re;    c.v1.im = a.v1.im - b.v1.im;
    c.v2.re = a.v2.re - b.v2.re;    c.v2.im = a.v2.im - b.v2.im;
    
    return c;
}

su_2 operator * (su_2 a, su_2 b){
    su_2 c;
    
        c.u1.re = -a.u1.im * b.u1.im + a.u1.re * b.u1.re - a.u2.im * b.v1.im + a.u2.re * b.v1.re;
        c.u1.im =  a.u1.re * b.u1.im + a.u1.im * b.u1.re + a.u2.re * b.v1.im + a.u2.im * b.v1.re;
        c.u2.re = -a.u1.im * b.u2.im + a.u1.re * b.u2.re - a.u2.im * b.v2.im + a.u2.re * b.v2.re;
        c.u2.im =  a.u1.re * b.u2.im + a.u1.im * b.u2.re + a.u2.re * b.v2.im + a.u2.im * b.v2.re;
        c.v1.re = -a.v1.im * b.u1.im + a.v1.re * b.u1.re - a.v2.im * b.v1.im + a.v2.re * b.v1.re;
        c.v1.im =  a.v1.re * b.u1.im + a.v1.im * b.u1.re + a.v2.re * b.v1.im + a.v2.im * b.v1.re;
        c.v2.re = -a.v1.im * b.u2.im + a.v1.re * b.u2.re - a.v2.im * b.v2.im + a.v2.re * b.v2.re;
        c.v2.im =  a.v1.re * b.u2.im + a.v1.im * b.u2.re + a.v2.re * b.v2.im + a.v2.im * b.v2.re;
    return c;
}

su_2 operator += (su_2 a, su_2 b){
    return a + b;
}

su_2 Herm(su_2 a){
    su_2 result;

        result.u1.re =  a.u1.re;
        result.u1.im = -a.u1.im;
        result.u2.re =  a.v1.re;
        result.u2.im = -a.v1.im;

        result.v1.re =  a.u2.re;
        result.v1.im = -a.u2.im;
        result.v2.re =  a.v2.re;
        result.v2.im = -a.v2.im;

    return result;
}

hgpu_complex Tr(su_2 a){
    hgpu_complex tr_a;
    
    tr_a.re = a.u1.re + a.v2.re;
    tr_a.im = a.u1.im + a.v2.im;
    
    return tr_a;
}

hgpu_double ReTr(su_2 a){
    return (hgpu_double)(a.u1.re + a.v2.re);
}

void matrix_reconstruct(su_2 *a){
    (*a).v1.re = - (*a).u2.re;
    (*a).v1.im = (*a).u2.im;
    (*a).v2.re = (*a).u1.re;
    (*a).v2.im = - (*a).u1.im;
}
