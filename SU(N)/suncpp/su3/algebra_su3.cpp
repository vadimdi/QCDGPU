/******************************************************************************
 * @file     algebra_su3.cpp
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.6
 *
 * @brief    [QCDGPU]
 *           Defines required matrix operations for the SU(3) group
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

#include "algebra_su3.h"

void        lattice_unity(su_3 *a){
        (*a).u1.re = 1.0; (*a).u1.im = 0.0;
        (*a).u2.re = 0.0; (*a).u2.im = 0.0;
        (*a).u3.re = 0.0; (*a).u3.im = 0.0;
        (*a).v1.re = 0.0; (*a).v1.im = 0.0;
        (*a).v2.re = 1.0; (*a).v2.im = 0.0;
        (*a).v3.re = 0.0; (*a).v3.im = 0.0;
        (*a).w1.re = 0.0; (*a).w1.im = 0.0;
        (*a).w2.re = 0.0; (*a).w2.im = 0.0;
        (*a).w3.re = 1.0; (*a).w3.im = 0.0;
}

void        lattice_zero(su_3 *a){
        (*a).u1.re = 0.0; (*a).u1.im = 0.0;
        (*a).u2.re = 0.0; (*a).u2.im = 0.0;
        (*a).u3.re = 0.0; (*a).u3.im = 0.0;
        (*a).v1.re = 0.0; (*a).v1.im = 0.0;
        (*a).v2.re = 0.0; (*a).v2.im = 0.0;
        (*a).v3.re = 0.0; (*a).v3.im = 0.0;
        (*a).w1.re = 0.0; (*a).w1.im = 0.0;
        (*a).w2.re = 0.0; (*a).w2.im = 0.0;
        (*a).w3.re = 0.0; (*a).w3.im = 0.0;
}

void lattice_matrixGID(su_3 *a, int gid, int dir, int sites){
        (*a).u1.re = sin(((dir + 1) / 1000.0 + 10.0/sites)*gid)/2.0; 
        (*a).u1.im = cos(((dir + 1) / 1000.0-10.0/sites)*gid)/3.0;
        (*a).u2.re = sin(((dir + 1) / 1000.0-18.0/sites)*gid)/4.0; 
        (*a).u2.im = cos(((dir + 1) / 1000.0+12.0/sites)*gid)/4.0;
        (*a).u3.re = sin(((dir + 1) / 1000.0-12.0/sites)*gid)/3.0; 
        (*a).u3.im = cos(((dir + 1) / 1000.0-15.0/sites)*gid)/2.0;
        
        (*a).v1.re = sin(((dir + 1) / 1000.0+15.0/sites)*gid)/3.0; 
        (*a).v1.im = cos(((dir + 1) / 1000.0-11.0/sites)*gid)/4.0;
        (*a).v2.re = sin(((dir + 1) / 1000.0-11.0/sites)*gid)/3.0; 
        (*a).v2.im = cos(((dir + 1) / 1000.0+16.0/sites)*gid)/2.0;
        (*a).v3.re = sin(((dir + 1) / 1000.0-16.0/sites)*gid)/3.0; 
        (*a).v3.im = cos(((dir + 1) / 1000.0-18.0/sites)*gid)/2.0;
        
        GramSchmidt(a);
}

su_3 operator + (su_3 a, su_3 b){
    su_3 c;
    
    c.u1.re = a.u1.re + b.u1.re;    c.u1.im = a.u1.im + b.u1.im;
    c.u2.re = a.u2.re + b.u2.re;    c.u2.im = a.u2.im + b.u2.im;
    c.u3.re = a.u3.re + b.u3.re;    c.u3.im = a.u3.im + b.u3.im;
    
    c.v1.re = a.v1.re + b.v1.re;    c.v1.im = a.v1.im + b.v1.im;
    c.v2.re = a.v2.re + b.v2.re;    c.v2.im = a.v2.im + b.v2.im;
    c.v3.re = a.v3.re + b.v3.re;    c.v3.im = a.v3.im + b.v3.im;
    
    c.w1.re = a.w1.re + b.w1.re;    c.w1.im = a.w1.im + b.w1.im;
    c.w2.re = a.w2.re + b.w2.re;    c.w2.im = a.w2.im + b.w2.im;
    c.w3.re = a.w3.re + b.w3.re;    c.w3.im = a.w3.im + b.w3.im;
    
    return c;
}

su_3 operator - (su_3 a, su_3 b){
    su_3 c;
    
    c.u1.re = a.u1.re - b.u1.re;    c.u1.im = a.u1.im - b.u1.im;
    c.u2.re = a.u2.re - b.u2.re;    c.u2.im = a.u2.im - b.u2.im;
    c.u3.re = a.u3.re - b.u3.re;    c.u3.im = a.u3.im - b.u3.im;
    
    c.v1.re = a.v1.re - b.v1.re;    c.v1.im = a.v1.im - b.v1.im;
    c.v2.re = a.v2.re - b.v2.re;    c.v2.im = a.v2.im - b.v2.im;
    c.v3.re = a.v3.re - b.v3.re;    c.v3.im = a.v3.im - b.v3.im;
    
    c.w1.re = a.w1.re - b.w1.re;    c.w1.im = a.w1.im - b.w1.im;
    c.w2.re = a.w2.re - b.w2.re;    c.w2.im = a.w2.im - b.w2.im;
    c.w3.re = a.w3.re - b.w3.re;    c.w3.im = a.w3.im - b.w3.im;
    
    return c;
}

su_3 operator * (su_3 a, su_3 b){
    su_3 c;
    
     c.u1.re = -a.u1.im * b.u1.im + a.u1.re * b.u1.re - a.u2.im * b.v1.im + a.u2.re * b.v1.re - a.u3.im * b.w1.im + a.u3.re * b.w1.re;
    c.u1.im =  a.u1.re * b.u1.im + a.u1.im * b.u1.re + a.u2.re * b.v1.im + a.u2.im * b.v1.re + a.u3.re * b.w1.im + a.u3.im * b.w1.re;
    c.u2.re = -a.u1.im * b.u2.im + a.u1.re * b.u2.re - a.u2.im * b.v2.im + a.u2.re * b.v2.re - a.u3.im * b.w2.im + a.u3.re * b.w2.re;
    c.u2.im =  a.u1.re * b.u2.im + a.u1.im * b.u2.re + a.u2.re * b.v2.im + a.u2.im * b.v2.re + a.u3.re * b.w2.im + a.u3.im * b.w2.re;
    c.u3.re = -a.u1.im * b.u3.im + a.u1.re * b.u3.re - a.u2.im * b.v3.im + a.u2.re * b.v3.re - a.u3.im * b.w3.im + a.u3.re * b.w3.re;
    c.u3.im =  a.u1.re * b.u3.im + a.u1.im * b.u3.re + a.u2.re * b.v3.im + a.u2.im * b.v3.re + a.u3.re * b.w3.im + a.u3.im * b.w3.re;

    c.v1.re = -a.v1.im * b.u1.im + a.v1.re * b.u1.re - a.v2.im * b.v1.im + a.v2.re * b.v1.re - a.v3.im * b.w1.im + a.v3.re * b.w1.re;
    c.v1.im =  a.v1.re * b.u1.im + a.v1.im * b.u1.re + a.v2.re * b.v1.im + a.v2.im * b.v1.re + a.v3.re * b.w1.im + a.v3.im * b.w1.re;
    c.v2.re = -a.v1.im * b.u2.im + a.v1.re * b.u2.re - a.v2.im * b.v2.im + a.v2.re * b.v2.re - a.v3.im * b.w2.im + a.v3.re * b.w2.re;
    c.v2.im =  a.v1.re * b.u2.im + a.v1.im * b.u2.re + a.v2.re * b.v2.im + a.v2.im * b.v2.re + a.v3.re * b.w2.im + a.v3.im * b.w2.re;
    c.v3.re = -a.v1.im * b.u3.im + a.v1.re * b.u3.re - a.v2.im * b.v3.im + a.v2.re * b.v3.re - a.v3.im * b.w3.im + a.v3.re * b.w3.re;
    c.v3.im =  a.v1.re * b.u3.im + a.v1.im * b.u3.re + a.v2.re * b.v3.im + a.v2.im * b.v3.re + a.v3.re * b.w3.im + a.v3.im * b.w3.re;

    c.w1.re = -a.w1.im * b.u1.im + a.w1.re * b.u1.re - a.w2.im * b.v1.im + a.w2.re * b.v1.re - a.w3.im * b.w1.im + a.w3.re * b.w1.re;
    c.w1.im =  a.w1.re * b.u1.im + a.w1.im * b.u1.re + a.w2.re * b.v1.im + a.w2.im * b.v1.re + a.w3.re * b.w1.im + a.w3.im * b.w1.re;
    c.w2.re = -a.w1.im * b.u2.im + a.w1.re * b.u2.re - a.w2.im * b.v2.im + a.w2.re * b.v2.re - a.w3.im * b.w2.im + a.w3.re * b.w2.re;
    c.w2.im =  a.w1.re * b.u2.im + a.w1.im * b.u2.re + a.w2.re * b.v2.im + a.w2.im * b.v2.re + a.w3.re * b.w2.im + a.w3.im * b.w2.re;
    c.w3.re = -a.w1.im * b.u3.im + a.w1.re * b.u3.re - a.w2.im * b.v3.im + a.w2.re * b.v3.re - a.w3.im * b.w3.im + a.w3.re * b.w3.re;
    c.w3.im =  a.w1.re * b.u3.im + a.w1.im * b.u3.re + a.w2.re * b.v3.im + a.w2.im * b.v3.re + a.w3.re * b.w3.im + a.w3.im * b.w3.re;
    
    return c;
}

su_3 Herm(su_3 a){
    su_3 result;
    
        result.u1.re =  a.u1.re;
        result.u1.im = -a.u1.im;
        result.u2.re =  a.v1.re;
        result.u2.im = -a.v1.im;
        result.u3.re =  a.w1.re;
        result.u3.im = -a.w1.im;

        result.v1.re =  a.u2.re;
        result.v1.im = -a.u2.im;
        result.v2.re =  a.v2.re;
        result.v2.im = -a.v2.im;
        result.v3.re =  a.w2.re;
        result.v3.im = -a.w2.im;

        result.w1.re =  a.u3.re;
        result.w1.im = -a.u3.im;
        result.w2.re =  a.v3.re;
        result.w2.im = -a.v3.im;
        result.w3.re =  a.w3.re;
        result.w3.im = -a.w3.im;
    
    return result;
}

hgpu_complex Tr(su_3 a){
    hgpu_complex tr_a;
    
    tr_a.re = a.u1.re + a.v2.re + a.w3.re;
    tr_a.im = a.u1.im + a.v2.im + a.w3.im;
    
    return tr_a;
}

hgpu_double ReTr(su_3 a){
    return (hgpu_double)(a.u1.re + a.v2.re + a.w3.re);
}

void GramSchmidt(su_3 *a){
    su_3 result;
    
    hgpu_complex sp;
    
    hgpu_double norm_u = sqrt((*a).u1.re * (*a).u1.re + (*a).u2.re * (*a).u2.re + (*a).u3.re * (*a).u3.re + (*a).u1.im * (*a).u1.im + (*a).u2.im * (*a).u2.im + (*a).u3.im * (*a).u3.im);
    hgpu_double norm_v;
    
    result.u1.re = (*a).u1.re / norm_u;
    result.u1.im = (*a).u1.im / norm_u;
    result.u2.re = (*a).u2.re / norm_u;
    result.u2.im = (*a).u2.im / norm_u;
    result.u3.re = (*a).u3.re / norm_u;
    result.u3.im = (*a).u3.im / norm_u;
    
    sp = hgpu_add(hgpu_add(hgpu_mul((*a).v1,hgpu_conjugate(result.u1)),hgpu_mul((*a).v2,hgpu_conjugate(result.u2))),hgpu_mul((*a).v3,hgpu_conjugate(result.u3)));
    result.v1 = hgpu_sub((*a).v1,hgpu_mul(result.u1,sp));
    result.v2 = hgpu_sub((*a).v2,hgpu_mul(result.u2,sp));
    result.v3 = hgpu_sub((*a).v3,hgpu_mul(result.u3,sp));
    
    norm_v = sqrt(result.v1.re * result.v1.re + result.v2.re * result.v2.re + result.v3.re * result.v3.re +
                      result.v1.im * result.v1.im + result.v2.im * result.v2.im + result.v3.im * result.v3.im);
    
    result.v1.re /= norm_v;
    result.v1.im /= norm_v;
    result.v2.re /= norm_v;
    result.v2.im /= norm_v;
    result.v3.re /= norm_v;
    result.v3.im /= norm_v;
    
    matrix_reconstruct(&result);
    *a = result;
}

void matrix_reconstruct(su_3 *a){
        (*a).w1.re =   (*a).u3.im * (*a).v2.im - (*a).u3.re * (*a).v2.re - (*a).u2.im * (*a).v3.im + (*a).u2.re * (*a).v3.re;
        (*a).w1.im =   (*a).u3.re * (*a).v2.im + (*a).u3.im * (*a).v2.re - (*a).u2.re * (*a).v3.im - (*a).u2.im * (*a).v3.re;
        (*a).w2.re =  -(*a).u3.im * (*a).v1.im + (*a).u3.re * (*a).v1.re + (*a).u1.im * (*a).v3.im - (*a).u1.re * (*a).v3.re;
        (*a).w2.im =  -(*a).u3.re * (*a).v1.im - (*a).u3.im * (*a).v1.re + (*a).u1.re * (*a).v3.im + (*a).u1.im * (*a).v3.re;
        (*a).w3.re =   (*a).u2.im * (*a).v1.im - (*a).u2.re * (*a).v1.re - (*a).u1.im * (*a).v2.im + (*a).u1.re * (*a).v2.re;
        (*a).w3.im =   (*a).u2.re * (*a).v1.im + (*a).u2.im * (*a).v1.re - (*a).u1.re * (*a).v2.im - (*a).u1.im * (*a).v2.re;
}
