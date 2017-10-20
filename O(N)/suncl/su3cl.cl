/******************************************************************************
 * @file     su3cl.cl
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.0
 *
 * @brief    [QCDGPU]
 *           Defines general procedures for lattice update (SU(3) gauge theory)
 *
 * @section  LICENSE
 *
 * Copyright (c) 2013-2017, Vadim Demchik, Natalia Kolomoyets
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

#ifndef SU3CL_CL
#define SU3CL_CL

#include "neighbours.cl"

                    HGPU_INLINE_PREFIX gpu_su_3
matrix_times3(gpu_su_3* u,gpu_su_3* v)
{
    gpu_su_3 tmp;
    hgpu_float4 u1,u2,u3;
    hgpu_float4 v1,v2,v3;

    hgpu_float u11,u12,u13;
    hgpu_float u21,u22,u23;
    hgpu_float v11,v12,v13;
    hgpu_float v21,v22,v23;

    hgpu_float w11,w12,w13;
    hgpu_float w21,w22,w23;
    hgpu_float t11,t12,t13;
    hgpu_float t21,t22,t23;

    u1 = (*u).uv1;   u2 = (*u).uv2;   u3 = (*u).uv3;
    v1 = (*v).uv1;   v2 = (*v).uv2;   v3 = (*v).uv3;

    u11 = u1.x;   u12 = u1.y;   u13 = u1.z;
    u21 = u3.x;   u22 = u3.y;   u23 = u1.w;

    v11 = u2.x;   v12 = u2.y;   v13 = u2.z;
    v21 = u3.z;   v22 = u3.w;   v23 = u2.w;

    w11 = v1.x;   w12 = v1.y;   w13 = v1.z;
    w21 = v3.x;   w22 = v3.y;   w23 = v1.w;

    t11 = v2.x;   t12 = v2.y;   t13 = v2.z;
    t21 = v3.z;   t22 = v3.w;   t23 = v2.w;

    (tmp.uv1).x = - t11*v11 - t21*v12 + u11*w11 + u12*w21 + (t13*t22 - t12*t23)*u13 - (t13*v13 + u13*w13)*w22 + (t12*v13 + u13*w12)*w23 + (t23*w12 - t22*w13)*v13;
    (tmp.uv1).y = - t12*v11 - t22*v12 + u11*w12 + u12*w22 + (t11*t23 - t13*t21)*u13 + (t21*w13 - t23*w11)*v13 + (t13*v13 + u13*w13)*w21 - (t11*v13 + u13*w11)*w23;
    (tmp.uv1).z = - t13*v11 - t23*v12 + u11*w13 + u12*w23 + (t12*t21 - t11*t22)*u13 + (t22*w11 - t21*w12)*v13 - (t12*v13 + u13*w12)*w21 + (t11*v13 + u13*w11)*w22;
    (tmp.uv1).w = - t13*v21 - t23*v22 + u21*w13 + u22*w23 + (t12*t21 - t11*t22)*u23 + (t22*w11 - t21*w12)*v23 - (t12*v23 + u23*w12)*w21 + (t11*v23 + u23*w11)*w22;

    (tmp.uv2).x =   t11*u11 + t21*u12 + v11*w11 + v12*w21 + (t22*w13 - t23*w12)*u13 + (t13*t22 - t12*t23)*v13 + (t13*u13 - v13*w13)*w22 + (v13*w12 - t12*u13)*w23;
    (tmp.uv2).y =   t12*u11 + t22*u12 + v11*w12 + v12*w22 + (t23*w11 - t21*w13)*u13 + (t11*t23 - t13*t21)*v13 + (v13*w13 - t13*u13)*w21 + (t11*u13 - v13*w11)*w23;
    (tmp.uv2).z =   t13*u11 + t23*u12 + v11*w13 + v12*w23 + (t21*w12 - t22*w11)*u13 + (t12*t21 - t11*t22)*v13 + (t12*u13 - v13*w12)*w21 + (v13*w11 - t11*u13)*w22;
    (tmp.uv2).w =   t13*u21 + t23*u22 + v21*w13 + v22*w23 + (t21*w12 - t22*w11)*u23 + (t12*t21 - t11*t22)*v23 + (t12*u23 - v23*w12)*w21 + (v23*w11 - t11*u23)*w22;

    (tmp.uv3).x = - t11*v21 - t21*v22 + u21*w11 + u22*w21 + (t13*t22 - t12*t23)*u23 + (t23*w12 - t22*w13)*v23 - (t13*v23 + u23*w13)*w22 + (t12*v23 + u23*w12)*w23;
    (tmp.uv3).y = - t12*v21 - t22*v22 + u21*w12 + u22*w22 + (t11*t23 - t13*t21)*u23 + (t21*w13 - t23*w11)*v23 + (t13*v23 + u23*w13)*w21 - (t11*v23 + u23*w11)*w23;
    (tmp.uv3).z =   t11*u21 + t21*u22 + v21*w11 + v22*w21 + (t22*w13 - t23*w12)*u23 + (t13*t22 - t12*t23)*v23 + (t13*u23 - v23*w13)*w22 + (v23*w12 - t12*u23)*w23;
    (tmp.uv3).w =   t12*u21 + t22*u22 + v21*w12 + v22*w22 + (t23*w11 - t21*w13)*u23 + (t11*t23 - t13*t21)*v23 + (v23*w13 - t13*u23)*w21 + (t11*u23 - v23*w11)*w23;

    return tmp;
}

                    HGPU_INLINE_PREFIX su_3
matrix_times_su3(su_3* u,su_3* v)
{
    su_3 tmp;

    tmp.u1.re = -(*u).u1.im * (*v).u1.im + (*u).u1.re * (*v).u1.re - (*u).u2.im * (*v).v1.im + (*u).u2.re * (*v).v1.re - (*u).u3.im * (*v).w1.im + (*u).u3.re * (*v).w1.re;
    tmp.u1.im =  (*u).u1.re * (*v).u1.im + (*u).u1.im * (*v).u1.re + (*u).u2.re * (*v).v1.im + (*u).u2.im * (*v).v1.re + (*u).u3.re * (*v).w1.im + (*u).u3.im * (*v).w1.re;

    tmp.u2.re = -(*u).u1.im * (*v).u2.im + (*u).u1.re * (*v).u2.re - (*u).u2.im * (*v).v2.im + (*u).u2.re * (*v).v2.re - (*u).u3.im * (*v).w2.im + (*u).u3.re * (*v).w2.re;
    tmp.u2.im =  (*u).u1.re * (*v).u2.im + (*u).u1.im * (*v).u2.re + (*u).u2.re * (*v).v2.im + (*u).u2.im * (*v).v2.re + (*u).u3.re * (*v).w2.im + (*u).u3.im * (*v).w2.re; 


    tmp.u3.re = -(*u).u1.im * (*v).u3.im + (*u).u1.re * (*v).u3.re - (*u).u2.im * (*v).v3.im + (*u).u2.re * (*v).v3.re - (*u).u3.im * (*v).w3.im + (*u).u3.re * (*v).w3.re;
    tmp.u3.im =  (*u).u1.re * (*v).u3.im + (*u).u1.im * (*v).u3.re + (*u).u2.re * (*v).v3.im + (*u).u2.im * (*v).v3.re + (*u).u3.re * (*v).w3.im + (*u).u3.im * (*v).w3.re;

    tmp.v1.re = -(*u).v1.im * (*v).u1.im + (*u).v1.re * (*v).u1.re - (*u).v2.im * (*v).v1.im + (*u).v2.re * (*v).v1.re - (*u).v3.im * (*v).w1.im + (*u).v3.re * (*v).w1.re;
    tmp.v1.im =  (*u).v1.re * (*v).u1.im + (*u).v1.im * (*v).u1.re + (*u).v2.re * (*v).v1.im + (*u).v2.im * (*v).v1.re + (*u).v3.re * (*v).w1.im + (*u).v3.im * (*v).w1.re;

    tmp.v2.re = -(*u).v1.im * (*v).u2.im + (*u).v1.re * (*v).u2.re - (*u).v2.im * (*v).v2.im + (*u).v2.re * (*v).v2.re - (*u).v3.im * (*v).w2.im + (*u).v3.re * (*v).w2.re;
    tmp.v2.im =  (*u).v1.re * (*v).u2.im + (*u).v1.im * (*v).u2.re + (*u).v2.re * (*v).v2.im + (*u).v2.im * (*v).v2.re + (*u).v3.re * (*v).w2.im + (*u).v3.im * (*v).w2.re;

    tmp.v3.re = -(*u).v1.im * (*v).u3.im + (*u).v1.re * (*v).u3.re - (*u).v2.im * (*v).v3.im + (*u).v2.re * (*v).v3.re - (*u).v3.im * (*v).w3.im + (*u).v3.re * (*v).w3.re;
    tmp.v3.im =  (*u).v1.re * (*v).u3.im + (*u).v1.im * (*v).u3.re + (*u).v2.re * (*v).v3.im + (*u).v2.im * (*v).v3.re + (*u).v3.re * (*v).w3.im + (*u).v3.im * (*v).w3.re;

    tmp.w1.re = -(*u).w1.im * (*v).u1.im + (*u).w1.re * (*v).u1.re - (*u).w2.im * (*v).v1.im + (*u).w2.re * (*v).v1.re - (*u).w3.im * (*v).w1.im + (*u).w3.re * (*v).w1.re;
    tmp.w1.im =  (*u).w1.re * (*v).u1.im + (*u).w1.im * (*v).u1.re + (*u).w2.re * (*v).v1.im + (*u).w2.im * (*v).v1.re + (*u).w3.re * (*v).w1.im + (*u).w3.im * (*v).w1.re;

    tmp.w2.re = -(*u).w1.im * (*v).u2.im + (*u).w1.re * (*v).u2.re - (*u).w2.im * (*v).v2.im + (*u).w2.re * (*v).v2.re - (*u).w3.im * (*v).w2.im + (*u).w3.re * (*v).w2.re;
    tmp.w2.im =  (*u).w1.re * (*v).u2.im + (*u).w1.im * (*v).u2.re + (*u).w2.re * (*v).v2.im + (*u).w2.im * (*v).v2.re + (*u).w3.re * (*v).w2.im + (*u).w3.im * (*v).w2.re;

    tmp.w3.re = -(*u).w1.im * (*v).u3.im + (*u).w1.re * (*v).u3.re - (*u).w2.im * (*v).v3.im + (*u).w2.re * (*v).v3.re - (*u).w3.im * (*v).w3.im + (*u).w3.re * (*v).w3.re;
    tmp.w3.im =  (*u).w1.re * (*v).u3.im + (*u).w1.im * (*v).u3.re + (*u).w2.re * (*v).v3.im + (*u).w2.im * (*v).v3.re + (*u).w3.re * (*v).w3.im + (*u).w3.im * (*v).w3.re;

    return tmp;
}

                    HGPU_INLINE_PREFIX su_3
lattice_reconstruct3(gpu_su_3* a)
{
    su_3 b;

    b.u1.re = (*a).uv1.x;
    b.u1.im = (*a).uv2.x;
    b.u2.re = (*a).uv1.y;
    b.u2.im = (*a).uv2.y;
    b.u3.re = (*a).uv1.z;
    b.u3.im = (*a).uv2.z;

    b.v1.re = (*a).uv3.x;
    b.v1.im = (*a).uv3.z;
    b.v2.re = (*a).uv3.y;
    b.v2.im = (*a).uv3.w;
    b.v3.re = (*a).uv1.w;
    b.v3.im = (*a).uv2.w;

    b.w1.re =   b.u2.re * b.v3.re - b.u2.im * b.v3.im - b.u3.re * b.v2.re + b.u3.im * b.v2.im;
    b.w1.im =  -b.u2.re * b.v3.im - b.u2.im * b.v3.re + b.u3.re * b.v2.im + b.u3.im * b.v2.re;
    b.w2.re =   b.u3.re * b.v1.re - b.u3.im * b.v1.im - b.u1.re * b.v3.re + b.u1.im * b.v3.im;
    b.w2.im =  -b.u3.re * b.v1.im - b.u3.im * b.v1.re + b.u1.re * b.v3.im + b.u1.im * b.v3.re;
    b.w3.re =   b.u1.re * b.v2.re - b.u1.im * b.v2.im - b.u2.re * b.v1.re + b.u2.im * b.v1.im;
    b.w3.im =  -b.u1.re * b.v2.im - b.u1.im * b.v2.re + b.u2.re * b.v1.im + b.u2.im * b.v1.re;

    return b;
}

#endif
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
