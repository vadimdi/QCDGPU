/******************************************************************************
 * @file     su3_measurements.cl
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.5
 *
 * @brief    [QCDGPU]
 *           Definition of general functions used in measurements, corresponding to the SU(3) gauge group
 *
 * @section  LICENSE
 *
 * Copyright (c) 2013, 2014 Vadim Demchik, Natalia Kolomoyets
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
#ifndef SU3MEASUREMENTSCL_CL
#define SU3MEASUREMENTSCL_CL

                     __attribute__((always_inline)) void
lattice_lambda1(double_su_3* matrix)
{
    (*matrix).u1.re = 0.0;
    (*matrix).u1.im = 0.0;

    (*matrix).u2.re = 1.0;
    (*matrix).u2.im = 0.0;

    (*matrix).u3.re = 0.0;
    (*matrix).u3.im = 0.0;

    (*matrix).v1.re = 1.0;
    (*matrix).v1.im = 0.0;

    (*matrix).v2.re = 0.0;
    (*matrix).v2.im = 0.0;

    (*matrix).v3.re = 0.0;
    (*matrix).v3.im = 0.0;

    (*matrix).w1.re = 0.0;
    (*matrix).w1.im = 0.0;

    (*matrix).w2.re = 0.0;
    (*matrix).w2.im = 0.0;

    (*matrix).w3.re = 0.0;
    (*matrix).w3.im = 0.0;
}

                     __attribute__((always_inline)) void
lattice_lambda2(double_su_3* matrix)
{
    (*matrix).u1.re = 0.0;
    (*matrix).u1.im = 0.0;

    (*matrix).u2.re = 0.0;
    (*matrix).u2.im =-1.0;

    (*matrix).u3.re = 0.0;
    (*matrix).u3.im = 0.0;

    (*matrix).v1.re = 0.0;
    (*matrix).v1.im = 1.0;

    (*matrix).v2.re = 0.0;
    (*matrix).v2.im = 0.0;

    (*matrix).v3.re = 0.0;
    (*matrix).v3.im = 0.0;

    (*matrix).w1.re = 0.0;
    (*matrix).w1.im = 0.0;

    (*matrix).w2.re = 0.0;
    (*matrix).w2.im = 0.0;

    (*matrix).w3.re = 0.0;
    (*matrix).w3.im = 0.0;
}

                     __attribute__((always_inline)) void
lattice_lambda3(double_su_3* matrix)
{
    (*matrix).u1.re = 1.0;
    (*matrix).u1.im = 0.0;

    (*matrix).u2.re = 0.0;
    (*matrix).u2.im = 0.0;

    (*matrix).u3.re = 0.0;
    (*matrix).u3.im = 0.0;

    (*matrix).v1.re = 0.0;
    (*matrix).v1.im = 0.0;

    (*matrix).v2.re =-1.0;
    (*matrix).v2.im = 0.0;

    (*matrix).v3.re = 0.0;
    (*matrix).v3.im = 0.0;

    (*matrix).w1.re = 0.0;
    (*matrix).w1.im = 0.0;

    (*matrix).w2.re = 0.0;
    (*matrix).w2.im = 0.0;

    (*matrix).w3.re = 0.0;
    (*matrix).w3.im = 0.0;
}

                     __attribute__((always_inline)) void
lattice_lambda4(double_su_3* matrix)
{
    (*matrix).u1.re = 0.0;
    (*matrix).u1.im = 0.0;

    (*matrix).u2.re = 0.0;
    (*matrix).u2.im = 0.0;

    (*matrix).u3.re = 1.0;
    (*matrix).u3.im = 0.0;

    (*matrix).v1.re = 0.0;
    (*matrix).v1.im = 0.0;

    (*matrix).v2.re = 0.0;
    (*matrix).v2.im = 0.0;

    (*matrix).v3.re = 0.0;
    (*matrix).v3.im = 0.0;

    (*matrix).w1.re = 1.0;
    (*matrix).w1.im = 0.0;

    (*matrix).w2.re = 0.0;
    (*matrix).w2.im = 0.0;

    (*matrix).w3.re = 0.0;
    (*matrix).w3.im = 0.0;
}

                     __attribute__((always_inline)) void
lattice_lambda5(double_su_3* matrix)
{
    (*matrix).u1.re = 0.0;
    (*matrix).u1.im = 0.0;

    (*matrix).u2.re = 0.0;
    (*matrix).u2.im = 0.0;

    (*matrix).u3.re = 0.0;
    (*matrix).u3.im =-1.0;

    (*matrix).v1.re = 0.0;
    (*matrix).v1.im = 0.0;

    (*matrix).v2.re = 0.0;
    (*matrix).v2.im = 0.0;

    (*matrix).v3.re = 0.0;
    (*matrix).v3.im = 0.0;

    (*matrix).w1.re = 0.0;
    (*matrix).w1.im = 1.0;

    (*matrix).w2.re = 0.0;
    (*matrix).w2.im = 0.0;

    (*matrix).w3.re = 0.0;
    (*matrix).w3.im = 0.0;
}

                     __attribute__((always_inline)) void
lattice_lambda6(double_su_3* matrix)
{
    (*matrix).u1.re = 0.0;
    (*matrix).u1.im = 0.0;

    (*matrix).u2.re = 0.0;
    (*matrix).u2.im = 0.0;

    (*matrix).u3.re = 0.0;
    (*matrix).u3.im = 0.0;

    (*matrix).v1.re = 0.0;
    (*matrix).v1.im = 0.0;

    (*matrix).v2.re = 0.0;
    (*matrix).v2.im = 0.0;

    (*matrix).v3.re = 1.0;
    (*matrix).v3.im = 0.0;

    (*matrix).w1.re = 0.0;
    (*matrix).w1.im = 0.0;

    (*matrix).w2.re = 1.0;
    (*matrix).w2.im = 0.0;

    (*matrix).w3.re = 0.0;
    (*matrix).w3.im = 0.0;
}

                     __attribute__((always_inline)) void
lattice_lambda7(double_su_3* matrix)
{
    (*matrix).u1.re = 0.0;
    (*matrix).u1.im = 0.0;

    (*matrix).u2.re = 0.0;
    (*matrix).u2.im = 0.0;

    (*matrix).u3.re = 0.0;
    (*matrix).u3.im = 0.0;

    (*matrix).v1.re = 0.0;
    (*matrix).v1.im = 0.0;

    (*matrix).v2.re = 0.0;
    (*matrix).v2.im = 0.0;

    (*matrix).v3.re = 0.0;
    (*matrix).v3.im =-1.0;

    (*matrix).w1.re = 0.0;
    (*matrix).w1.im = 0.0;

    (*matrix).w2.re = 0.0;
    (*matrix).w2.im = 1.0;

    (*matrix).w3.re = 0.0;
    (*matrix).w3.im = 0.0;
}

                     __attribute__((always_inline)) void
lattice_lambda8(double_su_3* matrix)
{
    (*matrix).u1.re = 0.57735026918962576450914878050196;   // 1/Sqrt(3)
    (*matrix).u1.im = 0.0;

    (*matrix).u2.re = 0.0;
    (*matrix).u2.im = 0.0;

    (*matrix).u3.re = 0.0;
    (*matrix).u3.im = 0.0;

    (*matrix).v1.re = 0.0;
    (*matrix).v1.im = 0.0;

    (*matrix).v2.re = 0.57735026918962576450914878050196;   // 1/Sqrt(3)
    (*matrix).v2.im = 0.0;

    (*matrix).v3.re = 0.0;
    (*matrix).v3.im = 0.0;

    (*matrix).w1.re = 0.0;
    (*matrix).w1.im = 0.0;

    (*matrix).w2.re = 0.0;
    (*matrix).w2.im = 0.0;

    (*matrix).w3.re = -1.1547005383792515290182975610039;   //-2/Sqrt(3)
    (*matrix).w3.im = 0.0;
}



// _________________ double precision ________________________
                    __attribute__((always_inline)) double_su_3
lattice_reconstruct3_double(gpu_su_3* a)
{
    double_su_3 b;

    b.u1.re = (hgpu_double) (*a).uv1.x;
    b.u1.im = (hgpu_double) (*a).uv2.x;
    b.u2.re = (hgpu_double) (*a).uv1.y;
    b.u2.im = (hgpu_double) (*a).uv2.y;
    b.u3.re = (hgpu_double) (*a).uv1.z;
    b.u3.im = (hgpu_double) (*a).uv2.z;

    b.v1.re = (hgpu_double) (*a).uv3.x;
    b.v1.im = (hgpu_double) (*a).uv3.z;
    b.v2.re = (hgpu_double) (*a).uv3.y;
    b.v2.im = (hgpu_double) (*a).uv3.w;
    b.v3.re = (hgpu_double) (*a).uv1.w;
    b.v3.im = (hgpu_double) (*a).uv2.w;

    b.w1.re =   b.u2.re * b.v3.re - b.u2.im * b.v3.im - b.u3.re * b.v2.re + b.u3.im * b.v2.im;
    b.w1.im =  -b.u2.re * b.v3.im - b.u2.im * b.v3.re + b.u3.re * b.v2.im + b.u3.im * b.v2.re;
    b.w2.re =   b.u3.re * b.v1.re - b.u3.im * b.v1.im - b.u1.re * b.v3.re + b.u1.im * b.v3.im;
    b.w2.im =  -b.u3.re * b.v1.im - b.u3.im * b.v1.re + b.u1.re * b.v3.im + b.u1.im * b.v3.re;
    b.w3.re =   b.u1.re * b.v2.re - b.u1.im * b.v2.im - b.u2.re * b.v1.re + b.u2.im * b.v1.im;
    b.w3.im =  -b.u1.re * b.v2.im - b.u1.im * b.v2.re + b.u2.re * b.v1.im + b.u2.im * b.v1.re;

    return b;
}

                    __attribute__((always_inline)) __private double_su_3
matrix_times_su3_double(double_su_3* u,double_su_3* v)
{
    double_su_3 tmp;

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

                    __attribute__((always_inline)) __private hgpu_complex_double
matrix_trace3_double(const double_su_3* u)
{
    hgpu_complex_double tmp;
       tmp.re = (*u).u1.re + (*u).v2.re + (*u).w3.re;
       tmp.im = (*u).u1.im + (*u).v2.im + (*u).w3.im;
    return tmp;
}

                    __attribute__((always_inline)) __private double_su_3
matrix_hermitian_su3_double(double_su_3* a)
{
    double_su_3 result;

        result.u1.re =  (*a).u1.re;
        result.u1.im = -(*a).u1.im;
        result.u2.re =  (*a).v1.re;
        result.u2.im = -(*a).v1.im;
        result.u3.re =  (*a).w1.re;
        result.u3.im = -(*a).w1.im;

        result.v1.re =  (*a).u2.re;
        result.v1.im = -(*a).u2.im;
        result.v2.re =  (*a).v2.re;
        result.v2.im = -(*a).v2.im;
        result.v3.re =  (*a).w2.re;
        result.v3.im = -(*a).w2.im;

        result.w1.re =  (*a).u3.re;
        result.w1.im = -(*a).u3.im;
        result.w2.re =  (*a).v3.re;
        result.w2.im = -(*a).v3.im;
        result.w3.re =  (*a).w3.re;
        result.w3.im = -(*a).w3.im;

    return result;
}

                    __attribute__((always_inline)) __private hgpu_double
matrix_retrace3_double(const double_su_3* u)
{
    hgpu_double tmp;
       tmp = (*u).u1.re + (*u).v2.re + (*u).w3.re;
    return tmp;
}

                    __attribute__((always_inline)) __private hgpu_double
lattice_retrace_plaquette3(gpu_su_3* u1, gpu_su_3* u2, gpu_su_3* u3, gpu_su_3* u4)
{
    double_su_3 m1, m2, m3, m4;
    double_su_3 w1, w2, w3;
    hgpu_double result;

    m1 = lattice_reconstruct3_double(u1);
    m2 = lattice_reconstruct3_double(u2);
    m3 = lattice_reconstruct3_double(u3);
    m4 = lattice_reconstruct3_double(u4);

    w1 = matrix_times_su3_double(&m1,&m2);
    w2 = matrix_hermitian_su3_double(&m3);
    w3 = matrix_times_su3_double(&w1,&w2);
    w2 = matrix_hermitian_su3_double(&m4);
    w1 = matrix_times_su3_double(&w3,&w2);

    result = matrix_retrace3_double(&w1);

    return result;
}

                    __attribute__((always_inline)) __private hgpu_double
lattice_retrace_plaquette3_F(gpu_su_3* u1, gpu_su_3* u2, gpu_su_3* u3, gpu_su_3* u4, hgpu_complex_double* F3, hgpu_complex_double* F8){
    double_su_3 m1, m2, m3, m4;
    double_su_3 w1, w2, w3, w4, w5;
    hgpu_double result;

    m1 = lattice_reconstruct3_double(u1);
    m2 = lattice_reconstruct3_double(u2);
    m3 = lattice_reconstruct3_double(u3);
    m4 = lattice_reconstruct3_double(u4);

    w1 = matrix_times_su3_double(&m1,&m2);
    w2 = matrix_hermitian_su3_double(&m3);
    w3 = matrix_times_su3_double(&w1,&w2);
    w2 = matrix_hermitian_su3_double(&m4);
    w1 = matrix_times_su3_double(&w3,&w2);

#ifdef FMUNU1
    lattice_lambda1(&w4);
#elif defined(FMUNU2)
    lattice_lambda2(&w4);
#elif defined(FMUNU4)
    lattice_lambda4(&w4);
#else
    lattice_lambda3(&w4);
#endif
    w5 = matrix_times_su3_double(&w1,&w4);
    (*F3) = matrix_trace3_double(&w5);

#ifdef FMUNU5
    lattice_lambda5(&w4);
#elif defined(FMUNU6)
    lattice_lambda6(&w4);
#elif defined(FMUNU7)
    lattice_lambda7(&w4);
#else
    lattice_lambda8(&w4);
#endif
    w5 = matrix_times_su3_double(&w1,&w4);
    (*F8) = matrix_trace3_double(&w5);

    result = matrix_retrace3_double(&w1);

    return result;
}

                    __attribute__((always_inline)) __private hgpu_double
matrix_retrace_su3(const su_3* u)
{
    hgpu_double tmp = (hgpu_double) (*u).u1.re + (hgpu_double) (*u).v2.re + (hgpu_double) (*u).w3.re;
    return tmp;
}

                    __attribute__((always_inline)) __private hgpu_double
matrix_imtrace_su3(const su_3* u)
{
    hgpu_double tmp = (hgpu_double) (*u).u1.im + (hgpu_double) (*u).v2.im + (hgpu_double) (*u).w3.im;
    return tmp;
}

                    __attribute__((always_inline)) __private hgpu_double
lattice_retrace_plaquette3_double(double_su_3* m1, double_su_3* m2, double_su_3* m3, double_su_3* m4)
{
    double_su_3 w1, w2, w3;
    hgpu_double result;

    w1 = matrix_times_su3_double(m1,m2);
    w2 = matrix_hermitian_su3_double(m3);
    w3 = matrix_times_su3_double(&w1,&w2);
    w2 = matrix_hermitian_su3_double(m4);
    w1 = matrix_times_su3_double(&w3,&w2);

    result = matrix_retrace3_double(&w1);

    return result;
}


#endif
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                   
                                                                                                                                                                   
