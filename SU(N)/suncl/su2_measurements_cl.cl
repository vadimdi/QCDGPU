/******************************************************************************
 * @file     su2_measurements.cl
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.6
 *
 * @brief    [QCDGPU]
 *           Definition of general functions used in measurements, corresponding to the SU(2) gauge group
 *
 * @section  LICENSE
 *
 * Copyright (c) 2013-2017 Vadim Demchik, Natalia Kolomoyets
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
#ifndef SU2MEASUREMENTSCL_CL
#define SU2MEASUREMENTSCL_CL

                     HGPU_INLINE_PREFIX_VOID void
lattice_sigma1(double_su_2* matrix)
{
    (*matrix).u1.re = 0.0;
    (*matrix).u1.im = 0.0;

    (*matrix).u2.re = 1.0;
    (*matrix).u2.im = 0.0;

    (*matrix).v1.re = 1.0;
    (*matrix).v1.im = 0.0;

    (*matrix).v2.re = 0.0;
    (*matrix).v2.im = 0.0;
}

                     HGPU_INLINE_PREFIX_VOID void
lattice_sigma2(double_su_2* matrix)
{
    (*matrix).u1.re = 0.0;
    (*matrix).u1.im = 0.0;

    (*matrix).u2.re = 0.0;
    (*matrix).u2.im =-1.0;

    (*matrix).v1.re = 0.0;
    (*matrix).v1.im = 1.0;

    (*matrix).v2.re = 0.0;
    (*matrix).v2.im = 0.0;
}

                     HGPU_INLINE_PREFIX_VOID void
lattice_sigma3(double_su_2* matrix)
{
    (*matrix).u1.re = 1.0;
    (*matrix).u1.im = 0.0;

    (*matrix).u2.re = 0.0;
    (*matrix).u2.im = 0.0;

    (*matrix).v1.re = 0.0;
    (*matrix).v1.im = 0.0;

    (*matrix).v2.re =-1.0;
    (*matrix).v2.im = 0.0;
}

                    HGPU_INLINE_PREFIX hgpu_complex_double
trace_sigma1_double(double_su_2* u)
{
    // -i Tr(u*sigma_1/2)
    hgpu_complex_double tr;
    
    tr.re = (*u).u2.im;
    tr.im = 0.0;
    
    return tr;
}

                    HGPU_INLINE_PREFIX hgpu_complex_double
trace_sigma2_double(double_su_2* u)
{
    // -i Tr(u*sigma_2/2)
    hgpu_complex_double tr;
    
    tr.re = (*u).u2.re;
    tr.im = 0.0;
    
    return tr;
}

                    HGPU_INLINE_PREFIX hgpu_complex_double
trace_sigma3_double(double_su_2* u)
{
    // -i Tr(u*sigma_3/2)
    hgpu_complex_double tr;
    
    tr.re = (*u).u1.im;
    tr.im = 0.0;
    
    return tr;
}

// _________________ double precision ________________________
                    HGPU_INLINE_PREFIX_VOID double_su_2
lattice_reconstruct2_double(gpu_su_2* a)
{
    double_su_2 b;

    b.u1.re = (hgpu_double) (*a).uv1.x;
    b.u1.im = (hgpu_double) (*a).uv1.z;
    b.u2.re = (hgpu_double) (*a).uv1.y;
    b.u2.im = (hgpu_double) (*a).uv1.w;

    b.v1.re = -(hgpu_double) (*a).uv1.y;
    b.v1.im =  (hgpu_double) (*a).uv1.w;
    b.v2.re =  (hgpu_double) (*a).uv1.x;
    b.v2.im = -(hgpu_double) (*a).uv1.z;

    return b;
}

                    HGPU_INLINE_PREFIX double_su_2
matrix_times_su2_double(double_su_2* u,double_su_2* v)
{
    double_su_2 tmp;

    tmp.u1.re = -(*u).u1.im * (*v).u1.im + (*u).u1.re * (*v).u1.re - (*u).u2.im * (*v).v1.im + (*u).u2.re * (*v).v1.re;
    tmp.u1.im =  (*u).u1.re * (*v).u1.im + (*u).u1.im * (*v).u1.re + (*u).u2.re * (*v).v1.im + (*u).u2.im * (*v).v1.re;

    tmp.u2.re = -(*u).u1.im * (*v).u2.im + (*u).u1.re * (*v).u2.re - (*u).u2.im * (*v).v2.im + (*u).u2.re * (*v).v2.re;
    tmp.u2.im =  (*u).u1.re * (*v).u2.im + (*u).u1.im * (*v).u2.re + (*u).u2.re * (*v).v2.im + (*u).u2.im * (*v).v2.re; 


    tmp.v1.re = -(*u).v1.im * (*v).u1.im + (*u).v1.re * (*v).u1.re - (*u).v2.im * (*v).v1.im + (*u).v2.re * (*v).v1.re;
    tmp.v1.im =  (*u).v1.re * (*v).u1.im + (*u).v1.im * (*v).u1.re + (*u).v2.re * (*v).v1.im + (*u).v2.im * (*v).v1.re;

    tmp.v2.re = -(*u).v1.im * (*v).u2.im + (*u).v1.re * (*v).u2.re - (*u).v2.im * (*v).v2.im + (*u).v2.re * (*v).v2.re;
    tmp.v2.im =  (*u).v1.re * (*v).u2.im + (*u).v1.im * (*v).u2.re + (*u).v2.re * (*v).v2.im + (*u).v2.im * (*v).v2.re;

    return tmp;
}

                    HGPU_INLINE_PREFIX hgpu_complex_double
matrix_trace2_double(const double_su_2* u)
{
    hgpu_complex_double tmp;
       tmp.re = (*u).u1.re + (*u).v2.re;
       tmp.im = (*u).u1.im + (*u).v2.im;
    return tmp;
}

#ifdef BIGLAT
                    HGPU_INLINE_PREFIX gpu_su_2
matrix_hermitian_gpu_su_2(gpu_su_2* a)
{
    gpu_su_2 m;
    
    m.uv1 = (hgpu_float4)((*a).uv1.x, -(*a).uv1.y, -(*a).uv1.z, -(*a).uv1.w);

    return m;
}
#endif

                    HGPU_INLINE_PREFIX double_su_2
matrix_hermitian_su2_double(double_su_2* a)
{
    double_su_2 result;

        result.u1.re =  (*a).u1.re;
        result.u1.im = -(*a).u1.im;
        result.u2.re =  (*a).v1.re;
        result.u2.im = -(*a).v1.im;

        result.v1.re =  (*a).u2.re;
        result.v1.im = -(*a).u2.im;
        result.v2.re =  (*a).v2.re;
        result.v2.im = -(*a).v2.im;

    return result;
}

                    HGPU_INLINE_PREFIX hgpu_double
matrix_retrace2_double(const double_su_2* u)
{
    hgpu_double tmp;
       tmp = (*u).u1.re + (*u).v2.re;
    return tmp;
}

                    HGPU_INLINE_PREFIX hgpu_double
lattice_retrace_plaquette2(gpu_su_2* u1, gpu_su_2* u2, gpu_su_2* u3, gpu_su_2* u4)
{
    double_su_2 m1, m2, m3, m4;
    double_su_2 w1, w2, w3;
    hgpu_double result;

    m1 = lattice_reconstruct2_double(u1);
    m2 = lattice_reconstruct2_double(u2);
    m3 = lattice_reconstruct2_double(u3);
    m4 = lattice_reconstruct2_double(u4);

    w1 = matrix_times_su2_double(&m1,&m2);
    w2 = matrix_hermitian_su2_double(&m3);
    w3 = matrix_times_su2_double(&w1,&w2);
    w2 = matrix_hermitian_su2_double(&m4);
    w1 = matrix_times_su2_double(&w3,&w2);
    result = matrix_retrace2_double(&w1);

    return result;
}

                    HGPU_INLINE_PREFIX hgpu_double
lattice_retrace_plaquette2_F(gpu_su_2* u1, gpu_su_2* u2, gpu_su_2* u3, gpu_su_2* u4, hgpu_complex_double* F3, hgpu_complex_double* F8){
    double_su_2 m1, m2, m3, m4;
    double_su_2 w1, w2, w3;
    hgpu_double result;

    m1 = lattice_reconstruct2_double(u1);
    m2 = lattice_reconstruct2_double(u2);
    m3 = lattice_reconstruct2_double(u3);
    m4 = lattice_reconstruct2_double(u4);

    w1 = matrix_times_su2_double(&m1,&m2);
    w2 = matrix_hermitian_su2_double(&m3);
    w3 = matrix_times_su2_double(&w1,&w2);
    w2 = matrix_hermitian_su2_double(&m4);
    w1 = matrix_times_su2_double(&w3,&w2);

#ifdef FMUNU1
    (*F3) = trace_sigma1_double(&w1);
#else
    (*F3) = trace_sigma3_double(&w1);
#endif

//FMUNU2
    (*F8) = trace_sigma2_double(&w1);

    result = matrix_retrace2_double(&w1);

    return result;
}

                    HGPU_INLINE_PREFIX hgpu_double
matrix_retrace_su2(const su_2* u)
{
    hgpu_double tmp = (hgpu_double) (*u).u1.re + (hgpu_double) (*u).v2.re;
    return tmp;
}

                    HGPU_INLINE_PREFIX hgpu_double
matrix_imtrace_su2(const su_2* u)
{
    hgpu_double tmp = (hgpu_double) (*u).u1.im + (hgpu_double) (*u).v2.im;
    return tmp;
}

                    HGPU_INLINE_PREFIX hgpu_double
lattice_retrace_plaquette2_double(double_su_2* m1, double_su_2* m2, double_su_2* m3, double_su_2* m4)
{
    double_su_2 w1, w2, w3;
    hgpu_double result;

    w1 = matrix_times_su2_double(m1,m2);
    w2 = matrix_hermitian_su2_double(m3);
    w3 = matrix_times_su2_double(&w1,&w2);
    w2 = matrix_hermitian_su2_double(m4);
    w1 = matrix_times_su2_double(&w3,&w2);
    result = matrix_retrace2_double(&w1);

    return result;
}


#endif
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                   
                                                                                                                                                                   
                                                                                                                                                                   
