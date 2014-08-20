/******************************************************************************
 * @file     complex.h
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.5
 *
 * @brief    [QCDGPU]
 *           Defines basic constants, types and algebra used in program
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
#ifndef COMPLEX_H
#define COMPLEX_H

// complex algebra
#if defined(cl_khr_fp64)

#define hgpu_double     double
#define hgpu_double2    double2
#define hgpu_double3    double3
#define hgpu_double4    double4

#if defined(cl_amd_fp64)
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
#elif defined(cl_khr_fp64)
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
#endif

#ifndef PRECISION
#define PRECISION   0       // single precision by default
#endif

#if (PRECISION==0)
#define PRECISION_SINGLE            // define single precision
#elif (PRECISION==1)
#define PRECISION_DOUBLE            // define double precision
#elif (PRECISION==2)
#define PRECISION_MIXED             // define mixed precision
#endif

#ifdef PRECISION_DOUBLE
#define hgpu_float   hgpu_double
#define hgpu_float2  hgpu_double2
#define hgpu_float3  hgpu_double3
#define hgpu_float4  hgpu_double4
#ifndef PI
#define PI    3.1415926535897932384626433832795    // pi
#endif
#ifndef PI2
#define PI2   6.2831853071795864769252867665590    // 2*pi
#endif
#else
#define hgpu_float  float
#define hgpu_float2 float2
#define hgpu_float3 float3
#define hgpu_float4 float4
#ifndef PI
#define PI    3.141592653589793f    // pi
#endif
#ifndef PI2
#define PI2   6.283185307179586f    // 2*pi
#endif
#endif

#define hgpu_single  float
#define hgpu_single4 float4
/**
 * Defines complex data type for ortogonalization of vectors(Gram-Schmidt procedure) with optional single or double precision.
 */
typedef struct hgpu_complex{
    hgpu_float re; /**< Real part of complex number (optional float/double)*/
    hgpu_float im; /**< Imaginary part of complex number (optional float/double)*/
} hgpu_complex;

/**
 * Defines complex data type for calculations with single precision.
 */
typedef struct hgpu_complex_single{
    float re; /**< Real part of complex number (float)*/
    float im; /**< Imaginary part of complex number (float)*/
} hgpu_complex_single;

/**
 * Defines complex data type for measurements with optional single/double precision.
 */
typedef struct hgpu_complex_double{
    hgpu_double re; /**< Real part of complex number (double)*/
    hgpu_double im; /**< Imaginary part of complex number (double)*/
} hgpu_complex_double;

#ifndef NATIVESPEEDUP
    #define native_sqrt(x)      sqrt(x)
    #define native_exp(x)       exp(x)
    #define native_divide(x,y)  (x/y)
    /**
    * The OpenCL sincos function adopted to the environment (returns sine of number and plases the cosine in the second argument)
    * @param x -- the argument of trigonometric function.
    * @param y -- pointer to the cosine result.
    */
    inline hgpu_float hgpu_sincos(hgpu_float x, hgpu_float* y){
        *y = cos((hgpu_float) x);
        return sin((hgpu_float) x);
    }
#else
    inline hgpu_float hgpu_sincos(hgpu_float x, hgpu_float* y){
        x = sincos(x,y);
        return ((hgpu_float) x);
    }
#endif

    hgpu_float hgpu_Re(hgpu_complex c);
    hgpu_float hgpu_Im(hgpu_complex c);
    bool hgpu_cmp(hgpu_complex a,hgpu_complex b);
    hgpu_float hgpu_abs(hgpu_complex c);
    hgpu_float hgpu_phase(hgpu_complex c);
    hgpu_complex hgpu_conjugate(hgpu_complex c);
    hgpu_complex hgpu_power(hgpu_complex c, hgpu_float y);
    hgpu_complex hgpu_sqrt(hgpu_complex c);
    hgpu_complex hgpu_exp(hgpu_complex c);
    hgpu_complex hgpu_sin(hgpu_complex c);
    hgpu_complex hgpu_cos(hgpu_complex c);
    hgpu_complex hgpu_minus(hgpu_complex c);
    hgpu_complex hgpu_add(hgpu_complex a,hgpu_complex b);
    hgpu_complex hgpu_sub(hgpu_complex a,hgpu_complex b);
    hgpu_complex hgpu_mul(hgpu_complex a,hgpu_complex b);
    hgpu_complex hgpu_div(hgpu_complex a,hgpu_complex b);
    hgpu_complex hgpu_float_to_complex(hgpu_float a);
    hgpu_complex hgpu_int_to_complex(int a);
    hgpu_complex hgpu_I();

    hgpu_complex_double hgpu_conjugate_double(hgpu_complex_double c);
    hgpu_complex_double hgpu_mul_double(hgpu_complex_double a,hgpu_complex_double b);
    hgpu_complex_double hgpu_add_double(hgpu_complex_double a,hgpu_complex_double b);
    hgpu_complex_double hgpu_sub_double(hgpu_complex_double a,hgpu_complex_double b);
#else

#pragma warning(disable:4505)

#define hgpu_double     double
#define hgpu_double2    double2
#define hgpu_double3    double3
#define hgpu_double4    double4

#define hgpu_float  hgpu_double
#define hgpu_single hgpu_double

#ifndef PI
#define PI    3.1415926535897932384626433832795    // pi
#endif
#ifndef PI2
#define PI2   6.2831853071795864769252867665590    // 2*pi
#endif

#define native_sqrt(x)      sqrt(x)
#define native_exp(x)       exp(x)
#define native_divide(x,y)  (x/y)

typedef struct hgpu_complex{
    hgpu_double re; /**< Real part of complex number (optional float/double)*/
    hgpu_double im; /**< Imaginary part of complex number (optional float/double)*/
} hgpu_complex;

typedef struct hgpu_complex_double{
    hgpu_double re; /**< Real part of complex number (optional float/double)*/
    hgpu_double im; /**< Imaginary part of complex number (optional float/double)*/
} hgpu_complex_double;

    inline hgpu_float sincos(hgpu_float x, hgpu_float* y){
        *y = cos((hgpu_float) x);
        return sin((hgpu_float) x);
    }

    inline hgpu_float hgpu_sincos(hgpu_float x, hgpu_float* y){
        *y = cos((hgpu_float) x);
        return sin((hgpu_float) x);
    }
    static hgpu_float hgpu_Re(hgpu_complex c);
    static hgpu_float hgpu_Im(hgpu_complex c);
    static bool hgpu_cmp(hgpu_complex a,hgpu_complex b);
    static hgpu_float hgpu_abs(hgpu_complex c);
    static hgpu_float hgpu_phase(hgpu_complex c);
    static hgpu_complex hgpu_conjugate(hgpu_complex c);
    static hgpu_complex hgpu_power(hgpu_complex c, hgpu_float y);
    static hgpu_complex hgpu_sqrt(hgpu_complex c);
    static hgpu_complex hgpu_exp(hgpu_complex c);
    static hgpu_complex hgpu_sin(hgpu_complex c);
    static hgpu_complex hgpu_cos(hgpu_complex c);
    static hgpu_complex hgpu_minus(hgpu_complex c);
    static hgpu_complex hgpu_add(hgpu_complex a,hgpu_complex b);
    static hgpu_complex hgpu_sub(hgpu_complex a,hgpu_complex b);
    static hgpu_complex hgpu_mul(hgpu_complex a,hgpu_complex b);
    static hgpu_complex hgpu_div(hgpu_complex a,hgpu_complex b);
    static hgpu_complex hgpu_float_to_complex(hgpu_float a);
    static hgpu_complex hgpu_int_to_complex(int a);
    static hgpu_complex hgpu_I();

    static hgpu_complex_double hgpu_conjugate_double(hgpu_complex_double c);
    static hgpu_complex_double hgpu_mul_double(hgpu_complex_double a,hgpu_complex_double b);
    static hgpu_complex_double hgpu_add_double(hgpu_complex_double a,hgpu_complex_double b);
    static hgpu_complex_double hgpu_sub_double(hgpu_complex_double a,hgpu_complex_double b);
#endif

hgpu_float hgpu_Re(hgpu_complex c){
    return c.re;
}

hgpu_float hgpu_Im(hgpu_complex c){
    return c.im;
}

bool hgpu_cmp(hgpu_complex a,hgpu_complex b) {
    return ((a.re==b.re) && (a.im==b.im));
}

hgpu_float hgpu_abs(hgpu_complex c) {
    hgpu_float d = c.re*c.re+c.im*c.im;
    return (hgpu_float) native_sqrt((hgpu_single) d);
}

hgpu_float hgpu_phase(hgpu_complex c) {
    return (hgpu_float) atan(native_divide((hgpu_single) c.im, (hgpu_single) c.re));
}

hgpu_complex hgpu_conjugate(hgpu_complex c) {
    hgpu_complex a;
        a.re =  c.re;
        a.im = -c.im;
    return a;
}


hgpu_complex hgpu_power(hgpu_complex c, hgpu_float y) {
    hgpu_complex a;
    hgpu_float bcos,bsin;
    hgpu_float r = (hgpu_float) pow((hgpu_float) hgpu_abs(c),(hgpu_float) y);
    hgpu_float p = y * hgpu_phase(c);
    bsin = (hgpu_float) hgpu_sincos(p,&bcos);
    a.re = r * bcos;
    a.im = r * bsin;
    return a;
}

hgpu_complex hgpu_sqrt(hgpu_complex c) {
    hgpu_complex a;
    hgpu_float bcos,bsin;
    hgpu_float r = native_sqrt((hgpu_single) hgpu_abs(c));
    hgpu_float p = 0.5 * hgpu_phase(c);
    bsin = hgpu_sincos(p,&bcos);
    a.re = r * bcos;
    a.im = r * bsin;
    return a;
}

hgpu_complex hgpu_exp(hgpu_complex c) {
    hgpu_complex a;
    hgpu_float bcos,bsin;
    hgpu_float r = native_exp((hgpu_single) hgpu_Re(c));
    bsin = hgpu_sincos(hgpu_Im(c),&bcos);
    a.re = r * bcos;
    a.im = r * bsin;
    return a;
}

hgpu_complex hgpu_sin(hgpu_complex c) {
    hgpu_complex a;
    hgpu_float bcos,bsin;
    hgpu_float d=hgpu_Im(c);
    bsin = hgpu_sincos(hgpu_Re(c),&bcos);
    a.re = cosh(d)*bsin;
    a.im = sinh(d)*bcos;
    return a;
}

hgpu_complex hgpu_cos(hgpu_complex c) {
    hgpu_complex a;
    hgpu_float bcos,bsin;
    hgpu_float d=hgpu_Im(c);
    bsin = hgpu_sincos(hgpu_Re(c),&bcos);
    a.re = cosh(d)*bcos;
    a.im = sinh(d)*bsin;
    return a;
}

hgpu_complex hgpu_minus(hgpu_complex c) {
    hgpu_complex a;
    a.re = -hgpu_Re(c);
    a.im = -hgpu_Im(c);
    return a;
}


/**
* Get a sum of two complex numbers with optional single/double precision
* @param a, b -- complex numbers to be summed
*/
hgpu_complex hgpu_add(hgpu_complex a,hgpu_complex b) {
    hgpu_complex c;
    c.re = a.re + b.re;
    c.im = a.im + b.im;
    return c;
}

/**
* Get a difference of two complex numbers with optional single/double precision
* @param a -- minuend
* @param b -- subtrahend
*/

hgpu_complex hgpu_sub(hgpu_complex a,hgpu_complex b) {
    hgpu_complex c;
    c.re = a.re - b.re;
    c.im = a.im - b.im;
    return c;
}

/**
* Get a product of two complex numbers with optional single/double precision
* @param a, b -- complex numbers to be multiplied
*/

hgpu_complex hgpu_mul(hgpu_complex a,hgpu_complex b) {
    hgpu_complex c;
    c.re = a.re * b.re - a.im * b.im;
    c.im = a.re * b.im + a.im * b.re;
    return c;
}

hgpu_complex hgpu_div(hgpu_complex a,hgpu_complex b) {
    hgpu_complex c;
    hgpu_float d = b.re * b.re + b.im * b.im;
    c.re = (a.re * b.re + a.im * b.im)/d;
    c.im = (a.im * b.re - a.re * b.im)/d;
    return c;
}
                    
hgpu_complex hgpu_float_to_complex(hgpu_float a) {
    hgpu_complex c;
    c.re = a;
    c.im = 0.0;
    return c;
}
                    
hgpu_complex hgpu_int_to_complex(int a) {
    hgpu_complex c;
    c.re = (hgpu_float) a;
    c.im = 0.0;
    return c;
}

hgpu_complex hgpu_I() {
    hgpu_complex c;
    c.re = 0.0;
    c.im = 1.0;
    return c;
}

//______________________________________________ double precision subroutine _____

hgpu_complex_double hgpu_conjugate_double(hgpu_complex_double c) {
    hgpu_complex_double a;
        a.re =  c.re;
        a.im = -c.im;
    return a;
}

hgpu_complex_double hgpu_mul_double(hgpu_complex_double a,hgpu_complex_double b) {
    hgpu_complex_double c;
    c.re = a.re * b.re - a.im * b.im;
    c.im = a.re * b.im + a.im * b.re;
    return c;
}

hgpu_complex_double hgpu_add_double(hgpu_complex_double a,hgpu_complex_double b) {
    hgpu_complex_double c;
    c.re = a.re + b.re;
    c.im = a.im + b.im;
    return c;
}

hgpu_complex_double hgpu_sub_double(hgpu_complex_double a,hgpu_complex_double b) {
    hgpu_complex_double c;
    c.re = a.re - b.re;
    c.im = a.im - b.im;
    return c;
}
#endif
                                                                                                                                          
                                                                                                                                          
                                                                                                                                          
                                                                                                                                          
                                                                                                                                          
                                                                                                                                          
                                                                                                                                          
                                                                                                                                          
                                                                                                                                          
                                                                                                                                          
                                                                                                                                          
                                                                                                                                          
                                                                                                                                          
                                                                                                                                           
                                                                                                                                           
