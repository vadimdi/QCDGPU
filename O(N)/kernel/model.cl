/******************************************************************************
 * @file     model.cl
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.0
 *
 * @brief    [QCDGPU]
 *           Description of lattice geometry, gauge groups, etc.
 *
 * @section  LICENSE
 *
 * Copyright (c) 2013, Vadim Demchik, Natalia Kolomoyets
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

#ifndef MODEL_CL
#define MODEL_CL

#if defined(cl_amd_fp64)  // AMD extension available?
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
#elif defined(cl_khr_fp64)  // Khronos extension available?
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

//#define ON_SUB_SCHEME   // for symmetric scheme for differences

// ________________ additional functions
#define MIN(a,b)    (((a) < (b)) ? (a) : (b))
#define MAX(a,b)    (((a) > (b)) ? (a) : (b))
#define CEIL(a)     ((a - (int)a)==0 ? (int)a : (int)a+1)

#define GID_SIZE    (get_global_size(0) * get_global_size(1) * get_global_size(2))
#define GID         (get_global_id(0) + get_global_id(1) * get_global_size(0) + get_global_id(2) * get_global_size(0) * get_global_size(1))
#define TID         (get_local_id(0))
#define BID         (get_group_id(0))
#define GROUP_SIZE  (get_local_size(0))

#ifndef N1N2
#define N1N2        (N1 * N2)
#endif
#ifndef N2N3
#define N2N3        (N2 * N3)
#endif
#ifndef N2N3N4
#define N2N3N4      (N2 * N3 * N4)
#endif
#ifndef N1N2N3
#define N1N2N3      (N1 * N2 * N3)
#endif
#ifndef N1EXACTN2N3
#define N1EXACTN2N3 (N1EXACT * N2 * N3)
#endif
#ifndef N1N2N3N4
#define N1N2N3N4    (N1 * N2 * N3 * N4)
#endif

#ifdef ON_MODEL
#define GROUPELEMENTS   1
#else
#if SUN==3
#define GROUPELEMENTS   12
#elif SUN==2
#define GROUPELEMENTS   4
#elif SUN==1
#define GROUPELEMENTS   1
#endif
#endif

#ifndef SKGROUP
#define SKGROUP         4
#endif
#define SKGROUPM1_2     (0.5*(SKGROUP-1))

#define X   0
#define Y   1
#define Z   2
#define T   3

#define SITES           (N1 * N2 * N3 * N4)
#define SITESEXACT      (N1EXACT * N2 * N3 * N4)
#define LINKS           (SITES * ND)
#define SITESHALF       (SITES / 2)
#define SITESHALFEXACT  (SITESEXACT / 2)

#ifndef ROWSIZE
#define ROWSIZE (SITES) // (LINKS / 4)
#endif

#ifndef ROWLINKS
#define ROWLINKS (ROWSIZE * ND)
#endif

#ifndef PLSIZE
#define PLSIZE  (N1N2N3)
#endif

#define EVEN    0
#define ODD     1


#ifndef PI
#define PI    3.1415926535897932384626433832795    // pi
#endif
#ifndef PI2
#define PI2   6.2831853071795864769252867665590    // 2*pi
#endif

#ifndef NHIT
#define NHIT    10
#endif

#ifndef NHIT_2
#define NHIT_2      (NHIT>>1)   // NHIT/2
#endif

#ifndef PRNGSTEP
#define PRNGSTEP    (SITES / 2)
#endif

typedef union _Float_and_Double{         // Float <---> Double converter
                float       flVal[2];
                hgpu_double dbVal;
} Float_and_Double;


#ifdef PRECISION_DOUBLE
typedef struct{
    hgpu_double uv1;
} gpu_o_1;

typedef struct{
    hgpu_double4 uv1;
} gpu_su_2;

typedef struct{
    hgpu_double4 uv1;
    hgpu_double4 uv2;
    hgpu_double4 uv3;
} gpu_su_3;
#else
typedef struct{
    float uv1;
} gpu_o_1;

typedef struct{
    float4 uv1;
} gpu_su_2;

typedef struct{
    float4 uv1;
    float4 uv2;
    float4 uv3;
} gpu_su_3;
#endif


typedef struct{
    hgpu_complex u1;
    hgpu_complex u2;
    hgpu_complex v1;
    hgpu_complex v2;
} su_2;

typedef struct {
    hgpu_complex u1;
    hgpu_complex u2;
    hgpu_complex u3;
    hgpu_complex v1;
    hgpu_complex v2;
    hgpu_complex v3;
    hgpu_complex w1;
    hgpu_complex w2;
    hgpu_complex w3;
} su_3;

typedef struct {
    hgpu_complex_double u1;
    hgpu_complex_double u2;
    hgpu_complex_double u3;
    hgpu_complex_double v1;
    hgpu_complex_double v2;
    hgpu_complex_double v3;
    hgpu_complex_double w1;
    hgpu_complex_double w2;
    hgpu_complex_double w3;
} double_su_3;


typedef struct{
    uint x;
    uint y;
    uint z;
    uint t;
    uint xy;
    uint xyz;
} coords_4;

typedef struct{
    hgpu_float phi;
    hgpu_float omega;
} su3_twist;

#endif
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
