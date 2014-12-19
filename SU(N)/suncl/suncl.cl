/******************************************************************************
 * @file     suncl.cl
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.5
 *
 * @brief    [QCDGPU]
 *           Contains lattice initialization and update procedures, matrix reunitarization Gram-Schmidt procedure
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
#ifndef SUNCL_CL
#define SUNCL_CL
 
#include "complex.h"
#include "model.cl"
#include "misc.cl"
#if SUN == 2
#include "su2cl.cl"
#include "su2_matrix_memory.cl"
#include "su2_update_cl.cl"
#endif
#if SUN == 3
#include "su3cl.cl"
#include "su3_matrix_memory.cl"
#include "su3_update_cl.cl"
#endif

                                        __kernel void
lattice_init_hot_X(__global hgpu_float4 * lattice_table,
                   __global const hgpu_prng_float4 * prns)
{

#if SUN == 2
    uint gidprn1 = GID;
    
    gpu_su_2 matrix;
    if (GID < SITES) {
        lattice_random2(&matrix,prns,gidprn1);
        lattice_store_2(lattice_table,&matrix,GID,X);
    }
#endif

#if SUN == 3
    uint gidprn1 = GID + 0 * ROWSIZE;
    uint gidprn2 = GID + 1 * ROWSIZE;
    uint gidprn3 = GID + 2 * ROWSIZE;
    gpu_su_3 matrix;
    if (GID < SITES) {
        lattice_random3(&matrix,prns,gidprn1,gidprn2,gidprn3);
        lattice_store_3(lattice_table,&matrix,GID,X);
    }
#endif
}

                                        __kernel void
lattice_init_hot_Y(__global hgpu_float4 * lattice_table,
                   __global const hgpu_prng_float4 * prns)
{

#if SUN == 2
    uint gidprn1 = GID;
    
    gpu_su_2 matrix;
    if (GID < SITES) {
        lattice_random2(&matrix,prns,gidprn1);
        lattice_store_2(lattice_table,&matrix,GID,Y);
    }
#endif

#if SUN == 3
    uint gidprn1 = GID + 0 * ROWSIZE;
    uint gidprn2 = GID + 1 * ROWSIZE;
    uint gidprn3 = GID + 2 * ROWSIZE;
    gpu_su_3 matrix;
    if (GID < SITES) {
        lattice_random3(&matrix,prns,gidprn1,gidprn2,gidprn3);
        lattice_store_3(lattice_table,&matrix,GID,Y);
    }
#endif
}

                                        __kernel void
lattice_init_hot_Z(__global hgpu_float4 * lattice_table,
                   __global const hgpu_prng_float4 * prns)
{

#if SUN == 2
    uint gidprn1 = GID;
    
    gpu_su_2 matrix;
    if (GID < SITES) {
        lattice_random2(&matrix,prns,gidprn1);
        lattice_store_2(lattice_table,&matrix,GID,Z);
    }
#endif

#if SUN == 3
    uint gidprn1 = GID + 0 * ROWSIZE;
    uint gidprn2 = GID + 1 * ROWSIZE;
    uint gidprn3 = GID + 2 * ROWSIZE;
    gpu_su_3 matrix;
    if (GID < SITES) {
        lattice_random3(&matrix,prns,gidprn1,gidprn2,gidprn3);
        lattice_store_3(lattice_table,&matrix,GID,Z);
    }
#endif
}

                                        __kernel void
lattice_init_hot_T(__global hgpu_float4 * lattice_table,
                   __global const hgpu_prng_float4 * prns)
{

#if SUN == 2
    uint gidprn1 = GID;
    
    gpu_su_2 matrix;
    if (GID < SITES) {
        lattice_random2(&matrix,prns,gidprn1);
        lattice_store_2(lattice_table,&matrix,GID,T);
    }
#endif

#if SUN == 3
    uint gidprn1 = GID + 0 * ROWSIZE;
    uint gidprn2 = GID + 1 * ROWSIZE;
    uint gidprn3 = GID + 2 * ROWSIZE;
    gpu_su_3 matrix;
    if (GID < SITES) {
        lattice_random3(&matrix,prns,gidprn1,gidprn2,gidprn3);
        lattice_store_3(lattice_table,&matrix,GID,T);
    }
#endif
}

                                        __kernel void
lattice_init_cold(__global hgpu_float4 * lattice_table)
{
#if SUN == 2
    gpu_su_2 matrix;
    lattice_unity2(&matrix);
    if (GID < SITES) {
        lattice_store_2(lattice_table,&matrix,GID,X);
        lattice_store_2(lattice_table,&matrix,GID,Y);
        lattice_store_2(lattice_table,&matrix,GID,Z);
        lattice_store_2(lattice_table,&matrix,GID,T);
    }
#endif

#if SUN == 3
    gpu_su_3 matrix;
    lattice_unity3(&matrix);
    if (GID < SITES) {
        lattice_store_3(lattice_table,&matrix,GID,X);
        lattice_store_3(lattice_table,&matrix,GID,Y);
        lattice_store_3(lattice_table,&matrix,GID,Z);
        lattice_store_3(lattice_table,&matrix,GID,T);
    }
#endif
}

                                        __kernel void
lattice_init_gid(__global hgpu_float4 * lattice_table)
{
#if SUN == 2
    gpu_su_2 matrix, matrix_y, matrix_z, matrix_t;

  if(GID < SITES)
  {
    hgpu_double f1 = sin((0.001+10.0/SITES)*GID)/2.0;
    hgpu_double f2 = cos((0.001-10.0/SITES)*GID)/3.0;
    hgpu_double f3 = sin((0.001-18.0/SITES)*GID)/4.0;
    hgpu_double f4 = sqrt(1.0-f1*f1-f2*f2-f3*f3);

    matrix.uv1 = (hgpu_float4) (f1, f3, f4, f2);

    hgpu_double f1y = sin((0.002+10.0/SITES)*GID)/2.0;
    hgpu_double f2y = cos((0.002-10.0/SITES)*GID)/3.0;
    hgpu_double f3y = sin((0.002-18.0/SITES)*GID)/4.0;

    hgpu_double f4y = sqrt(1.0-f1y*f1y-f2y*f2y-f3y*f3y);

    matrix_y.uv1 = (hgpu_float4) (f1y, f3y, f4y, f2y);

    hgpu_double f1z = sin((0.003+10.0/SITES)*GID)/2.0;
    hgpu_double f2z = cos((0.003-10.0/SITES)*GID)/3.0;
    hgpu_double f3z = sin((0.003-18.0/SITES)*GID)/4.0;
    
    hgpu_double f4z = sqrt(1.0-f1z*f1z-f2z*f2z-f3z*f3z);
    matrix_z.uv1 = (hgpu_float4) (f1z, f3z, f4z, f2z);

    hgpu_double f1t = sin((0.004+10.0/SITES)*GID)/2.0;
    hgpu_double f2t = cos((0.004-10.0/SITES)*GID)/3.0;
    hgpu_double f3t = sin((0.004-18.0/SITES)*GID)/4.0;

    hgpu_double f4t = sqrt(1.0-f1t*f1t-f2t*f2t-f3t*f3t);

    matrix_t.uv1 = (hgpu_float4) (f1t, f3t, f4t, f2t);
    
    lattice_store_2(lattice_table,&matrix,GID,X);
    lattice_store_2(lattice_table,&matrix_y,GID,Y);
    lattice_store_2(lattice_table,&matrix_z,GID,Z);
    lattice_store_2(lattice_table,&matrix_t,GID,T);
  }
#endif

#if SUN == 3
    gpu_su_3 matrix, matrix_y, matrix_z, matrix_t;

  if(GID < SITES)
  {
    hgpu_double f1  = sin((0.001+10.0/SITES)*GID)/2.0;
    hgpu_double f2  = cos((0.001-10.0/SITES)*GID)/3.0;
    hgpu_double f3  = sin((0.001-18.0/SITES)*GID)/4.0;
    hgpu_double f4  = cos((0.001+12.0/SITES)*GID)/4.0;
    hgpu_double f5  = sin((0.001-12.0/SITES)*GID)/3.0;
    hgpu_double f6  = cos((0.001-15.0/SITES)*GID)/2.0;
    hgpu_double f7  = sin((0.001+15.0/SITES)*GID)/3.0;
    hgpu_double f8  = cos((0.001-11.0/SITES)*GID)/4.0;
    hgpu_double f9  = sin((0.001-11.0/SITES)*GID)/3.0;
    hgpu_double f10 = cos((0.001+16.0/SITES)*GID)/2.0;
    hgpu_double f11 = sin((0.001-16.0/SITES)*GID)/3.0;
    hgpu_double f12 = cos((0.001-18.0/SITES)*GID)/2.0;
    

    matrix.uv1 = (hgpu_float4) (f1, f3, f5, f11);
    matrix.uv2 = (hgpu_float4) (f2, f4, f6, f12);
    matrix.uv3 = (hgpu_float4) (f7, f9, f8, f10);
     lattice_GramSchmidt3(&matrix);

    hgpu_double f1y  = sin((0.002+10.0/SITES)*GID)/2.0;
    hgpu_double f2y  = cos((0.002-10.0/SITES)*GID)/3.0;
    hgpu_double f3y  = sin((0.002-18.0/SITES)*GID)/4.0;
    hgpu_double f4y  = cos((0.002+12.0/SITES)*GID)/4.0;
    hgpu_double f5y  = sin((0.002-12.0/SITES)*GID)/3.0;
    hgpu_double f6y  = cos((0.002-15.0/SITES)*GID)/2.0;
    hgpu_double f7y  = sin((0.002+15.0/SITES)*GID)/3.0;
    hgpu_double f8y  = cos((0.002-11.0/SITES)*GID)/4.0;
    hgpu_double f9y  = sin((0.002-11.0/SITES)*GID)/3.0;
    hgpu_double f10y = cos((0.002+16.0/SITES)*GID)/2.0;
    hgpu_double f11y = sin((0.002-16.0/SITES)*GID)/3.0;
    hgpu_double f12y = cos((0.002-18.0/SITES)*GID)/2.0;
    

    matrix_y.uv1 = (hgpu_float4) (f1y, f3y, f5y, f11y);
    matrix_y.uv2 = (hgpu_float4) (f2y, f4y, f6y, f12y);
    matrix_y.uv3 = (hgpu_float4) (f7y, f9y, f8y, f10y);
     lattice_GramSchmidt3(&matrix_y);

    hgpu_double f1z  = sin((0.003+10.0/SITES)*GID)/2.0;
    hgpu_double f2z  = cos((0.003-10.0/SITES)*GID)/3.0;
    hgpu_double f3z  = sin((0.003-18.0/SITES)*GID)/4.0;
    hgpu_double f4z  = cos((0.003+12.0/SITES)*GID)/4.0;
    hgpu_double f5z  = sin((0.003-12.0/SITES)*GID)/3.0;
    hgpu_double f6z  = cos((0.003-15.0/SITES)*GID)/2.0;
    hgpu_double f7z  = sin((0.003+15.0/SITES)*GID)/3.0;
    hgpu_double f8z  = cos((0.003-11.0/SITES)*GID)/4.0;
    hgpu_double f9z  = sin((0.003-11.0/SITES)*GID)/3.0;
    hgpu_double f10z = cos((0.003+16.0/SITES)*GID)/2.0;
    hgpu_double f11z = sin((0.003-16.0/SITES)*GID)/3.0;
    hgpu_double f12z = cos((0.003-18.0/SITES)*GID)/2.0;
    

    matrix_z.uv1 = (hgpu_float4) (f1z, f3z, f5z, f11z);
    matrix_z.uv2 = (hgpu_float4) (f2z, f4z, f6z, f12z);
    matrix_z.uv3 = (hgpu_float4) (f7z, f9z, f8z, f10z);
     lattice_GramSchmidt3(&matrix_z);

    hgpu_double f1t  = sin((0.004+10.0/SITES)*GID)/2.0;
    hgpu_double f2t  = cos((0.004-10.0/SITES)*GID)/3.0;
    hgpu_double f3t  = sin((0.004-18.0/SITES)*GID)/4.0;
    hgpu_double f4t  = cos((0.004+12.0/SITES)*GID)/4.0;
    hgpu_double f5t  = sin((0.004-12.0/SITES)*GID)/3.0;
    hgpu_double f6t  = cos((0.004-15.0/SITES)*GID)/2.0;
    hgpu_double f7t  = sin((0.004+15.0/SITES)*GID)/3.0;
    hgpu_double f8t  = cos((0.004-11.0/SITES)*GID)/4.0;
    hgpu_double f9t  = sin((0.004-11.0/SITES)*GID)/3.0;
    hgpu_double f10t = cos((0.004+16.0/SITES)*GID)/2.0;
    hgpu_double f11t = sin((0.004-16.0/SITES)*GID)/3.0;
    hgpu_double f12t = cos((0.004-18.0/SITES)*GID)/2.0;
    

    matrix_t.uv1 = (hgpu_float4) (f1t, f3t, f5t, f11t);
    matrix_t.uv2 = (hgpu_float4) (f2t, f4t, f6t, f12t);
    matrix_t.uv3 = (hgpu_float4) (f7t, f9t, f8t, f10t);
     lattice_GramSchmidt3(&matrix_t);
    
    lattice_store_3(lattice_table,&matrix,GID,X);
    lattice_store_3(lattice_table,&matrix_y,GID,Y);
    lattice_store_3(lattice_table,&matrix_z,GID,Z);
    lattice_store_3(lattice_table,&matrix_t,GID,T);
  }
#endif
}

                                        __kernel void
lattice_GramSchmidt(__global hgpu_float4 * lattice_table,
                    __global hgpu_float *  lattice_parameters)
{
#if SUN == 2
    gpu_su_2 matrix;
    coords_4 coord;
    uint gindex = GID;
    lattice_gid_to_coords(&gindex,&coord);

    if (GID < SITES) {
        matrix = lattice_table_notwist_2(lattice_table,GID,X);
        lattice_su2_Normalize(&matrix);
        lattice_store_2(lattice_table,&matrix,GID,X);
        matrix = lattice_table_notwist_2(lattice_table,GID,Y);
        lattice_su2_Normalize(&matrix);
        lattice_store_2(lattice_table,&matrix,GID,Y);
        matrix = lattice_table_notwist_2(lattice_table,GID,Z);
        lattice_su2_Normalize(&matrix);
        lattice_store_2(lattice_table,&matrix,GID,Z);
        matrix = lattice_table_notwist_2(lattice_table,GID,T);
        lattice_su2_Normalize(&matrix);
        lattice_store_2(lattice_table,&matrix,GID,T);
    }
#endif

#if SUN == 3
    gpu_su_3 matrix;
    coords_4 coord;
    uint gindex = GID;
    lattice_gid_to_coords(&gindex,&coord);

    if (GID < SITES) {
        matrix = lattice_table_notwist_3(lattice_table,GID,X);
        lattice_GramSchmidt3(&matrix);
        lattice_store_3(lattice_table,&matrix,GID,X);
        matrix = lattice_table_notwist_3(lattice_table,GID,Y);
        lattice_GramSchmidt3(&matrix);
        lattice_store_3(lattice_table,&matrix,GID,Y);
        matrix = lattice_table_notwist_3(lattice_table,GID,Z);
        lattice_GramSchmidt3(&matrix);
        lattice_store_3(lattice_table,&matrix,GID,Z);
        matrix = lattice_table_notwist_3(lattice_table,GID,T);
        lattice_GramSchmidt3(&matrix);
        lattice_store_3(lattice_table,&matrix,GID,T);
    }
#endif
}

                                        __kernel void
update_even_X(__global hgpu_float4 * lattice_table,
              __global hgpu_float * lattice_parameters,
              __global const hgpu_prng_float4 * prns)
{
    coords_4 coord;
    uint gindex = lattice_even_gid();  // x_+/-_y,z,t
    hgpu_float bet      = lattice_parameters[0];
#if SUN == 2
    su2_twist twist;
    twist.phi   = lattice_parameters[1];

    if (GID < SITESHALF) {
        gpu_su_2 m0,mU;
        su_2 staple;
        lattice_gid_to_coords(&gindex,&coord);

        m0     = lattice_table_2(lattice_table,&coord,gindex,X,&twist);    // [p,X]
        staple = lattice_staple_2(lattice_table,gindex,X,&twist);
#ifdef GID_UPD
prns[GID].x = (float) gindex;
#endif
        mU     = lattice_heatbath_2(&staple,&m0,&bet,prns);

#ifndef BULK_UPDATES
           lattice_store_2(lattice_table,&mU,gindex,X);    // update lattice
#endif
#endif

#if SUN == 3
    su3_twist twist;
    twist.phi   = lattice_parameters[1];
    twist.omega = lattice_parameters[2];

    if (GID < SITESHALF) {
        gpu_su_3 m0,mU;
        su_3 staple;
        lattice_gid_to_coords(&gindex,&coord);

        m0     = lattice_table_3(lattice_table,&coord,gindex,X,&twist);    // [p,X]
        staple = lattice_staple_3(lattice_table,gindex,X,&twist);
#ifdef GID_UPD
prns[GID].x = (float) gindex;
#endif
        mU     = lattice_heatbath3(&staple,&m0,&bet,prns);

#ifndef BULK_UPDATES
           lattice_store_3(lattice_table,&mU,gindex,X);    // update lattice
#endif
#endif
    }
}
                                        __kernel void
update_even_Y(__global hgpu_float4 * lattice_table,
              __global hgpu_float * lattice_parameters,
              __global const hgpu_prng_float4 * prns)
{
    coords_4 coord;
    uint gindex = lattice_even_gid();  // y_+/-_x,z,t
    hgpu_float bet   = lattice_parameters[0];
#if SUN == 2
    su2_twist twist;
    twist.phi   = lattice_parameters[1];

    if (GID < SITESHALF) {
        gpu_su_2 m0,mU;
        su_2 staple;
        lattice_gid_to_coords(&gindex,&coord);

        m0     = lattice_table_2(lattice_table,&coord,gindex,Y,&twist);    // [p,Y]
        staple = lattice_staple_2(lattice_table,gindex,Y,&twist);
#ifdef GID_UPD
prns[GID].x = (float) gindex;
#endif
        mU     = lattice_heatbath_2(&staple,&m0,&bet,prns);

#ifndef BULK_UPDATES
           lattice_store_2(lattice_table,&mU,gindex,Y);    // update lattice
#endif
#endif

#if SUN == 3
    su3_twist twist;
    twist.phi   = lattice_parameters[1];
    twist.omega = lattice_parameters[2];

    if (GID < SITESHALF) {
        gpu_su_3 m0,mU;
        su_3 staple;
        lattice_gid_to_coords(&gindex,&coord);

        m0     = lattice_table_3(lattice_table,&coord,gindex,Y,&twist);    // [p,Y]
        staple = lattice_staple_3(lattice_table,gindex,Y,&twist);
#ifdef GID_UPD
prns[GID].x = (float) gindex;
#endif
        mU     = lattice_heatbath3(&staple,&m0,&bet,prns);

#ifndef BULK_UPDATES
           lattice_store_3(lattice_table,&mU,gindex,Y);    // update lattice
#endif
#endif
    }
}
                                        __kernel void
update_even_Z(__global hgpu_float4 * lattice_table,
              __global hgpu_float * lattice_parameters,
              __global const hgpu_prng_float4 * prns)
{
    coords_4 coord;
    uint gindex = lattice_even_gid();  // z_+/-_x,y,t
    hgpu_float bet   = lattice_parameters[0];
#if SUN == 2
    su2_twist twist;
    twist.phi   = lattice_parameters[1];

    if (GID < SITESHALF) {
        gpu_su_2 m0,mU;
        su_2 staple;
        lattice_gid_to_coords(&gindex,&coord);

        m0     = lattice_table_2(lattice_table,&coord,gindex,Z,&twist);    // [p,Z]
        staple = lattice_staple_2(lattice_table,gindex,Z,&twist);
#ifdef GID_UPD
prns[GID].x = (float) gindex;
#endif
        mU     = lattice_heatbath_2(&staple,&m0,&bet,prns);

#ifndef BULK_UPDATES
           lattice_store_2(lattice_table,&mU,gindex,Z);    // update lattice
#endif
#endif

#if SUN == 3
    su3_twist twist;
    twist.phi   = lattice_parameters[1];
    twist.omega = lattice_parameters[2];

    if (GID < SITESHALF) {
        gpu_su_3 m0,mU;
        su_3 staple;
        lattice_gid_to_coords(&gindex,&coord);

        m0     = lattice_table_3(lattice_table,&coord,gindex,Z,&twist);    // [p,Z]
        staple = lattice_staple_3(lattice_table,gindex,Z,&twist);
#ifdef GID_UPD
prns[GID].x = (float) gindex;
#endif
        mU     = lattice_heatbath3(&staple,&m0,&bet,prns);

#ifndef BULK_UPDATES
           lattice_store_3(lattice_table,&mU,gindex,Z);    // update lattice
#endif
#endif
    }
}
                                        __kernel void
update_even_T(__global hgpu_float4 * lattice_table,
              __global hgpu_float * lattice_parameters,
              __global const hgpu_prng_float4 * prns)
{
    coords_4 coord;
    uint gindex = lattice_even_gid();  // t_+/-_x,y,z
    hgpu_float bet   = lattice_parameters[0];
#if SUN == 2
    su2_twist twist;
    twist.phi   = lattice_parameters[1];

    if (GID < SITESHALF) {
        gpu_su_2 m0,mU;
        su_2 staple;
        lattice_gid_to_coords(&gindex,&coord);

        m0     = lattice_table_2(lattice_table,&coord,gindex,T,&twist);    // [p,T]
        staple = lattice_staple_2(lattice_table,gindex,T,&twist);
#ifdef GID_UPD
prns[GID].x = (float) gindex;
#endif
        mU     = lattice_heatbath_2(&staple,&m0,&bet,prns);

#ifndef BULK_UPDATES
           lattice_store_2(lattice_table,&mU,gindex,T);    // update lattice
#endif
#endif

#if SUN == 3
    su3_twist twist;
    twist.phi   = lattice_parameters[1];
    twist.omega = lattice_parameters[2];

    if (GID < SITESHALF) {
        gpu_su_3 m0,mU;
        su_3 staple;
        lattice_gid_to_coords(&gindex,&coord);

        m0     = lattice_table_3(lattice_table,&coord,gindex,T,&twist);    // [p,T]
        staple = lattice_staple_3(lattice_table,gindex,T,&twist);
#ifdef GID_UPD
prns[GID].x = (float) gindex;
#endif
        mU     = lattice_heatbath3(&staple,&m0,&bet,prns);

#ifndef BULK_UPDATES
           lattice_store_3(lattice_table,&mU,gindex,T);    // update lattice
#endif
#endif
    }
}
                                        __kernel void
update_odd_X(__global hgpu_float4 * lattice_table,
             __global hgpu_float * lattice_parameters,
             __global const hgpu_prng_float4 * prns)
{
    coords_4 coord;
    uint gindex = lattice_odd_gid();  // x_+/-_y,z,t
    hgpu_float bet   = lattice_parameters[0];
#if SUN == 2
    su2_twist twist;
    twist.phi   = lattice_parameters[1];

    if (GID < SITESHALF) {
        gpu_su_2 m0,mU;
        su_2 staple;
        lattice_gid_to_coords(&gindex,&coord);

        m0     = lattice_table_2(lattice_table,&coord,gindex,X,&twist);    // [p,X]
        staple = lattice_staple_2(lattice_table,gindex,X,&twist);
#ifdef GID_UPD
prns[GID].x = (float) gindex;
#endif
        mU     = lattice_heatbath_2(&staple,&m0,&bet,prns);

#ifndef BULK_UPDATES
           lattice_store_2(lattice_table,&mU,gindex,X);    // update lattice
#endif
#endif

#if SUN == 3
    su3_twist twist;
    twist.phi   = lattice_parameters[1];
    twist.omega = lattice_parameters[2];

    if (GID < SITESHALF) {
        gpu_su_3 m0,mU;
        su_3 staple;
        lattice_gid_to_coords(&gindex,&coord);

        m0     = lattice_table_3(lattice_table,&coord,gindex,X,&twist);    // [p,X]
        staple = lattice_staple_3(lattice_table,gindex,X,&twist);
#ifdef GID_UPD
prns[GID].x = (float) gindex;
#endif
        mU     = lattice_heatbath3(&staple,&m0,&bet,prns);

#ifndef BULK_UPDATES
           lattice_store_3(lattice_table,&mU,gindex,X);    // update lattice
#endif
#endif
    }
}
                                        __kernel void
update_odd_Y(__global hgpu_float4 * lattice_table,
             __global hgpu_float * lattice_parameters,
             __global const hgpu_prng_float4 * prns)
{
    coords_4 coord;
    uint gindex = lattice_odd_gid();  // y_+/-_x,z,t
    hgpu_float bet   = lattice_parameters[0];
#if SUN == 2
    su2_twist twist;
    twist.phi   = lattice_parameters[1];

    if (GID < SITESHALF) {
        gpu_su_2 m0,mU;
        su_2 staple;
        lattice_gid_to_coords(&gindex,&coord);

        m0     = lattice_table_2(lattice_table,&coord,gindex,Y,&twist);    // [p,Y]
        staple = lattice_staple_2(lattice_table,gindex,Y,&twist);
#ifdef GID_UPD
prns[GID].x = (float) gindex;
#endif
        mU     = lattice_heatbath_2(&staple,&m0,&bet,prns);

#ifndef BULK_UPDATES
           lattice_store_2(lattice_table,&mU,gindex,Y);    // update lattice
#endif
#endif

#if SUN == 3
    su3_twist twist;
    twist.phi   = lattice_parameters[1];
    twist.omega = lattice_parameters[2];

    if (GID < SITESHALF) {
        gpu_su_3 m0,mU;
        su_3 staple;
        lattice_gid_to_coords(&gindex,&coord);

        m0     = lattice_table_3(lattice_table,&coord,gindex,Y,&twist);    // [p,Y]
        staple = lattice_staple_3(lattice_table,gindex,Y,&twist);
#ifdef GID_UPD
prns[GID].x = (float) gindex;
#endif
        mU     = lattice_heatbath3(&staple,&m0,&bet,prns);

#ifndef BULK_UPDATES
           lattice_store_3(lattice_table,&mU,gindex,Y);    // update lattice
#endif
#endif
    }
}
                                        __kernel void
update_odd_Z(__global hgpu_float4 * lattice_table,
             __global hgpu_float * lattice_parameters,
             __global const hgpu_prng_float4 * prns)
{
    coords_4 coord;
    uint gindex = lattice_odd_gid();  // z_+/-_x,y,t
    hgpu_float bet   = lattice_parameters[0];
#if SUN == 2
    su2_twist twist;
    twist.phi   = lattice_parameters[1];

    if (GID < SITESHALF) {
        gpu_su_2 m0,mU;
        su_2 staple;
        lattice_gid_to_coords(&gindex,&coord);

        m0     = lattice_table_2(lattice_table,&coord,gindex,Z,&twist);    // [p,Z]
        staple = lattice_staple_2(lattice_table,gindex,Z,&twist);
#ifdef GID_UPD
prns[GID].x = (float) gindex;
#endif
        mU     = lattice_heatbath_2(&staple,&m0,&bet,prns);

#ifndef BULK_UPDATES
           lattice_store_2(lattice_table,&mU,gindex,Z);    // update lattice
#endif
#endif

#if SUN == 3
    su3_twist twist;
    twist.phi   = lattice_parameters[1];
    twist.omega = lattice_parameters[2];

    if (GID < SITESHALF) {
        gpu_su_3 m0,mU;
        su_3 staple;
        lattice_gid_to_coords(&gindex,&coord);

        m0     = lattice_table_3(lattice_table,&coord,gindex,Z,&twist);    // [p,Z]
        staple = lattice_staple_3(lattice_table,gindex,Z,&twist);
#ifdef GID_UPD
prns[GID].x = (float) gindex;
#endif
        mU     = lattice_heatbath3(&staple,&m0,&bet,prns);

#ifndef BULK_UPDATES
           lattice_store_3(lattice_table,&mU,gindex,Z);    // update lattice
#endif
#endif

    }
}
                                        __kernel void
update_odd_T(__global hgpu_float4 * lattice_table,
             __global hgpu_float * lattice_parameters,
             __global const hgpu_prng_float4 * prns)
{
    coords_4 coord;
    uint gindex = lattice_odd_gid();  // t_+/-_x,y,z
    hgpu_float bet   = lattice_parameters[0];
#if SUN == 2
    su2_twist twist;
    twist.phi   = lattice_parameters[1];

    if (GID < SITESHALF) {
        gpu_su_2 m0,mU;
        su_2 staple;
        lattice_gid_to_coords(&gindex,&coord);

        m0     = lattice_table_2(lattice_table,&coord,gindex,T,&twist);    // [p,T]
        staple = lattice_staple_2(lattice_table,gindex,T,&twist);
#ifdef GID_UPD
prns[GID].x = (float) gindex;
#endif
        mU     = lattice_heatbath_2(&staple,&m0,&bet,prns);

#ifndef BULK_UPDATES
           lattice_store_2(lattice_table,&mU,gindex,T);    // update lattice
#endif
#endif

#if SUN == 3
    su3_twist twist;
    twist.phi   = lattice_parameters[1];
    twist.omega = lattice_parameters[2];

    if (GID < SITESHALF) {
        gpu_su_3 m0,mU;
        su_3 staple;
        lattice_gid_to_coords(&gindex,&coord);

        m0     = lattice_table_3(lattice_table,&coord,gindex,T,&twist);    // [p,T]
        staple = lattice_staple_3(lattice_table,gindex,T,&twist);
#ifdef GID_UPD
prns[GID].x = (float) gindex;
#endif
        mU     = lattice_heatbath3(&staple,&m0,&bet,prns);

#ifndef BULK_UPDATES
           lattice_store_3(lattice_table,&mU,gindex,T);    // update lattice
#endif
#endif
    }
}

#endif
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                  
                                                                                                                                                                  
