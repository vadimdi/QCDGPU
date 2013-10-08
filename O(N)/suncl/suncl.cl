/******************************************************************************
 * @file     suncl.cl
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.0
 *
 * @brief    [QCDGPU]
 *           Contains lattice initialization and update procedures, matrix reunitarization Gram-Schmidt procedure
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

#ifndef SUNCL_CL
#define SUNCL_CL

#include "complex.h"
#include "model.cl"
#include "misc.cl"
#include "su3cl.cl"
#include "su3_matrix_memory.cl"
#include "su3_update_cl.cl"

                                        __kernel void
lattice_init_hot_X(__global hgpu_float4 * lattice_table,
                   __global const hgpu_single4 * prns)
{

#if SUN == 3
    uint gidprn1 = GID + 0 * ROWSIZE;
    uint gidprn2 = GID + 1 * ROWSIZE;
    uint gidprn3 = GID + 2 * ROWSIZE;
    gpu_su_3 matrix;
    if (GID < SITESEXACT) {
        lattice_random3(&matrix,prns,gidprn1,gidprn2,gidprn3);
        lattice_store_3(lattice_table,&matrix,GID,X);
    }
#endif
}

                                        __kernel void
lattice_init_hot_Y(__global hgpu_float4 * lattice_table,
                   __global const hgpu_single4 * prns)
{

#if SUN == 3
    uint gidprn1 = GID + 0 * ROWSIZE;
    uint gidprn2 = GID + 1 * ROWSIZE;
    uint gidprn3 = GID + 2 * ROWSIZE;
    gpu_su_3 matrix;
    if (GID < SITESEXACT) {
        lattice_random3(&matrix,prns,gidprn1,gidprn2,gidprn3);
        lattice_store_3(lattice_table,&matrix,GID,Y);
    }
#endif
}

                                        __kernel void
lattice_init_hot_Z(__global hgpu_float4 * lattice_table,
                   __global const hgpu_single4 * prns)
{

#if SUN == 3
    uint gidprn1 = GID + 0 * ROWSIZE;
    uint gidprn2 = GID + 1 * ROWSIZE;
    uint gidprn3 = GID + 2 * ROWSIZE;
    gpu_su_3 matrix;
    if (GID < SITESEXACT) {
        lattice_random3(&matrix,prns,gidprn1,gidprn2,gidprn3);
        lattice_store_3(lattice_table,&matrix,GID,Z);
    }
#endif
}

                                        __kernel void
lattice_init_hot_T(__global hgpu_float4 * lattice_table,
                   __global const hgpu_single4 * prns)
{

#if SUN == 3
    uint gidprn1 = GID + 0 * ROWSIZE;
    uint gidprn2 = GID + 1 * ROWSIZE;
    uint gidprn3 = GID + 2 * ROWSIZE;
    gpu_su_3 matrix;
    if (GID < SITESEXACT) {
        lattice_random3(&matrix,prns,gidprn1,gidprn2,gidprn3);
        lattice_store_3(lattice_table,&matrix,GID,T);
    }
#endif
}

                                        __kernel void
lattice_init_cold(__global hgpu_float4 * lattice_table)
{
#if SUN == 3
    gpu_su_3 matrix;
    lattice_unity3(&matrix);
    if (GID < SITESEXACT) {
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
#if SUN == 3
    gpu_su_3 matrix;
    matrix.uv1 = (hgpu_float4) (GID, 0.0, 0.0, 0.0);
    matrix.uv2 = (hgpu_float4) (0.0, 0.0, 0.0, 0.0);
    matrix.uv3 = (hgpu_float4) (0.0, 1.0/GID, 0.0, 0.0);
    lattice_store_3(lattice_table,&matrix,GID,X);
    lattice_store_3(lattice_table,&matrix,GID,Y);
    lattice_store_3(lattice_table,&matrix,GID,Z);
    lattice_store_3(lattice_table,&matrix,GID,T);
#endif
}

                                        __kernel void
lattice_GramSchmidt(__global hgpu_float4 * lattice_table,
                    __global hgpu_float *  lattice_parameters)
{
#if SUN == 3
    gpu_su_3 matrix;
    coords_4 coord;
    uint gindex = GID;
    lattice_gid_to_coords(&gindex,&coord);

    su3_twist twist;
    twist.phi   = lattice_parameters[1];
    twist.omega = lattice_parameters[2];

    if (GID < SITESEXACT) {
        matrix = lattice_table_3(lattice_table,&coord,GID,X,&twist);
        lattice_GramSchmidt3(&matrix);
        lattice_store_3(lattice_table,&matrix,GID,X);
        matrix = lattice_table_3(lattice_table,&coord,GID,Y,&twist);
        lattice_GramSchmidt3(&matrix);
        lattice_store_3(lattice_table,&matrix,GID,Y);
        matrix = lattice_table_3(lattice_table,&coord,GID,Z,&twist);
        lattice_GramSchmidt3(&matrix);
        lattice_store_3(lattice_table,&matrix,GID,Z);
        matrix = lattice_table_3(lattice_table,&coord,GID,T,&twist);
        lattice_GramSchmidt3(&matrix);
        lattice_store_3(lattice_table,&matrix,GID,T);
    }
#endif
}

                                        __kernel void
update_even_X(__global hgpu_float4 * lattice_table,
              __global hgpu_float * lattice_parameters,
              __global const hgpu_single4 * prns
#ifdef ACC_RATE
             ,__local hgpu_double2  * lattice_lds
             ,__global hgpu_double2 * lattice_measurement
#endif
              )
{
#ifdef ACC_RATE
    hgpu_double out  = 0.0;
    hgpu_double out2 = 0.0;
    lattice_lds[TID]  = (hgpu_double2) 0.0; barrier(CLK_LOCAL_MEM_FENCE);
#endif
    coords_4 coord;
    uint gindex = lattice_even_gid();  // x_+/-_y,z,t
    hgpu_float bet      = lattice_parameters[0];
    su3_twist twist;
    twist.phi   = lattice_parameters[1];
    twist.omega = lattice_parameters[2];

    if (GID < SITESHALFEXACT) {
#ifdef ACC_RATE
        out  = (double) (NHIT+NHIT+NHIT);
#endif
        gpu_su_3 m0,mU;
        su_3 staple;
        lattice_gid_to_coords(&gindex,&coord);

        m0     = lattice_table_3(lattice_table,&coord,gindex,X,&twist);    // [p,X]
        staple = lattice_staple_3(lattice_table,gindex,X,&twist);
#ifdef ACC_RATE
        mU     = lattice_heatbath3(&staple,&m0,&bet,prns,&out);
#else
        mU     = lattice_heatbath3(&staple,&m0,&bet,prns);
#endif

#ifndef BULK_UPDATES
        lattice_store_3(lattice_table,&mU,gindex,X);    // update lattice
#endif
    }

#ifdef ACC_RATE
    reduce_first_step_val_double(lattice_lds,&out, &out2);

    if(TID == 0) lattice_measurement[BID].x += out2;
#endif
}
                                        __kernel void
update_even_Y(__global hgpu_float4 * lattice_table,
              __global hgpu_float * lattice_parameters,
              __global const hgpu_single4 * prns
#ifdef ACC_RATE
             ,__local hgpu_double2  * lattice_lds
             ,__global hgpu_double2 * lattice_measurement
#endif
              )
{
#ifdef ACC_RATE
    hgpu_double out  = 0.0;
    hgpu_double out2 = 0.0;
    lattice_lds[TID]  = (hgpu_double2) 0.0; barrier(CLK_LOCAL_MEM_FENCE);
#endif
    coords_4 coord;
    uint gindex = lattice_even_gid();  // y_+/-_x,z,t
    hgpu_float bet   = lattice_parameters[0];
    su3_twist twist;
    twist.phi   = lattice_parameters[1];
    twist.omega = lattice_parameters[2];

    if (GID < SITESHALFEXACT) {
#ifdef ACC_RATE
        out  = (double) (NHIT+NHIT+NHIT);
#endif
        gpu_su_3 m0,mU;
        su_3 staple;
        lattice_gid_to_coords(&gindex,&coord);

        m0     = lattice_table_3(lattice_table,&coord,gindex,Y,&twist);    // [p,Y]
        staple = lattice_staple_3(lattice_table,gindex,Y,&twist);
#ifdef ACC_RATE
        mU     = lattice_heatbath3(&staple,&m0,&bet,prns,&out);
#else
        mU     = lattice_heatbath3(&staple,&m0,&bet,prns);
#endif

#ifndef BULK_UPDATES
        lattice_store_3(lattice_table,&mU,gindex,Y);    // update lattice
#endif
    }

#ifdef ACC_RATE
    reduce_first_step_val_double(lattice_lds,&out, &out2);

    if(TID == 0) lattice_measurement[BID].x += out2;
#endif
}
                                        __kernel void
update_even_Z(__global hgpu_float4 * lattice_table,
              __global hgpu_float * lattice_parameters,
              __global const hgpu_single4 * prns
#ifdef ACC_RATE
             ,__local hgpu_double2  * lattice_lds
             ,__global hgpu_double2 * lattice_measurement
#endif
              )
{
#ifdef ACC_RATE
    hgpu_double out  = 0.0;
    hgpu_double out2 = 0.0;
    lattice_lds[TID]  = (hgpu_double2) 0.0; barrier(CLK_LOCAL_MEM_FENCE);
#endif
    coords_4 coord;
    uint gindex = lattice_even_gid();  // z_+/-_x,y,t
    hgpu_float bet   = lattice_parameters[0];
    su3_twist twist;
    twist.phi   = lattice_parameters[1];
    twist.omega = lattice_parameters[2];

    if (GID < SITESHALFEXACT) {
#ifdef ACC_RATE
        out  = (double) (NHIT+NHIT+NHIT);
#endif
        gpu_su_3 m0,mU;
        su_3 staple;
        lattice_gid_to_coords(&gindex,&coord);

        m0     = lattice_table_3(lattice_table,&coord,gindex,Z,&twist);    // [p,Z]
        staple = lattice_staple_3(lattice_table,gindex,Z,&twist);
#ifdef ACC_RATE
        mU     = lattice_heatbath3(&staple,&m0,&bet,prns,&out);
#else
        mU     = lattice_heatbath3(&staple,&m0,&bet,prns);
#endif

#ifndef BULK_UPDATES
        lattice_store_3(lattice_table,&mU,gindex,Z);    // update lattice
#endif
    }

#ifdef ACC_RATE
    reduce_first_step_val_double(lattice_lds,&out, &out2);

    if(TID == 0) lattice_measurement[BID].x += out2;
#endif
}
                                        __kernel void
update_even_T(__global hgpu_float4 * lattice_table,
              __global hgpu_float * lattice_parameters,
              __global const hgpu_single4 * prns
#ifdef ACC_RATE
             ,__local hgpu_double2  * lattice_lds
             ,__global hgpu_double2 * lattice_measurement
#endif
              )
{
#ifdef ACC_RATE
    hgpu_double out  = 0.0;
    hgpu_double out2 = 0.0;
    lattice_lds[TID]  = (hgpu_double2) 0.0; barrier(CLK_LOCAL_MEM_FENCE);
#endif
    coords_4 coord;
    uint gindex = lattice_even_gid();  // t_+/-_x,y,z
    hgpu_float bet   = lattice_parameters[0];
    su3_twist twist;
    twist.phi   = lattice_parameters[1];
    twist.omega = lattice_parameters[2];

    if (GID < SITESHALFEXACT) {
#ifdef ACC_RATE
        out  = (double) (NHIT+NHIT+NHIT);
#endif
        gpu_su_3 m0,mU;
        su_3 staple;
        lattice_gid_to_coords(&gindex,&coord);

        m0     = lattice_table_3(lattice_table,&coord,gindex,T,&twist);    // [p,T]
        staple = lattice_staple_3(lattice_table,gindex,T,&twist);
#ifdef ACC_RATE
        mU     = lattice_heatbath3(&staple,&m0,&bet,prns,&out);
#else
        mU     = lattice_heatbath3(&staple,&m0,&bet,prns);
#endif

#ifndef BULK_UPDATES
        lattice_store_3(lattice_table,&mU,gindex,T);    // update lattice
#endif
    }

#ifdef ACC_RATE
    reduce_first_step_val_double(lattice_lds,&out, &out2);

    if(TID == 0) lattice_measurement[BID].x += out2;
#endif
}
                                        __kernel void
update_odd_X(__global hgpu_float4 * lattice_table,
             __global hgpu_float * lattice_parameters,
             __global const hgpu_single4 * prns
#ifdef ACC_RATE
             ,__local hgpu_double2  * lattice_lds
             ,__global hgpu_double2 * lattice_measurement
#endif
             )
{
#ifdef ACC_RATE
    hgpu_double out  = 0.0;
    hgpu_double out2 = 0.0;
    lattice_lds[TID]  = (hgpu_double2) 0.0; barrier(CLK_LOCAL_MEM_FENCE);
#endif
    coords_4 coord;
    uint gindex = lattice_odd_gid();  // x_+/-_y,z,t
    hgpu_float bet   = lattice_parameters[0];
    su3_twist twist;
    twist.phi   = lattice_parameters[1];
    twist.omega = lattice_parameters[2];

    if (GID < SITESHALFEXACT) {
#ifdef ACC_RATE
        out  = (double) (NHIT+NHIT+NHIT);
#endif
        gpu_su_3 m0,mU;
        su_3 staple;
        lattice_gid_to_coords(&gindex,&coord);

        m0     = lattice_table_3(lattice_table,&coord,gindex,X,&twist);    // [p,X]
        staple = lattice_staple_3(lattice_table,gindex,X,&twist);
#ifdef ACC_RATE
        mU     = lattice_heatbath3(&staple,&m0,&bet,prns,&out);
#else
        mU     = lattice_heatbath3(&staple,&m0,&bet,prns);
#endif

#ifndef BULK_UPDATES
        lattice_store_3(lattice_table,&mU,gindex,X);    // update lattice
#endif
    }

#ifdef ACC_RATE
    reduce_first_step_val_double(lattice_lds,&out, &out2);

    if(TID == 0) lattice_measurement[BID].y += out2;
#endif
}
                                        __kernel void
update_odd_Y(__global hgpu_float4 * lattice_table,
             __global hgpu_float * lattice_parameters,
             __global const hgpu_single4 * prns
#ifdef ACC_RATE
             ,__local hgpu_double2  * lattice_lds
             ,__global hgpu_double2 * lattice_measurement
#endif
             )
{
#ifdef ACC_RATE
    hgpu_double out  = 0.0;
    hgpu_double out2 = 0.0;
    lattice_lds[TID]  = (hgpu_double2) 0.0; barrier(CLK_LOCAL_MEM_FENCE);
#endif
    coords_4 coord;
    uint gindex = lattice_odd_gid();  // y_+/-_x,z,t
    hgpu_float bet   = lattice_parameters[0];
    su3_twist twist;
    twist.phi   = lattice_parameters[1];
    twist.omega = lattice_parameters[2];

    if (GID < SITESHALFEXACT) {
#ifdef ACC_RATE
        out  = (double) (NHIT+NHIT+NHIT);
#endif
        gpu_su_3 m0,mU;
        su_3 staple;
        lattice_gid_to_coords(&gindex,&coord);

        m0     = lattice_table_3(lattice_table,&coord,gindex,Y,&twist);    // [p,Y]
        staple = lattice_staple_3(lattice_table,gindex,Y,&twist);
#ifdef ACC_RATE
        mU     = lattice_heatbath3(&staple,&m0,&bet,prns,&out);
#else
        mU     = lattice_heatbath3(&staple,&m0,&bet,prns);
#endif

#ifndef BULK_UPDATES
        lattice_store_3(lattice_table,&mU,gindex,Y);    // update lattice
#endif
    }

#ifdef ACC_RATE
    reduce_first_step_val_double(lattice_lds,&out, &out2);

    if(TID == 0) lattice_measurement[BID].y += out2;
#endif
}
                                        __kernel void
update_odd_Z(__global hgpu_float4 * lattice_table,
             __global hgpu_float * lattice_parameters,
             __global const hgpu_single4 * prns
#ifdef ACC_RATE
             ,__local hgpu_double2  * lattice_lds
             ,__global hgpu_double2 * lattice_measurement
#endif
             )
{
#ifdef ACC_RATE
    hgpu_double out  = 0.0;
    hgpu_double out2 = 0.0;
    lattice_lds[TID]  = (hgpu_double2) 0.0; barrier(CLK_LOCAL_MEM_FENCE);
#endif
    coords_4 coord;
    uint gindex = lattice_odd_gid();  // z_+/-_x,y,t
    hgpu_float bet   = lattice_parameters[0];
    su3_twist twist;
    twist.phi   = lattice_parameters[1];
    twist.omega = lattice_parameters[2];

    if (GID < SITESHALFEXACT) {
#ifdef ACC_RATE
        out  = (double) (NHIT+NHIT+NHIT);
#endif
        gpu_su_3 m0,mU;
        su_3 staple;
        lattice_gid_to_coords(&gindex,&coord);

        m0     = lattice_table_3(lattice_table,&coord,gindex,Z,&twist);    // [p,Z]
        staple = lattice_staple_3(lattice_table,gindex,Z,&twist);
#ifdef ACC_RATE
        mU     = lattice_heatbath3(&staple,&m0,&bet,prns,&out);
#else
        mU     = lattice_heatbath3(&staple,&m0,&bet,prns);
#endif

#ifndef BULK_UPDATES
        lattice_store_3(lattice_table,&mU,gindex,Z);    // update lattice
#endif
    }

#ifdef ACC_RATE
    reduce_first_step_val_double(lattice_lds,&out, &out2);

    if(TID == 0) lattice_measurement[BID].y += out2;
#endif
}
                                        __kernel void
update_odd_T(__global hgpu_float4 * lattice_table,
             __global hgpu_float * lattice_parameters,
             __global const hgpu_single4 * prns
#ifdef ACC_RATE
             ,__local hgpu_double2  * lattice_lds
             ,__global hgpu_double2 * lattice_measurement
#endif
             )
{
#ifdef ACC_RATE
    hgpu_double out  = 0.0;
    hgpu_double out2 = 0.0;
    lattice_lds[TID]  = (hgpu_double2) 0.0; barrier(CLK_LOCAL_MEM_FENCE);
#endif
    coords_4 coord;
    uint gindex = lattice_odd_gid();  // t_+/-_x,y,z
    hgpu_float bet   = lattice_parameters[0];
    su3_twist twist;
    twist.phi   = lattice_parameters[1];
    twist.omega = lattice_parameters[2];

    if (GID < SITESHALFEXACT) {
#ifdef ACC_RATE
        out  = (double) (NHIT+NHIT+NHIT);
#endif
        gpu_su_3 m0,mU;
        su_3 staple;
        lattice_gid_to_coords(&gindex,&coord);

        m0     = lattice_table_3(lattice_table,&coord,gindex,T,&twist);    // [p,T]
        staple = lattice_staple_3(lattice_table,gindex,T,&twist);
#ifdef ACC_RATE
        mU     = lattice_heatbath3(&staple,&m0,&bet,prns,&out);
#else
        mU     = lattice_heatbath3(&staple,&m0,&bet,prns);
#endif

#ifndef BULK_UPDATES
        lattice_store_3(lattice_table,&mU,gindex,T);    // update lattice
#endif
    }

#ifdef ACC_RATE
    reduce_first_step_val_double(lattice_lds,&out, &out2);

    if(TID == 0) lattice_measurement[BID].y += out2;
#endif
}


#ifdef ACC_RATE
                                        __kernel void
reduce_acceptance_rate(__global hgpu_double2 * lattice_measurement,
                       __global hgpu_double2 * lattice_acceptance_rate,
                       __local  hgpu_double2 * lattice_lds,
                           uint size,
                           uint index)
{
    reduce_final_step_double2(lattice_lds,lattice_measurement,size);
    hgpu_double2 out = lattice_lds[TID];
    if (GID==0) lattice_acceptance_rate[index] += out;
}
#endif

#endif
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
