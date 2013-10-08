/******************************************************************************
 * @file     oncl.cl
 * @author   Vadim Demchik <vadimdi@yahoo.com>
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

#ifndef ONCL_CL
#define ONCL_CL

#include "complex.h"
#include "model.cl"
#include "misc.cl"
#include "neighbours.cl"
#include "o1cl.cl"
#include "o1_matrix_memory.cl"
#include "o1_update_cl.cl"

                                       __kernel void
lattice_init_hot(__global hgpu_float * lattice_table,
                 __global const hgpu_single4 * prns,
                 __global const hgpu_float * lattice_parameters
                 )
{
    uint gidprn1 = GID;
    hgpu_float max_U  = lattice_parameters[22];
    gpu_o_1 matrix;
    if (GID < SITESEXACT) {
        lattice_random_o_1(&matrix,prns,gidprn1,&max_U);
        lattice_store_o_1(lattice_table,&matrix,gidprn1);
    }
}

                                        __kernel void
lattice_init_cold(__global hgpu_float * lattice_table)
{
    gpu_o_1 matrix;
    lattice_ground_o_1(&matrix);
    if (GID < SITESEXACT) {
        lattice_store_o_1(lattice_table,&matrix,GID);
    }
}

                                        __kernel void
update_even(__global hgpu_float * lattice_table,
            __global const hgpu_float * lattice_parameters,
            __global const hgpu_single4 * prns
#ifdef ACC_RATE
           ,__local hgpu_double2  * lattice_lds
           ,__global hgpu_double2 * lattice_measurement
#endif
            )
{
    uint gindex = lattice_even_gid();

#ifdef ACC_RATE
    hgpu_double out  = 0.0;
    hgpu_double out2 = 0.0;
    lattice_lds[TID]  = (hgpu_double2) 0.0; barrier(CLK_LOCAL_MEM_FENCE);
#endif
    if (GID < SITESHALFEXACT) {
#ifdef ACC_RATE
        out  = (double) NHIT;
#endif
        coords_4 coord;
        coords_4 coordX,coordY,coordZ,coordT;
        gpu_o_1     Ux,Ux1,Ux2,Ux3,Ux4,Ux_new;
        hgpu_float  Uxmu0,Uxmu_new;
        hgpu_float4 Uxmu;
#ifdef ON_SUB_SCHEME
        hgpu_float4 Uxmu_minus;
#endif
        hgpu_float  action_old,action_new,deltaS;
        uint gdiX,gdiY,gdiZ,gdiT;

        hgpu_float4 rnd;
        uint indprng = GID;

        hgpu_float max_U  = lattice_parameters[22];

        lattice_gid_to_coords(&gindex,&coord);

        // prepare neighbours
        lattice_neighbours_gid(&coord,&coordX,&gdiX,X);
        lattice_neighbours_gid(&coord,&coordY,&gdiY,Y);
        lattice_neighbours_gid(&coord,&coordZ,&gdiZ,Z);
        lattice_neighbours_gid(&coord,&coordT,&gdiT,T);

            Ux  = lattice_table_o_1(lattice_table,gindex);  // [p]
            Ux1 = lattice_table_o_1(lattice_table,gdiX);    // [p+X]
            Ux2 = lattice_table_o_1(lattice_table,gdiY);    // [p+Y]
            Ux3 = lattice_table_o_1(lattice_table,gdiZ);    // [p+Z]
            Ux4 = lattice_table_o_1(lattice_table,gdiT);    // [p+T]
        
        Uxmu0 = Ux.uv1;
        Uxmu  = (hgpu_float4) (Ux1.uv1,Ux2.uv1,Ux3.uv1,Ux4.uv1);

#ifdef ON_SUB_SCHEME
        // minus neihgbours        
        lattice_neighbours_gid_minus(&coord,&coordX,&gdiX,X);
        lattice_neighbours_gid_minus(&coord,&coordY,&gdiY,Y);
        lattice_neighbours_gid_minus(&coord,&coordZ,&gdiZ,Z);
        lattice_neighbours_gid_minus(&coord,&coordT,&gdiT,T);

            Ux1 = lattice_table_o_1(lattice_table,gdiX);    // [p+X]
            Ux2 = lattice_table_o_1(lattice_table,gdiY);    // [p+Y]
            Ux3 = lattice_table_o_1(lattice_table,gdiZ);    // [p+Z]
            Ux4 = lattice_table_o_1(lattice_table,gdiT);    // [p+T]

        Uxmu_minus = (hgpu_float4) (Ux1.uv1,Ux2.uv1,Ux3.uv1,Ux4.uv1);
        o1_action(&action_old,&Uxmu0,&Uxmu,&Uxmu_minus,lattice_parameters); // for symmetric difference
#else
        o1_action(&action_old,&Uxmu0,&Uxmu,lattice_parameters);
#endif

        bool flag = false;
        for (int i=0;i<NHIT_2;i++){

            rnd.x = (hgpu_float) prns[indprng].x;
            rnd.y = (hgpu_float) prns[indprng].y;
            rnd.z = (hgpu_float) prns[indprng].z;
            rnd.w = (hgpu_float) prns[indprng].w;
            indprng += PRNGSTEP;

            Uxmu_new   = rnd.x * max_U;
#ifdef ON_SUB_SCHEME
            o1_action(&action_new,&Uxmu_new,&Uxmu,&Uxmu_minus,lattice_parameters); // for symmetric difference
#else
            o1_action(&action_new,&Uxmu_new,&Uxmu,lattice_parameters);
#endif
            deltaS = exp(action_old-action_new);
            if (rnd.y<deltaS) {
                Ux_new.uv1 = Uxmu_new;
                flag = true;
                break;
            }
#ifdef ACC_RATE
            out -=1.0;
#endif

            Uxmu_new   = rnd.z * max_U;
#ifdef ON_SUB_SCHEME
            o1_action(&action_new,&Uxmu_new,&Uxmu,&Uxmu_minus,lattice_parameters); // for symmetric difference
#else
            o1_action(&action_new,&Uxmu_new,&Uxmu,lattice_parameters);
#endif
            deltaS = exp(action_old-action_new);
            if (rnd.w<deltaS) {
                Ux_new.uv1 = Uxmu_new;
                flag = true;
                break;
            }
#ifdef ACC_RATE
            out -=1.0;
#endif
        }
#if ((NHIT-NHIT_2-NHIT_2)==1)
if (!flag) {
            rnd.x = (hgpu_float) prns[indprng].x;
            rnd.y = (hgpu_float) prns[indprng].y;
            rnd.z = (hgpu_float) prns[indprng].z;
            rnd.w = (hgpu_float) prns[indprng].w;
            indprng += PRNGSTEP;

            Uxmu_new   = rnd.x * max_U;
#ifdef ON_SUB_SCHEME
            o1_action(&action_new,&Uxmu_new,&Uxmu,&Uxmu_minus,lattice_parameters); // for symmetric difference
#else
            o1_action(&action_new,&Uxmu_new,&Uxmu,lattice_parameters);
#endif
            deltaS = exp(action_old-action_new);
            if (rnd.y<deltaS) {
                Ux_new.uv1 = Uxmu_new;
                flag = true;
            }
#ifdef ACC_RATE
            else
                out -=1.0;
#endif
}
#endif

        if (flag) lattice_store_o_1(lattice_table,&Ux_new,gindex);
    }
#ifdef ACC_RATE
    reduce_first_step_val_double(lattice_lds,&out, &out2);
    if(TID == 0) lattice_measurement[BID].x += out2;
#endif
}

                                        __kernel void
update_odd(__global hgpu_float * lattice_table,
           __global const hgpu_float * lattice_parameters,
           __global const hgpu_single4 * prns
#ifdef ACC_RATE
           ,__local hgpu_double2  * lattice_lds
           ,__global hgpu_double2 * lattice_measurement
#endif
           )
{
    uint gindex = lattice_odd_gid();

#ifdef ACC_RATE
    hgpu_double out  = 0.0;
    hgpu_double out2 = 0.0;
    lattice_lds[TID]  = (hgpu_double2) 0.0; barrier(CLK_LOCAL_MEM_FENCE);
#endif
    if (GID < SITESHALFEXACT) {
#ifdef ACC_RATE
        out  = (double) NHIT;
#endif
        coords_4 coord;
        coords_4 coordX,coordY,coordZ,coordT;
        gpu_o_1     Ux,Ux1,Ux2,Ux3,Ux4,Ux_new;
        hgpu_float  Uxmu0,Uxmu_new;
        hgpu_float4 Uxmu;
#ifdef ON_SUB_SCHEME
        hgpu_float4 Uxmu_minus;
#endif
        hgpu_float  action_old,action_new,deltaS;
        uint gdiX,gdiY,gdiZ,gdiT;

        hgpu_float4 rnd;
        uint indprng = GID;

        hgpu_float max_U  = lattice_parameters[22];

        lattice_gid_to_coords(&gindex,&coord);

        // prepare neighbours
        lattice_neighbours_gid(&coord,&coordX,&gdiX,X);
        lattice_neighbours_gid(&coord,&coordY,&gdiY,Y);
        lattice_neighbours_gid(&coord,&coordZ,&gdiZ,Z);
        lattice_neighbours_gid(&coord,&coordT,&gdiT,T);

            Ux  = lattice_table_o_1(lattice_table,gindex);  // [p]
            Ux1 = lattice_table_o_1(lattice_table,gdiX);    // [p+X]
            Ux2 = lattice_table_o_1(lattice_table,gdiY);    // [p+Y]
            Ux3 = lattice_table_o_1(lattice_table,gdiZ);    // [p+Z]
            Ux4 = lattice_table_o_1(lattice_table,gdiT);    // [p+T]
        
        Uxmu0 = Ux.uv1;
        Uxmu  = (hgpu_float4) (Ux1.uv1,Ux2.uv1,Ux3.uv1,Ux4.uv1);

#ifdef ON_SUB_SCHEME
        // minus neihgbours        
        lattice_neighbours_gid_minus(&coord,&coordX,&gdiX,X);
        lattice_neighbours_gid_minus(&coord,&coordY,&gdiY,Y);
        lattice_neighbours_gid_minus(&coord,&coordZ,&gdiZ,Z);
        lattice_neighbours_gid_minus(&coord,&coordT,&gdiT,T);

            Ux1 = lattice_table_o_1(lattice_table,gdiX);    // [p+X]
            Ux2 = lattice_table_o_1(lattice_table,gdiY);    // [p+Y]
            Ux3 = lattice_table_o_1(lattice_table,gdiZ);    // [p+Z]
            Ux4 = lattice_table_o_1(lattice_table,gdiT);    // [p+T]

        Uxmu_minus = (hgpu_float4) (Ux1.uv1,Ux2.uv1,Ux3.uv1,Ux4.uv1);
        o1_action(&action_old,&Uxmu0,&Uxmu,&Uxmu_minus,lattice_parameters); // for symmetric difference
#else
        o1_action(&action_old,&Uxmu0,&Uxmu,lattice_parameters);
#endif

        bool flag = false;
        for (int i=0;i<NHIT_2;i++){
            rnd.x = (hgpu_float) prns[indprng].x;
            rnd.y = (hgpu_float) prns[indprng].y;
            rnd.z = (hgpu_float) prns[indprng].z;
            rnd.w = (hgpu_float) prns[indprng].w;
            indprng += PRNGSTEP;

            Uxmu_new   = rnd.x * max_U;
#ifdef ON_SUB_SCHEME
            o1_action(&action_new,&Uxmu_new,&Uxmu,&Uxmu_minus,lattice_parameters); // for symmetric difference
#else
            o1_action(&action_new,&Uxmu_new,&Uxmu,lattice_parameters);
#endif
            deltaS = exp(action_old-action_new);
            if (rnd.y<deltaS) {
                Ux_new.uv1 = Uxmu_new;
                flag = true;
                break;
            }
#ifdef ACC_RATE
            out -=1.0;
#endif

            Uxmu_new   = rnd.z * max_U;
#ifdef ON_SUB_SCHEME
            o1_action(&action_new,&Uxmu_new,&Uxmu,&Uxmu_minus,lattice_parameters); // for symmetric difference
#else
            o1_action(&action_new,&Uxmu_new,&Uxmu,lattice_parameters);
#endif
            deltaS = exp(action_old-action_new);
            if (rnd.w<deltaS) {
                Ux_new.uv1 = Uxmu_new;
                flag = true;
                break;
            }
#ifdef ACC_RATE
            out -=1.0;
#endif
        }
#if ((NHIT-NHIT_2-NHIT_2)==1)
if (!flag) {
            rnd.x = (hgpu_float) prns[indprng].x;
            rnd.y = (hgpu_float) prns[indprng].y;
            rnd.z = (hgpu_float) prns[indprng].z;
            rnd.w = (hgpu_float) prns[indprng].w;
            indprng += PRNGSTEP;

            Uxmu_new   = rnd.x * max_U;
#ifdef ON_SUB_SCHEME
            o1_action(&action_new,&Uxmu_new,&Uxmu,&Uxmu_minus,lattice_parameters); // for symmetric difference
#else
            o1_action(&action_new,&Uxmu_new,&Uxmu,lattice_parameters);
#endif
            deltaS = exp(action_old-action_new);
            if (rnd.y<deltaS) {
                Ux_new.uv1 = Uxmu_new;
                flag = true;
            }
#ifdef ACC_RATE
            else
                out -=1.0;
#endif
}
#endif

        if (flag) lattice_store_o_1(lattice_table,&Ux_new,gindex);
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
