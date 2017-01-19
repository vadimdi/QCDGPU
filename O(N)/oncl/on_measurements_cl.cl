/******************************************************************************
 * @file     on_measurements.cl
 * @author   Vadim Demchik <vadimdi@yahoo.com>
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.0
 *
 * @brief    [QCDGPU]
 *           Definition of general functions used in measurements, corresponding to the O(N) group
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

#ifndef ONMEASUREMENTSCL_CL
#define ONMEASUREMENTSCL_CL

#include "complex.h"
#include "model.cl"
#include "misc.cl"
#include "neighbours.cl"
#include "o1cl.cl"
#include "o1_matrix_memory.cl"

                                        __kernel void
lattice_measurement(__global hgpu_float   * lattice_table,
                    __global hgpu_double2 * lattice_measurement,
                    __global const hgpu_float   * lattice_parameters,
                    __local hgpu_double2  * lattice_lds)
{
    hgpu_double2 out  = (hgpu_double2) 0.0;
    hgpu_double2 out2 = (hgpu_double2) 0.0;
#if ON_MODEL == 1
    uint gindex = GID;

    lattice_lds[TID] = (hgpu_double2) 0.0; barrier(CLK_LOCAL_MEM_FENCE);
    if (GID<SITESEXACT) {
        coords_4 coord;
        coords_4 coordX,coordY,coordZ,coordT;
        uint gdiX,gdiY,gdiZ,gdiT;

        gpu_o_1     Ux,Ux1,Ux2,Ux3,Ux4;
        hgpu_float  Uxmu0;
        hgpu_float4 Uxmu;
#ifdef ON_SUB_SCHEME
        hgpu_float4 Uxmu_minus;
#endif

        lattice_gid_to_coords(&gindex,&coord);

        // prepare neighbours
        lattice_neighbours_gid(&coord,&coordX,&gdiX,X);
        lattice_neighbours_gid(&coord,&coordY,&gdiY,Y);
        lattice_neighbours_gid(&coord,&coordZ,&gdiZ,Z);
        lattice_neighbours_gid(&coord,&coordT,&gdiT,T);

            Ux  = lattice_table_o_1(lattice_table,gindex);// [p]
            Ux1 = lattice_table_o_1(lattice_table,gdiX);  // [p+X]
            Ux2 = lattice_table_o_1(lattice_table,gdiY);  // [p+Y]
            Ux3 = lattice_table_o_1(lattice_table,gdiZ);  // [p+Z]
            Ux4 = lattice_table_o_1(lattice_table,gdiT);  // [p+T]
        
        Uxmu0 = Ux.uv1;
        Uxmu  = (hgpu_float4) (Ux1.uv1,Ux2.uv1,Ux3.uv1,Ux4.uv1);

        hgpu_float S;
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
        o1_action(&S,&Uxmu0,&Uxmu,&Uxmu_minus,lattice_parameters); // for symmetric difference
#else
        o1_action(&S,&Uxmu0,&Uxmu,lattice_parameters);
#endif
        out.x = S;
    }
    reduce_first_step_val_double2(lattice_lds,&out, &out2);
    if(TID == 0) lattice_measurement[BID] = out2;

#endif
}
                                        __kernel void
reduce_measurement_double2(__global hgpu_double2 * lattice_measurement,
                           __global hgpu_double2 * lattice_energies,
                           __local hgpu_double2  * lattice_lds,
                           uint size,
                           uint index)
{
    reduce_final_step_double2(lattice_lds,lattice_measurement,size);
    hgpu_double2 out = lattice_lds[TID];
    if (GID==0) lattice_energies[index] += out;
}


                                        __kernel void
lattice_measurement_plq(__global hgpu_float   * lattice_table,
                        __global hgpu_double2 * lattice_measurement,
                        __global hgpu_float   * lattice_parameters,
                        __local hgpu_double2  * lattice_lds)
{
#if ON_MODEL == 1
    hgpu_double2 out  = (hgpu_double2) 0.0;
    hgpu_double2 out2 = (hgpu_double2) 0.0;
    uint gindex = GID;
    hgpu_float eta    = lattice_parameters[19];
    hgpu_float b      = lattice_parameters[20];


    lattice_lds[TID] = (hgpu_double2) 0.0; barrier(CLK_LOCAL_MEM_FENCE);
    if (GID<SITESEXACT) {
        coords_4 coord;

        gpu_o_1     Ux;

        lattice_gid_to_coords(&gindex,&coord);
            Ux  = lattice_table_o_1(lattice_table,gindex);// [p]
                
        hgpu_float varnce = o1_to_physical_field(&Ux,&b,&eta);

        out.x += varnce;
        out.y += varnce*varnce;

    }
    reduce_first_step_val_double2(lattice_lds,&out, &out2);
    if(TID == 0) lattice_measurement[BID] = out2;

#endif
}

                               __kernel void
reduce_measurement_plq_double2(__global hgpu_double2 * lattice_measurement,
                               __global hgpu_double2 * lattice_energies_plq,
                               __local  hgpu_double2 * lattice_lds,
                               uint size,
                               uint index)
{
    reduce_final_step_double2(lattice_lds,lattice_measurement,size);
    if (GID==0) lattice_energies_plq[index] += lattice_lds[0];
}

                                        __kernel void
lattice_measurement_correlator(__global hgpu_float   * lattice_table,
                               __global hgpu_double2 * lattice_measurement,
                               __global hgpu_float   * lattice_parameters,
                               __local  hgpu_double2 * lattice_lds,
                                        uint4          lattice_stepz)
{
#if ON_MODEL == 1
    hgpu_double2 out  = (hgpu_double2) 0.0;
    hgpu_double2 out2 = (hgpu_double2) 0.0;
    uint gindex = GID;
    hgpu_float eta    = lattice_parameters[19];
    hgpu_float b      = lattice_parameters[20];


    lattice_lds[TID] = (hgpu_double2) 0.0; barrier(CLK_LOCAL_MEM_FENCE);
    if (GID<SITESEXACT) {
        coords_4 coord,coord2/*,coord3*/,coord_new;

        gpu_o_1     Ux,Ux2,Ux3;
        uint        gdi2,gdi3;

        lattice_gid_to_coords(&gindex,&coord);
            Ux  = lattice_table_o_1(lattice_table,gindex);

        // prepare neighbours
        // correlator1 - step (lattice_stepz)
        coord2.x = lattice_stepz.x;
        coord2.y = lattice_stepz.y;
        coord2.z = lattice_stepz.z;
        coord2.t = lattice_stepz.w;

        lattice_neighbours_step2_gid(&coord,&coord_new,&gdi2,&coord2);
            Ux2 = lattice_table_o_1(lattice_table,gdi2);
    
        out.x = o1_correlator(&Ux,&Ux2,&b,&eta);

        // correlator2 - diagonal (+1,+1,+1,+1)
        lattice_neighbours_diagonal_gid(&coord,&coord_new,&gdi3);
            Ux3 = lattice_table_o_1(lattice_table,gdi3);

        out.y = o1_correlator(&Ux,&Ux3,&b,&eta);
    }
    reduce_first_step_val_double2(lattice_lds,&out, &out2);
    if(TID == 0) lattice_measurement[BID] = out2;

#endif
}

                                        __kernel void
clear_measurement(__global hgpu_double2 * lattice_measurement)
{
    lattice_measurement[GID] = (hgpu_double2) 0.0;
}


#endif
