/******************************************************************************
 * @file     polyakov.cl
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.4
 *
 * @brief    [QCDGPU]
 *           Measurement of the Polyakov loop
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
#ifndef POLYAKOV_CL
#define POLYAKOV_CL

#include "complex.h"
#include "model.cl"
#include "misc.cl"
#if SUN == 2
#include "su2cl.cl"
#include "su2_matrix_memory.cl"
#include "su2_measurements_cl.cl"
#endif
#if SUN == 3
#include "su3cl.cl"
#include "su3_matrix_memory.cl"
#include "su3_measurements_cl.cl"
#endif

                 __kernel void
lattice_polyakov(__global hgpu_float4  * lattice_table,
                 __global hgpu_double2 * lattice_measurement,
                 __global hgpu_float   * lattice_parameters,
                 __local hgpu_double2  * lattice_lds,
                 uint offset)
{
#if SUN == 2
    hgpu_double2 out  = (hgpu_double2) 0.0;
    hgpu_double2 out2 = (hgpu_double2) 0.0;
    hgpu_double2 out3 = (hgpu_double2) 0.0;
    hgpu_double2 out4 = (hgpu_double2) 0.0;
    hgpu_double polyakov_loop_re = 0.0;
    hgpu_double polyakov_loop_im = 0.0;
    hgpu_double polyakov_loop_p2 = 0.0;
    hgpu_double polyakov_loop_p4 = 0.0;

    uint gindex;
    uint gdi = GID;             //gdi = y + z * N2 + x * N2N3 + t * N1N2N3
    lattice_gid_to_gid_xyz(&gdi,&gindex); // <-- to use function lattice_gid_to_coords(&gindex, ...)

    gpu_su_2 m0,m1;
    su_2 v0,v1,v2;
    coords_4 coord;
    coords_4 coord10;
    uint gdiT;
    su2_twist twist;
    twist.phi   = lattice_parameters[1];

    lattice_lds[TID] = (hgpu_double2) 0.0;
    if (GID<N1N2N3) {
        lattice_gid_to_coords(&gindex,&coord);

        m0 = lattice_table_2(lattice_table,&coord,gindex,T,&twist);       // [p,T]
        v0 = lattice_reconstruct2(&m0);
        for (int i = 1; i < N4; i++){
            lattice_neighbours_gid(&coord,&coord10,&gdiT,T);
            m1 = lattice_table_2(lattice_table,&coord10,gdiT,T,&twist);   // [p,T]
            v1 = lattice_reconstruct2(&m1);

            v2 = matrix_times_su2(&v0,&v1);
            v0 = v2;
            coord = coord10;
        }
        polyakov_loop_re = matrix_retrace_su2(&v0);
        polyakov_loop_im = matrix_imtrace_su2(&v0);
        out = (hgpu_double2) (polyakov_loop_re,polyakov_loop_im);
        polyakov_loop_p2 = polyakov_loop_re * polyakov_loop_re + polyakov_loop_im * polyakov_loop_im;
        polyakov_loop_p4 = polyakov_loop_p2 * polyakov_loop_p2;
        out3 = (hgpu_double2) (polyakov_loop_p2,polyakov_loop_p4);
    }

    // first reduction
    reduce_first_step_val_double2(lattice_lds,&out, &out2);

#if (PL > 1)
     reduce_first_step_val_double2(lattice_lds,&out3,&out4);
#endif

    if(TID == 0) {
        lattice_measurement[BID] = out2;
#if (PL > 1)
        lattice_measurement[BID + offset] = out4;
#endif
    }

#endif

#if SUN == 3
    hgpu_double2 out  = (hgpu_double2) 0.0;
    hgpu_double2 out2 = (hgpu_double2) 0.0;
    hgpu_double2 out3 = (hgpu_double2) 0.0;
    hgpu_double2 out4 = (hgpu_double2) 0.0;
    hgpu_double polyakov_loop_re = 0.0;
    hgpu_double polyakov_loop_im = 0.0;
    hgpu_double polyakov_loop_p2 = 0.0;
    hgpu_double polyakov_loop_p4 = 0.0;

    uint gindex;
    uint gdi = GID;
    lattice_gid_to_gid_xyz(&gdi,&gindex);

    gpu_su_3 m0,m1;
    su_3 v0,v1,v2;
    coords_4 coord;
    coords_4 coord10;
    uint gdiT;
    su3_twist twist;
    twist.phi   = lattice_parameters[1];
    twist.omega = lattice_parameters[2];

    lattice_lds[TID] = (hgpu_double2) 0.0;
    if (GID<N1N2N3) {
        lattice_gid_to_coords(&gindex,&coord);

        m0 = lattice_table_3(lattice_table,&coord,gindex,T,&twist);       // [p,T]
        v0 = lattice_reconstruct3(&m0);
        for (int i = 1; i < N4; i++){
            lattice_neighbours_gid(&coord,&coord10,&gdiT,T);
            m1 = lattice_table_3(lattice_table,&coord10,gdiT,T,&twist);   // [p,T]
            v1 = lattice_reconstruct3(&m1);

            v2 = matrix_times_su3(&v0,&v1);
            v0 = v2;
            coord = coord10;
        }
        polyakov_loop_re = matrix_retrace_su3(&v0);
        polyakov_loop_im = matrix_imtrace_su3(&v0);
        out = (hgpu_double2) (polyakov_loop_re,polyakov_loop_im);
        polyakov_loop_p2 = polyakov_loop_re * polyakov_loop_re + polyakov_loop_im * polyakov_loop_im;
        polyakov_loop_p4 = polyakov_loop_p2 * polyakov_loop_p2;
        out3 = (hgpu_double2) (polyakov_loop_p2,polyakov_loop_p4);
    }

    // first reduction
    reduce_first_step_val_double2(lattice_lds,&out, &out2);

#if (PL > 1)
     reduce_first_step_val_double2(lattice_lds,&out3,&out4);
#endif

    if(TID == 0) {
        lattice_measurement[BID] = out2;
#if (PL > 1)
        lattice_measurement[BID + offset] = out4;
#endif
    }

#endif
}

                        __kernel void
reduce_polyakov_double2(__global hgpu_double2 * lattice_measurement,
                        __global hgpu_double2 * lattice_polyakov_loop,
                        __local  hgpu_double2 * lattice_lds,
                        uint4 param,
                        uint index)
{
    hgpu_double2 out  = (hgpu_double2) 0.0;
    hgpu_double2 out2 = (hgpu_double2) 0.0;
    uint size    = param.x;

    reduce_final_step_double2(lattice_lds,lattice_measurement,size);
    out = lattice_lds[TID];

#if (PL == 2)
    uint offset  = param.y;
    uint offset2 = param.z;
    reduce_final_step_double2_offset(lattice_lds,lattice_measurement,size,offset);
    out2 = lattice_lds[TID];
#endif

    if (GID==0) {
        lattice_polyakov_loop[index] = out;
#if (PL == 2)
        lattice_polyakov_loop[index + offset2] = out2;
#endif
    }

}

#endif
