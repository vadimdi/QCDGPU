/******************************************************************************
 * @file     wilson_loop.cl
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.4
 *
 * @brief    [QCDGPU]
 *           Measurement of the Wilson loop
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
#ifndef WILSON_LOOP_CL
#define WILSON_LOOP_CL

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
lattice_measurement_wilson(__global hgpu_float4  * lattice_table,
                           __global hgpu_double2 * lattice_measurement,
                           __global hgpu_float   * lattice_parameters,
                           __local hgpu_double2  * lattice_lds)
{
#if SUN == 2
    hgpu_double out  = 0.0;
    hgpu_double out2 = 0.0;
    uint gdi = GID;
    hgpu_double wilson_loop;
    coords_4 coord,  coord2;
    coords_4 coordX, coordX2;
    coords_4 coord_1,coord_2,coord_4;
    uint gdi_1,gdi_2,gdi_4,gdiT;
    gpu_su_2 u1;
    double_su_2 m1, m2, m3, m4, m5, w1;

    su2_twist twist;
        twist.phi   = lattice_parameters[1];
    hgpu_float wilson_R,wilson_T,tmp_R,tmp_T;
        wilson_R    = lattice_parameters[3];
        wilson_T    = lattice_parameters[4];
    wilson_loop = 0.0;

    lattice_lds[TID] = (hgpu_double2) 0.0;
    if (GID<SITES) {
        // _______________________________ bottom link
        lattice_gid_to_coords(&gdi,&coord_1);
        gdi_1 = gdi;
        coord_2 = coord_1;
        u1 = lattice_table_2(lattice_table,&coord_2, gdi,T,&twist);         // [p,T]
        m1 = lattice_reconstruct2_double(&u1);
        tmp_T = 1.1;

        while (tmp_T<wilson_T){
            coord = coord_2;
            // prepare neighbours
            lattice_neighbours_gid(&coord,&coord_2,&gdiT,T);
            u1 = lattice_table_2(lattice_table,&coord_2,gdiT,T,&twist);     // [p,T]
            m5 = lattice_reconstruct2_double(&u1);
            w1 = matrix_times_su2_double(&m1,&m5);
            m1 = w1;

            tmp_T += 1.0;
        }
        lattice_neighbours_gid(&coord_2,&coord,&gdi_2,T);
        coord_2 = coord;

        // _______________________________ left and right X links
        // coord_1 - (p,X)
        // coord_2 - (p+T,X)
        u1 = lattice_table_2(lattice_table,&coord_1,gdi_1,X,&twist);        // [p,X]
        m4 = lattice_reconstruct2_double(&u1);
        u1 = lattice_table_2(lattice_table,&coord_2,gdi_2,X,&twist);        // [p+T,X]
        m2 = lattice_reconstruct2_double(&u1);
        tmp_R  = 1.1;
        coord  = coord_1;
        coord2 = coord_2;

        while (tmp_R<wilson_R){
            coordX  = coord;
            coordX2 = coord2;
            // prepare neighbours
            lattice_neighbours_gid(&coordX,  &coord,  &gdi, X);
            lattice_neighbours_gid(&coordX2, &coord2, &gdiT,X);

            u1 = lattice_table_2(lattice_table,&coord, gdi, X,&twist);      // [p,X]
            m5 = lattice_reconstruct2_double(&u1);
            w1 = matrix_times_su2_double(&m4,&m5);
            m4 = w1;

            u1 = lattice_table_2(lattice_table,&coord2,gdiT,X,&twist);      // [p+T,X]
            m5 = lattice_reconstruct2_double(&u1);
            w1 = matrix_times_su2_double(&m2,&m5);
            m2 = w1;

            tmp_R += 1.0;
        }
        lattice_neighbours_gid(&coord, &coord_4, &gdi_4,X);

        // _______________________________ top link
        u1 = lattice_table_2(lattice_table,&coord_4,gdi_4,T,&twist);        // [p+X,T]
        m3 = lattice_reconstruct2_double(&u1);
        tmp_T = 1.1;

        while (tmp_T<wilson_T){
            coord = coord_4;
            // prepare neighbours
            lattice_neighbours_gid(&coord,&coord_4, &gdi_4,T);
            u1 = lattice_table_2(lattice_table,&coord_4,gdi_4,T,&twist);    // [p,T]
            m5 = lattice_reconstruct2_double(&u1);
            w1 = matrix_times_su2_double(&m3,&m5);
            m3 = w1;

            tmp_T += 1.0;
        }
        wilson_loop += lattice_retrace_plaquette2_double(&m1,&m2,&m3,&m4); 

        // _______________________________ left and right Y links
        // coord_1 - (p,Y)
        // coord_2 - (p+T,Y)
        u1 = lattice_table_2(lattice_table,&coord_1,gdi_1,Y,&twist);        // [p,Y]
        m4 = lattice_reconstruct2_double(&u1);
        u1 = lattice_table_2(lattice_table,&coord_2,gdi_2,Y,&twist);        // [p+T,Y]
        m2 = lattice_reconstruct2_double(&u1);
        tmp_R  = 1.1;
        coord  = coord_1;
        coord2 = coord_2;

        while (tmp_R<wilson_R){
            coordX  = coord;
            coordX2 = coord2;
            // prepare neighbours
            lattice_neighbours_gid(&coordX,  &coord,  &gdi, Y);
            lattice_neighbours_gid(&coordX2, &coord2, &gdiT,Y);

            u1 = lattice_table_2(lattice_table,&coord, gdi, Y,&twist);      // [p,Y]
            m5 = lattice_reconstruct2_double(&u1);
            w1 = matrix_times_su2_double(&m4,&m5);
            m4 = w1;

            u1 = lattice_table_2(lattice_table,&coord2,gdiT,Y,&twist);      // [p+T,Y]
            m5 = lattice_reconstruct2_double(&u1);
            w1 = matrix_times_su2_double(&m2,&m5);
            m2 = w1;

            tmp_R += 1.0;
        }
        lattice_neighbours_gid(&coord, &coord_4, &gdi_4,Y);

        // _______________________________ top link
        u1 = lattice_table_2(lattice_table,&coord_4,gdi_4,T,&twist);        // [p+Y,T]
        m3 = lattice_reconstruct2_double(&u1);
        tmp_T = 1.1;

        while (tmp_T<wilson_T){
            coord = coord_4;
            // prepare neighbours
            lattice_neighbours_gid(&coord,&coord_4, &gdi_4,T);
            u1 = lattice_table_2(lattice_table,&coord_4,gdi_4,T,&twist);    // [p,T]
            m5 = lattice_reconstruct2_double(&u1);
            w1 = matrix_times_su2_double(&m3,&m5);
            m3 = w1;

            tmp_T += 1.0;
        }
        wilson_loop += lattice_retrace_plaquette2_double(&m1,&m2,&m3,&m4); 

        // _______________________________ left and right Z links
        // coord_1 - (p,Z)
        // coord_2 - (p+T,Z)
        u1 = lattice_table_2(lattice_table,&coord_1,gdi_1,Z,&twist);        // [p,Z]
        m4 = lattice_reconstruct2_double(&u1);
        u1 = lattice_table_2(lattice_table,&coord_2,gdi_2,Z,&twist);        // [p+T,Z]
        m2 = lattice_reconstruct2_double(&u1);
        tmp_R  = 1.1;
        coord  = coord_1;
        coord2 = coord_2;

        while (tmp_R<wilson_R){
            coordX  = coord;
            coordX2 = coord2;
            // prepare neighbours
            lattice_neighbours_gid(&coordX,  &coord,  &gdi, Z);
            lattice_neighbours_gid(&coordX2, &coord2, &gdiT,Z);

            u1 = lattice_table_2(lattice_table,&coord, gdi, Z,&twist);      // [p,Z]
            m5 = lattice_reconstruct2_double(&u1);
            w1 = matrix_times_su2_double(&m4,&m5);
            m4 = w1;

            u1 = lattice_table_2(lattice_table,&coord2,gdiT,Z,&twist);      // [p+T,Z]
            m5 = lattice_reconstruct2_double(&u1);
            w1 = matrix_times_su2_double(&m2,&m5);
            m2 = w1;

            tmp_R += 1.0;
        }
        lattice_neighbours_gid(&coord, &coord_4, &gdi_4,Z);

        // _______________________________ top link
        u1 = lattice_table_2(lattice_table,&coord_4,gdi_4,T,&twist);        // [p+Z,T]
        m3 = lattice_reconstruct2_double(&u1);
        tmp_T = 1.1;

        while (tmp_T<wilson_T){
            coord = coord_4;
            // prepare neighbours
            lattice_neighbours_gid(&coord,&coord_4, &gdi_4,T);
            u1 = lattice_table_2(lattice_table,&coord_4,gdi_4,T,&twist);    // [p,T]
            m5 = lattice_reconstruct2_double(&u1);
            w1 = matrix_times_su2_double(&m3,&m5);
            m3 = w1;

            tmp_T += 1.0;
        }
        wilson_loop += lattice_retrace_plaquette2_double(&m1,&m2,&m3,&m4); 

        out = wilson_loop;
    }
    // first reduction
    reduce_first_step_val_double(lattice_lds,&out, &out2);
    if((TID == 0)) lattice_measurement[BID].x = out2;
#endif

#if SUN == 3
    hgpu_double out  = 0.0;
    hgpu_double out2 = 0.0;
    uint gdi = GID;
    hgpu_double wilson_loop;
    coords_4 coord,  coord2;
    coords_4 coordX, coordX2;
    coords_4 coord_1,coord_2,coord_4;
    uint gdi_1,gdi_2,gdi_4,gdiT;
    gpu_su_3 u1;
    double_su_3 m1, m2, m3, m4, m5, w1;

    su3_twist twist;
        twist.phi   = lattice_parameters[1];
        twist.omega = lattice_parameters[2];
    hgpu_float wilson_R,wilson_T,tmp_R,tmp_T;
        wilson_R    = lattice_parameters[3];
        wilson_T    = lattice_parameters[4];
    wilson_loop = 0.0;

    lattice_lds[TID] = (hgpu_double2) 0.0;
    if (GID<SITES) {
        // _______________________________ bottom link
        lattice_gid_to_coords(&gdi,&coord_1);
        gdi_1 = gdi;
        coord_2 = coord_1;
        u1 = lattice_table_3(lattice_table,&coord_2, gdi,T,&twist);         // [p,T]
        m1 = lattice_reconstruct3_double(&u1);
        tmp_T = 1.1;

        while (tmp_T<wilson_T){
            coord = coord_2;
            // prepare neighbours
            lattice_neighbours_gid(&coord,&coord_2,&gdiT,T);
            u1 = lattice_table_3(lattice_table,&coord_2,gdiT,T,&twist);     // [p,T]
            m5 = lattice_reconstruct3_double(&u1);
            w1 = matrix_times_su3_double(&m1,&m5);
            m1 = w1;

            tmp_T += 1.0;
        }
        lattice_neighbours_gid(&coord_2,&coord,&gdi_2,T);
        coord_2 = coord;

        // _______________________________ left and right X links
        u1 = lattice_table_3(lattice_table,&coord_1,gdi_1,X,&twist);        // [p,X]
        m4 = lattice_reconstruct3_double(&u1);
        u1 = lattice_table_3(lattice_table,&coord_2,gdi_2,X,&twist);        // [p+T,X]
        m2 = lattice_reconstruct3_double(&u1);
        tmp_R  = 1.1;
        coord  = coord_1;
        coord2 = coord_2;

        while (tmp_R<wilson_R){
            coordX  = coord;
            coordX2 = coord2;
            // prepare neighbours
            lattice_neighbours_gid(&coordX,  &coord,  &gdi, X);
            lattice_neighbours_gid(&coordX2, &coord2, &gdiT,X);

            u1 = lattice_table_3(lattice_table,&coord, gdi, X,&twist);      // [p,X]
            m5 = lattice_reconstruct3_double(&u1);
            w1 = matrix_times_su3_double(&m4,&m5);
            m4 = w1;

            u1 = lattice_table_3(lattice_table,&coord2,gdiT,X,&twist);      // [p+T,X]
            m5 = lattice_reconstruct3_double(&u1);
            w1 = matrix_times_su3_double(&m2,&m5);
            m2 = w1;

            tmp_R += 1.0;
        }
        lattice_neighbours_gid(&coord, &coord_4, &gdi_4,X);

        // _______________________________ top link
        u1 = lattice_table_3(lattice_table,&coord_4,gdi_4,T,&twist);        // [p+X,T]
        m3 = lattice_reconstruct3_double(&u1);
        tmp_T = 1.1;

        while (tmp_T<wilson_T){
            coord = coord_4;
            // prepare neighbours
            lattice_neighbours_gid(&coord,&coord_4, &gdi_4,T);
            u1 = lattice_table_3(lattice_table,&coord_4,gdi_4,T,&twist);    // [p,T]
            m5 = lattice_reconstruct3_double(&u1);
            w1 = matrix_times_su3_double(&m3,&m5);
            m3 = w1;

            tmp_T += 1.0;
        }
        wilson_loop += lattice_retrace_plaquette3_double(&m1,&m2,&m3,&m4); 

        // _______________________________ left and right Y links
        u1 = lattice_table_3(lattice_table,&coord_1,gdi_1,Y,&twist);        // [p,Y]
        m4 = lattice_reconstruct3_double(&u1);
        u1 = lattice_table_3(lattice_table,&coord_2,gdi_2,Y,&twist);        // [p+T,Y]
        m2 = lattice_reconstruct3_double(&u1);
        tmp_R  = 1.1;
        coord  = coord_1;
        coord2 = coord_2;

        while (tmp_R<wilson_R){
            coordX  = coord;
            coordX2 = coord2;
            // prepare neighbours
            lattice_neighbours_gid(&coordX,  &coord,  &gdi, Y);
            lattice_neighbours_gid(&coordX2, &coord2, &gdiT,Y);

            u1 = lattice_table_3(lattice_table,&coord, gdi, Y,&twist);      // [p,Y]
            m5 = lattice_reconstruct3_double(&u1);
            w1 = matrix_times_su3_double(&m4,&m5);
            m4 = w1;

            u1 = lattice_table_3(lattice_table,&coord2,gdiT,Y,&twist);      // [p+T,Y]
            m5 = lattice_reconstruct3_double(&u1);
            w1 = matrix_times_su3_double(&m2,&m5);
            m2 = w1;

            tmp_R += 1.0;
        }
        lattice_neighbours_gid(&coord, &coord_4, &gdi_4,Y);

        // _______________________________ top link
        u1 = lattice_table_3(lattice_table,&coord_4,gdi_4,T,&twist);        // [p+Y,T]
        m3 = lattice_reconstruct3_double(&u1);
        tmp_T = 1.1;

        while (tmp_T<wilson_T){
            coord = coord_4;
            // prepare neighbours
            lattice_neighbours_gid(&coord,&coord_4, &gdi_4,T);
            u1 = lattice_table_3(lattice_table,&coord_4,gdi_4,T,&twist);    // [p,T]
            m5 = lattice_reconstruct3_double(&u1);
            w1 = matrix_times_su3_double(&m3,&m5);
            m3 = w1;

            tmp_T += 1.0;
        }
        wilson_loop += lattice_retrace_plaquette3_double(&m1,&m2,&m3,&m4); 

        // _______________________________ left and right Z links
        u1 = lattice_table_3(lattice_table,&coord_1,gdi_1,Z,&twist);        // [p,Z]
        m4 = lattice_reconstruct3_double(&u1);
        u1 = lattice_table_3(lattice_table,&coord_2,gdi_2,Z,&twist);        // [p+T,Z]
        m2 = lattice_reconstruct3_double(&u1);
        tmp_R  = 1.1;
        coord  = coord_1;
        coord2 = coord_2;

        while (tmp_R<wilson_R){
            coordX  = coord;
            coordX2 = coord2;
            // prepare neighbours
            lattice_neighbours_gid(&coordX,  &coord,  &gdi, Z);
            lattice_neighbours_gid(&coordX2, &coord2, &gdiT,Z);

            u1 = lattice_table_3(lattice_table,&coord, gdi, Z,&twist);      // [p,Z]
            m5 = lattice_reconstruct3_double(&u1);
            w1 = matrix_times_su3_double(&m4,&m5);
            m4 = w1;

            u1 = lattice_table_3(lattice_table,&coord2,gdiT,Z,&twist);      // [p+T,Z]
            m5 = lattice_reconstruct3_double(&u1);
            w1 = matrix_times_su3_double(&m2,&m5);
            m2 = w1;

            tmp_R += 1.0;
        }
        lattice_neighbours_gid(&coord, &coord_4, &gdi_4,Z);

        // _______________________________ top link
        u1 = lattice_table_3(lattice_table,&coord_4,gdi_4,T,&twist);        // [p+Z,T]
        m3 = lattice_reconstruct3_double(&u1);
        tmp_T = 1.1;

        while (tmp_T<wilson_T){
            coord = coord_4;
            // prepare neighbours
            lattice_neighbours_gid(&coord,&coord_4, &gdi_4,T);
            u1 = lattice_table_3(lattice_table,&coord_4,gdi_4,T,&twist);    // [p,T]
            m5 = lattice_reconstruct3_double(&u1);
            w1 = matrix_times_su3_double(&m3,&m5);
            m3 = w1;

            tmp_T += 1.0;
        }
        wilson_loop += lattice_retrace_plaquette3_double(&m1,&m2,&m3,&m4); 

        out = wilson_loop;
    }
    // first reduction
    reduce_first_step_val_double(lattice_lds,&out, &out2);
    if((TID == 0)) lattice_measurement[BID].x = out2;
#endif
}

                                        __kernel void
reduce_wilson_double2(__global hgpu_double2 * lattice_measurement,
                      __global hgpu_double  * lattice_wilson_loop,
                      __local hgpu_double2  * lattice_lds,
                      uint size,
                      uint index)
{
    reduce_final_step_double2(lattice_lds,lattice_measurement,size);
    if (GID==0) lattice_wilson_loop[index] = lattice_lds[0].x;
}





#endif
