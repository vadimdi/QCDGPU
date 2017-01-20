/******************************************************************************
 * @file     sun_measurements_cl.cl
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.6
 *
 * @brief    [QCDGPU]
 *           Measurements of the Wilson action, plaquette average and components of the chromoelectromagnetic field tensor
 *
 * @section  LICENSE
 *
 * Copyright (c) 2013-2016 Vadim Demchik, Natalia Kolomoyets
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
#ifndef SUNMEASUREMENTSCL_CL
#define SUNMEASUREMENTSCL_CL

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
lattice_measurement(__global hgpu_float4  * lattice_table,
                    __global hgpu_double2 * lattice_measurement,
                    __global hgpu_float   * lattice_parameters,
                    __local hgpu_double2  * lattice_lds)
{
#if SUN == 2
    hgpu_double2 out  = (hgpu_double2) 0.0;
    hgpu_double2 out2 = (hgpu_double2) 0.0;
    hgpu_double retrac_spat = 0.0;
    hgpu_double retrac_temp = 0.0;
    uint gindex = GID;
    hgpu_float bet = lattice_parameters[0];
    su2_twist twist;
    twist.phi   = lattice_parameters[1];

    if (GID<SITES) {
        coords_4 coord;
        coords_4 coordX,coordY,coordZ,coordT;
        uint gdiX,gdiY,gdiZ,gdiT;

        gpu_su_2 m1,m2,m3,m4,m5,m6;

        lattice_gid_to_coords(&gindex,&coord);

        // prepare neighbours
        lattice_neighbours_gid(&coord,&coordX,&gdiX,X);
        lattice_neighbours_gid(&coord,&coordY,&gdiY,Y);
        lattice_neighbours_gid(&coord,&coordZ,&gdiZ,Z);
        lattice_neighbours_gid(&coord,&coordT,&gdiT,T);

            m1 = lattice_table_2(lattice_table,&coord, GID,X,&twist);    // [p,x]
            m2 = lattice_table_2(lattice_table,&coordX,gdiX,Y,&twist);   // [p+X,Y]
            m3 = lattice_table_2(lattice_table,&coordY,gdiY,X,&twist);   // [p+Y,X]
            m4 = lattice_table_2(lattice_table,&coord, GID,Y,&twist);    // [p,Y]
        retrac_spat = lattice_retrace_plaquette2(&m1,&m2,&m3,&m4); // x-y: [p,X]-[p+X,Y]-[p+Y,X]*-[p,Y]

            m2 = lattice_table_2(lattice_table,&coordX,gdiX,Z,&twist);   // [p+X,Z]
            m3 = lattice_table_2(lattice_table,&coordZ,gdiZ,X,&twist);   // [p+Z,X]
            m5 = lattice_table_2(lattice_table,&coord, GID,Z,&twist);    // [p,Z]
        retrac_spat += lattice_retrace_plaquette2(&m1,&m2,&m3,&m5); // x-z: [p,X]-[p+X,Z]-[p+Z,X]*-[p,Z]*

            m2 = lattice_table_2(lattice_table,&coordX,gdiX,T,&twist);   // [p+X,T]
            m3 = lattice_table_2(lattice_table,&coordT,gdiT,X,&twist);   // [p+T,X]
            m6 = lattice_table_2(lattice_table,&coord, GID,T,&twist);    // [p,T]
        retrac_temp = lattice_retrace_plaquette2(&m1,&m2,&m3,&m6); // x-t: [p,X]-[p+X,T]-[p+T,X]*-[p,T]*

            m2 = lattice_table_2(lattice_table,&coordY,gdiY,Z,&twist);   // [p+Y,Z]
            m3 = lattice_table_2(lattice_table,&coordZ,gdiZ,Y,&twist);   // [p+Z,Y]
        retrac_spat += lattice_retrace_plaquette2(&m4,&m2,&m3,&m5); // y-z: [p,Y]-[p+Y,Z]-[p+Z,Y]*-[p,Z]*

            m2 = lattice_table_2(lattice_table,&coordY,gdiY,T,&twist);   // [p+Y,T]
            m3 = lattice_table_2(lattice_table,&coordT,gdiT,Y,&twist);   // [p+T,Y]
        retrac_temp += lattice_retrace_plaquette2(&m4,&m2,&m3,&m6); // y-t: [p,Y]-[p+Y,T]-[p+T,Y]*-[p,T]*

            m2 = lattice_table_2(lattice_table,&coordZ,gdiZ,T,&twist);   // [p+Z,T]
            m3 = lattice_table_2(lattice_table,&coordT,gdiT,Z,&twist);   // [p+T,Z]
        retrac_temp += lattice_retrace_plaquette2(&m5,&m2,&m3,&m6); // z-t: [p,Z]-[p+Z,T]-[p+T,Z]*-[p,T]*

        // first reduction
        out.x = bet * (6.0 - retrac_spat);
        out.y = bet * (6.0 - retrac_temp);
    }
    reduce_first_step_val_double2(lattice_lds,&out, &out2);
    if(TID == 0) lattice_measurement[BID] = out2;
#endif

#if SUN == 3
    hgpu_double2 out  = (hgpu_double2) 0.0;
    hgpu_double2 out2 = (hgpu_double2) 0.0;
    hgpu_double retrac_spat = 0.0;
    hgpu_double retrac_temp = 0.0;
    uint gindex = GID;
    hgpu_float bet = lattice_parameters[0];
    su3_twist twist;
    twist.phi   = lattice_parameters[1];
    twist.omega = lattice_parameters[2];

    if (GID<SITES) {
        coords_4 coord;
        coords_4 coordX,coordY,coordZ,coordT;
        uint gdiX,gdiY,gdiZ,gdiT;

        gpu_su_3 m1,m2,m3,m4,m5,m6;

        lattice_gid_to_coords(&gindex,&coord);

        // prepare neighbours
        lattice_neighbours_gid(&coord,&coordX,&gdiX,X);
        lattice_neighbours_gid(&coord,&coordY,&gdiY,Y);
        lattice_neighbours_gid(&coord,&coordZ,&gdiZ,Z);
        lattice_neighbours_gid(&coord,&coordT,&gdiT,T);

            m1 = lattice_table_3(lattice_table,&coord, GID,X,&twist);    // [p,x]
            m2 = lattice_table_3(lattice_table,&coordX,gdiX,Y,&twist);   // [p+X,Y]
            m3 = lattice_table_3(lattice_table,&coordY,gdiY,X,&twist);   // [p+Y,X]
            m4 = lattice_table_3(lattice_table,&coord, GID,Y,&twist);    // [p,Y]
        retrac_spat = lattice_retrace_plaquette3(&m1,&m2,&m3,&m4); // x-y: [p,X]-[p+X,Y]-[p+Y,X]*-[p,Y]

            m2 = lattice_table_3(lattice_table,&coordX,gdiX,Z,&twist);   // [p+X,Z]
            m3 = lattice_table_3(lattice_table,&coordZ,gdiZ,X,&twist);   // [p+Z,X]
            m5 = lattice_table_3(lattice_table,&coord, GID,Z,&twist);    // [p,Z]
        retrac_spat += lattice_retrace_plaquette3(&m1,&m2,&m3,&m5); // x-z: [p,X]-[p+X,Z]-[p+Z,X]*-[p,Z]*

            m2 = lattice_table_3(lattice_table,&coordX,gdiX,T,&twist);   // [p+X,T]
            m3 = lattice_table_3(lattice_table,&coordT,gdiT,X,&twist);   // [p+T,X]
            m6 = lattice_table_3(lattice_table,&coord, GID,T,&twist);    // [p,T]
        retrac_temp = lattice_retrace_plaquette3(&m1,&m2,&m3,&m6); // x-t: [p,X]-[p+X,T]-[p+T,X]*-[p,T]*

            m2 = lattice_table_3(lattice_table,&coordY,gdiY,Z,&twist);   // [p+Y,Z]
            m3 = lattice_table_3(lattice_table,&coordZ,gdiZ,Y,&twist);   // [p+Z,Y]
        retrac_spat += lattice_retrace_plaquette3(&m4,&m2,&m3,&m5); // y-z: [p,Y]-[p+Y,Z]-[p+Z,Y]*-[p,Z]*

            m2 = lattice_table_3(lattice_table,&coordY,gdiY,T,&twist);   // [p+Y,T]
            m3 = lattice_table_3(lattice_table,&coordT,gdiT,Y,&twist);   // [p+T,Y]
        retrac_temp += lattice_retrace_plaquette3(&m4,&m2,&m3,&m6); // y-t: [p,Y]-[p+Y,T]-[p+T,Y]*-[p,T]*

            m2 = lattice_table_3(lattice_table,&coordZ,gdiZ,T,&twist);   // [p+Z,T]
            m3 = lattice_table_3(lattice_table,&coordT,gdiT,Z,&twist);   // [p+T,Z]
        retrac_temp += lattice_retrace_plaquette3(&m5,&m2,&m3,&m6); // z-t: [p,Z]-[p+Z,T]-[p+T,Z]*-[p,T]*

        // first reduction
		out.x = bet * (9.0 - retrac_spat);
		out.y = bet * (9.0 - retrac_temp);
    }
    
    reduce_first_step_val_double2(lattice_lds,&out, &out2);

    if(TID == 0) lattice_measurement[BID] = out2;
#endif
}

#if (defined PLK) || (defined PLKx)
                                        __kernel void
lattice_action_diff_x(__global hgpu_float4  * lattice_table,
                    __global hgpu_double2 * lattice_measurement,
                    __global hgpu_float   * lattice_parameters,
                    __local hgpu_double2  * lattice_lds)
{
#if SUN == 2
    hgpu_double2 out  = (hgpu_double2) 0.0;
    hgpu_double2 out2 = (hgpu_double2) 0.0;
    hgpu_double retrac_spat = 0.0;
    hgpu_double retrac_temp = 0.0;
    
    uint gindex;// -- index on the ordinar lattice
    uint gdiK = GID;// -- index on the enlarged lattice
    
    hgpu_float bet = lattice_parameters[0];
    su2_twist twist;
    twist.phi   = lattice_parameters[1];
    
    lattice_gidK_x_to_gid(&gdiK,&gindex); // => 1 workgroup corresponds to the same x for N2 = N1 {z1 {=gid} = z3 * N2*PLK; N2 * PLK ~ workgroup_size}.
    
#ifndef BIGLAT
    if (gindex<SITES) {
#else
    if (gindex<(N1+1)*N2N3N4) {
#endif
        coords_4 coord;
        coords_4 coordX,coordY,coordZ,coordT;
        uint gdiX,gdiY,gdiZ,gdiT;

        gpu_su_2 m1,m2,m3,m4,m5,m6;

        lattice_gid_to_coords(&gindex,&coord);

        // prepare neighbours
        lattice_neighbours_gid(&coord,&coordX,&gdiX,X);
        lattice_neighbours_gid(&coord,&coordY,&gdiY,Y);
        lattice_neighbours_gid(&coord,&coordZ,&gdiZ,Z);
        lattice_neighbours_gid(&coord,&coordT,&gdiT,T);

            m1 = lattice_table_2(lattice_table,&coord, gindex,X,&twist);    // [p,x]
            m2 = lattice_table_2(lattice_table,&coordX,gdiX,Y,&twist);   // [p+X,Y]
            m3 = lattice_table_2(lattice_table,&coordY,gdiY,X,&twist);   // [p+Y,X]
            m4 = lattice_table_2(lattice_table,&coord, gindex,Y,&twist);    // [p,Y]
        retrac_spat = lattice_retrace_plaquette2(&m1,&m2,&m3,&m4); // x-y: [p,X]-[p+X,Y]-[p+Y,X]*-[p,Y]

            m2 = lattice_table_2(lattice_table,&coordX,gdiX,Z,&twist);   // [p+X,Z]
            m3 = lattice_table_2(lattice_table,&coordZ,gdiZ,X,&twist);   // [p+Z,X]
            m5 = lattice_table_2(lattice_table,&coord, gindex,Z,&twist);    // [p,Z]
        retrac_spat += lattice_retrace_plaquette2(&m1,&m2,&m3,&m5); // x-z: [p,X]-[p+X,Z]-[p+Z,X]*-[p,Z]*

            m2 = lattice_table_2(lattice_table,&coordX,gdiX,T,&twist);   // [p+X,T]
            m3 = lattice_table_2(lattice_table,&coordT,gdiT,X,&twist);   // [p+T,X]
            m6 = lattice_table_2(lattice_table,&coord, gindex,T,&twist);    // [p,T]
        retrac_temp = lattice_retrace_plaquette2(&m1,&m2,&m3,&m6); // x-t: [p,X]-[p+X,T]-[p+T,X]*-[p,T]*

            m2 = lattice_table_2(lattice_table,&coordY,gdiY,Z,&twist);   // [p+Y,Z]
            m3 = lattice_table_2(lattice_table,&coordZ,gdiZ,Y,&twist);   // [p+Z,Y]
        retrac_spat += lattice_retrace_plaquette2(&m4,&m2,&m3,&m5); // y-z: [p,Y]-[p+Y,Z]-[p+Z,Y]*-[p,Z]*

            m2 = lattice_table_2(lattice_table,&coordY,gdiY,T,&twist);   // [p+Y,T]
            m3 = lattice_table_2(lattice_table,&coordT,gdiT,Y,&twist);   // [p+T,Y]
        retrac_temp += lattice_retrace_plaquette2(&m4,&m2,&m3,&m6); // y-t: [p,Y]-[p+Y,T]-[p+T,Y]*-[p,T]*

            m2 = lattice_table_2(lattice_table,&coordZ,gdiZ,T,&twist);   // [p+Z,T]
            m3 = lattice_table_2(lattice_table,&coordT,gdiT,Z,&twist);   // [p+T,Z]
        retrac_temp += lattice_retrace_plaquette2(&m5,&m2,&m3,&m6); // z-t: [p,Z]-[p+Z,T]-[p+T,Z]*-[p,T]*

        out.x = bet * (6.0 - retrac_spat);
        out.y = bet * (6.0 - retrac_temp);
    }
    
    // first reduction
    reduce_first_step_val_double2(lattice_lds,&out, &out2);
    if(TID == 0) lattice_measurement[BID] = out2;
#elif SUN == 3
    hgpu_double2 out  = (hgpu_double2) 0.0;
    hgpu_double2 out2 = (hgpu_double2) 0.0;
    hgpu_double retrac_spat = 0.0;
    hgpu_double retrac_temp = 0.0;
    
    uint gindex;// -- index on the ordinar lattice
    uint gdiK = GID;// -- index on the enlarged lattice
    
    hgpu_float bet = lattice_parameters[0];
    su3_twist twist;
    twist.phi   = lattice_parameters[1];
    twist.omega   = lattice_parameters[2];
    
    lattice_gidK_x_to_gid(&gdiK,&gindex);
    
#ifndef BIGLAT
    if (gindex<SITES) {
#else
    if (gindex<(N1+1)*N2N3N4) {
#endif
        coords_4 coord;
        coords_4 coordX,coordY,coordZ,coordT;
        uint gdiX,gdiY,gdiZ,gdiT;

        gpu_su_3 m1,m2,m3,m4,m5,m6;

        lattice_gid_to_coords(&gindex,&coord);

        // prepare neighbours
        lattice_neighbours_gid(&coord,&coordX,&gdiX,X);
        lattice_neighbours_gid(&coord,&coordY,&gdiY,Y);
        lattice_neighbours_gid(&coord,&coordZ,&gdiZ,Z);
        lattice_neighbours_gid(&coord,&coordT,&gdiT,T);

            m1 = lattice_table_3(lattice_table,&coord, gindex,X,&twist);    // [p,x]
            m2 = lattice_table_3(lattice_table,&coordX,gdiX,Y,&twist);   // [p+X,Y]
            m3 = lattice_table_3(lattice_table,&coordY,gdiY,X,&twist);   // [p+Y,X]
            m4 = lattice_table_3(lattice_table,&coord, gindex,Y,&twist);    // [p,Y]
        retrac_spat = lattice_retrace_plaquette3(&m1,&m2,&m3,&m4); // x-y: [p,X]-[p+X,Y]-[p+Y,X]*-[p,Y]

            m2 = lattice_table_3(lattice_table,&coordX,gdiX,Z,&twist);   // [p+X,Z]
            m3 = lattice_table_3(lattice_table,&coordZ,gdiZ,X,&twist);   // [p+Z,X]
            m5 = lattice_table_3(lattice_table,&coord, gindex,Z,&twist);    // [p,Z]
        retrac_spat += lattice_retrace_plaquette3(&m1,&m2,&m3,&m5); // x-z: [p,X]-[p+X,Z]-[p+Z,X]*-[p,Z]*

            m2 = lattice_table_3(lattice_table,&coordX,gdiX,T,&twist);   // [p+X,T]
            m3 = lattice_table_3(lattice_table,&coordT,gdiT,X,&twist);   // [p+T,X]
            m6 = lattice_table_3(lattice_table,&coord, gindex,T,&twist);    // [p,T]
        retrac_temp = lattice_retrace_plaquette3(&m1,&m2,&m3,&m6); // x-t: [p,X]-[p+X,T]-[p+T,X]*-[p,T]*

            m2 = lattice_table_3(lattice_table,&coordY,gdiY,Z,&twist);   // [p+Y,Z]
            m3 = lattice_table_3(lattice_table,&coordZ,gdiZ,Y,&twist);   // [p+Z,Y]
        retrac_spat += lattice_retrace_plaquette3(&m4,&m2,&m3,&m5); // y-z: [p,Y]-[p+Y,Z]-[p+Z,Y]*-[p,Z]*

            m2 = lattice_table_3(lattice_table,&coordY,gdiY,T,&twist);   // [p+Y,T]
            m3 = lattice_table_3(lattice_table,&coordT,gdiT,Y,&twist);   // [p+T,Y]
        retrac_temp += lattice_retrace_plaquette3(&m4,&m2,&m3,&m6); // y-t: [p,Y]-[p+Y,T]-[p+T,Y]*-[p,T]*

            m2 = lattice_table_3(lattice_table,&coordZ,gdiZ,T,&twist);   // [p+Z,T]
            m3 = lattice_table_3(lattice_table,&coordT,gdiT,Z,&twist);   // [p+T,Z]
        retrac_temp += lattice_retrace_plaquette3(&m5,&m2,&m3,&m6); // z-t: [p,Z]-[p+Z,T]-[p+T,Z]*-[p,T]*

        out.x = bet * (9.0 - retrac_spat);
        out.y = bet * (9.0 - retrac_temp);
    }
    
    // first reduction
    reduce_first_step_val_double2(lattice_lds,&out, &out2);
    if(TID == 0) lattice_measurement[BID] = out2;
    
#endif
}

                                        __kernel void
lattice_action_diff_y(__global hgpu_float4  * lattice_table,
                    __global hgpu_double2 * lattice_measurement,
                    __global hgpu_float   * lattice_parameters,
                    __local hgpu_double2  * lattice_lds)
{
#if SUN == 2
    hgpu_double2 out  = (hgpu_double2) 0.0;
    hgpu_double2 out2 = (hgpu_double2) 0.0;
    hgpu_double retrac_spat = 0.0;
    hgpu_double retrac_temp = 0.0;
    
    uint gindex;// -- index on the ordinar lattice
    uint gdiK = GID;// -- index on the enlarged lattice
    
    hgpu_float bet = lattice_parameters[0];
    su2_twist twist;
    twist.phi   = lattice_parameters[1];
    
    lattice_gidK_y_to_gid(&gdiK,&gindex);
    
    if (gindex<SITES) {
        coords_4 coord;
        coords_4 coordX,coordY,coordZ,coordT;
        uint gdiX,gdiY,gdiZ,gdiT;

        gpu_su_2 m1,m2,m3,m4,m5,m6;

        lattice_gid_to_coords(&gindex,&coord);

        // prepare neighbours
        lattice_neighbours_gid(&coord,&coordX,&gdiX,X);
        lattice_neighbours_gid(&coord,&coordY,&gdiY,Y);
        lattice_neighbours_gid(&coord,&coordZ,&gdiZ,Z);
        lattice_neighbours_gid(&coord,&coordT,&gdiT,T);

            m1 = lattice_table_2(lattice_table,&coord, gindex,X,&twist);    // [p,x]
            m2 = lattice_table_2(lattice_table,&coordX,gdiX,Y,&twist);   // [p+X,Y]
            m3 = lattice_table_2(lattice_table,&coordY,gdiY,X,&twist);   // [p+Y,X]
            m4 = lattice_table_2(lattice_table,&coord, gindex,Y,&twist);    // [p,Y]
        retrac_spat = lattice_retrace_plaquette2(&m1,&m2,&m3,&m4); // x-y: [p,X]-[p+X,Y]-[p+Y,X]*-[p,Y]

            m2 = lattice_table_2(lattice_table,&coordX,gdiX,Z,&twist);   // [p+X,Z]
            m3 = lattice_table_2(lattice_table,&coordZ,gdiZ,X,&twist);   // [p+Z,X]
            m5 = lattice_table_2(lattice_table,&coord, gindex,Z,&twist);    // [p,Z]
        retrac_spat += lattice_retrace_plaquette2(&m1,&m2,&m3,&m5); // x-z: [p,X]-[p+X,Z]-[p+Z,X]*-[p,Z]*

            m2 = lattice_table_2(lattice_table,&coordX,gdiX,T,&twist);   // [p+X,T]
            m3 = lattice_table_2(lattice_table,&coordT,gdiT,X,&twist);   // [p+T,X]
            m6 = lattice_table_2(lattice_table,&coord, gindex,T,&twist);    // [p,T]
        retrac_temp = lattice_retrace_plaquette2(&m1,&m2,&m3,&m6); // x-t: [p,X]-[p+X,T]-[p+T,X]*-[p,T]*

            m2 = lattice_table_2(lattice_table,&coordY,gdiY,Z,&twist);   // [p+Y,Z]
            m3 = lattice_table_2(lattice_table,&coordZ,gdiZ,Y,&twist);   // [p+Z,Y]
        retrac_spat += lattice_retrace_plaquette2(&m4,&m2,&m3,&m5); // y-z: [p,Y]-[p+Y,Z]-[p+Z,Y]*-[p,Z]*

            m2 = lattice_table_2(lattice_table,&coordY,gdiY,T,&twist);   // [p+Y,T]
            m3 = lattice_table_2(lattice_table,&coordT,gdiT,Y,&twist);   // [p+T,Y]
        retrac_temp += lattice_retrace_plaquette2(&m4,&m2,&m3,&m6); // y-t: [p,Y]-[p+Y,T]-[p+T,Y]*-[p,T]*

            m2 = lattice_table_2(lattice_table,&coordZ,gdiZ,T,&twist);   // [p+Z,T]
            m3 = lattice_table_2(lattice_table,&coordT,gdiT,Z,&twist);   // [p+T,Z]
        retrac_temp += lattice_retrace_plaquette2(&m5,&m2,&m3,&m6); // z-t: [p,Z]-[p+Z,T]-[p+T,Z]*-[p,T]*

        out.x = bet * (6.0 - retrac_spat);
        out.y = bet * (6.0 - retrac_temp);
    }
    
    // first reduction
    reduce_first_step_val_double2(lattice_lds,&out, &out2);
    if(TID == 0) lattice_measurement[BID] = out2;
#elif SUN == 3
    hgpu_double2 out  = (hgpu_double2) 0.0;
    hgpu_double2 out2 = (hgpu_double2) 0.0;
    hgpu_double retrac_spat = 0.0;
    hgpu_double retrac_temp = 0.0;
    
    uint gindex;// -- index on the ordinar lattice
    uint gdiK = GID;// -- index on the enlarged lattice
    
    hgpu_float bet = lattice_parameters[0];
    su3_twist twist;
    twist.phi   = lattice_parameters[1];
    twist.omega   = lattice_parameters[2];
    
    lattice_gidK_y_to_gid(&gdiK,&gindex);
    
    if (gindex<SITES) {
        coords_4 coord;
        coords_4 coordX,coordY,coordZ,coordT;
        uint gdiX,gdiY,gdiZ,gdiT;

        gpu_su_3 m1,m2,m3,m4,m5,m6;

        lattice_gid_to_coords(&gindex,&coord);

        // prepare neighbours
        lattice_neighbours_gid(&coord,&coordX,&gdiX,X);
        lattice_neighbours_gid(&coord,&coordY,&gdiY,Y);
        lattice_neighbours_gid(&coord,&coordZ,&gdiZ,Z);
        lattice_neighbours_gid(&coord,&coordT,&gdiT,T);

            m1 = lattice_table_3(lattice_table,&coord, gindex,X,&twist);    // [p,x]
            m2 = lattice_table_3(lattice_table,&coordX,gdiX,Y,&twist);   // [p+X,Y]
            m3 = lattice_table_3(lattice_table,&coordY,gdiY,X,&twist);   // [p+Y,X]
            m4 = lattice_table_3(lattice_table,&coord, gindex,Y,&twist);    // [p,Y]
        retrac_spat = lattice_retrace_plaquette3(&m1,&m2,&m3,&m4); // x-y: [p,X]-[p+X,Y]-[p+Y,X]*-[p,Y]

            m2 = lattice_table_3(lattice_table,&coordX,gdiX,Z,&twist);   // [p+X,Z]
            m3 = lattice_table_3(lattice_table,&coordZ,gdiZ,X,&twist);   // [p+Z,X]
            m5 = lattice_table_3(lattice_table,&coord, gindex,Z,&twist);    // [p,Z]
        retrac_spat += lattice_retrace_plaquette3(&m1,&m2,&m3,&m5); // x-z: [p,X]-[p+X,Z]-[p+Z,X]*-[p,Z]*

            m2 = lattice_table_3(lattice_table,&coordX,gdiX,T,&twist);   // [p+X,T]
            m3 = lattice_table_3(lattice_table,&coordT,gdiT,X,&twist);   // [p+T,X]
            m6 = lattice_table_3(lattice_table,&coord, gindex,T,&twist);    // [p,T]
        retrac_temp = lattice_retrace_plaquette3(&m1,&m2,&m3,&m6); // x-t: [p,X]-[p+X,T]-[p+T,X]*-[p,T]*

            m2 = lattice_table_3(lattice_table,&coordY,gdiY,Z,&twist);   // [p+Y,Z]
            m3 = lattice_table_3(lattice_table,&coordZ,gdiZ,Y,&twist);   // [p+Z,Y]
        retrac_spat += lattice_retrace_plaquette3(&m4,&m2,&m3,&m5); // y-z: [p,Y]-[p+Y,Z]-[p+Z,Y]*-[p,Z]*

            m2 = lattice_table_3(lattice_table,&coordY,gdiY,T,&twist);   // [p+Y,T]
            m3 = lattice_table_3(lattice_table,&coordT,gdiT,Y,&twist);   // [p+T,Y]
        retrac_temp += lattice_retrace_plaquette3(&m4,&m2,&m3,&m6); // y-t: [p,Y]-[p+Y,T]-[p+T,Y]*-[p,T]*

            m2 = lattice_table_3(lattice_table,&coordZ,gdiZ,T,&twist);   // [p+Z,T]
            m3 = lattice_table_3(lattice_table,&coordT,gdiT,Z,&twist);   // [p+T,Z]
        retrac_temp += lattice_retrace_plaquette3(&m5,&m2,&m3,&m6); // z-t: [p,Z]-[p+Z,T]-[p+T,Z]*-[p,T]*

        out.x = bet * (9.0 - retrac_spat);
        out.y = bet * (9.0 - retrac_temp);
    }
    
    // first reduction
    reduce_first_step_val_double2(lattice_lds,&out, &out2);
    if(TID == 0) lattice_measurement[BID] = out2;
#endif
}

                                        __kernel void
lattice_action_diff_z(__global hgpu_float4  * lattice_table,
                    __global hgpu_double2 * lattice_measurement,
                    __global hgpu_float   * lattice_parameters,
                    __local hgpu_double2  * lattice_lds)
{
#if SUN == 2
    hgpu_double2 out  = (hgpu_double2) 0.0;
    hgpu_double2 out2 = (hgpu_double2) 0.0;
    hgpu_double retrac_spat = 0.0;
    hgpu_double retrac_temp = 0.0;
    
    uint gindex;// -- index on the ordinar lattice
    uint gdiK = GID;// -- index on the enlarged lattice
    
    hgpu_float bet = lattice_parameters[0];
    su2_twist twist;
    twist.phi   = lattice_parameters[1];
    
    lattice_gidK_z_to_gid(&gdiK,&gindex);
    
    if (gindex<SITES) {
        coords_4 coord;
        coords_4 coordX,coordY,coordZ,coordT;
        uint gdiX,gdiY,gdiZ,gdiT;

        gpu_su_2 m1,m2,m3,m4,m5,m6;

        lattice_gid_to_coords(&gindex,&coord);

        // prepare neighbours
        lattice_neighbours_gid(&coord,&coordX,&gdiX,X);
        lattice_neighbours_gid(&coord,&coordY,&gdiY,Y);
        lattice_neighbours_gid(&coord,&coordZ,&gdiZ,Z);
        lattice_neighbours_gid(&coord,&coordT,&gdiT,T);

            m1 = lattice_table_2(lattice_table,&coord, gindex,X,&twist);    // [p,x]
            m2 = lattice_table_2(lattice_table,&coordX,gdiX,Y,&twist);   // [p+X,Y]
            m3 = lattice_table_2(lattice_table,&coordY,gdiY,X,&twist);   // [p+Y,X]
            m4 = lattice_table_2(lattice_table,&coord, gindex,Y,&twist);    // [p,Y]
        retrac_spat = lattice_retrace_plaquette2(&m1,&m2,&m3,&m4); // x-y: [p,X]-[p+X,Y]-[p+Y,X]*-[p,Y]

            m2 = lattice_table_2(lattice_table,&coordX,gdiX,Z,&twist);   // [p+X,Z]
            m3 = lattice_table_2(lattice_table,&coordZ,gdiZ,X,&twist);   // [p+Z,X]
            m5 = lattice_table_2(lattice_table,&coord, gindex,Z,&twist);    // [p,Z]
        retrac_spat += lattice_retrace_plaquette2(&m1,&m2,&m3,&m5); // x-z: [p,X]-[p+X,Z]-[p+Z,X]*-[p,Z]*

            m2 = lattice_table_2(lattice_table,&coordX,gdiX,T,&twist);   // [p+X,T]
            m3 = lattice_table_2(lattice_table,&coordT,gdiT,X,&twist);   // [p+T,X]
            m6 = lattice_table_2(lattice_table,&coord, gindex,T,&twist);    // [p,T]
        retrac_temp = lattice_retrace_plaquette2(&m1,&m2,&m3,&m6); // x-t: [p,X]-[p+X,T]-[p+T,X]*-[p,T]*

            m2 = lattice_table_2(lattice_table,&coordY,gdiY,Z,&twist);   // [p+Y,Z]
            m3 = lattice_table_2(lattice_table,&coordZ,gdiZ,Y,&twist);   // [p+Z,Y]
        retrac_spat += lattice_retrace_plaquette2(&m4,&m2,&m3,&m5); // y-z: [p,Y]-[p+Y,Z]-[p+Z,Y]*-[p,Z]*

            m2 = lattice_table_2(lattice_table,&coordY,gdiY,T,&twist);   // [p+Y,T]
            m3 = lattice_table_2(lattice_table,&coordT,gdiT,Y,&twist);   // [p+T,Y]
        retrac_temp += lattice_retrace_plaquette2(&m4,&m2,&m3,&m6); // y-t: [p,Y]-[p+Y,T]-[p+T,Y]*-[p,T]*

            m2 = lattice_table_2(lattice_table,&coordZ,gdiZ,T,&twist);   // [p+Z,T]
            m3 = lattice_table_2(lattice_table,&coordT,gdiT,Z,&twist);   // [p+T,Z]
        retrac_temp += lattice_retrace_plaquette2(&m5,&m2,&m3,&m6); // z-t: [p,Z]-[p+Z,T]-[p+T,Z]*-[p,T]*

        out.x = bet * (6.0 - retrac_spat);
        out.y = bet * (6.0 - retrac_temp);
    }
    
    // first reduction
    reduce_first_step_val_double2(lattice_lds,&out, &out2);
    if(TID == 0) lattice_measurement[BID] = out2;
#elif SUN == 3
    hgpu_double2 out  = (hgpu_double2) 0.0;
    hgpu_double2 out2 = (hgpu_double2) 0.0;
    hgpu_double retrac_spat = 0.0;
    hgpu_double retrac_temp = 0.0;
    
    uint gindex;// -- index on the ordinar lattice
    uint gdiK = GID;// -- index on the enlarged lattice
    
    hgpu_float bet = lattice_parameters[0];
    su3_twist twist;
    twist.phi   = lattice_parameters[1];
    twist.omega   = lattice_parameters[2];
    
    lattice_gidK_z_to_gid(&gdiK,&gindex);
    
    if (gindex<SITES) {
        coords_4 coord;
        coords_4 coordX,coordY,coordZ,coordT;
        uint gdiX,gdiY,gdiZ,gdiT;

        gpu_su_3 m1,m2,m3,m4,m5,m6;

        lattice_gid_to_coords(&gindex,&coord);

        // prepare neighbours
        lattice_neighbours_gid(&coord,&coordX,&gdiX,X);
        lattice_neighbours_gid(&coord,&coordY,&gdiY,Y);
        lattice_neighbours_gid(&coord,&coordZ,&gdiZ,Z);
        lattice_neighbours_gid(&coord,&coordT,&gdiT,T);

            m1 = lattice_table_3(lattice_table,&coord, gindex,X,&twist);    // [p,x]
            m2 = lattice_table_3(lattice_table,&coordX,gdiX,Y,&twist);   // [p+X,Y]
            m3 = lattice_table_3(lattice_table,&coordY,gdiY,X,&twist);   // [p+Y,X]
            m4 = lattice_table_3(lattice_table,&coord, gindex,Y,&twist);    // [p,Y]
        retrac_spat = lattice_retrace_plaquette3(&m1,&m2,&m3,&m4); // x-y: [p,X]-[p+X,Y]-[p+Y,X]*-[p,Y]

            m2 = lattice_table_3(lattice_table,&coordX,gdiX,Z,&twist);   // [p+X,Z]
            m3 = lattice_table_3(lattice_table,&coordZ,gdiZ,X,&twist);   // [p+Z,X]
            m5 = lattice_table_3(lattice_table,&coord, gindex,Z,&twist);    // [p,Z]
        retrac_spat += lattice_retrace_plaquette3(&m1,&m2,&m3,&m5); // x-z: [p,X]-[p+X,Z]-[p+Z,X]*-[p,Z]*

            m2 = lattice_table_3(lattice_table,&coordX,gdiX,T,&twist);   // [p+X,T]
            m3 = lattice_table_3(lattice_table,&coordT,gdiT,X,&twist);   // [p+T,X]
            m6 = lattice_table_3(lattice_table,&coord, gindex,T,&twist);    // [p,T]
        retrac_temp = lattice_retrace_plaquette3(&m1,&m2,&m3,&m6); // x-t: [p,X]-[p+X,T]-[p+T,X]*-[p,T]*

            m2 = lattice_table_3(lattice_table,&coordY,gdiY,Z,&twist);   // [p+Y,Z]
            m3 = lattice_table_3(lattice_table,&coordZ,gdiZ,Y,&twist);   // [p+Z,Y]
        retrac_spat += lattice_retrace_plaquette3(&m4,&m2,&m3,&m5); // y-z: [p,Y]-[p+Y,Z]-[p+Z,Y]*-[p,Z]*

            m2 = lattice_table_3(lattice_table,&coordY,gdiY,T,&twist);   // [p+Y,T]
            m3 = lattice_table_3(lattice_table,&coordT,gdiT,Y,&twist);   // [p+T,Y]
        retrac_temp += lattice_retrace_plaquette3(&m4,&m2,&m3,&m6); // y-t: [p,Y]-[p+Y,T]-[p+T,Y]*-[p,T]*

            m2 = lattice_table_3(lattice_table,&coordZ,gdiZ,T,&twist);   // [p+Z,T]
            m3 = lattice_table_3(lattice_table,&coordT,gdiT,Z,&twist);   // [p+T,Z]
        retrac_temp += lattice_retrace_plaquette3(&m5,&m2,&m3,&m6); // z-t: [p,Z]-[p+Z,T]-[p+T,Z]*-[p,T]*

        out.x = bet * (9.0 - retrac_spat);
        out.y = bet * (9.0 - retrac_temp);
    }
    
    // first reduction
    reduce_first_step_val_double2(lattice_lds,&out, &out2);
    if(TID == 0) lattice_measurement[BID] = out2;
#endif
}
#endif

                                        __kernel void
reduce_measurement_double2(__global hgpu_double2 * lattice_measurement,
                           __global hgpu_double2 * lattice_energies,
                           __local hgpu_double2  * lattice_lds,
                           uint size,
                           uint index)
{
    reduce_final_step_double2(lattice_lds,lattice_measurement,size);
    hgpu_double2 out = lattice_lds[TID];
    if (GID==0) lattice_energies[index] = out;
}

                                        __kernel void
reduce_action_diff_x_double2(__global hgpu_double2 * lattice_measurement,
                           __global hgpu_double2 * lattice_energies,
                           __local hgpu_double2  * lattice_lds,
                           uint4 param,
                           uint index)
{
    hgpu_double2 out[N1];
    uint size1 = param.x / N1;
    int i;
    
    reduce_final_step_double2(lattice_lds, lattice_measurement, size1);
    out[0] = lattice_lds[TID];
    
    uint offset  = size1;
    for(i = 1; i < N1; i++)
    {
       reduce_final_step_double2_offset(lattice_lds, lattice_measurement, size1, i * offset);
       out[i] = lattice_lds[TID];
    }
    
    uint offset2 = param.y;
    if (GID==0)
      for(i = 0; i < N1; i++)
	lattice_energies[index + offset2 * i] = out[i];
}

                                        __kernel void
reduce_action_diff_y_double2(__global hgpu_double2 * lattice_measurement,
                           __global hgpu_double2 * lattice_energies,
                           __local hgpu_double2  * lattice_lds,
                           uint4 param,
                           uint index)
{
    hgpu_double2 out[N2];
    uint size1 = param.x / N2;
    int i;
    
    reduce_final_step_double2(lattice_lds, lattice_measurement, size1);
    out[0] = lattice_lds[TID];
    
    uint offset  = size1;
    for(i = 1; i < N2; i++)
    {
       reduce_final_step_double2_offset(lattice_lds, lattice_measurement, size1, i * offset);
       out[i] = lattice_lds[TID];
    }
    
    uint offset2 = param.y;
    if (GID==0)
      for(i = 0; i < N2; i++)
	lattice_energies[index + offset2 * i] = out[i];
}

                                        __kernel void
reduce_action_diff_z_double2(__global hgpu_double2 * lattice_measurement,
                           __global hgpu_double2 * lattice_energies,
                           __local hgpu_double2  * lattice_lds,
                           uint4 param,
                           uint index)
{
    hgpu_double2 out[N3];
    uint size1 = param.x / N3;
    int i;
    
    reduce_final_step_double2(lattice_lds, lattice_measurement, size1);
    out[0] = lattice_lds[TID];
    
    uint offset  = size1;
    for(i = 1; i < N3; i++)
    {
       reduce_final_step_double2_offset(lattice_lds, lattice_measurement, size1, i * offset);
       out[i] = lattice_lds[TID];
    }
    
    uint offset2 = param.y;
    if (GID==0)
      for(i = 0; i < N3; i++)
	lattice_energies[index + offset2 * i] = out[i];
}

#ifdef FMUNU
// Fmunu measurement for H field
                                        __kernel void
lattice_measurement_plq(__global hgpu_float4  * lattice_table,
                        __global hgpu_double2 * lattice_measurement,
                        __global hgpu_float   * lattice_parameters,
                        __local  hgpu_double2 * lattice_lds,
                        uint offset)
{
#if SUN == 2
    hgpu_double2 out  = (hgpu_double2) 0.0;
    hgpu_double2 out2 = (hgpu_double2) 0.0;

    hgpu_double2 Fxy3_reduced = (hgpu_double2) 0.0;
    hgpu_double2 Fxz3_reduced = (hgpu_double2) 0.0;
    hgpu_double2 Fyz3_reduced = (hgpu_double2) 0.0;
    hgpu_double2 Fxy8_reduced = (hgpu_double2) 0.0;
    hgpu_double2 Fxz8_reduced = (hgpu_double2) 0.0;
    hgpu_double2 Fyz8_reduced = (hgpu_double2) 0.0;


    hgpu_double2 out_Fxy3 = (hgpu_double2) 0.0;
    hgpu_double2 out_Fxz3 = (hgpu_double2) 0.0;
    hgpu_double2 out_Fyz3 = (hgpu_double2) 0.0;
    hgpu_double2 out_Fxy8 = (hgpu_double2) 0.0;
    hgpu_double2 out_Fxz8 = (hgpu_double2) 0.0;
    hgpu_double2 out_Fyz8 = (hgpu_double2) 0.0;

    hgpu_double retrac_spat = 0.0;
    hgpu_double retrac_temp = 0.0;

    uint gindex = GID;
    su2_twist twist;
    twist.phi   = lattice_parameters[1];
    hgpu_complex_double Fxy3, Fxz3, Fyz3;
    hgpu_complex_double Fxy8, Fxz8, Fyz8;

    lattice_lds[TID] = (hgpu_double2) 0.0;
    if (GID<SITES) {
        coords_4 coord;
        coords_4 coordX,coordY,coordZ,coordT;
        uint gdiX,gdiY,gdiZ,gdiT;

        gpu_su_2 m1,m2,m3,m4,m5,m6;

        lattice_gid_to_coords(&gindex,&coord);

        // prepare neighbours
        lattice_neighbours_gid(&coord,&coordX,&gdiX,X);
        lattice_neighbours_gid(&coord,&coordY,&gdiY,Y);
        lattice_neighbours_gid(&coord,&coordZ,&gdiZ,Z);
        lattice_neighbours_gid(&coord,&coordT,&gdiT,T);

            m1 = lattice_table_2(lattice_table,&coord, GID,X,&twist);    // [p,x]
            m2 = lattice_table_2(lattice_table,&coordX,gdiX,Y,&twist);   // [p+X,Y]
            m3 = lattice_table_2(lattice_table,&coordY,gdiY,X,&twist);   // [p+Y,X]
            m4 = lattice_table_2(lattice_table,&coord, GID,Y,&twist);    // [p,Y]
        retrac_spat = lattice_retrace_plaquette2_F(&m1,&m2,&m3,&m4,&Fxy3,&Fxy8); // x-y: [p,X]-[p+X,Y]-[p+Y,X]*-[p,Y]

            m2 = lattice_table_2(lattice_table,&coordX,gdiX,Z,&twist);   // [p+X,Z]
            m3 = lattice_table_2(lattice_table,&coordZ,gdiZ,X,&twist);   // [p+Z,X]
            m5 = lattice_table_2(lattice_table,&coord, GID,Z,&twist);    // [p,Z]
        retrac_spat += lattice_retrace_plaquette2_F(&m1,&m2,&m3,&m5,&Fxz3,&Fxz8); // x-z: [p,X]-[p+X,Z]-[p+Z,X]*-[p,Z]*

            m2 = lattice_table_2(lattice_table,&coordX,gdiX,T,&twist);   // [p+X,T]
            m3 = lattice_table_2(lattice_table,&coordT,gdiT,X,&twist);   // [p+T,X]
            m6 = lattice_table_2(lattice_table,&coord, GID,T,&twist);    // [p,T]
        retrac_temp = lattice_retrace_plaquette2(&m1,&m2,&m3,&m6); // x-t: [p,X]-[p+X,T]-[p+T,X]*-[p,T]*

            m2 = lattice_table_2(lattice_table,&coordY,gdiY,Z,&twist);   // [p+Y,Z]
            m3 = lattice_table_2(lattice_table,&coordZ,gdiZ,Y,&twist);   // [p+Z,Y]
        retrac_spat += lattice_retrace_plaquette2_F(&m4,&m2,&m3,&m5,&Fyz3,&Fyz8); // y-z: [p,Y]-[p+Y,Z]-[p+Z,Y]*-[p,Z]*

            m2 = lattice_table_2(lattice_table,&coordY,gdiY,T,&twist);   // [p+Y,T]
            m3 = lattice_table_2(lattice_table,&coordT,gdiT,Y,&twist);   // [p+T,Y]
        retrac_temp += lattice_retrace_plaquette2(&m4,&m2,&m3,&m6); // y-t: [p,Y]-[p+Y,T]-[p+T,Y]*-[p,T]*

            m2 = lattice_table_2(lattice_table,&coordZ,gdiZ,T,&twist);   // [p+Z,T]
            m3 = lattice_table_2(lattice_table,&coordT,gdiT,Z,&twist);   // [p+T,Z]
        retrac_temp += lattice_retrace_plaquette2(&m5,&m2,&m3,&m6); // z-t: [p,Z]-[p+Z,T]-[p+T,Z]*-[p,T]*

        out      = (hgpu_double2) {retrac_spat,retrac_temp};
        out_Fxy3 = (hgpu_double2) {Fxy3.re,Fxy3.im};
        out_Fxz3 = (hgpu_double2) {Fxz3.re,Fxz3.im};
        out_Fyz3 = (hgpu_double2) {Fyz3.re,Fyz3.im};
        out_Fxy8 = (hgpu_double2) {Fxy8.re,Fxy8.im};
        out_Fxz8 = (hgpu_double2) {Fxz8.re,Fxz8.im};
        out_Fyz8 = (hgpu_double2) {Fyz8.re,Fyz8.im};
    }

    // first reduction
    reduce_first_step_val_double2(lattice_lds,&out, &out2);

    // first reduction - Fxy3
    reduce_first_step_val_double2(lattice_lds,&out_Fxy3, &Fxy3_reduced);

    // first reduction - Fxz3
    reduce_first_step_val_double2(lattice_lds,&out_Fxz3, &Fxz3_reduced);

    // first reduction - Fyz3
    reduce_first_step_val_double2(lattice_lds,&out_Fyz3, &Fyz3_reduced);

    // first reduction - Fxy8
    reduce_first_step_val_double2(lattice_lds,&out_Fxy8, &Fxy8_reduced);

    // first reduction - Fxz8
    reduce_first_step_val_double2(lattice_lds,&out_Fxz8, &Fxz8_reduced);

    // first reduction - Fyz8
    reduce_first_step_val_double2(lattice_lds,&out_Fyz8, &Fyz8_reduced);

    if(TID == 0) {
        lattice_measurement[BID             ] = out2;
        lattice_measurement[BID +     offset] = Fxy3_reduced;
        lattice_measurement[BID + 2 * offset] = Fxz3_reduced;
        lattice_measurement[BID + 3 * offset] = Fyz3_reduced;
        lattice_measurement[BID + 4 * offset] = Fxy8_reduced;
        lattice_measurement[BID + 5 * offset] = Fxz8_reduced;
        lattice_measurement[BID + 6 * offset] = Fyz8_reduced;
    }
#endif

#if SUN == 3
    hgpu_double2 out  = (hgpu_double2) 0.0;
    hgpu_double2 out2 = (hgpu_double2) 0.0;

    hgpu_double2 Fxy3_reduced = (hgpu_double2) 0.0;
    hgpu_double2 Fxz3_reduced = (hgpu_double2) 0.0;
    hgpu_double2 Fyz3_reduced = (hgpu_double2) 0.0;
    hgpu_double2 Fxy8_reduced = (hgpu_double2) 0.0;
    hgpu_double2 Fxz8_reduced = (hgpu_double2) 0.0;
    hgpu_double2 Fyz8_reduced = (hgpu_double2) 0.0;


    hgpu_double2 out_Fxy3 = (hgpu_double2) 0.0;
    hgpu_double2 out_Fxz3 = (hgpu_double2) 0.0;
    hgpu_double2 out_Fyz3 = (hgpu_double2) 0.0;
    hgpu_double2 out_Fxy8 = (hgpu_double2) 0.0;
    hgpu_double2 out_Fxz8 = (hgpu_double2) 0.0;
    hgpu_double2 out_Fyz8 = (hgpu_double2) 0.0;

    hgpu_double retrac_spat = 0.0;
    hgpu_double retrac_temp = 0.0;

    uint gindex = GID;
    su3_twist twist;
    twist.phi   = lattice_parameters[1];
    twist.omega = lattice_parameters[2];
    hgpu_complex_double Fxy3, Fxz3, Fyz3;
    hgpu_complex_double Fxy8, Fxz8, Fyz8;

    lattice_lds[TID] = (hgpu_double2) 0.0;
    if (GID<SITES) {
        coords_4 coord;
        coords_4 coordX,coordY,coordZ,coordT;
        uint gdiX,gdiY,gdiZ,gdiT;

        gpu_su_3 m1,m2,m3,m4,m5,m6;

        lattice_gid_to_coords(&gindex,&coord);

        // prepare neighbours
        lattice_neighbours_gid(&coord,&coordX,&gdiX,X);
        lattice_neighbours_gid(&coord,&coordY,&gdiY,Y);
        lattice_neighbours_gid(&coord,&coordZ,&gdiZ,Z);
        lattice_neighbours_gid(&coord,&coordT,&gdiT,T);

            m1 = lattice_table_3(lattice_table,&coord, GID,X,&twist);    // [p,x]
            m2 = lattice_table_3(lattice_table,&coordX,gdiX,Y,&twist);   // [p+X,Y]
            m3 = lattice_table_3(lattice_table,&coordY,gdiY,X,&twist);   // [p+Y,X]
            m4 = lattice_table_3(lattice_table,&coord, GID,Y,&twist);    // [p,Y]
        retrac_spat = lattice_retrace_plaquette3_F(&m1,&m2,&m3,&m4,&Fxy3,&Fxy8); // x-y: [p,X]-[p+X,Y]-[p+Y,X]*-[p,Y]

            m2 = lattice_table_3(lattice_table,&coordX,gdiX,Z,&twist);   // [p+X,Z]
            m3 = lattice_table_3(lattice_table,&coordZ,gdiZ,X,&twist);   // [p+Z,X]
            m5 = lattice_table_3(lattice_table,&coord, GID,Z,&twist);    // [p,Z]
        retrac_spat += lattice_retrace_plaquette3_F(&m1,&m2,&m3,&m5,&Fxz3,&Fxz8); // x-z: [p,X]-[p+X,Z]-[p+Z,X]*-[p,Z]*

            m2 = lattice_table_3(lattice_table,&coordX,gdiX,T,&twist);   // [p+X,T]
            m3 = lattice_table_3(lattice_table,&coordT,gdiT,X,&twist);   // [p+T,X]
            m6 = lattice_table_3(lattice_table,&coord, GID,T,&twist);    // [p,T]
        retrac_temp = lattice_retrace_plaquette3(&m1,&m2,&m3,&m6); // x-t: [p,X]-[p+X,T]-[p+T,X]*-[p,T]*

            m2 = lattice_table_3(lattice_table,&coordY,gdiY,Z,&twist);   // [p+Y,Z]
            m3 = lattice_table_3(lattice_table,&coordZ,gdiZ,Y,&twist);   // [p+Z,Y]
        retrac_spat += lattice_retrace_plaquette3_F(&m4,&m2,&m3,&m5,&Fyz3,&Fyz8); // y-z: [p,Y]-[p+Y,Z]-[p+Z,Y]*-[p,Z]*

            m2 = lattice_table_3(lattice_table,&coordY,gdiY,T,&twist);   // [p+Y,T]
            m3 = lattice_table_3(lattice_table,&coordT,gdiT,Y,&twist);   // [p+T,Y]
        retrac_temp += lattice_retrace_plaquette3(&m4,&m2,&m3,&m6); // y-t: [p,Y]-[p+Y,T]-[p+T,Y]*-[p,T]*

            m2 = lattice_table_3(lattice_table,&coordZ,gdiZ,T,&twist);   // [p+Z,T]
            m3 = lattice_table_3(lattice_table,&coordT,gdiT,Z,&twist);   // [p+T,Z]
        retrac_temp += lattice_retrace_plaquette3(&m5,&m2,&m3,&m6); // z-t: [p,Z]-[p+Z,T]-[p+T,Z]*-[p,T]*

        out      = (hgpu_double2) {retrac_spat,retrac_temp};
        out_Fxy3 = (hgpu_double2) {Fxy3.re,Fxy3.im};
        out_Fxz3 = (hgpu_double2) {Fxz3.re,Fxz3.im};
        out_Fyz3 = (hgpu_double2) {Fyz3.re,Fyz3.im};
        out_Fxy8 = (hgpu_double2) {Fxy8.re,Fxy8.im};
        out_Fxz8 = (hgpu_double2) {Fxz8.re,Fxz8.im};
        out_Fyz8 = (hgpu_double2) {Fyz8.re,Fyz8.im};
    }

    // first reduction
    reduce_first_step_val_double2(lattice_lds,&out, &out2);

    // first reduction - Fxy3
    reduce_first_step_val_double2(lattice_lds,&out_Fxy3, &Fxy3_reduced);

    // first reduction - Fxz3
    reduce_first_step_val_double2(lattice_lds,&out_Fxz3, &Fxz3_reduced);

    // first reduction - Fyz3
    reduce_first_step_val_double2(lattice_lds,&out_Fyz3, &Fyz3_reduced);

    // first reduction - Fxy8
    reduce_first_step_val_double2(lattice_lds,&out_Fxy8, &Fxy8_reduced);

    // first reduction - Fxz8
    reduce_first_step_val_double2(lattice_lds,&out_Fxz8, &Fxz8_reduced);

    // first reduction - Fyz8
    reduce_first_step_val_double2(lattice_lds,&out_Fyz8, &Fyz8_reduced);

    if(TID == 0) {
        lattice_measurement[BID             ] = out2;
        lattice_measurement[BID +     offset] = Fxy3_reduced;
        lattice_measurement[BID + 2 * offset] = Fxz3_reduced;
        lattice_measurement[BID + 3 * offset] = Fyz3_reduced;
        lattice_measurement[BID + 4 * offset] = Fxy8_reduced;
        lattice_measurement[BID + 5 * offset] = Fxz8_reduced;
        lattice_measurement[BID + 6 * offset] = Fyz8_reduced;
    }
#endif
}
#elif defined(F0MU)
// Fmunu measurement for E field
                                        __kernel void
lattice_measurement_plq(__global hgpu_float4  * lattice_table,
                        __global hgpu_double2 * lattice_measurement,
                        __global hgpu_float   * lattice_parameters,
                        __local hgpu_double2  * lattice_lds,
                        uint offset)
{
#if SUN == 2
    hgpu_double2 out  = (hgpu_double2) 0.0;
    hgpu_double2 out2 = (hgpu_double2) 0.0;

    hgpu_double2 Fxt3_reduced = (hgpu_double2) 0.0;
    hgpu_double2 Fyt3_reduced = (hgpu_double2) 0.0;
    hgpu_double2 Fzt3_reduced = (hgpu_double2) 0.0;
    hgpu_double2 Fxt8_reduced = (hgpu_double2) 0.0;
    hgpu_double2 Fyt8_reduced = (hgpu_double2) 0.0;
    hgpu_double2 Fzt8_reduced = (hgpu_double2) 0.0;

    hgpu_double2 out_Fxt3 = (hgpu_double2) 0.0;
    hgpu_double2 out_Fyt3 = (hgpu_double2) 0.0;
    hgpu_double2 out_Fzt3 = (hgpu_double2) 0.0;
    hgpu_double2 out_Fxt8 = (hgpu_double2) 0.0;
    hgpu_double2 out_Fyt8 = (hgpu_double2) 0.0;
    hgpu_double2 out_Fzt8 = (hgpu_double2) 0.0;

    hgpu_double retrac_spat = 0.0;
    hgpu_double retrac_temp = 0.0;

    uint gindex = GID;
    su2_twist twist;
    twist.phi   = lattice_parameters[1];
    hgpu_complex_double Fxt3, Fyt3, Fzt3;
    hgpu_complex_double Fxt8, Fyt8, Fzt8;

    lattice_lds[TID] = (hgpu_double2) 0.0;
    if (GID<SITES) {
        coords_4 coord;
        coords_4 coordX,coordY,coordZ,coordT;
        uint gdiX,gdiY,gdiZ,gdiT;

        gpu_su_2 m1,m2,m3,m4,m5,m6;

        lattice_gid_to_coords(&gindex,&coord);

        // prepare neighbours
        lattice_neighbours_gid(&coord,&coordX,&gdiX,X);
        lattice_neighbours_gid(&coord,&coordY,&gdiY,Y);
        lattice_neighbours_gid(&coord,&coordZ,&gdiZ,Z);
        lattice_neighbours_gid(&coord,&coordT,&gdiT,T);

            m1 = lattice_table_2(lattice_table,&coord, GID,X,&twist);    // [p,x]
            m2 = lattice_table_2(lattice_table,&coordX,gdiX,Y,&twist);   // [p+X,Y]
            m3 = lattice_table_2(lattice_table,&coordY,gdiY,X,&twist);   // [p+Y,X]
            m4 = lattice_table_2(lattice_table,&coord, GID,Y,&twist);    // [p,Y]
        retrac_spat += lattice_retrace_plaquette2(&m1,&m2,&m3,&m4);       // x-y: [p,X]-[p+X,Y]-[p+Y,X]*-[p,Y]

            m2 = lattice_table_2(lattice_table,&coordX,gdiX,Z,&twist);   // [p+X,Z]
            m3 = lattice_table_2(lattice_table,&coordZ,gdiZ,X,&twist);   // [p+Z,X]
            m5 = lattice_table_2(lattice_table,&coord, GID,Z,&twist);    // [p,Z]
        retrac_spat += lattice_retrace_plaquette2(&m1,&m2,&m3,&m5);      // x-z: [p,X]-[p+X,Z]-[p+Z,X]*-[p,Z]*

            m2 = lattice_table_2(lattice_table,&coordX,gdiX,T,&twist);   // [p+X,T]
            m3 = lattice_table_2(lattice_table,&coordT,gdiT,X,&twist);   // [p+T,X]
            m6 = lattice_table_2(lattice_table,&coord, GID,T,&twist);    // [p,T]
        retrac_temp += lattice_retrace_plaquette2_F(&m1,&m2,&m3,&m6,&Fxt3,&Fxt8); // x-t: [p,X]-[p+X,T]-[p+T,X]*-[p,T]*

            m2 = lattice_table_2(lattice_table,&coordY,gdiY,Z,&twist);   // [p+Y,Z]
            m3 = lattice_table_2(lattice_table,&coordZ,gdiZ,Y,&twist);   // [p+Z,Y]
        retrac_spat += lattice_retrace_plaquette2(&m4,&m2,&m3,&m5);      // y-z: [p,Y]-[p+Y,Z]-[p+Z,Y]*-[p,Z]*

            m2 = lattice_table_2(lattice_table,&coordY,gdiY,T,&twist);   // [p+Y,T]
            m3 = lattice_table_2(lattice_table,&coordT,gdiT,Y,&twist);   // [p+T,Y]
        retrac_temp += lattice_retrace_plaquette2_F(&m4,&m2,&m3,&m6,&Fyt3,&Fyt8); // y-t: [p,Y]-[p+Y,T]-[p+T,Y]*-[p,T]*

            m2 = lattice_table_2(lattice_table,&coordZ,gdiZ,T,&twist);   // [p+Z,T]
            m3 = lattice_table_2(lattice_table,&coordT,gdiT,Z,&twist);   // [p+T,Z]
        retrac_temp += lattice_retrace_plaquette2_F(&m5,&m2,&m3,&m6,&Fzt3,&Fzt8); // z-t: [p,Z]-[p+Z,T]-[p+T,Z]*-[p,T]*

        out      = (hgpu_double2) {retrac_spat,retrac_temp};
        out_Fxt3 = (hgpu_double2) {Fxt3.re,Fxt3.im};
        out_Fyt3 = (hgpu_double2) {Fyt3.re,Fyt3.im};
        out_Fzt3 = (hgpu_double2) {Fzt3.re,Fzt3.im};
        out_Fxt8 = (hgpu_double2) {Fxt8.re,Fxt8.im};
        out_Fyt8 = (hgpu_double2) {Fyt8.re,Fyt8.im};
        out_Fzt8 = (hgpu_double2) {Fzt8.re,Fzt8.im};
    }

    // first reduction
    reduce_first_step_val_double2(lattice_lds,&out, &out2);
    if(TID == 0) lattice_measurement[BID             ] = out2;

    // first reduction - Fxt3
    reduce_first_step_val_double2(lattice_lds,&out_Fxt3, &Fxt3_reduced);
    if(TID == 0) lattice_measurement[BID +     offset] = Fxt3_reduced;

    // first reduction - Fyt3
    reduce_first_step_val_double2(lattice_lds,&out_Fyt3, &Fyt3_reduced);
    if(TID == 0) lattice_measurement[BID + 2 * offset] = Fyt3_reduced;

    // first reduction - Fzt3
    reduce_first_step_val_double2(lattice_lds,&out_Fzt3, &Fzt3_reduced);
    if(TID == 0) lattice_measurement[BID + 3 * offset] = Fzt3_reduced;

    // first reduction - Fxt8
    reduce_first_step_val_double2(lattice_lds,&out_Fxt8, &Fxt8_reduced);
    if(TID == 0) lattice_measurement[BID + 4 * offset] = Fxt8_reduced;

    // first reduction - Fyt8
    reduce_first_step_val_double2(lattice_lds,&out_Fyt8, &Fyt8_reduced);
    if(TID == 0) lattice_measurement[BID + 5 * offset] = Fyt8_reduced;

    // first reduction - Fzt8
    reduce_first_step_val_double2(lattice_lds,&out_Fzt8, &Fzt8_reduced);
    if(TID == 0) lattice_measurement[BID + 6 * offset] = Fzt8_reduced;
#endif

#if SUN == 3
    hgpu_double2 out  = (hgpu_double2) 0.0;
    hgpu_double2 out2 = (hgpu_double2) 0.0;

    hgpu_double2 Fxt3_reduced = (hgpu_double2) 0.0;
    hgpu_double2 Fyt3_reduced = (hgpu_double2) 0.0;
    hgpu_double2 Fzt3_reduced = (hgpu_double2) 0.0;
    hgpu_double2 Fxt8_reduced = (hgpu_double2) 0.0;
    hgpu_double2 Fyt8_reduced = (hgpu_double2) 0.0;
    hgpu_double2 Fzt8_reduced = (hgpu_double2) 0.0;

    hgpu_double2 out_Fxt3 = (hgpu_double2) 0.0;
    hgpu_double2 out_Fyt3 = (hgpu_double2) 0.0;
    hgpu_double2 out_Fzt3 = (hgpu_double2) 0.0;
    hgpu_double2 out_Fxt8 = (hgpu_double2) 0.0;
    hgpu_double2 out_Fyt8 = (hgpu_double2) 0.0;
    hgpu_double2 out_Fzt8 = (hgpu_double2) 0.0;

    hgpu_double retrac_spat = 0.0;
    hgpu_double retrac_temp = 0.0;

    uint gindex = GID;
    su3_twist twist;
    twist.phi   = lattice_parameters[1];
    twist.omega = lattice_parameters[2];
    hgpu_complex_double Fxt3, Fyt3, Fzt3;
    hgpu_complex_double Fxt8, Fyt8, Fzt8;

    lattice_lds[TID] = (hgpu_double2) 0.0;
    if (GID<SITES) {
        coords_4 coord;
        coords_4 coordX,coordY,coordZ,coordT;
        uint gdiX,gdiY,gdiZ,gdiT;

        gpu_su_3 m1,m2,m3,m4,m5,m6;

        lattice_gid_to_coords(&gindex,&coord);

        // prepare neighbours
        lattice_neighbours_gid(&coord,&coordX,&gdiX,X);
        lattice_neighbours_gid(&coord,&coordY,&gdiY,Y);
        lattice_neighbours_gid(&coord,&coordZ,&gdiZ,Z);
        lattice_neighbours_gid(&coord,&coordT,&gdiT,T);

            m1 = lattice_table_3(lattice_table,&coord, GID,X,&twist);    // [p,x]
            m2 = lattice_table_3(lattice_table,&coordX,gdiX,Y,&twist);   // [p+X,Y]
            m3 = lattice_table_3(lattice_table,&coordY,gdiY,X,&twist);   // [p+Y,X]
            m4 = lattice_table_3(lattice_table,&coord, GID,Y,&twist);    // [p,Y]
        retrac_spat += lattice_retrace_plaquette3(&m1,&m2,&m3,&m4);       // x-y: [p,X]-[p+X,Y]-[p+Y,X]*-[p,Y]

            m2 = lattice_table_3(lattice_table,&coordX,gdiX,Z,&twist);   // [p+X,Z]
            m3 = lattice_table_3(lattice_table,&coordZ,gdiZ,X,&twist);   // [p+Z,X]
            m5 = lattice_table_3(lattice_table,&coord, GID,Z,&twist);    // [p,Z]
        retrac_spat += lattice_retrace_plaquette3(&m1,&m2,&m3,&m5);      // x-z: [p,X]-[p+X,Z]-[p+Z,X]*-[p,Z]*

            m2 = lattice_table_3(lattice_table,&coordX,gdiX,T,&twist);   // [p+X,T]
            m3 = lattice_table_3(lattice_table,&coordT,gdiT,X,&twist);   // [p+T,X]
            m6 = lattice_table_3(lattice_table,&coord, GID,T,&twist);    // [p,T]
        retrac_temp += lattice_retrace_plaquette3_F(&m1,&m2,&m3,&m6,&Fxt3,&Fxt8); // x-t: [p,X]-[p+X,T]-[p+T,X]*-[p,T]*

            m2 = lattice_table_3(lattice_table,&coordY,gdiY,Z,&twist);   // [p+Y,Z]
            m3 = lattice_table_3(lattice_table,&coordZ,gdiZ,Y,&twist);   // [p+Z,Y]
        retrac_spat += lattice_retrace_plaquette3(&m4,&m2,&m3,&m5);      // y-z: [p,Y]-[p+Y,Z]-[p+Z,Y]*-[p,Z]*

            m2 = lattice_table_3(lattice_table,&coordY,gdiY,T,&twist);   // [p+Y,T]
            m3 = lattice_table_3(lattice_table,&coordT,gdiT,Y,&twist);   // [p+T,Y]
        retrac_temp += lattice_retrace_plaquette3_F(&m4,&m2,&m3,&m6,&Fyt3,&Fyt8); // y-t: [p,Y]-[p+Y,T]-[p+T,Y]*-[p,T]*

            m2 = lattice_table_3(lattice_table,&coordZ,gdiZ,T,&twist);   // [p+Z,T]
            m3 = lattice_table_3(lattice_table,&coordT,gdiT,Z,&twist);   // [p+T,Z]
        retrac_temp += lattice_retrace_plaquette3_F(&m5,&m2,&m3,&m6,&Fzt3,&Fzt8); // z-t: [p,Z]-[p+Z,T]-[p+T,Z]*-[p,T]*

        out      = (hgpu_double2) {retrac_spat,retrac_temp};
        out_Fxt3 = (hgpu_double2) {Fxt3.re,Fxt3.im};
        out_Fyt3 = (hgpu_double2) {Fyt3.re,Fyt3.im};
        out_Fzt3 = (hgpu_double2) {Fzt3.re,Fzt3.im};
        out_Fxt8 = (hgpu_double2) {Fxt8.re,Fxt8.im};
        out_Fyt8 = (hgpu_double2) {Fyt8.re,Fyt8.im};
        out_Fzt8 = (hgpu_double2) {Fzt8.re,Fzt8.im};
    }

    // first reduction
    reduce_first_step_val_double2(lattice_lds,&out, &out2);
    if(TID == 0) lattice_measurement[BID             ] = out2;

    // first reduction - Fxt3
    reduce_first_step_val_double2(lattice_lds,&out_Fxt3, &Fxt3_reduced);
    if(TID == 0) lattice_measurement[BID +     offset] = Fxt3_reduced;

    // first reduction - Fyt3
    reduce_first_step_val_double2(lattice_lds,&out_Fyt3, &Fyt3_reduced);
    if(TID == 0) lattice_measurement[BID + 2 * offset] = Fyt3_reduced;

    // first reduction - Fzt3
    reduce_first_step_val_double2(lattice_lds,&out_Fzt3, &Fzt3_reduced);
    if(TID == 0) lattice_measurement[BID + 3 * offset] = Fzt3_reduced;

    // first reduction - Fxt8
    reduce_first_step_val_double2(lattice_lds,&out_Fxt8, &Fxt8_reduced);
    if(TID == 0) lattice_measurement[BID + 4 * offset] = Fxt8_reduced;

    // first reduction - Fyt8
    reduce_first_step_val_double2(lattice_lds,&out_Fyt8, &Fyt8_reduced);
    if(TID == 0) lattice_measurement[BID + 5 * offset] = Fyt8_reduced;

    // first reduction - Fzt8
    reduce_first_step_val_double2(lattice_lds,&out_Fzt8, &Fzt8_reduced);
    if(TID == 0) lattice_measurement[BID + 6 * offset] = Fzt8_reduced;
#endif
}

#else
// no Fmunu measurement
                                        __kernel void
lattice_measurement_plq(__global hgpu_float4  * lattice_table,
                        __global hgpu_double2 * lattice_measurement,
                        __global hgpu_float   * lattice_parameters,
                        __local hgpu_double2  * lattice_lds,
                        uint offset)
{
#if SUN == 2
    hgpu_double2 out  = (hgpu_double2) 0.0;
    hgpu_double2 out2 = (hgpu_double2) 0.0;

    uint gindex = GID;
    su2_twist twist;
    twist.phi   = lattice_parameters[1];

    hgpu_double retrac_spat = 0.0;
    hgpu_double retrac_temp = 0.0;

    lattice_lds[TID] = (hgpu_double2) 0.0;
    if (GID<SITES) {
        coords_4 coord;
        coords_4 coordX,coordY,coordZ,coordT;
        uint gdiX,gdiY,gdiZ,gdiT;

        gpu_su_2 m1,m2,m3,m4,m5,m6;

        lattice_gid_to_coords(&gindex,&coord);

        // prepare neighbours
        lattice_neighbours_gid(&coord,&coordX,&gdiX,X);
        lattice_neighbours_gid(&coord,&coordY,&gdiY,Y);
        lattice_neighbours_gid(&coord,&coordZ,&gdiZ,Z);
        lattice_neighbours_gid(&coord,&coordT,&gdiT,T);

            m1 = lattice_table_2(lattice_table,&coord, GID,X,&twist);    // [p,x]
            m2 = lattice_table_2(lattice_table,&coordX,gdiX,Y,&twist);   // [p+X,Y]
            m3 = lattice_table_2(lattice_table,&coordY,gdiY,X,&twist);   // [p+Y,X]
            m4 = lattice_table_2(lattice_table,&coord, GID,Y,&twist);    // [p,Y]
        retrac_spat = lattice_retrace_plaquette2(&m1,&m2,&m3,&m4); // x-y: [p,X]-[p+X,Y]-[p+Y,X]*-[p,Y]

            m2 = lattice_table_2(lattice_table,&coordX,gdiX,Z,&twist);   // [p+X,Z]
            m3 = lattice_table_2(lattice_table,&coordZ,gdiZ,X,&twist);   // [p+Z,X]
            m5 = lattice_table_2(lattice_table,&coord, GID,Z,&twist);    // [p,Z]
        retrac_spat += lattice_retrace_plaquette2(&m1,&m2,&m3,&m5); // x-z: [p,X]-[p+X,Z]-[p+Z,X]*-[p,Z]*

            m2 = lattice_table_2(lattice_table,&coordX,gdiX,T,&twist);   // [p+X,T]
            m3 = lattice_table_2(lattice_table,&coordT,gdiT,X,&twist);   // [p+T,X]
            m6 = lattice_table_2(lattice_table,&coord, GID,T,&twist);    // [p,T]
        retrac_temp = lattice_retrace_plaquette2(&m1,&m2,&m3,&m6); // x-t: [p,X]-[p+X,T]-[p+T,X]*-[p,T]*

            m2 = lattice_table_2(lattice_table,&coordY,gdiY,Z,&twist);   // [p+Y,Z]
            m3 = lattice_table_2(lattice_table,&coordZ,gdiZ,Y,&twist);   // [p+Z,Y]
        retrac_spat += lattice_retrace_plaquette2(&m4,&m2,&m3,&m5); // y-z: [p,Y]-[p+Y,Z]-[p+Z,Y]*-[p,Z]*

            m2 = lattice_table_2(lattice_table,&coordY,gdiY,T,&twist);   // [p+Y,T]
            m3 = lattice_table_2(lattice_table,&coordT,gdiT,Y,&twist);   // [p+T,Y]
        retrac_temp += lattice_retrace_plaquette2(&m4,&m2,&m3,&m6); // y-t: [p,Y]-[p+Y,T]-[p+T,Y]*-[p,T]*

            m2 = lattice_table_2(lattice_table,&coordZ,gdiZ,T,&twist);   // [p+Z,T]
            m3 = lattice_table_2(lattice_table,&coordT,gdiT,Z,&twist);   // [p+T,Z]
        retrac_temp += lattice_retrace_plaquette2(&m5,&m2,&m3,&m6); // z-t: [p,Z]-[p+Z,T]-[p+T,Z]*-[p,T]*

        // first reduction
        out      = (hgpu_double2) (retrac_spat,retrac_temp);
    }
    reduce_first_step_val_double2(lattice_lds,&out, &out2);
    if(TID == 0) lattice_measurement[BID             ] = out2;
#endif

#if SUN == 3
    hgpu_double2 out  = (hgpu_double2) 0.0;
    hgpu_double2 out2 = (hgpu_double2) 0.0;

    uint gindex = GID;
    su3_twist twist;
    twist.phi   = lattice_parameters[1];
    twist.omega = lattice_parameters[2];

    hgpu_double retrac_spat = 0.0;
    hgpu_double retrac_temp = 0.0;

    lattice_lds[TID] = (hgpu_double2) 0.0;
    if (GID<SITES) {
        coords_4 coord;
        coords_4 coordX,coordY,coordZ,coordT;
        uint gdiX,gdiY,gdiZ,gdiT;

        gpu_su_3 m1,m2,m3,m4,m5,m6;

        lattice_gid_to_coords(&gindex,&coord);

        // prepare neighbours
        lattice_neighbours_gid(&coord,&coordX,&gdiX,X);
        lattice_neighbours_gid(&coord,&coordY,&gdiY,Y);
        lattice_neighbours_gid(&coord,&coordZ,&gdiZ,Z);
        lattice_neighbours_gid(&coord,&coordT,&gdiT,T);

            m1 = lattice_table_3(lattice_table,&coord, GID,X,&twist);    // [p,x]
            m2 = lattice_table_3(lattice_table,&coordX,gdiX,Y,&twist);   // [p+X,Y]
            m3 = lattice_table_3(lattice_table,&coordY,gdiY,X,&twist);   // [p+Y,X]
            m4 = lattice_table_3(lattice_table,&coord, GID,Y,&twist);    // [p,Y]
        retrac_spat = lattice_retrace_plaquette3(&m1,&m2,&m3,&m4); // x-y: [p,X]-[p+X,Y]-[p+Y,X]*-[p,Y]

            m2 = lattice_table_3(lattice_table,&coordX,gdiX,Z,&twist);   // [p+X,Z]
            m3 = lattice_table_3(lattice_table,&coordZ,gdiZ,X,&twist);   // [p+Z,X]
            m5 = lattice_table_3(lattice_table,&coord, GID,Z,&twist);    // [p,Z]
        retrac_spat += lattice_retrace_plaquette3(&m1,&m2,&m3,&m5); // x-z: [p,X]-[p+X,Z]-[p+Z,X]*-[p,Z]*

            m2 = lattice_table_3(lattice_table,&coordX,gdiX,T,&twist);   // [p+X,T]
            m3 = lattice_table_3(lattice_table,&coordT,gdiT,X,&twist);   // [p+T,X]
            m6 = lattice_table_3(lattice_table,&coord, GID,T,&twist);    // [p,T]
        retrac_temp = lattice_retrace_plaquette3(&m1,&m2,&m3,&m6); // x-t: [p,X]-[p+X,T]-[p+T,X]*-[p,T]*

            m2 = lattice_table_3(lattice_table,&coordY,gdiY,Z,&twist);   // [p+Y,Z]
            m3 = lattice_table_3(lattice_table,&coordZ,gdiZ,Y,&twist);   // [p+Z,Y]
        retrac_spat += lattice_retrace_plaquette3(&m4,&m2,&m3,&m5); // y-z: [p,Y]-[p+Y,Z]-[p+Z,Y]*-[p,Z]*

            m2 = lattice_table_3(lattice_table,&coordY,gdiY,T,&twist);   // [p+Y,T]
            m3 = lattice_table_3(lattice_table,&coordT,gdiT,Y,&twist);   // [p+T,Y]
        retrac_temp += lattice_retrace_plaquette3(&m4,&m2,&m3,&m6); // y-t: [p,Y]-[p+Y,T]-[p+T,Y]*-[p,T]*

            m2 = lattice_table_3(lattice_table,&coordZ,gdiZ,T,&twist);   // [p+Z,T]
            m3 = lattice_table_3(lattice_table,&coordT,gdiT,Z,&twist);   // [p+T,Z]
        retrac_temp += lattice_retrace_plaquette3(&m5,&m2,&m3,&m6); // z-t: [p,Z]-[p+Z,T]-[p+T,Z]*-[p,T]*

        // first reduction
        out      = (hgpu_double2) {retrac_spat,retrac_temp};
    }

    reduce_first_step_val_double2(lattice_lds,&out, &out2);
    if(TID == 0) lattice_measurement[BID             ] = out2;
#endif    
}

#endif

                               __kernel void
reduce_measurement_plq_double2(__global hgpu_double2 * lattice_measurement,
                               __global hgpu_double2 * lattice_energies_plq,
                               __local  hgpu_double2 * lattice_lds,
                               uint4 param,
                               uint index)
{
    uint size    = param.x;

    reduce_final_step_double2(lattice_lds,lattice_measurement,size);
    if (GID==0) lattice_energies_plq[index] = lattice_lds[0];
#if (defined(FMUNU) || defined(F0MU))
    uint offset  = param.y;
    uint offset2 = param.z;
    //______________________________ Fxy3
    uint idx = offset;
    reduce_final_step_double2_offset(lattice_lds,lattice_measurement,size,idx);
    if (GID==0) lattice_energies_plq[index +     offset2] = lattice_lds[0];

    //______________________________ Fxz3
    idx = 2 * offset;
    reduce_final_step_double2_offset(lattice_lds,lattice_measurement,size,idx);
    if (GID==0) lattice_energies_plq[index + 2 * offset2] = lattice_lds[0];

    //______________________________ Fyz3
    idx = 3 * offset;
    reduce_final_step_double2_offset(lattice_lds,lattice_measurement,size,idx);
    if (GID==0) lattice_energies_plq[index + 3 * offset2] = lattice_lds[0];

    //______________________________ Fxy8
    idx = 4 * offset;
    reduce_final_step_double2_offset(lattice_lds,lattice_measurement,size,idx);
    if (GID==0) lattice_energies_plq[index + 4 * offset2] = lattice_lds[0];

    //______________________________ Fxz8
    idx = 5 * offset;
    reduce_final_step_double2_offset(lattice_lds,lattice_measurement,size,idx);
    if (GID==0) lattice_energies_plq[index + 5 * offset2] = lattice_lds[0];

    //______________________________ Fyz8
    idx = 6 * offset;
    reduce_final_step_double2_offset(lattice_lds,lattice_measurement,size,idx);
    if (GID==0) lattice_energies_plq[index + 6 * offset2] = lattice_lds[0];
#endif
}

                                        __kernel void
clear_measurement(__global hgpu_double2 * lattice_measurement)
{
    lattice_measurement[GID] = (hgpu_double2) 0.0;
}

#endif
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                  
                                                                                                                                                                  
