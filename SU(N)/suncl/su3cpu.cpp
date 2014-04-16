/******************************************************************************
 * @file     suncpu.cpp
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.4
 *
 * @brief    [QCDGPU]
 *           Procedures for host simulations (SU(3) gauge group)
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

#include "suncl.h"
#include "su3cpu.h"

namespace SU3_CPU{
using SU3_CPU::SU;
using model_CL::model;

      char direction_X[] = "X";
      char direction_Y[] = "Y";
      char direction_Z[] = "Z";
      char direction_T[] = "T";
      
      char**    SU::directions   = NULL;
SU::SU(void)
{
        SU::directions = (char**) calloc(4,sizeof(char*));

    SU::directions[0] = (char*) calloc(strlen(direction_X) + 1,sizeof(char));
    SU::directions[1] = (char*) calloc(strlen(direction_Y) + 1,sizeof(char));
    SU::directions[2] = (char*) calloc(strlen(direction_Z) + 1,sizeof(char));
    SU::directions[3] = (char*) calloc(strlen(direction_T) + 1,sizeof(char));

        strcpy_s(SU::directions[0],(strlen(direction_X) + 1),direction_X);
        strcpy_s(SU::directions[1],(strlen(direction_Y) + 1),direction_Y);
        strcpy_s(SU::directions[2],(strlen(direction_Z) + 1),direction_Z);
        strcpy_s(SU::directions[3],(strlen(direction_T) + 1),direction_T);

    lattice_data = NULL;    // lattice data for CPU simulation
}

SU::~SU(void)
{
    for(int i=0; i<4; i++) free(directions[i]);
    free(directions);
}

void            SU::lattice_coords_print(SU::coords_4 coords){
    printf("(%u;%u;%u;%u)\n",coords.x,coords.y,coords.z,coords.t);
}

unsigned int    SU::lattice_coords_to_gid(model* lat,SU::coords_4 coords){
    unsigned int result = coords.y + coords.z * lat->lattice_domain_size[1] + coords.t * lat->lattice_domain_n2n3 + coords.x * lat->lattice_domain_n2n3n4;
    return result;
}

unsigned int    SU::lattice_coords_to_gid_half(model* lat,SU::coords_4 coords){
    unsigned int result = (coords.y + coords.z * lat->lattice_domain_size[1] + coords.t * lat->lattice_domain_n2n3 + coords.x * lat->lattice_domain_n2n3n4) / 2;
    return result;
}

SU::coords_4    SU::lattice_neighbours_coords(model* lat,const SU::coords_4 coord,int dir){
    SU::coords_4 coord2 = coord;
    switch (dir){
        case 0:
            coord2.x++; if (coord2.x>=lat->lattice_domain_n1) coord2.x = 0;
            break;
        case 1:
            coord2.y++; if (coord2.y>=lat->lattice_domain_size[1]) coord2.y = 0;
            break;
        case 2:
            coord2.z++; if (coord2.z>=lat->lattice_domain_size[2]) coord2.z = 0;
            break;
        case 3:
            coord2.t++; if (coord2.t>=lat->lattice_domain_size[3]) coord2.t = 0;
            break;

        default:
            break;
    }
    return coord2;
}

SU::coords_4    SU::lattice_neighbours_coords_backward(model* lat,const SU::coords_4 coord,int dir){
    coords_4 coord2 = coord;
    switch (dir){
        case 0:
            if (coord2.x<=0) {coord2.x = lat->lattice_domain_n1;} coord2.x--;
            break;
        case 1:
            if (coord2.y<=0) {coord2.y = lat->lattice_domain_size[1];} coord2.y--;
            break;
        case 2:
            if (coord2.z<=0) {coord2.z = lat->lattice_domain_size[2];} coord2.z--;
            break;
        case 3:
            if (coord2.t<=0) {coord2.t = lat->lattice_domain_size[3];} coord2.t--;
            break;

        default:
            break;
    }
    return coord2;
}

SU::coords_4    SU::lattice_gid_to_coords(model* lat,unsigned int gindex){
    coords_4 result;
        unsigned int z1,z2,z3,z4;

        z4 = gindex / lat->lattice_domain_n2n3n4;
        z1 = gindex - z4 * lat->lattice_domain_n2n3n4;
        z3 = z1 / lat->lattice_domain_n2n3;
        z1 = z1 - z3 * lat->lattice_domain_n2n3;
        z2 = z1 / lat->lattice_domain_size[1];
        z1 = z1 - z2 * lat->lattice_domain_size[1];

        result.x = z4;
        result.y = z1;
        result.z = z2;
        result.t = z3;

    return result;
}

SU::su_3        SU::lattice_matrix_reconstruct3(SU::su_3 a){
    su_3 result;

        result = a;

        result.w1.re =   a.u3.im * a.v2.im - a.u3.re * a.v2.re - a.u2.im * a.v3.im + a.u2.re * a.v3.re;
        result.w1.im =   a.u3.re * a.v2.im + a.u3.im * a.v2.re - a.u2.re * a.v3.im - a.u2.im * a.v3.re;
        result.w2.re =  -a.u3.im * a.v1.im + a.u3.re * a.v1.re + a.u1.im * a.v3.im - a.u1.re * a.v3.re;
        result.w2.im =  -a.u3.re * a.v1.im - a.u3.im * a.v1.re + a.u1.re * a.v3.im + a.u1.im * a.v3.re;
        result.w3.re =   a.u2.im * a.v1.im - a.u2.re * a.v1.re - a.u1.im * a.v2.im + a.u1.re * a.v2.re;
        result.w3.im =   a.u2.re * a.v1.im + a.u2.im * a.v1.re - a.u1.re * a.v2.im - a.u1.im * a.v2.re;

    return result;
}

SU::su_3        SU::lattice_table_3(model* lat,coords_4 coords,unsigned int gindex,int dir){
    SU::su_3 result;
    SU::su_3 m_omega;
    result = lattice_data[gindex + lat->lattice_table_row_size * dir];
    if ((dir==1) && (coords.z == (lat->lattice_domain_size[2]-1))) {
        double phi_p_omega_2   = 0.5 * (lat->PHI   + lat->OMEGA);
        double omega_m_phi_2   = 0.5 * (lat->OMEGA - lat->PHI);

        m_omega.u1.re = cos(phi_p_omega_2);
        m_omega.u1.im = sin(phi_p_omega_2);
        m_omega.u2.re = 0.0;
        m_omega.u2.im = 0.0;
        m_omega.u3.re = 0.0;
        m_omega.u3.im = 0.0;
        m_omega.v1.re = 0.0;
        m_omega.v1.im = 0.0;
        m_omega.v2.re = cos(omega_m_phi_2);
        m_omega.v2.im = sin(omega_m_phi_2);
        m_omega.v3.re = 0.0;
        m_omega.v3.im = 0.0;
        m_omega.w1.re = 0.0;
        m_omega.w1.im = 0.0;
        m_omega.w2.re = 0.0;
        m_omega.w2.im = 0.0;
        m_omega.w3.re = cos(-lat->OMEGA);
        m_omega.w3.im = sin(-lat->OMEGA);

        result = lattice_matrix_times3(m_omega,result);
    }
    return result;
}

void            SU::lattice_store_3(model* lat,SU::su_3 matrix,unsigned int gindex,int dir){
    lattice_data[gindex + lat->lattice_table_row_size * dir] = matrix;
}

SU::su_3        SU::lattice_get_3(model* lat,unsigned int * lattice_table,unsigned int gindex,int dir){
    SU::su_3 result,result2;

  if (lat->precision == model::model_precision_single) {
    switch (dir){
        case 0:
            result.u1.re = GPU0->convert_to_double(lattice_table[gindex*4     + 0 * lat->rowsize4]);
            result.u2.re = GPU0->convert_to_double(lattice_table[gindex*4 + 1 + 0 * lat->rowsize4]);
            result.u3.re = GPU0->convert_to_double(lattice_table[gindex*4 + 2 + 0 * lat->rowsize4]);
            result.v3.re = GPU0->convert_to_double(lattice_table[gindex*4 + 3 + 0 * lat->rowsize4]);

            result.u1.im = GPU0->convert_to_double(lattice_table[gindex*4     + 4 * lat->rowsize4]);
            result.u2.im = GPU0->convert_to_double(lattice_table[gindex*4 + 1 + 4 * lat->rowsize4]);
            result.u3.im = GPU0->convert_to_double(lattice_table[gindex*4 + 2 + 4 * lat->rowsize4]);
            result.v3.im = GPU0->convert_to_double(lattice_table[gindex*4 + 3 + 4 * lat->rowsize4]);

            result.v1.re = GPU0->convert_to_double(lattice_table[gindex*4     + 8 * lat->rowsize4]);
            result.v2.re = GPU0->convert_to_double(lattice_table[gindex*4 + 1 + 8 * lat->rowsize4]);
            result.v1.im = GPU0->convert_to_double(lattice_table[gindex*4 + 2 + 8 * lat->rowsize4]);
            result.v2.im = GPU0->convert_to_double(lattice_table[gindex*4 + 3 + 8 * lat->rowsize4]);
            break;

        case 1:
            result.u1.re = GPU0->convert_to_double(lattice_table[gindex*4     + 1 * lat->rowsize4]);
            result.u2.re = GPU0->convert_to_double(lattice_table[gindex*4 + 1 + 1 * lat->rowsize4]);
            result.u3.re = GPU0->convert_to_double(lattice_table[gindex*4 + 2 + 1 * lat->rowsize4]);
            result.v3.re = GPU0->convert_to_double(lattice_table[gindex*4 + 3 + 1 * lat->rowsize4]);

            result.u1.im = GPU0->convert_to_double(lattice_table[gindex*4     + 5 * lat->rowsize4]);
            result.u2.im = GPU0->convert_to_double(lattice_table[gindex*4 + 1 + 5 * lat->rowsize4]);
            result.u3.im = GPU0->convert_to_double(lattice_table[gindex*4 + 2 + 5 * lat->rowsize4]);
            result.v3.im = GPU0->convert_to_double(lattice_table[gindex*4 + 3 + 5 * lat->rowsize4]);

            result.v1.re = GPU0->convert_to_double(lattice_table[gindex*4     + 9 * lat->rowsize4]);
            result.v2.re = GPU0->convert_to_double(lattice_table[gindex*4 + 1 + 9 * lat->rowsize4]);
            result.v1.im = GPU0->convert_to_double(lattice_table[gindex*4 + 2 + 9 * lat->rowsize4]);
            result.v2.im = GPU0->convert_to_double(lattice_table[gindex*4 + 3 + 9 * lat->rowsize4]);
            break;

        case 2:
            result.u1.re = GPU0->convert_to_double(lattice_table[gindex*4     + 2 * lat->rowsize4]);
            result.u2.re = GPU0->convert_to_double(lattice_table[gindex*4 + 1 + 2 * lat->rowsize4]);
            result.u3.re = GPU0->convert_to_double(lattice_table[gindex*4 + 2 + 2 * lat->rowsize4]);
            result.v3.re = GPU0->convert_to_double(lattice_table[gindex*4 + 3 + 2 * lat->rowsize4]);

            result.u1.im = GPU0->convert_to_double(lattice_table[gindex*4     + 6 * lat->rowsize4]);
            result.u2.im = GPU0->convert_to_double(lattice_table[gindex*4 + 1 + 6 * lat->rowsize4]);
            result.u3.im = GPU0->convert_to_double(lattice_table[gindex*4 + 2 + 6 * lat->rowsize4]);
            result.v3.im = GPU0->convert_to_double(lattice_table[gindex*4 + 3 + 6 * lat->rowsize4]);

            result.v1.re = GPU0->convert_to_double(lattice_table[gindex*4     + 10 * lat->rowsize4]);
            result.v2.re = GPU0->convert_to_double(lattice_table[gindex*4 + 1 + 10 * lat->rowsize4]);
            result.v1.im = GPU0->convert_to_double(lattice_table[gindex*4 + 2 + 10 * lat->rowsize4]);
            result.v2.im = GPU0->convert_to_double(lattice_table[gindex*4 + 3 + 10 * lat->rowsize4]);
            break;

        case 3:
            result.u1.re = GPU0->convert_to_double(lattice_table[gindex*4     + 3 * lat->rowsize4]);
            result.u2.re = GPU0->convert_to_double(lattice_table[gindex*4 + 1 + 3 * lat->rowsize4]);
            result.u3.re = GPU0->convert_to_double(lattice_table[gindex*4 + 2 + 3 * lat->rowsize4]);
            result.v3.re = GPU0->convert_to_double(lattice_table[gindex*4 + 3 + 3 * lat->rowsize4]);

            result.u1.im = GPU0->convert_to_double(lattice_table[gindex*4     + 7 * lat->rowsize4]);
            result.u2.im = GPU0->convert_to_double(lattice_table[gindex*4 + 1 + 7 * lat->rowsize4]);
            result.u3.im = GPU0->convert_to_double(lattice_table[gindex*4 + 2 + 7 * lat->rowsize4]);
            result.v3.im = GPU0->convert_to_double(lattice_table[gindex*4 + 3 + 7 * lat->rowsize4]);

            result.v1.re = GPU0->convert_to_double(lattice_table[gindex*4     + 11 * lat->rowsize4]);
            result.v2.re = GPU0->convert_to_double(lattice_table[gindex*4 + 1 + 11 * lat->rowsize4]);
            result.v1.im = GPU0->convert_to_double(lattice_table[gindex*4 + 2 + 11 * lat->rowsize4]);
            result.v2.im = GPU0->convert_to_double(lattice_table[gindex*4 + 3 + 11 * lat->rowsize4]);
            break;

        default:
            result = lattice_unity3();
            break;
    }
  } else {
    switch (dir){
        case 0:
            result.u1.re = GPU0->convert_to_double(lattice_table[gindex*8 + 0 +  0 * lat->rowsize4],lattice_table[gindex*8 + 1 +  0 * lat->rowsize4]);
            result.u2.re = GPU0->convert_to_double(lattice_table[gindex*8 + 2 +  0 * lat->rowsize4],lattice_table[gindex*8 + 3 +  0 * lat->rowsize4]);
            result.u3.re = GPU0->convert_to_double(lattice_table[gindex*8 + 4 +  0 * lat->rowsize4],lattice_table[gindex*8 + 5 +  0 * lat->rowsize4]);
            result.v3.re = GPU0->convert_to_double(lattice_table[gindex*8 + 6 +  0 * lat->rowsize4],lattice_table[gindex*8 + 7 +  0 * lat->rowsize4]);

            result.u1.im = GPU0->convert_to_double(lattice_table[gindex*8 + 0 +  8 * lat->rowsize4],lattice_table[gindex*8 + 1 +  8 * lat->rowsize4]);
            result.u2.im = GPU0->convert_to_double(lattice_table[gindex*8 + 2 +  8 * lat->rowsize4],lattice_table[gindex*8 + 3 +  8 * lat->rowsize4]);
            result.u3.im = GPU0->convert_to_double(lattice_table[gindex*8 + 4 +  8 * lat->rowsize4],lattice_table[gindex*8 + 5 +  8 * lat->rowsize4]);
            result.v3.im = GPU0->convert_to_double(lattice_table[gindex*8 + 6 +  8 * lat->rowsize4],lattice_table[gindex*8 + 7 +  8 * lat->rowsize4]);

            result.v1.re = GPU0->convert_to_double(lattice_table[gindex*8 + 0 + 16 * lat->rowsize4],lattice_table[gindex*8 + 1 + 16 * lat->rowsize4]);
            result.v2.re = GPU0->convert_to_double(lattice_table[gindex*8 + 2 + 16 * lat->rowsize4],lattice_table[gindex*8 + 3 + 16 * lat->rowsize4]);
            result.v1.im = GPU0->convert_to_double(lattice_table[gindex*8 + 4 + 16 * lat->rowsize4],lattice_table[gindex*8 + 5 + 16 * lat->rowsize4]);
            result.v2.im = GPU0->convert_to_double(lattice_table[gindex*8 + 6 + 16 * lat->rowsize4],lattice_table[gindex*8 + 7 + 16 * lat->rowsize4]);
            break;

        case 1:
            result.u1.re = GPU0->convert_to_double(lattice_table[gindex*8 + 0 +  2 * lat->rowsize4],lattice_table[gindex*8 + 1 +  2 * lat->rowsize4]);
            result.u2.re = GPU0->convert_to_double(lattice_table[gindex*8 + 2 +  2 * lat->rowsize4],lattice_table[gindex*8 + 3 +  2 * lat->rowsize4]);
            result.u3.re = GPU0->convert_to_double(lattice_table[gindex*8 + 4 +  2 * lat->rowsize4],lattice_table[gindex*8 + 5 +  2 * lat->rowsize4]);
            result.v3.re = GPU0->convert_to_double(lattice_table[gindex*8 + 6 +  2 * lat->rowsize4],lattice_table[gindex*8 + 7 +  2 * lat->rowsize4]);

            result.u1.im = GPU0->convert_to_double(lattice_table[gindex*8 + 0 + 10 * lat->rowsize4],lattice_table[gindex*8 + 1 + 10 * lat->rowsize4]);
            result.u2.im = GPU0->convert_to_double(lattice_table[gindex*8 + 2 + 10 * lat->rowsize4],lattice_table[gindex*8 + 3 + 10 * lat->rowsize4]);
            result.u3.im = GPU0->convert_to_double(lattice_table[gindex*8 + 4 + 10 * lat->rowsize4],lattice_table[gindex*8 + 5 + 10 * lat->rowsize4]);
            result.v3.im = GPU0->convert_to_double(lattice_table[gindex*8 + 6 + 10 * lat->rowsize4],lattice_table[gindex*8 + 7 + 10 * lat->rowsize4]);

            result.v1.re = GPU0->convert_to_double(lattice_table[gindex*8 + 0 + 18 * lat->rowsize4],lattice_table[gindex*8 + 1 + 18 * lat->rowsize4]);
            result.v2.re = GPU0->convert_to_double(lattice_table[gindex*8 + 2 + 18 * lat->rowsize4],lattice_table[gindex*8 + 3 + 18 * lat->rowsize4]);
            result.v1.im = GPU0->convert_to_double(lattice_table[gindex*8 + 4 + 18 * lat->rowsize4],lattice_table[gindex*8 + 5 + 18 * lat->rowsize4]);
            result.v2.im = GPU0->convert_to_double(lattice_table[gindex*8 + 6 + 18 * lat->rowsize4],lattice_table[gindex*8 + 7 + 18 * lat->rowsize4]);
            break;

        case 2:
            result.u1.re = GPU0->convert_to_double(lattice_table[gindex*8 + 0 +  4 * lat->rowsize4],lattice_table[gindex*8 + 1 +  4 * lat->rowsize4]);
            result.u2.re = GPU0->convert_to_double(lattice_table[gindex*8 + 2 +  4 * lat->rowsize4],lattice_table[gindex*8 + 3 +  4 * lat->rowsize4]);
            result.u3.re = GPU0->convert_to_double(lattice_table[gindex*8 + 4 +  4 * lat->rowsize4],lattice_table[gindex*8 + 5 +  4 * lat->rowsize4]);
            result.v3.re = GPU0->convert_to_double(lattice_table[gindex*8 + 6 +  4 * lat->rowsize4],lattice_table[gindex*8 + 7 +  4 * lat->rowsize4]);

            result.u1.im = GPU0->convert_to_double(lattice_table[gindex*8 + 0 + 12 * lat->rowsize4],lattice_table[gindex*8 + 1 + 12 * lat->rowsize4]);
            result.u2.im = GPU0->convert_to_double(lattice_table[gindex*8 + 2 + 12 * lat->rowsize4],lattice_table[gindex*8 + 3 + 12 * lat->rowsize4]);
            result.u3.im = GPU0->convert_to_double(lattice_table[gindex*8 + 4 + 12 * lat->rowsize4],lattice_table[gindex*8 + 5 + 12 * lat->rowsize4]);
            result.v3.im = GPU0->convert_to_double(lattice_table[gindex*8 + 6 + 12 * lat->rowsize4],lattice_table[gindex*8 + 7 + 12 * lat->rowsize4]);

            result.v1.re = GPU0->convert_to_double(lattice_table[gindex*8 + 0 + 20 * lat->rowsize4],lattice_table[gindex*8 + 1 + 20 * lat->rowsize4]);
            result.v2.re = GPU0->convert_to_double(lattice_table[gindex*8 + 2 + 20 * lat->rowsize4],lattice_table[gindex*8 + 3 + 20 * lat->rowsize4]);
            result.v1.im = GPU0->convert_to_double(lattice_table[gindex*8 + 4 + 20 * lat->rowsize4],lattice_table[gindex*8 + 5 + 20 * lat->rowsize4]);
            result.v2.im = GPU0->convert_to_double(lattice_table[gindex*8 + 6 + 20 * lat->rowsize4],lattice_table[gindex*8 + 7 + 20 * lat->rowsize4]);
            break;

        case 3:
            result.u1.re = GPU0->convert_to_double(lattice_table[gindex*8 + 0 +  6 * lat->rowsize4],lattice_table[gindex*8 + 1 +  6 * lat->rowsize4]);
            result.u2.re = GPU0->convert_to_double(lattice_table[gindex*8 + 2 +  6 * lat->rowsize4],lattice_table[gindex*8 + 3 +  6 * lat->rowsize4]);
            result.u3.re = GPU0->convert_to_double(lattice_table[gindex*8 + 4 +  6 * lat->rowsize4],lattice_table[gindex*8 + 5 +  6 * lat->rowsize4]);
            result.v3.re = GPU0->convert_to_double(lattice_table[gindex*8 + 6 +  6 * lat->rowsize4],lattice_table[gindex*8 + 7 +  6 * lat->rowsize4]);

            result.u1.im = GPU0->convert_to_double(lattice_table[gindex*8 + 0 + 14 * lat->rowsize4],lattice_table[gindex*8 + 1 + 14 * lat->rowsize4]);
            result.u2.im = GPU0->convert_to_double(lattice_table[gindex*8 + 2 + 14 * lat->rowsize4],lattice_table[gindex*8 + 3 + 14 * lat->rowsize4]);
            result.u3.im = GPU0->convert_to_double(lattice_table[gindex*8 + 4 + 14 * lat->rowsize4],lattice_table[gindex*8 + 5 + 14 * lat->rowsize4]);
            result.v3.im = GPU0->convert_to_double(lattice_table[gindex*8 + 6 + 14 * lat->rowsize4],lattice_table[gindex*8 + 7 + 14 * lat->rowsize4]);

            result.v1.re = GPU0->convert_to_double(lattice_table[gindex*8 + 0 + 22 * lat->rowsize4],lattice_table[gindex*8 + 1 + 22 * lat->rowsize4]);
            result.v2.re = GPU0->convert_to_double(lattice_table[gindex*8 + 2 + 22 * lat->rowsize4],lattice_table[gindex*8 + 3 + 22 * lat->rowsize4]);
            result.v1.im = GPU0->convert_to_double(lattice_table[gindex*8 + 4 + 22 * lat->rowsize4],lattice_table[gindex*8 + 5 + 22 * lat->rowsize4]);
            result.v2.im = GPU0->convert_to_double(lattice_table[gindex*8 + 6 + 22 * lat->rowsize4],lattice_table[gindex*8 + 7 + 22 * lat->rowsize4]);
            break;

        default:
            result = lattice_unity3();
            break;
    }
  }
    result2 = lattice_GramSchmidt_3(result);

    return result2;
}

SU::su_3        SU::lattice_get_raw_3(model* lat,unsigned int * lattice_table,unsigned int gindex,int dir){
    SU::su_3 result;

  if (lat->precision == model::model_precision_single) {
    switch (dir){
        case 0:
            result.u1.re = GPU0->convert_to_double(lattice_table[gindex*4     + 0 * lat->rowsize4]);
            result.u2.re = GPU0->convert_to_double(lattice_table[gindex*4 + 1 + 0 * lat->rowsize4]);
            result.u3.re = GPU0->convert_to_double(lattice_table[gindex*4 + 2 + 0 * lat->rowsize4]);
            result.v3.re = GPU0->convert_to_double(lattice_table[gindex*4 + 3 + 0 * lat->rowsize4]);

            result.u1.im = GPU0->convert_to_double(lattice_table[gindex*4     + 4 * lat->rowsize4]);
            result.u2.im = GPU0->convert_to_double(lattice_table[gindex*4 + 1 + 4 * lat->rowsize4]);
            result.u3.im = GPU0->convert_to_double(lattice_table[gindex*4 + 2 + 4 * lat->rowsize4]);
            result.v3.im = GPU0->convert_to_double(lattice_table[gindex*4 + 3 + 4 * lat->rowsize4]);

            result.v1.re = GPU0->convert_to_double(lattice_table[gindex*4     + 8 * lat->rowsize4]);
            result.v2.re = GPU0->convert_to_double(lattice_table[gindex*4 + 1 + 8 * lat->rowsize4]);
            result.v1.im = GPU0->convert_to_double(lattice_table[gindex*4 + 2 + 8 * lat->rowsize4]);
            result.v2.im = GPU0->convert_to_double(lattice_table[gindex*4 + 3 + 8 * lat->rowsize4]);
            break;

        case 1:
            result.u1.re = GPU0->convert_to_double(lattice_table[gindex*4     + 1 * lat->rowsize4]);
            result.u2.re = GPU0->convert_to_double(lattice_table[gindex*4 + 1 + 1 * lat->rowsize4]);
            result.u3.re = GPU0->convert_to_double(lattice_table[gindex*4 + 2 + 1 * lat->rowsize4]);
            result.v3.re = GPU0->convert_to_double(lattice_table[gindex*4 + 3 + 1 * lat->rowsize4]);

            result.u1.im = GPU0->convert_to_double(lattice_table[gindex*4     + 5 * lat->rowsize4]);
            result.u2.im = GPU0->convert_to_double(lattice_table[gindex*4 + 1 + 5 * lat->rowsize4]);
            result.u3.im = GPU0->convert_to_double(lattice_table[gindex*4 + 2 + 5 * lat->rowsize4]);
            result.v3.im = GPU0->convert_to_double(lattice_table[gindex*4 + 3 + 5 * lat->rowsize4]);

            result.v1.re = GPU0->convert_to_double(lattice_table[gindex*4     + 9 * lat->rowsize4]);
            result.v2.re = GPU0->convert_to_double(lattice_table[gindex*4 + 1 + 9 * lat->rowsize4]);
            result.v1.im = GPU0->convert_to_double(lattice_table[gindex*4 + 2 + 9 * lat->rowsize4]);
            result.v2.im = GPU0->convert_to_double(lattice_table[gindex*4 + 3 + 9 * lat->rowsize4]);
            break;

        case 2:
            result.u1.re = GPU0->convert_to_double(lattice_table[gindex*4     + 2 * lat->rowsize4]);
            result.u2.re = GPU0->convert_to_double(lattice_table[gindex*4 + 1 + 2 * lat->rowsize4]);
            result.u3.re = GPU0->convert_to_double(lattice_table[gindex*4 + 2 + 2 * lat->rowsize4]);
            result.v3.re = GPU0->convert_to_double(lattice_table[gindex*4 + 3 + 2 * lat->rowsize4]);

            result.u1.im = GPU0->convert_to_double(lattice_table[gindex*4     + 6 * lat->rowsize4]);
            result.u2.im = GPU0->convert_to_double(lattice_table[gindex*4 + 1 + 6 * lat->rowsize4]);
            result.u3.im = GPU0->convert_to_double(lattice_table[gindex*4 + 2 + 6 * lat->rowsize4]);
            result.v3.im = GPU0->convert_to_double(lattice_table[gindex*4 + 3 + 6 * lat->rowsize4]);

            result.v1.re = GPU0->convert_to_double(lattice_table[gindex*4     + 10 * lat->rowsize4]);
            result.v2.re = GPU0->convert_to_double(lattice_table[gindex*4 + 1 + 10 * lat->rowsize4]);
            result.v1.im = GPU0->convert_to_double(lattice_table[gindex*4 + 2 + 10 * lat->rowsize4]);
            result.v2.im = GPU0->convert_to_double(lattice_table[gindex*4 + 3 + 10 * lat->rowsize4]);
            break;

        case 3:
            result.u1.re = GPU0->convert_to_double(lattice_table[gindex*4     + 3 * lat->rowsize4]);
            result.u2.re = GPU0->convert_to_double(lattice_table[gindex*4 + 1 + 3 * lat->rowsize4]);
            result.u3.re = GPU0->convert_to_double(lattice_table[gindex*4 + 2 + 3 * lat->rowsize4]);
            result.v3.re = GPU0->convert_to_double(lattice_table[gindex*4 + 3 + 3 * lat->rowsize4]);

            result.u1.im = GPU0->convert_to_double(lattice_table[gindex*4     + 7 * lat->rowsize4]);
            result.u2.im = GPU0->convert_to_double(lattice_table[gindex*4 + 1 + 7 * lat->rowsize4]);
            result.u3.im = GPU0->convert_to_double(lattice_table[gindex*4 + 2 + 7 * lat->rowsize4]);
            result.v3.im = GPU0->convert_to_double(lattice_table[gindex*4 + 3 + 7 * lat->rowsize4]);

            result.v1.re = GPU0->convert_to_double(lattice_table[gindex*4     + 11 * lat->rowsize4]);
            result.v2.re = GPU0->convert_to_double(lattice_table[gindex*4 + 1 + 11 * lat->rowsize4]);
            result.v1.im = GPU0->convert_to_double(lattice_table[gindex*4 + 2 + 11 * lat->rowsize4]);
            result.v2.im = GPU0->convert_to_double(lattice_table[gindex*4 + 3 + 11 * lat->rowsize4]);
            break;

        default:
            break;
    }
  } else {
    switch (dir){
        case 0:
            result.u1.re = GPU0->convert_to_double(lattice_table[gindex*8 + 0 +  0 * lat->rowsize4],lattice_table[gindex*8 + 1 +  0 * lat->rowsize4]);
            result.u2.re = GPU0->convert_to_double(lattice_table[gindex*8 + 2 +  0 * lat->rowsize4],lattice_table[gindex*8 + 3 +  0 * lat->rowsize4]);
            result.u3.re = GPU0->convert_to_double(lattice_table[gindex*8 + 4 +  0 * lat->rowsize4],lattice_table[gindex*8 + 5 +  0 * lat->rowsize4]);
            result.v3.re = GPU0->convert_to_double(lattice_table[gindex*8 + 6 +  0 * lat->rowsize4],lattice_table[gindex*8 + 7 +  0 * lat->rowsize4]);

            result.u1.im = GPU0->convert_to_double(lattice_table[gindex*8 + 0 +  8 * lat->rowsize4],lattice_table[gindex*8 + 1 +  8 * lat->rowsize4]);
            result.u2.im = GPU0->convert_to_double(lattice_table[gindex*8 + 2 +  8 * lat->rowsize4],lattice_table[gindex*8 + 3 +  8 * lat->rowsize4]);
            result.u3.im = GPU0->convert_to_double(lattice_table[gindex*8 + 4 +  8 * lat->rowsize4],lattice_table[gindex*8 + 5 +  8 * lat->rowsize4]);
            result.v3.im = GPU0->convert_to_double(lattice_table[gindex*8 + 6 +  8 * lat->rowsize4],lattice_table[gindex*8 + 7 +  8 * lat->rowsize4]);

            result.v1.re = GPU0->convert_to_double(lattice_table[gindex*8 + 0 + 16 * lat->rowsize4],lattice_table[gindex*8 + 1 + 16 * lat->rowsize4]);
            result.v2.re = GPU0->convert_to_double(lattice_table[gindex*8 + 2 + 16 * lat->rowsize4],lattice_table[gindex*8 + 3 + 16 * lat->rowsize4]);
            result.v1.im = GPU0->convert_to_double(lattice_table[gindex*8 + 4 + 16 * lat->rowsize4],lattice_table[gindex*8 + 5 + 16 * lat->rowsize4]);
            result.v2.im = GPU0->convert_to_double(lattice_table[gindex*8 + 6 + 16 * lat->rowsize4],lattice_table[gindex*8 + 7 + 16 * lat->rowsize4]);
            break;

        case 1:
            result.u1.re = GPU0->convert_to_double(lattice_table[gindex*8 + 0 +  2 * lat->rowsize4],lattice_table[gindex*8 + 1 +  2 * lat->rowsize4]);
            result.u2.re = GPU0->convert_to_double(lattice_table[gindex*8 + 2 +  2 * lat->rowsize4],lattice_table[gindex*8 + 3 +  2 * lat->rowsize4]);
            result.u3.re = GPU0->convert_to_double(lattice_table[gindex*8 + 4 +  2 * lat->rowsize4],lattice_table[gindex*8 + 5 +  2 * lat->rowsize4]);
            result.v3.re = GPU0->convert_to_double(lattice_table[gindex*8 + 6 +  2 * lat->rowsize4],lattice_table[gindex*8 + 7 +  2 * lat->rowsize4]);

            result.u1.im = GPU0->convert_to_double(lattice_table[gindex*8 + 0 + 10 * lat->rowsize4],lattice_table[gindex*8 + 1 + 10 * lat->rowsize4]);
            result.u2.im = GPU0->convert_to_double(lattice_table[gindex*8 + 2 + 10 * lat->rowsize4],lattice_table[gindex*8 + 3 + 10 * lat->rowsize4]);
            result.u3.im = GPU0->convert_to_double(lattice_table[gindex*8 + 4 + 10 * lat->rowsize4],lattice_table[gindex*8 + 5 + 10 * lat->rowsize4]);
            result.v3.im = GPU0->convert_to_double(lattice_table[gindex*8 + 6 + 10 * lat->rowsize4],lattice_table[gindex*8 + 7 + 10 * lat->rowsize4]);

            result.v1.re = GPU0->convert_to_double(lattice_table[gindex*8 + 0 + 18 * lat->rowsize4],lattice_table[gindex*8 + 1 + 18 * lat->rowsize4]);
            result.v2.re = GPU0->convert_to_double(lattice_table[gindex*8 + 2 + 18 * lat->rowsize4],lattice_table[gindex*8 + 3 + 18 * lat->rowsize4]);
            result.v1.im = GPU0->convert_to_double(lattice_table[gindex*8 + 4 + 18 * lat->rowsize4],lattice_table[gindex*8 + 5 + 18 * lat->rowsize4]);
            result.v2.im = GPU0->convert_to_double(lattice_table[gindex*8 + 6 + 18 * lat->rowsize4],lattice_table[gindex*8 + 7 + 18 * lat->rowsize4]);
            break;

        case 2:
            result.u1.re = GPU0->convert_to_double(lattice_table[gindex*8 + 0 +  4 * lat->rowsize4],lattice_table[gindex*8 + 1 +  4 * lat->rowsize4]);
            result.u2.re = GPU0->convert_to_double(lattice_table[gindex*8 + 2 +  4 * lat->rowsize4],lattice_table[gindex*8 + 3 +  4 * lat->rowsize4]);
            result.u3.re = GPU0->convert_to_double(lattice_table[gindex*8 + 4 +  4 * lat->rowsize4],lattice_table[gindex*8 + 5 +  4 * lat->rowsize4]);
            result.v3.re = GPU0->convert_to_double(lattice_table[gindex*8 + 6 +  4 * lat->rowsize4],lattice_table[gindex*8 + 7 +  4 * lat->rowsize4]);

            result.u1.im = GPU0->convert_to_double(lattice_table[gindex*8 + 0 + 12 * lat->rowsize4],lattice_table[gindex*8 + 1 + 12 * lat->rowsize4]);
            result.u2.im = GPU0->convert_to_double(lattice_table[gindex*8 + 2 + 12 * lat->rowsize4],lattice_table[gindex*8 + 3 + 12 * lat->rowsize4]);
            result.u3.im = GPU0->convert_to_double(lattice_table[gindex*8 + 4 + 12 * lat->rowsize4],lattice_table[gindex*8 + 5 + 12 * lat->rowsize4]);
            result.v3.im = GPU0->convert_to_double(lattice_table[gindex*8 + 6 + 12 * lat->rowsize4],lattice_table[gindex*8 + 7 + 12 * lat->rowsize4]);

            result.v1.re = GPU0->convert_to_double(lattice_table[gindex*8 + 0 + 20 * lat->rowsize4],lattice_table[gindex*8 + 1 + 20 * lat->rowsize4]);
            result.v2.re = GPU0->convert_to_double(lattice_table[gindex*8 + 2 + 20 * lat->rowsize4],lattice_table[gindex*8 + 3 + 20 * lat->rowsize4]);
            result.v1.im = GPU0->convert_to_double(lattice_table[gindex*8 + 4 + 20 * lat->rowsize4],lattice_table[gindex*8 + 5 + 20 * lat->rowsize4]);
            result.v2.im = GPU0->convert_to_double(lattice_table[gindex*8 + 6 + 20 * lat->rowsize4],lattice_table[gindex*8 + 7 + 20 * lat->rowsize4]);
            break;

        case 3:
            result.u1.re = GPU0->convert_to_double(lattice_table[gindex*8 + 0 +  6 * lat->rowsize4],lattice_table[gindex*8 + 1 +  6 * lat->rowsize4]);
            result.u2.re = GPU0->convert_to_double(lattice_table[gindex*8 + 2 +  6 * lat->rowsize4],lattice_table[gindex*8 + 3 +  6 * lat->rowsize4]);
            result.u3.re = GPU0->convert_to_double(lattice_table[gindex*8 + 4 +  6 * lat->rowsize4],lattice_table[gindex*8 + 5 +  6 * lat->rowsize4]);
            result.v3.re = GPU0->convert_to_double(lattice_table[gindex*8 + 6 +  6 * lat->rowsize4],lattice_table[gindex*8 + 7 +  6 * lat->rowsize4]);

            result.u1.im = GPU0->convert_to_double(lattice_table[gindex*8 + 0 + 14 * lat->rowsize4],lattice_table[gindex*8 + 1 + 14 * lat->rowsize4]);
            result.u2.im = GPU0->convert_to_double(lattice_table[gindex*8 + 2 + 14 * lat->rowsize4],lattice_table[gindex*8 + 3 + 14 * lat->rowsize4]);
            result.u3.im = GPU0->convert_to_double(lattice_table[gindex*8 + 4 + 14 * lat->rowsize4],lattice_table[gindex*8 + 5 + 14 * lat->rowsize4]);
            result.v3.im = GPU0->convert_to_double(lattice_table[gindex*8 + 6 + 14 * lat->rowsize4],lattice_table[gindex*8 + 7 + 14 * lat->rowsize4]);

            result.v1.re = GPU0->convert_to_double(lattice_table[gindex*8 + 0 + 22 * lat->rowsize4],lattice_table[gindex*8 + 1 + 22 * lat->rowsize4]);
            result.v2.re = GPU0->convert_to_double(lattice_table[gindex*8 + 2 + 22 * lat->rowsize4],lattice_table[gindex*8 + 3 + 22 * lat->rowsize4]);
            result.v1.im = GPU0->convert_to_double(lattice_table[gindex*8 + 4 + 22 * lat->rowsize4],lattice_table[gindex*8 + 5 + 22 * lat->rowsize4]);
            result.v2.im = GPU0->convert_to_double(lattice_table[gindex*8 + 6 + 22 * lat->rowsize4],lattice_table[gindex*8 + 7 + 22 * lat->rowsize4]);
            break;

        default:
            break;
    }
  }

    result.w1.re = 0.0;
    result.w1.im = 0.0;
    result.w2.re = 0.0;
    result.w2.im = 0.0;
    result.w3.re = 0.0;
    result.w3.im = 0.0;

    return result;
}

SU::su_3        SU::lattice_staple_3(model* lat,unsigned int * lattice_staple,unsigned int gindex){
        SU::su_3 result;

  if (lat->precision == model::model_precision_single) {
            result.u1.re = GPU0->convert_to_double(lattice_staple[gindex*4     + 0 * lat->halfrowsize4]);
            result.u2.re = GPU0->convert_to_double(lattice_staple[gindex*4 + 1 + 0 * lat->halfrowsize4]);
            result.u3.re = GPU0->convert_to_double(lattice_staple[gindex*4 + 2 + 0 * lat->halfrowsize4]);
            result.v3.re = GPU0->convert_to_double(lattice_staple[gindex*4 + 3 + 0 * lat->halfrowsize4]);

            result.u1.im = GPU0->convert_to_double(lattice_staple[gindex*4     + 1 * lat->halfrowsize4]);
            result.u2.im = GPU0->convert_to_double(lattice_staple[gindex*4 + 1 + 1 * lat->halfrowsize4]);
            result.u3.im = GPU0->convert_to_double(lattice_staple[gindex*4 + 2 + 1 * lat->halfrowsize4]);
            result.v3.im = GPU0->convert_to_double(lattice_staple[gindex*4 + 3 + 1 * lat->halfrowsize4]);

            result.v1.re = GPU0->convert_to_double(lattice_staple[gindex*4     + 2 * lat->halfrowsize4]);
            result.v2.re = GPU0->convert_to_double(lattice_staple[gindex*4 + 1 + 2 * lat->halfrowsize4]);
            result.v1.im = GPU0->convert_to_double(lattice_staple[gindex*4 + 2 + 2 * lat->halfrowsize4]);
            result.v2.im = GPU0->convert_to_double(lattice_staple[gindex*4 + 3 + 2 * lat->halfrowsize4]);

            result.w1.re = GPU0->convert_to_double(lattice_staple[gindex*4     + 3 * lat->halfrowsize4]);
            result.w2.re = GPU0->convert_to_double(lattice_staple[gindex*4 + 1 + 3 * lat->halfrowsize4]);
            result.w3.re = GPU0->convert_to_double(lattice_staple[gindex*4 + 2 + 3 * lat->halfrowsize4]);

            result.w1.im = GPU0->convert_to_double(lattice_staple[gindex*4     + 4 * lat->halfrowsize4]);
            result.w2.im = GPU0->convert_to_double(lattice_staple[gindex*4 + 1 + 4 * lat->halfrowsize4]);
            result.w3.im = GPU0->convert_to_double(lattice_staple[gindex*4 + 2 + 4 * lat->halfrowsize4]);
  } else {
            result.u1.re = GPU0->convert_to_double(lattice_staple[gindex*8 + 0 + 0 * lat->halfrowsize4],lattice_staple[gindex*8 + 1 + 0 * lat->halfrowsize4]);
            result.u2.re = GPU0->convert_to_double(lattice_staple[gindex*8 + 2 + 0 * lat->halfrowsize4],lattice_staple[gindex*8 + 3 + 0 * lat->halfrowsize4]);
            result.u3.re = GPU0->convert_to_double(lattice_staple[gindex*8 + 4 + 0 * lat->halfrowsize4],lattice_staple[gindex*8 + 5 + 0 * lat->halfrowsize4]);
            result.v3.re = GPU0->convert_to_double(lattice_staple[gindex*8 + 6 + 0 * lat->halfrowsize4],lattice_staple[gindex*8 + 7 + 0 * lat->halfrowsize4]);

            result.u1.im = GPU0->convert_to_double(lattice_staple[gindex*8 + 0 + 2 * lat->halfrowsize4],lattice_staple[gindex*8 + 1 + 2 * lat->halfrowsize4]);
            result.u2.im = GPU0->convert_to_double(lattice_staple[gindex*8 + 2 + 2 * lat->halfrowsize4],lattice_staple[gindex*8 + 3 + 2 * lat->halfrowsize4]);
            result.u3.im = GPU0->convert_to_double(lattice_staple[gindex*8 + 4 + 2 * lat->halfrowsize4],lattice_staple[gindex*8 + 5 + 2 * lat->halfrowsize4]);
            result.v3.im = GPU0->convert_to_double(lattice_staple[gindex*8 + 6 + 2 * lat->halfrowsize4],lattice_staple[gindex*8 + 7 + 2 * lat->halfrowsize4]);

            result.v1.re = GPU0->convert_to_double(lattice_staple[gindex*8 + 0 + 4 * lat->halfrowsize4],lattice_staple[gindex*8 + 1 + 4 * lat->halfrowsize4]);
            result.v2.re = GPU0->convert_to_double(lattice_staple[gindex*8 + 2 + 4 * lat->halfrowsize4],lattice_staple[gindex*8 + 3 + 4 * lat->halfrowsize4]);
            result.v1.im = GPU0->convert_to_double(lattice_staple[gindex*8 + 4 + 4 * lat->halfrowsize4],lattice_staple[gindex*8 + 5 + 4 * lat->halfrowsize4]);
            result.v2.im = GPU0->convert_to_double(lattice_staple[gindex*8 + 6 + 4 * lat->halfrowsize4],lattice_staple[gindex*8 + 7 + 4 * lat->halfrowsize4]);

            result.w1.re = GPU0->convert_to_double(lattice_staple[gindex*8 + 0 + 6 * lat->halfrowsize4],lattice_staple[gindex*8 + 1 + 6 * lat->halfrowsize4]);
            result.w2.re = GPU0->convert_to_double(lattice_staple[gindex*8 + 2 + 6 * lat->halfrowsize4],lattice_staple[gindex*8 + 3 + 6 * lat->halfrowsize4]);
            result.w3.re = GPU0->convert_to_double(lattice_staple[gindex*8 + 4 + 6 * lat->halfrowsize4],lattice_staple[gindex*8 + 5 + 6 * lat->halfrowsize4]);

            result.w1.im = GPU0->convert_to_double(lattice_staple[gindex*8 + 0 + 8 * lat->halfrowsize4],lattice_staple[gindex*8 + 1 + 8 * lat->halfrowsize4]);
            result.w2.im = GPU0->convert_to_double(lattice_staple[gindex*8 + 2 + 8 * lat->halfrowsize4],lattice_staple[gindex*8 + 3 + 8 * lat->halfrowsize4]);
            result.w3.im = GPU0->convert_to_double(lattice_staple[gindex*8 + 4 + 8 * lat->halfrowsize4],lattice_staple[gindex*8 + 5 + 8 * lat->halfrowsize4]);
  }
        return result;
}

SU::su_3        SU::lattice_matrix_hermitian(SU::su_3 a){
    su_3 result;

        result.u1.re =  a.u1.re;
        result.u1.im = -a.u1.im;
        result.u2.re =  a.v1.re;
        result.u2.im = -a.v1.im;
        result.u3.re =  a.w1.re;
        result.u3.im = -a.w1.im;

        result.v1.re =  a.u2.re;
        result.v1.im = -a.u2.im;
        result.v2.re =  a.v2.re;
        result.v2.im = -a.v2.im;
        result.v3.re =  a.w2.re;
        result.v3.im = -a.w2.im;

        result.w1.re =  a.u3.re;
        result.w1.im = -a.u3.im;
        result.w2.re =  a.v3.re;
        result.w2.im = -a.v3.im;
        result.w3.re =  a.w3.re;
        result.w3.im = -a.w3.im;

    return result;
}

SU::su_3        SU::lattice_GramSchmidt_3(SU::su_3 a){
    su_3 result;
        hgpu_complex sp;
        double norm_u = sqrt(a.u1.re * a.u1.re + a.u2.re * a.u2.re + a.u3.re * a.u3.re + a.u1.im * a.u1.im + a.u2.im * a.u2.im + a.u3.im * a.u3.im);
        double norm_v;

        result.u1.re = a.u1.re / norm_u;
        result.u1.im = a.u1.im / norm_u;
        result.u2.re = a.u2.re / norm_u;
        result.u2.im = a.u2.im / norm_u;
        result.u3.re = a.u3.re / norm_u;
        result.u3.im = a.u3.im / norm_u;

        sp = hgpu_add(hgpu_add(hgpu_mul(a.v1,hgpu_conjugate(result.u1)),hgpu_mul(a.v2,hgpu_conjugate(result.u2))),hgpu_mul(a.v3,hgpu_conjugate(result.u3)));
        result.v1 = hgpu_sub(a.v1,hgpu_mul(result.u1,sp));
        result.v2 = hgpu_sub(a.v2,hgpu_mul(result.u2,sp));
        result.v3 = hgpu_sub(a.v3,hgpu_mul(result.u3,sp));
        norm_v = sqrt(result.v1.re * result.v1.re + result.v2.re * result.v2.re + result.v3.re * result.v3.re +
                      result.v1.im * result.v1.im + result.v2.im * result.v2.im + result.v3.im * result.v3.im);

        result.v1.re /= norm_v;
        result.v1.im /= norm_v;
        result.v2.re /= norm_v;
        result.v2.im /= norm_v;
        result.v3.re /= norm_v;
        result.v3.im /= norm_v;

    return lattice_matrix_reconstruct3(result);
}

SU::su_2        SU::lattice_matrix_times2(SU::su_2 a,SU::su_2 b){
    su_2 result;

        result.u1.re = -a.u1.im * b.u1.im + a.u1.re * b.u1.re - a.u2.im * b.v1.im + a.u2.re * b.v1.re;
        result.u1.im =  a.u1.re * b.u1.im + a.u1.im * b.u1.re + a.u2.re * b.v1.im + a.u2.im * b.v1.re;
        result.u2.re = -a.u1.im * b.u2.im + a.u1.re * b.u2.re - a.u2.im * b.v2.im + a.u2.re * b.v2.re;
        result.u2.im =  a.u1.re * b.u2.im + a.u1.im * b.u2.re + a.u2.re * b.v2.im + a.u2.im * b.v2.re;
        result.v1.re = -a.v1.im * b.u1.im + a.v1.re * b.u1.re - a.v2.im * b.v1.im + a.v2.re * b.v1.re;
        result.v1.im =  a.v1.re * b.u1.im + a.v1.im * b.u1.re + a.v2.re * b.v1.im + a.v2.im * b.v1.re;
        result.v2.re = -a.v1.im * b.u2.im + a.v1.re * b.u2.re - a.v2.im * b.v2.im + a.v2.re * b.v2.re;
        result.v2.im =  a.v1.re * b.u2.im + a.v1.im * b.u2.re + a.v2.re * b.v2.im + a.v2.im * b.v2.re;

    return result;
}

SU::su_3        SU::lattice_matrix_times3(SU::su_3 a,SU::su_3 b){
    su_3 result;
    result.u1.re = -a.u1.im * b.u1.im + a.u1.re * b.u1.re - a.u2.im * b.v1.im + a.u2.re * b.v1.re - a.u3.im * b.w1.im + a.u3.re * b.w1.re;
    result.u1.im =  a.u1.re * b.u1.im + a.u1.im * b.u1.re + a.u2.re * b.v1.im + a.u2.im * b.v1.re + a.u3.re * b.w1.im + a.u3.im * b.w1.re;
    result.u2.re = -a.u1.im * b.u2.im + a.u1.re * b.u2.re - a.u2.im * b.v2.im + a.u2.re * b.v2.re - a.u3.im * b.w2.im + a.u3.re * b.w2.re;
    result.u2.im =  a.u1.re * b.u2.im + a.u1.im * b.u2.re + a.u2.re * b.v2.im + a.u2.im * b.v2.re + a.u3.re * b.w2.im + a.u3.im * b.w2.re;
    result.u3.re = -a.u1.im * b.u3.im + a.u1.re * b.u3.re - a.u2.im * b.v3.im + a.u2.re * b.v3.re - a.u3.im * b.w3.im + a.u3.re * b.w3.re;
    result.u3.im =  a.u1.re * b.u3.im + a.u1.im * b.u3.re + a.u2.re * b.v3.im + a.u2.im * b.v3.re + a.u3.re * b.w3.im + a.u3.im * b.w3.re;

    result.v1.re = -a.v1.im * b.u1.im + a.v1.re * b.u1.re - a.v2.im * b.v1.im + a.v2.re * b.v1.re - a.v3.im * b.w1.im + a.v3.re * b.w1.re;
    result.v1.im =  a.v1.re * b.u1.im + a.v1.im * b.u1.re + a.v2.re * b.v1.im + a.v2.im * b.v1.re + a.v3.re * b.w1.im + a.v3.im * b.w1.re;
    result.v2.re = -a.v1.im * b.u2.im + a.v1.re * b.u2.re - a.v2.im * b.v2.im + a.v2.re * b.v2.re - a.v3.im * b.w2.im + a.v3.re * b.w2.re;
    result.v2.im =  a.v1.re * b.u2.im + a.v1.im * b.u2.re + a.v2.re * b.v2.im + a.v2.im * b.v2.re + a.v3.re * b.w2.im + a.v3.im * b.w2.re;
    result.v3.re = -a.v1.im * b.u3.im + a.v1.re * b.u3.re - a.v2.im * b.v3.im + a.v2.re * b.v3.re - a.v3.im * b.w3.im + a.v3.re * b.w3.re;
    result.v3.im =  a.v1.re * b.u3.im + a.v1.im * b.u3.re + a.v2.re * b.v3.im + a.v2.im * b.v3.re + a.v3.re * b.w3.im + a.v3.im * b.w3.re;

    result.w1.re = -a.w1.im * b.u1.im + a.w1.re * b.u1.re - a.w2.im * b.v1.im + a.w2.re * b.v1.re - a.w3.im * b.w1.im + a.w3.re * b.w1.re;
    result.w1.im =  a.w1.re * b.u1.im + a.w1.im * b.u1.re + a.w2.re * b.v1.im + a.w2.im * b.v1.re + a.w3.re * b.w1.im + a.w3.im * b.w1.re;
    result.w2.re = -a.w1.im * b.u2.im + a.w1.re * b.u2.re - a.w2.im * b.v2.im + a.w2.re * b.v2.re - a.w3.im * b.w2.im + a.w3.re * b.w2.re;
    result.w2.im =  a.w1.re * b.u2.im + a.w1.im * b.u2.re + a.w2.re * b.v2.im + a.w2.im * b.v2.re + a.w3.re * b.w2.im + a.w3.im * b.w2.re;
    result.w3.re = -a.w1.im * b.u3.im + a.w1.re * b.u3.re - a.w2.im * b.v3.im + a.w2.re * b.v3.re - a.w3.im * b.w3.im + a.w3.re * b.w3.re;
    result.w3.im =  a.w1.re * b.u3.im + a.w1.im * b.u3.re + a.w2.re * b.v3.im + a.w2.im * b.v3.re + a.w3.re * b.w3.im + a.w3.im * b.w3.re;

    return result;
}

SU::su_3        SU::lattice_matrix_add3(SU::su_3 a,SU::su_3 b){
    su_3 result;
        result.u1.re = a.u1.re + b.u1.re;  result.u1.im = a.u1.im + b.u1.im;
        result.u2.re = a.u2.re + b.u2.re;  result.u2.im = a.u2.im + b.u2.im;
        result.u3.re = a.u3.re + b.u3.re;  result.u3.im = a.u3.im + b.u3.im;

        result.v1.re = a.v1.re + b.v1.re;  result.v1.im = a.v1.im + b.v1.im;
        result.v2.re = a.v2.re + b.v2.re;  result.v2.im = a.v2.im + b.v2.im;
        result.v3.re = a.v3.re + b.v3.re;  result.v3.im = a.v3.im + b.v3.im;

        result.w1.re = a.w1.re + b.w1.re;  result.w1.im = a.w1.im + b.w1.im;
        result.w2.re = a.w2.re + b.w2.re;  result.w2.im = a.w2.im + b.w2.im;
        result.w3.re = a.w3.re + b.w3.re;  result.w3.im = a.w3.im + b.w3.im;

    return result;
}

SU::su_3        SU::lattice_staple_hermitian3(SU::su_3 m1,SU::su_3 m2,SU::su_3 m3){
        su_3 u1, u2, u3;

        u1 = lattice_matrix_hermitian(m2);
        u2 = lattice_matrix_times3(m1,u1);
        u1 = lattice_matrix_hermitian(m3);
        u3 = lattice_matrix_times3(u2,u1);


    return u3;
}

SU::su_3        SU::lattice_staple_hermitian_backward3(SU::su_3 m1,SU::su_3 m2,SU::su_3 m3){
        su_3 u1, u2, u3;

        u1 = lattice_matrix_hermitian(m1);
        u2 = lattice_matrix_hermitian(m2);
        u3 = lattice_matrix_times3(u1,u2);
        u1 = lattice_matrix_times3(u3,m3);

    return u1;
}

void            SU::lattice_matrix_print(SU::su_3 matrix){
    printf("| % e +I % e \t % e +I % e \t % e +I % e |\n",matrix.u1.re,matrix.u1.im,matrix.u2.re,matrix.u2.im,matrix.u3.re,matrix.u3.im);
    printf("| % e +I % e \t % e +I % e \t % e +I % e |\n",matrix.v1.re,matrix.v1.im,matrix.v2.re,matrix.v2.im,matrix.v3.re,matrix.v3.im);
    printf("| % e +I % e \t % e +I % e \t % e +I % e |\n",matrix.w1.re,matrix.w1.im,matrix.w2.re,matrix.w2.im,matrix.w3.re,matrix.w3.im);
}

void            SU::lattice_matrix_print_double(SU::su_3 matrix){
    printf("| % 16.14f +% 16.14f I\t% 16.14f +% 16.14f I\t% 16.14f +% 16.14f I |\n",matrix.u1.re,matrix.u1.im,matrix.u2.re,matrix.u2.im,matrix.u3.re,matrix.u3.im);
    printf("| % 16.14f +% 16.14f I\t% 16.14f +% 16.14f I\t% 16.14f +% 16.14f I |\n",matrix.v1.re,matrix.v1.im,matrix.v2.re,matrix.v2.im,matrix.v3.re,matrix.v3.im);
    printf("| % 16.14f +% 16.14f I\t% 16.14f +% 16.14f I\t% 16.14f +% 16.14f I |\n",matrix.w1.re,matrix.w1.im,matrix.w2.re,matrix.w2.im,matrix.w3.re,matrix.w3.im);
}

SU::su_3        SU::lattice_plaquette3(SU::su_3 a,SU::su_3 b,SU::su_3 c,SU::su_3 d){
    su_3 m1, m2, m3;

    m1 = lattice_matrix_times3(a,b);
    m2 = lattice_matrix_hermitian(c);
    m3 = lattice_matrix_times3(m1,m2);
    m2 = lattice_matrix_hermitian(d);
    m1 = lattice_matrix_times3(m3,m2);

    return m1;
}

double          SU::lattice_retrace(SU::su_3 a){
    double result;

           result = a.u1.re + a.v2.re + a.w3.re;

    return result;
}

double          SU::lattice_imtrace(SU::su_3 a){
    double result;

           result = a.u1.im + a.v2.im + a.w3.im;

    return result;
}

SU::su_2        SU::lattice_unity2(void){
    su_2 a;
        a.u1.re = 1.0;  a.u1.im = 0.0;
        a.u2.re = 0.0;  a.u2.im = 0.0;
        a.v1.re = 0.0;  a.v1.im = 0.0;
        a.v2.re = 1.0;  a.v2.im = 0.0;
    return a;
}

SU::su_3        SU::lattice_unity3(void){
    su_3 a;
        a.u1.re = 1.0; a.u1.im = 0.0;
        a.u2.re = 0.0; a.u2.im = 0.0;
        a.u3.re = 0.0; a.u3.im = 0.0;
        a.v1.re = 0.0; a.v1.im = 0.0;
        a.v2.re = 1.0; a.v2.im = 0.0;
        a.v3.re = 0.0; a.v3.im = 0.0;
        a.w1.re = 0.0; a.w1.im = 0.0;
        a.w2.re = 0.0; a.w2.im = 0.0;
        a.w3.re = 1.0; a.w3.im = 0.0;
    return a;
}

SU::su_3        SU::lattice_lambda(int index){
    su_3 a;
    switch (index){
        case 1:
            a = lattice_lambda1();
            break;
        case 2:
            a = lattice_lambda2();
            break;
        case 3:
            a = lattice_lambda3();
            break;
        case 4:
            a = lattice_lambda4();
            break;
        case 5:
            a = lattice_lambda5();
            break;
        case 6:
            a = lattice_lambda6();
            break;
        case 7:
            a = lattice_lambda7();
            break;
        case 8:
            a = lattice_lambda8();
            break;
        default:
            a = lattice_unity3();
            break;
    }
    return a;
}

SU::su_3        SU::lattice_lambda1(void){
    su_3 a;
        a.u1.re = 0.0; a.u1.im = 0.0;
        a.u2.re = 1.0; a.u2.im = 0.0;
        a.u3.re = 0.0; a.u3.im = 0.0;
        a.v1.re = 1.0; a.v1.im = 0.0;
        a.v2.re = 0.0; a.v2.im = 0.0;
        a.v3.re = 0.0; a.v3.im = 0.0;
        a.w1.re = 0.0; a.w1.im = 0.0;
        a.w2.re = 0.0; a.w2.im = 0.0;
        a.w3.re = 0.0; a.w3.im = 0.0;
    return a;
}

SU::su_3        SU::lattice_lambda2(void){
    su_3 a;
        a.u1.re = 0.0; a.u1.im = 0.0;
        a.u2.re = 0.0; a.u2.im =-1.0;
        a.u3.re = 0.0; a.u3.im = 0.0;
        a.v1.re = 0.0; a.v1.im = 1.0;
        a.v2.re = 0.0; a.v2.im = 0.0;
        a.v3.re = 0.0; a.v3.im = 0.0;
        a.w1.re = 0.0; a.w1.im = 0.0;
        a.w2.re = 0.0; a.w2.im = 0.0;
        a.w3.re = 0.0; a.w3.im = 0.0;
    return a;
}

SU::su_3        SU::lattice_lambda3(void){
    su_3 a;
        a.u1.re = 1.0; a.u1.im = 0.0;
        a.u2.re = 0.0; a.u2.im = 0.0;
        a.u3.re = 0.0; a.u3.im = 0.0;
        a.v1.re = 0.0; a.v1.im = 0.0;
        a.v2.re =-1.0; a.v2.im = 0.0;
        a.v3.re = 0.0; a.v3.im = 0.0;
        a.w1.re = 0.0; a.w1.im = 0.0;
        a.w2.re = 0.0; a.w2.im = 0.0;
        a.w3.re = 0.0; a.w3.im = 0.0;
    return a;
}

SU::su_3        SU::lattice_lambda4(void){
    su_3 a;
        a.u1.re = 0.0; a.u1.im = 0.0;
        a.u2.re = 0.0; a.u2.im = 0.0;
        a.u3.re = 1.0; a.u3.im = 0.0;
        a.v1.re = 0.0; a.v1.im = 0.0;
        a.v2.re = 0.0; a.v2.im = 0.0;
        a.v3.re = 0.0; a.v3.im = 0.0;
        a.w1.re = 1.0; a.w1.im = 0.0;
        a.w2.re = 0.0; a.w2.im = 0.0;
        a.w3.re = 0.0; a.w3.im = 0.0;
    return a;
}

SU::su_3        SU::lattice_lambda5(void){
    su_3 a;
        a.u1.re = 0.0; a.u1.im = 0.0;
        a.u2.re = 0.0; a.u2.im = 0.0;
        a.u3.re = 0.0; a.u3.im =-1.0;
        a.v1.re = 0.0; a.v1.im = 0.0;
        a.v2.re = 0.0; a.v2.im = 0.0;
        a.v3.re = 0.0; a.v3.im = 0.0;
        a.w1.re = 0.0; a.w1.im = 1.0;
        a.w2.re = 0.0; a.w2.im = 0.0;
        a.w3.re = 0.0; a.w3.im = 0.0;
    return a;
}

SU::su_3        SU::lattice_lambda6(void){
    su_3 a;
        a.u1.re = 0.0; a.u1.im = 0.0;
        a.u2.re = 0.0; a.u2.im = 0.0;
        a.u3.re = 0.0; a.u3.im = 0.0;
        a.v1.re = 0.0; a.v1.im = 0.0;
        a.v2.re = 0.0; a.v2.im = 0.0;
        a.v3.re = 1.0; a.v3.im = 0.0;
        a.w1.re = 0.0; a.w1.im = 0.0;
        a.w2.re = 1.0; a.w2.im = 0.0;
        a.w3.re = 0.0; a.w3.im = 0.0;
    return a;
}

SU::su_3        SU::lattice_lambda7(void){
    su_3 a;
        a.u1.re = 0.0; a.u1.im = 0.0;
        a.u2.re = 0.0; a.u2.im = 0.0;
        a.u3.re = 0.0; a.u3.im = 0.0;
        a.v1.re = 0.0; a.v1.im = 0.0;
        a.v2.re = 0.0; a.v2.im = 0.0;
        a.v3.re = 0.0; a.v3.im =-1.0;
        a.w1.re = 0.0; a.w1.im = 0.0;
        a.w2.re = 0.0; a.w2.im = 1.0;
        a.w3.re = 0.0; a.w3.im = 0.0;
    return a;
}

SU::su_3        SU::lattice_lambda8(void){
    su_3 a;
        a.u1.re = 0.57735026918962576450914878050196; a.u1.im = 0.0;
        a.u2.re = 0.0;                                a.u2.im = 0.0;
        a.u3.re = 0.0;                                a.u3.im = 0.0;
        a.v1.re = 0.0;                                a.v1.im = 0.0;
        a.v2.re = 0.57735026918962576450914878050196; a.v2.im = 0.0;
        a.v3.re = 0.0;                                a.v3.im = 0.0;
        a.w1.re = 0.0;                                a.w1.im = 0.0;
        a.w2.re = 0.0;                                a.w2.im = 0.0;
        a.w3.re =-1.1547005383792515290182975610039;  a.w3.im = 0.0;
    return a;
}

SU::su_3        SU::lattice_zero3(void){
    su_3 a;
        a.u1.re = 0.0; a.u1.im = 0.0;
        a.u2.re = 0.0; a.u2.im = 0.0;
        a.u3.re = 0.0; a.u3.im = 0.0;
        a.v1.re = 0.0; a.v1.im = 0.0;
        a.v2.re = 0.0; a.v2.im = 0.0;
        a.v3.re = 0.0; a.v3.im = 0.0;
        a.w1.re = 0.0; a.w1.im = 0.0;
        a.w2.re = 0.0; a.w2.im = 0.0;
        a.w3.re = 0.0; a.w3.im = 0.0;
    return a;
}

SU::su_2        SU::lattice_heatbath2(model* lat,SU::su_2 a,double beta,cl_float4* prns,unsigned int gindex){
    su_2 b,aH,c;

    bool flag = false;
    cl_float4 rnd;
    double Mx,My,Mz,Mw;
    double det,bdet,cosrnd;
    double delta = 0.0;
    double costh,sinth,cosal,sinal,phi,sinphi,cosphi;
    int i = 0;

    unsigned int indprng = gindex;

    Mx = (a.u1.re + a.v2.re);
    My = (a.u1.im - a.v2.im);
    Mz = (a.u2.re - a.v1.re);
    Mw = (a.u2.im + a.v1.im);

    det = sqrt(Mx * Mx + My * My + Mz * Mz + Mw * Mw);

    aH.u1.re =  Mx / det;
    aH.u1.im = -My / det;
    aH.u2.re = -Mz / det;
    aH.u2.im = -Mw / det;
    aH.v2.re =  aH.u1.re;
    aH.v2.im = -aH.u1.im;
    aH.v1.re = -aH.u2.re;
    aH.v1.im =  aH.u2.im;

    bdet = beta * det;

    while ((i < lat->NHIT) && (flag == false)){
            rnd.s[0] = (cl_float) (((double) rand()) / ((double) RAND_MAX + 1.0));
            rnd.s[1] = (cl_float) (((double) rand()) / ((double) RAND_MAX + 1.0));
            rnd.s[2] = (cl_float) (((double) rand()) / ((double) RAND_MAX + 1.0));
            rnd.s[3] = (cl_float) (((double) rand()) / ((double) RAND_MAX + 1.0));
        cosrnd = cos(PI2 * rnd.s[1]);
        delta = -(log(1.0 - rnd.s[0]) + cosrnd * cosrnd * log(1.0 - rnd.s[2])) / bdet;
        if ((rnd.s[3] * rnd.s[3])<=(1-0.5*delta)) flag=true;
        indprng += lat->prngstep;
        i++;
    }
        if (flag) {
            rnd.s[0] = (cl_float) (((double) rand()) / ((double) RAND_MAX + 1.0));
            rnd.s[1] = (cl_float) (((double) rand()) / ((double) RAND_MAX + 1.0));
            cosal = 1.0 - delta;
            costh = 2.0 * rnd.s[0] - 1.0;
            sinth = sqrt(1.0 - costh * costh);
            sinal = sqrt(1.0 - cosal * cosal);
            phi   = PI2 * rnd.s[1];

            sinphi = sin(phi);
            cosphi = cos(phi);

            c.u1.re =  cosal;
            c.u1.im =  sinal * costh;
            c.u2.re =  sinal * sinth * sinphi;
            c.u2.im =  sinal * sinth * cosphi;
            c.v2.re =  c.u1.re;
            c.v2.im = -c.u1.im;
            c.v1.re = -c.u2.re;
            c.v1.im =  c.u2.im;
            
            b = lattice_matrix_times2(c,aH);

        } else {
            b = lattice_unity2();
        }

    return b;
}

SU::su_3        SU::lattice_heatbath3(model* lat,SU::su_3 staple,SU::su_3 U0,double beta,cl_float4* prns,unsigned int gindex){
    su_3 result,x0,x1,U1,U2,r00;
       x0 = lattice_matrix_times3(U0,staple);

        su_2 r0,r;

        r.u1.re = x0.u1.re;
        r.u1.im = x0.u1.im;
        r.u2.re = x0.u2.re;
        r.u2.im = x0.u2.im;
        r.v1.re = x0.v1.re;
        r.v1.im = x0.v1.im;
        r.v2.re = x0.v2.re;
        r.v2.im = x0.v2.im;

        r0 = lattice_heatbath2(lat,r,beta,prns,gindex);  // r->(heatbath)->r0

        r00.u1.re = r0.u1.re;
        r00.u1.im = r0.u1.im;
        r00.u2.re = r0.u2.re;
        r00.u2.im = r0.u2.im;
        r00.u3.re = 0.0;
        r00.u3.im = 0.0;

        r00.v1.re = r0.v1.re;
        r00.v1.im = r0.v1.im;
        r00.v2.re = r0.v2.re;
        r00.v2.im = r0.v2.im;
        r00.v3.re = 0.0;
        r00.v3.im = 0.0;

        r00.w1.re = 0.0;
        r00.w1.im = 0.0;
        r00.w2.re = 0.0;
        r00.w2.im = 0.0;
        r00.w3.re = 1.0;
        r00.w3.im = 0.0;

        U1 = lattice_matrix_times3(r00,U0);

        x1 = lattice_matrix_times3(U1,staple);

        r.u1.re = x1.u1.re;
        r.u1.im = x1.u1.im;
        r.u2.re = x1.u3.re;
        r.u2.im = x1.u3.im;
        r.v1.re = x1.w1.re;
        r.v1.im = x1.w1.im;
        r.v2.re = x1.w3.re;
        r.v2.im = x1.w3.im;

        r0 = lattice_heatbath2(lat,r,beta,prns,gindex);  // r->(heatbath)->r0

        r00.u1.re = r0.u1.re;
        r00.u1.im = r0.u1.im;
        r00.u2.re = 0.0;
        r00.u2.im = 0.0;
        r00.u3.re = r0.u2.re;
        r00.u3.im = r0.u2.im;

        r00.v1.re = 0.0;
        r00.v1.im = 0.0;
        r00.v2.re = 1.0;
        r00.v2.im = 0.0;
        r00.v3.re = 0.0;
        r00.v3.im = 0.0;

        r00.w1.re = r0.v1.re;
        r00.w1.im = r0.v1.im;
        r00.w2.re = 0.0;
        r00.w2.im = 0.0;
        r00.w3.re = r0.v2.re;
        r00.w3.im = r0.v2.im;

        U2 = lattice_matrix_times3(r00,U1);

        x0 = lattice_matrix_times3(U2,staple);

        r.u1.re = x0.v2.re;
        r.u1.im = x0.v2.im;
        r.u2.re = x0.v3.re;
        r.u2.im = x0.v3.im;
        r.v1.re = x0.w2.re;
        r.v1.im = x0.w2.im;
        r.v2.re = x0.w3.re;
        r.v2.im = x0.w3.im;

        r0 = lattice_heatbath2(lat,r,beta,prns,gindex);  // r->(heatbath)->r0

        r00.u1.re = 1.0;
        r00.u1.im = 0.0;
        r00.u2.re = 0.0;
        r00.u2.im = 0.0;
        r00.u3.re = 0.0;
        r00.u3.im = 0.0;

        r00.v1.re = 0.0;
        r00.v1.im = 0.0;
        r00.v2.re = r0.u1.re;
        r00.v2.im = r0.u1.im;
        r00.v3.re = r0.u2.re;
        r00.v3.im = r0.u2.im;

        r00.w1.re = 0.0;
        r00.w1.im = 0.0;
        r00.w2.re = r0.v1.re;
        r00.w2.im = r0.v1.im;
        r00.w3.re = r0.v2.re;
        r00.w3.im = r0.v2.im;

        result = lattice_matrix_times3(r00,U2);

    return result;
}

unsigned int    SU::lattice_odd_gid(model* lat,unsigned int gid){
    unsigned int even_check,gindex,gde;
    gde = 2 * gid + 1;
    coords_4 coord;
    coord = lattice_gid_to_coords(lat,gde);
    even_check = ((coord.x + coord.y + coord.z + coord.t) & 1) ^ 1;
    gindex = gde - even_check;

    return gindex;
}

void            SU::lattice_load_cpu(model* lat,unsigned int* lattice_pointer){
        coords_4 coords;
        unsigned int gdi;
        su_3 matrix;
        for (int x1 = 0; x1 < lat->lattice_domain_n1; x1++)
            for (int x2 = 0; x2 < lat->lattice_domain_size[1]; x2++)
            for (int x3 = 0; x3 < lat->lattice_domain_size[2]; x3++)
            for (int x4 = 0; x4 < lat->lattice_domain_size[3]; x4++)
            for (int dir1 = 0; dir1<lat->lattice_nd; dir1++)
            {
                        coords.x = x1;
                        coords.y = x2;
                        coords.z = x3;
                        coords.t = x4;

                        gdi = lattice_coords_to_gid(lat,coords);
                        matrix = lattice_get_3(lat,lattice_pointer,gdi,dir1);
                        lattice_store_3(lat,matrix,gdi,dir1);
            }

}

double*         SU::lattice_avr_plaquette_cpu(model* lat){
    double* result = new double[2];
    coords_4 coords,coords2,coords3;

    double plq_spat = 0.0;
    double plq_temp = 0.0;

    double tmp_spat,tmp_temp;
    unsigned int gdi,gdi2,gdi3;
    su_3 matrix_1,matrix_2,matrix_3,matrix_4,plaquette;
    double mult = lat->BETA / ((double) lat->lattice_group);

    for (int x1 = 0; x1 < lat->lattice_domain_size[0]; x1++)
        for (int x2 = 0; x2 < lat->lattice_domain_size[1]; x2++)
            for (int x3 = 0; x3 < lat->lattice_domain_size[2]; x3++)
                for (int x4 = 0; x4 < lat->lattice_domain_size[3]; x4++)
                    for (int dir1 = 0; dir1<(lat->lattice_nd-1); dir1++)
                        for (int dir2 = (dir1+1); dir2<lat->lattice_nd; dir2++)
                {
                    coords.x = x1;
                    coords.y = x2;
                    coords.z = x3;
                    coords.t = x4;

                    tmp_spat = 0.0;
                    tmp_temp = 0.0;

                    gdi = lattice_coords_to_gid(lat,coords);
                    matrix_1 = lattice_table_3(lat,coords,gdi,dir1);

                    coords2 = lattice_neighbours_coords(lat,coords,dir1);
                    gdi2 = lattice_coords_to_gid(lat,coords2);
                    matrix_2 = lattice_table_3(lat,coords2,gdi2,dir2);

                    coords3 = lattice_neighbours_coords(lat,coords,dir2);
                    gdi3 = lattice_coords_to_gid(lat,coords3);
                    matrix_3 = lattice_table_3(lat,coords3,gdi3,dir1);

                    matrix_4 = lattice_table_3(lat,coords,gdi,dir2);

                    plaquette = lattice_plaquette3(matrix_1,matrix_2,matrix_3,matrix_4);

                    if (dir2==(lat->lattice_nd-1)) {
                        tmp_temp = lattice_retrace(plaquette);
                        plq_temp += mult * (3.0 - tmp_temp);
                    } else {
                        tmp_spat = lattice_retrace(plaquette);
                        plq_spat += mult * (3.0 - tmp_spat);
                    }
                }

        result[0] = plq_spat / ((double) (lat->lattice_full_site * (lat->lattice_nd - 1)));
        result[1] = plq_temp / ((double) (lat->lattice_full_site * (lat->lattice_nd - 1)));

        return result;
}

double*         SU::lattice_avr_plaquette_plq_cpu(model* lat){
    double* result = new double[26];
    coords_4 coords,coords2,coords3;

    double plq_spat = 0.0;
    double plq_temp = 0.0;

    double tmp_spat,tmp_temp;
    double F_xy_3_re_temp = 0.0;
    double F_xy_3_im_temp = 0.0;
    double F_xz_3_re_temp = 0.0;
    double F_xz_3_im_temp = 0.0;
    double F_yz_3_re_temp = 0.0;
    double F_yz_3_im_temp = 0.0;
    double F_xy_8_re_temp = 0.0;
    double F_xy_8_im_temp = 0.0;
    double F_xz_8_re_temp = 0.0;
    double F_xz_8_im_temp = 0.0;
    double F_yz_8_re_temp = 0.0;
    double F_yz_8_im_temp = 0.0;

    double re_trace = 0.0;
    double im_trace = 0.0;
    double F_xy_3_re_variance = 0.0;
    double F_xy_3_im_variance = 0.0;
    double F_xz_3_re_variance = 0.0;
    double F_xz_3_im_variance = 0.0;
    double F_yz_3_re_variance = 0.0;
    double F_yz_3_im_variance = 0.0;
    double F_xy_8_re_variance = 0.0;
    double F_xy_8_im_variance = 0.0;
    double F_xz_8_re_variance = 0.0;
    double F_xz_8_im_variance = 0.0;
    double F_yz_8_re_variance = 0.0;
    double F_yz_8_im_variance = 0.0;

    double denominator1 = ((double) (lat->lattice_full_site * (lat->lattice_nd - 1)));
    double denominator2 = ((double) (lat->lattice_full_site));

    unsigned int gdi,gdi2,gdi3;
    su_3 matrix_1,matrix_2,matrix_3,matrix_4,plaquette;
    su_3 matrix_5;

    for (int x1 = 0; x1 < lat->lattice_full_size[0]; x1++)
        for (int x2 = 0; x2 < lat->lattice_full_size[1]; x2++)
            for (int x3 = 0; x3 < lat->lattice_full_size[2]; x3++)
                for (int x4 = 0; x4 < lat->lattice_full_size[3]; x4++)
                    for (int dir1 = 0; dir1<(lat->lattice_nd-1); dir1++)
                        for (int dir2 = (dir1+1); dir2<lat->lattice_nd; dir2++)
                {
                    coords.x = x1;
                    coords.y = x2;
                    coords.z = x3;
                    coords.t = x4;

                    tmp_spat = 0.0;
                    tmp_temp = 0.0;

                    gdi = lattice_coords_to_gid(lat,coords);
                    matrix_1 = lattice_table_3(lat,coords,gdi,dir1);

                    coords2 = lattice_neighbours_coords(lat,coords,dir1);
                    gdi2 = lattice_coords_to_gid(lat,coords2);
                    matrix_2 = lattice_table_3(lat,coords2,gdi2,dir2);

                    coords3 = lattice_neighbours_coords(lat,coords,dir2);
                    gdi3 = lattice_coords_to_gid(lat,coords3);
                    matrix_3 = lattice_table_3(lat,coords3,gdi3,dir1);

                    matrix_4 = lattice_table_3(lat,coords,gdi,dir2);

                    plaquette = lattice_plaquette3(matrix_1,matrix_2,matrix_3,matrix_4);

                    if (dir2==(lat->lattice_nd-1)) {
                        tmp_temp = lattice_retrace(plaquette);
                        plq_temp += tmp_temp;
                    } else {
                        tmp_spat = lattice_retrace(plaquette);
                        plq_spat += tmp_spat;
                    }

                    if (((dir1==0) && (dir2==1) && (lat->get_Fmunu)) || // _____XY for get_Fmunu
                        ((dir1==0) && (dir2==3) && (lat->get_F0mu)))    // _____XT for get_F0mu
                    {
                        matrix_5 = lattice_matrix_times3(plaquette,lattice_lambda(lat->Fmunu_index1));
                        re_trace = lattice_retrace(matrix_5);
                        im_trace = lattice_imtrace(matrix_5);
                        F_xy_3_re_temp += re_trace;
                        F_xy_3_im_temp += im_trace;
                        F_xy_3_re_variance += re_trace * re_trace;
                        F_xy_3_im_variance += im_trace * im_trace;

                        matrix_5 = lattice_matrix_times3(plaquette,lattice_lambda(lat->Fmunu_index2));
                        re_trace = lattice_retrace(matrix_5);
                        im_trace = lattice_imtrace(matrix_5);
                        F_xy_8_re_temp += re_trace;
                        F_xy_8_im_temp += im_trace;
                        F_xy_8_re_variance += re_trace * re_trace;
                        F_xy_8_im_variance += im_trace * im_trace;
                    }
                    if (((dir1==0) && (dir2==2) && (lat->get_Fmunu)) || // _____XZ for get_Fmunu
                        ((dir1==1) && (dir2==3) && (lat->get_F0mu)))    // _____YT for get_F0mu
                    {
                        matrix_5 = lattice_matrix_times3(plaquette,lattice_lambda(lat->Fmunu_index1));
                        re_trace = lattice_retrace(matrix_5);
                        im_trace = lattice_imtrace(matrix_5);
                        F_xz_3_re_temp += re_trace;
                        F_xz_3_im_temp += im_trace;
                        F_xz_3_re_variance += re_trace * re_trace;
                        F_xz_3_im_variance += im_trace * im_trace;

                        matrix_5 = lattice_matrix_times3(plaquette,lattice_lambda(lat->Fmunu_index2));
                        re_trace = lattice_retrace(matrix_5);
                        im_trace = lattice_imtrace(matrix_5);
                        F_xz_8_re_temp += re_trace;
                        F_xz_8_im_temp += im_trace;
                        F_xz_8_re_variance += re_trace * re_trace;
                        F_xz_8_im_variance += im_trace * im_trace;
                    }
                    if (((dir1==1) && (dir2==2) && (lat->get_Fmunu)) || // _____YZ for get_Fmunu
                        ((dir1==2) && (dir2==3) && (lat->get_F0mu)))    // _____ZT for get_F0mu
                    {
                        matrix_5 = lattice_matrix_times3(plaquette,lattice_lambda(lat->Fmunu_index1));
                        re_trace = lattice_retrace(matrix_5);
                        im_trace = lattice_imtrace(matrix_5);
                        F_yz_3_re_temp += re_trace;
                        F_yz_3_im_temp += im_trace;
                        F_yz_3_re_variance += re_trace * re_trace;
                        F_yz_3_im_variance += im_trace * im_trace;

                        matrix_5 = lattice_matrix_times3(plaquette,lattice_lambda(lat->Fmunu_index2));
                        re_trace = lattice_retrace(matrix_5);
                        im_trace = lattice_imtrace(matrix_5);
                        F_yz_8_re_temp += re_trace;
                        F_yz_8_im_temp += im_trace;
                        F_yz_8_re_variance += re_trace * re_trace;
                        F_yz_8_im_variance += im_trace * im_trace;
                    }
                }

        result[ 0] = plq_spat       / denominator1;
        result[ 1] = plq_temp       / denominator1;
        result[ 2] = F_xy_3_re_temp / denominator2;
        result[ 3] = F_xy_3_im_temp / denominator2;
        result[ 4] = F_xz_3_re_temp / denominator2;
        result[ 5] = F_xz_3_im_temp / denominator2;
        result[ 6] = F_yz_3_re_temp / denominator2;
        result[ 7] = F_yz_3_im_temp / denominator2;
        result[ 8] = F_xy_8_re_temp / denominator2;
        result[ 9] = F_xy_8_im_temp / denominator2;
        result[10] = F_xz_8_re_temp / denominator2;
        result[11] = F_xz_8_im_temp / denominator2;
        result[12] = F_yz_8_re_temp / denominator2;
        result[13] = F_yz_8_im_temp / denominator2;

        result[14] = F_xy_3_re_variance / denominator2 - result[ 2] * result[ 2];
        result[15] = F_xy_3_im_variance / denominator2 - result[ 3] * result[ 3];
        result[16] = F_xz_3_re_variance / denominator2 - result[ 4] * result[ 4];
        result[17] = F_xz_3_im_variance / denominator2 - result[ 5] * result[ 5];
        result[18] = F_yz_3_re_variance / denominator2 - result[ 6] * result[ 6];
        result[19] = F_yz_3_im_variance / denominator2 - result[ 7] * result[ 7];
        result[20] = F_xy_8_re_variance / denominator2 - result[ 8] * result[ 8];
        result[21] = F_xy_8_im_variance / denominator2 - result[ 9] * result[ 9];
        result[22] = F_xz_8_re_variance / denominator2 - result[10] * result[10];
        result[23] = F_xz_8_im_variance / denominator2 - result[11] * result[11];
        result[24] = F_yz_8_re_variance / denominator2 - result[12] * result[12];
        result[25] = F_yz_8_im_variance / denominator2 - result[13] * result[13];

        return result;
}

double*         SU::lattice_avr_Polyakov_loop_cpu(model* lat){
    double* result = new double[4];
    coords_4 coords;
    double polyakov_loop    = 0.0;
    double polyakov_loop_im = 0.0;
    double polyakov_loop_P2 = 0.0;
    double polyakov_loop_P4 = 0.0;
    double pl_loop    = 0.0;
    double pl_loop_im = 0.0;
    double pl_loop_P2 = 0.0;
    unsigned int gdi;
    su_3 matrix_1,matrix_2,matrix_3;

    int dir1 = 3;
    for (int x1 = 0; x1 < lat->lattice_domain_size[0]; x1++)
        for (int x2 = 0; x2 < lat->lattice_domain_size[1]; x2++)
            for (int x3 = 0; x3 < lat->lattice_domain_size[2]; x3++) {
                coords.x = x1;
                coords.y = x2;
                coords.z = x3;
                coords.t = 0;
                gdi = lattice_coords_to_gid(lat,coords);
                matrix_2 = lattice_table_3(lat,coords,gdi,dir1);
                for (int x4 = 1; x4 < lat->lattice_domain_size[3]; x4++)
                {
                    coords.x = x1;
                    coords.y = x2;
                    coords.z = x3;
                    coords.t = x4;

                    gdi = lattice_coords_to_gid(lat,coords);
                    matrix_1 = lattice_table_3(lat,coords,gdi,dir1);

                    matrix_3 = lattice_matrix_times3(matrix_2,matrix_1);
                    matrix_2 = matrix_3;
                }
                pl_loop         = lattice_retrace(matrix_2);
                pl_loop_im      = lattice_imtrace(matrix_2);
                pl_loop_P2      = pl_loop * pl_loop + pl_loop_im * pl_loop_im;
                polyakov_loop    += pl_loop;
                polyakov_loop_im += pl_loop_im;
                polyakov_loop_P2 += pl_loop_P2;
                polyakov_loop_P4 += pl_loop_P2 * pl_loop_P2;
            }
    result[0] = (polyakov_loop    / (lat->lattice_full_n1n2n3 * lat->lattice_group));
    result[1] = (polyakov_loop_im / (lat->lattice_full_n1n2n3 * lat->lattice_group));
    result[2] = (polyakov_loop_P2 / (lat->lattice_full_n1n2n3 * lat->lattice_group * lat->lattice_group));
    result[3] = (polyakov_loop_P4 / (lat->lattice_full_n1n2n3 * lat->lattice_group * lat->lattice_group * lat->lattice_group * lat->lattice_group));

    return result;
}

double          SU::lattice_avr_Wilson_loop_cpu(model* lat){
    double result;
    coords_4 coords,coords2,coords3,coords_l,coords_r;
    coords_4 coords0,coords1;
    double wilson_loop    = 0.0;
    unsigned int gdi0,gdi1,gdi,gdi_r,gdi_l;
    su_3 matrix_2,matrix_3,matrix_l,matrix_r,matrix_t,matrix_b,plaquette;

    int dir1 = 3;
    for (int x1 = 0; x1 < lat->lattice_domain_size[0]; x1++)
        for (int x2 = 0; x2 < lat->lattice_domain_size[1]; x2++)
            for (int x3 = 0; x3 < lat->lattice_domain_size[2]; x3++)
                for (int x4 = 0; x4 < lat->lattice_domain_size[3]; x4++) {
                    coords.x = x1;
                    coords.y = x2;
                    coords.z = x3;
                    coords.t = x4;
                    coords0 = coords;
                    gdi = lattice_coords_to_gid(lat,coords);
                    gdi0 = gdi;
                    matrix_b = lattice_table_3(lat,coords,gdi,dir1);
                    for (int x_t = 1; x_t < lat->wilson_T; x_t++)    // bottom link
                    {
                        coords.x = x1;
                        coords.y = x2;
                        coords.z = x3;
                        coords.t = (x4 + x_t) % lat->lattice_domain_size[3];

                        gdi = lattice_coords_to_gid(lat,coords);
                        matrix_2 = lattice_table_3(lat,coords,gdi,dir1);

                        matrix_3 = lattice_matrix_times3(matrix_b,matrix_2);
                        matrix_b = matrix_3;
                    }
                    coords1.x = x1;
                    coords1.y = x2;
                    coords1.z = x3;
                    coords1.t = (x4 + lat->wilson_T) % lat->lattice_domain_size[3];
                    gdi1 = lattice_coords_to_gid(lat,coords1);
                    

                    for (int dir2 = 0; dir2 < (lat->lattice_nd - 1); dir2++){
                        coords_r = coords1;
                        gdi_r = gdi1;
                        matrix_r = lattice_table_3(lat,coords_r,gdi_r,dir2);

                        coords_l = coords0;
                        gdi_l = gdi0;
                        matrix_l = lattice_table_3(lat,coords_l,gdi_l,dir2);

                        for (int x_r = 1; x_r < lat->wilson_R; x_r++){  // left and right links
                            coords2 = lattice_neighbours_coords(lat,coords_l,dir2);
                            coords3 = lattice_neighbours_coords(lat,coords_r,dir2);
                            coords_l = coords2;
                            coords_r = coords3;

                            gdi_l = lattice_coords_to_gid(lat,coords_l);
                            matrix_2 = lattice_table_3(lat,coords_l,gdi_l,dir2);
                            matrix_3 = lattice_matrix_times3(matrix_l,matrix_2);
                            matrix_l = matrix_3;

                            gdi_r = lattice_coords_to_gid(lat,coords_r);
                            matrix_2 = lattice_table_3(lat,coords_r,gdi_r,dir2);
                            matrix_3 = lattice_matrix_times3(matrix_r,matrix_2);
                            matrix_r = matrix_3;
                        }
                        coords2 = lattice_neighbours_coords(lat,coords_l,dir2);
                        coords_l = coords2;
                        gdi = lattice_coords_to_gid(lat,coords_l);
                        matrix_t = lattice_table_3(lat,coords,gdi,dir1);
                        for (int x_t = 1; x_t < lat->wilson_T; x_t++)    // top link
                        {
                            coords2 = lattice_neighbours_coords(lat,coords_l,dir1);
                            coords_l = coords2;

                            gdi = lattice_coords_to_gid(lat,coords_l);
                            matrix_2 = lattice_table_3(lat,coords_l,gdi,dir1);

                            matrix_3 = lattice_matrix_times3(matrix_t,matrix_2);
                            matrix_t = matrix_3;
                        }
                        plaquette = lattice_plaquette3(matrix_b,matrix_r,matrix_t,matrix_l);
                        wilson_loop += lattice_retrace(plaquette);
                    }
            }
            result = wilson_loop / ((double) (lat->lattice_full_site * lat->lattice_group));
    return result;
}

double*         SU::lattice_plaquette_cpu(model* lat,coords_4 coords){
    double* result = new double[2];
    coords_4 coords2,coords3;

    double plq_spat = 0.0;
    double plq_temp = 0.0;

    double tmp_spat,tmp_temp;
    unsigned int gdi,gdi2,gdi3;
    su_3 matrix_1,matrix_2,matrix_3,matrix_4,plaquette;

    for (int dir1 = 0; dir1<(lat->lattice_nd-1); dir1++)
        for (int dir2 = dir1+1; dir2<lat->lattice_nd; dir2++)
                {
                    tmp_spat = 0.0;
                    tmp_temp = 0.0;

                    gdi = lattice_coords_to_gid(lat,coords);
                    matrix_1 = lattice_table_3(lat,coords,gdi,dir1);

                    coords2 = lattice_neighbours_coords(lat,coords,dir1);
                    gdi2 = lattice_coords_to_gid(lat,coords2);
                    matrix_2 = lattice_table_3(lat,coords2,gdi2,dir2);

                    coords3 = lattice_neighbours_coords(lat,coords,dir2);
                    gdi3 = lattice_coords_to_gid(lat,coords3);
                    matrix_3 = lattice_table_3(lat,coords3,gdi3,dir1);

                    matrix_4 = lattice_table_3(lat,coords,gdi,dir2);

                    plaquette = lattice_plaquette3(matrix_1,matrix_2,matrix_3,matrix_4);

                    if (dir2==(lat->lattice_nd-1))
                        tmp_temp = lattice_retrace(plaquette);
                    else 
                        tmp_spat = lattice_retrace(plaquette);

                    plq_temp += tmp_temp;
                    plq_spat += tmp_spat;

                }

        result[0] = plq_spat / (lat->lattice_group);
        result[1] = plq_temp / (lat->lattice_group);

        return result;
}

SU::su_3        SU::lattice_get_staple_cpu(model* lat,SU::coords_4 coords,int dir1){
        unsigned int gdi,gdi2,gdi3,gdi4;
        coords_4 coords2,coords3,coords4;
        su_3 matrix_1,matrix_2,matrix_3,matrix_4,matrix_5,matrix_6;

        // _____ staples check__________________________________________________________________
        matrix_5 = lattice_zero3();
        for (int dir2 = 0; dir2<lat->lattice_nd; dir2++)
        {
            if (dir2!=dir1) {
                gdi = lattice_coords_to_gid(lat,coords);

                coords2  = lattice_neighbours_coords(lat,coords,dir1);
                gdi2     = lattice_coords_to_gid(lat,coords2);                       // gdi2 = x + dir1
                matrix_1 = lattice_table_3(lat,coords2,gdi2,dir2);

                coords3  = lattice_neighbours_coords(lat,coords,dir2);
                gdi3     = lattice_coords_to_gid(lat,coords3);                      // gdi3 = x + dir2
                matrix_2 = lattice_table_3(lat,coords3,gdi3,dir1);

                matrix_3 = lattice_table_3(lat,coords,gdi,dir2);

                matrix_4 = lattice_staple_hermitian3(matrix_1,matrix_2,matrix_3);

                matrix_6 = lattice_matrix_add3(matrix_5,matrix_4);

                coords3  = lattice_neighbours_coords_backward(lat,coords2,dir2);
                gdi3     = lattice_coords_to_gid(lat,coords3);                      // gdi3 = x + dir1 - dir2
                matrix_1 = lattice_table_3(lat,coords3,gdi3,dir2);

                coords4  = lattice_neighbours_coords_backward(lat,coords,dir2);
                gdi4     = lattice_coords_to_gid(lat,coords4);                      // gdi4 = x - dir2
                matrix_2 = lattice_table_3(lat,coords4,gdi4,dir1);

                matrix_3 = lattice_table_3(lat,coords4,gdi4,dir2);
                matrix_4 = lattice_staple_hermitian_backward3(matrix_1,matrix_2,matrix_3);

                matrix_5 = lattice_matrix_add3(matrix_6,matrix_4);
            }
        }
        // _____ staples check (END)
        return matrix_5;
}

void            SU::lattice_check_cpu(model* lat){
    lat->Analysis[DM_S_spat].CPU_last_value      = 0.0;
    lat->Analysis[DM_S_temp].CPU_last_value      = 0.0;
    lat->Analysis[DM_S_total].CPU_last_value     = 0.0;
    lat->Analysis[DM_Plq_spat].CPU_last_value    = 0.0;
    lat->Analysis[DM_Plq_temp].CPU_last_value    = 0.0;
    lat->Analysis[DM_Plq_total].CPU_last_value   = 0.0;
    lat->Analysis[DM_Wilson_loop].CPU_last_value = 0.0;
    unsigned int lattice_measurement_s = lat->lattice_measurement_size;
    if ((lat->get_Fmunu)||((lat->get_F0mu))) lattice_measurement_s = lat->lattice_measurement_size_F;

    if (lat->lattice_pointer_last != NULL){
        lattice_data = (su_3*)  calloc(lat->lattice_table_row_size * lat->lattice_nd, sizeof(su_3));

        lattice_load_cpu(lat,lat->lattice_pointer_last);    // Load lattice state to CPU
        double* plaq = lattice_avr_plaquette_cpu(lat);
            lat->Analysis[DM_S_spat].CPU_last_value  = plaq[0];
            lat->Analysis[DM_S_temp].CPU_last_value  = plaq[1];
            lat->Analysis[DM_S_total].CPU_last_value = 0.5 * (plaq[0] + plaq[1]);
            delete[] plaq;

        if ((lat->get_plaquettes_avr)||(lat->get_Fmunu)||(lat->get_F0mu)) {
            double* plaq_plq = lattice_avr_plaquette_plq_cpu(lat);
                lat->Analysis[DM_Plq_spat].CPU_last_value  = plaq_plq[0];
                lat->Analysis[DM_Plq_temp].CPU_last_value  = plaq_plq[1];
                lat->Analysis[DM_Plq_total].CPU_last_value = 0.5 * (plaq_plq[0] + plaq_plq[1]);
                lat->Analysis[DM_Fmunu_xy_3_re].CPU_last_value = plaq_plq[ 2];
                lat->Analysis[DM_Fmunu_xy_3_im].CPU_last_value = plaq_plq[ 3];
                lat->Analysis[DM_Fmunu_xz_3_re].CPU_last_value = plaq_plq[ 4];
                lat->Analysis[DM_Fmunu_xz_3_im].CPU_last_value = plaq_plq[ 5];
                lat->Analysis[DM_Fmunu_yz_3_re].CPU_last_value = plaq_plq[ 6];
                lat->Analysis[DM_Fmunu_yz_3_im].CPU_last_value = plaq_plq[ 7];

                lat->Analysis[DM_Fmunu_xy_8_re].CPU_last_value = plaq_plq[ 8];
                lat->Analysis[DM_Fmunu_xy_8_im].CPU_last_value = plaq_plq[ 9];
                lat->Analysis[DM_Fmunu_xz_8_re].CPU_last_value = plaq_plq[10];
                lat->Analysis[DM_Fmunu_xz_8_im].CPU_last_value = plaq_plq[11];
                lat->Analysis[DM_Fmunu_yz_8_re].CPU_last_value = plaq_plq[12];
                lat->Analysis[DM_Fmunu_yz_8_im].CPU_last_value = plaq_plq[13];

                lat->Analysis[DM_Fmunu_xy_3_re].CPU_last_variance = plaq_plq[14];
                lat->Analysis[DM_Fmunu_xy_3_im].CPU_last_variance = plaq_plq[15];
                lat->Analysis[DM_Fmunu_xz_3_re].CPU_last_variance = plaq_plq[16];
                lat->Analysis[DM_Fmunu_xz_3_im].CPU_last_variance = plaq_plq[17];
                lat->Analysis[DM_Fmunu_yz_3_re].CPU_last_variance = plaq_plq[18];
                lat->Analysis[DM_Fmunu_yz_3_im].CPU_last_variance = plaq_plq[19];

                lat->Analysis[DM_Fmunu_xy_8_re].CPU_last_variance = plaq_plq[20];
                lat->Analysis[DM_Fmunu_xy_8_im].CPU_last_variance = plaq_plq[21];
                lat->Analysis[DM_Fmunu_xz_8_re].CPU_last_variance = plaq_plq[22];
                lat->Analysis[DM_Fmunu_xz_8_im].CPU_last_variance = plaq_plq[23];
                lat->Analysis[DM_Fmunu_yz_8_re].CPU_last_variance = plaq_plq[24];
                lat->Analysis[DM_Fmunu_yz_8_im].CPU_last_variance = plaq_plq[25];

                lat->Analysis[DM_Fmunu_abs_3_re].CPU_last_value = sqrt(lat->Analysis[DM_Fmunu_xy_3_re].CPU_last_value * lat->Analysis[DM_Fmunu_xy_3_re].CPU_last_value +
                                                                       lat->Analysis[DM_Fmunu_xz_3_re].CPU_last_value * lat->Analysis[DM_Fmunu_xz_3_re].CPU_last_value +
                                                                       lat->Analysis[DM_Fmunu_yz_3_re].CPU_last_value * lat->Analysis[DM_Fmunu_yz_3_re].CPU_last_value);

                lat->Analysis[DM_Fmunu_abs_3_im].CPU_last_value = sqrt(lat->Analysis[DM_Fmunu_xy_3_im].CPU_last_value * lat->Analysis[DM_Fmunu_xy_3_im].CPU_last_value +
                                                                       lat->Analysis[DM_Fmunu_xz_3_im].CPU_last_value * lat->Analysis[DM_Fmunu_xz_3_im].CPU_last_value +
                                                                       lat->Analysis[DM_Fmunu_yz_3_im].CPU_last_value * lat->Analysis[DM_Fmunu_yz_3_im].CPU_last_value);

                lat->Analysis[DM_Fmunu_abs_8_re].CPU_last_value = sqrt(lat->Analysis[DM_Fmunu_xy_8_re].CPU_last_value * lat->Analysis[DM_Fmunu_xy_8_re].CPU_last_value +
                                                                       lat->Analysis[DM_Fmunu_xz_8_re].CPU_last_value * lat->Analysis[DM_Fmunu_xz_8_re].CPU_last_value +
                                                                       lat->Analysis[DM_Fmunu_yz_8_re].CPU_last_value * lat->Analysis[DM_Fmunu_yz_8_re].CPU_last_value);

                lat->Analysis[DM_Fmunu_abs_8_im].CPU_last_value = sqrt(lat->Analysis[DM_Fmunu_xy_8_im].CPU_last_value * lat->Analysis[DM_Fmunu_xy_8_im].CPU_last_value +
                                                                       lat->Analysis[DM_Fmunu_xz_8_im].CPU_last_value * lat->Analysis[DM_Fmunu_xz_8_im].CPU_last_value +
                                                                       lat->Analysis[DM_Fmunu_yz_8_im].CPU_last_value * lat->Analysis[DM_Fmunu_yz_8_im].CPU_last_value);

                delete[] plaq_plq;
        }

        if (lat->get_wilson_loop) {
            double wilson_loop = lattice_avr_Wilson_loop_cpu(lat);
            lat->Analysis[DM_Wilson_loop].CPU_last_value = wilson_loop;
        }


        double* pl_avr = lattice_avr_Polyakov_loop_cpu(lat);
            lat->Analysis[DM_Polyakov_loop].CPU_last_value    = pl_avr[0];
            lat->Analysis[DM_Polyakov_loop_im].CPU_last_value = pl_avr[1];
            lat->Analysis[DM_Polyakov_loop_P2].CPU_last_value = pl_avr[2];
            lat->Analysis[DM_Polyakov_loop_P4].CPU_last_value = pl_avr[3];

            delete[] pl_avr;
        free(lattice_data);
    }
}

void            SU::lattice_analysis_cpu(model* lat){
    double pl            = 0.0;
    double pl_im         = 0.0;

    for (int i=0; i<lat->ITER; i++) {
        if (lat->get_plaquettes_avr){
            if ((i<5)||(i==(lat->ITER - 1))) printf("[%5u]: % 10.8f \t % 10.8f \t % 10.8f \t % 10.8f \t % 10.8f \t % 10.8f\n",
                    (i),lat->Analysis[DM_Plq_spat].CPU_data[i],lat->Analysis[DM_Plq_temp].CPU_data[i],lat->Analysis[DM_S_spat].CPU_data[i],lat->Analysis[DM_S_temp].CPU_data[i],lat->Analysis[DM_Polyakov_loop].CPU_data[i],lat->Analysis[DM_Polyakov_loop_im].CPU_data[i]);
        } else {
            if ((i<5)||(i==(lat->ITER - 1))) printf("[%5u]: % 10.8f \t % 10.8f \t % 10.8f \t % 10.8f\n",
                    (i),lat->Analysis[DM_S_spat].CPU_data[i],lat->Analysis[DM_S_temp].CPU_data[i],lat->Analysis[DM_Polyakov_loop].CPU_data[i],lat->Analysis[DM_Polyakov_loop_im].CPU_data[i]);
        }
    }

    pl            /= (lat->ITER - 1);
    pl_im         /= (lat->ITER - 1);

    printf("CPU Mean %s:  % 10.8f\n",lat->Analysis[DM_S_spat].data_name,lat->Analysis[DM_S_spat].CPU_mean_value);
    printf("CPU Mean %s:  % 10.8f\n",lat->Analysis[DM_S_temp].data_name,lat->Analysis[DM_S_temp].CPU_mean_value);
    printf("CPU Mean %s: % 10.8f\n",lat->Analysis[DM_S_total].data_name,lat->Analysis[DM_S_total].CPU_mean_value);

    printf("CPU Mean %s:  % 10.8f\n",lat->Analysis[DM_Plq_spat].data_name,lat->Analysis[DM_Plq_spat].CPU_mean_value);
    printf("CPU Mean %s:  % 10.8f\n",lat->Analysis[DM_Plq_temp].data_name,lat->Analysis[DM_Plq_temp].CPU_mean_value);
    printf("CPU Mean %s: % 10.8f\n",lat->Analysis[DM_Plq_total].data_name,lat->Analysis[DM_Plq_total].CPU_mean_value);

    printf("CPU Mean Polyakov loop:   % 10.8f\n",pl);
    printf("CPU Mean Polyakov loop im:% 10.8f\n",pl_im);


}

void            SU::lattice_simulate(model* lat,unsigned int* lattice_pointer){
    if (lattice_pointer != NULL){
        double* plaq;
        double* plaq_plq;
        double* pl_avr;
        double beta_effective = lat->BETA / lat->lattice_group;
        int prng_index = 0;

        lat->Analysis[DM_S_spat].CPU_data            = (double*) calloc((lat->ITER),sizeof(double));
        lat->Analysis[DM_S_temp].CPU_data            = (double*) calloc((lat->ITER),sizeof(double));
        lat->Analysis[DM_Plq_spat].CPU_data          = (double*) calloc((lat->ITER),sizeof(double));
        lat->Analysis[DM_Plq_temp].CPU_data          = (double*) calloc((lat->ITER),sizeof(double));
        lat->Analysis[DM_Polyakov_loop].CPU_data     = (double*) calloc((lat->ITER),sizeof(double));
        lat->Analysis[DM_Polyakov_loop_im].CPU_data  = (double*) calloc((lat->ITER),sizeof(double));
        lat->Analysis[DM_Polyakov_loop_P2].CPU_data  = (double*) calloc((lat->ITER),sizeof(double));
        lat->Analysis[DM_Polyakov_loop_P4].CPU_data  = (double*) calloc((lat->ITER),sizeof(double));

        lattice_data = (su_3*)  calloc(lat->lattice_table_row_size * lat->lattice_nd, sizeof(su_3));
        lattice_load_cpu(lat,lattice_pointer);    // Load lattice state to CPU

        int idx = 0;
            plaq     = lattice_avr_plaquette_cpu(lat);        // Lattice measurement
                lat->Analysis[DM_S_spat].CPU_data[idx] = plaq[0];
                lat->Analysis[DM_S_temp].CPU_data[idx] = plaq[1];

        if (lat->get_plaquettes_avr){
            plaq_plq = lattice_avr_plaquette_plq_cpu(lat);    // Lattice measurement (plaquettes)
                lat->Analysis[DM_Plq_spat].CPU_data[idx] = plaq_plq[0];
                lat->Analysis[DM_Plq_temp].CPU_data[idx] = plaq_plq[1];
        }
            pl_avr   = lattice_avr_Polyakov_loop_cpu(lat);    // Lattice Polyakov loop measurement
                lat->Analysis[DM_Polyakov_loop].CPU_data[idx]    = pl_avr[0];
                lat->Analysis[DM_Polyakov_loop_im].CPU_data[idx] = pl_avr[1];
                lat->Analysis[DM_Polyakov_loop_P2].CPU_data[idx] = pl_avr[2];
                lat->Analysis[DM_Polyakov_loop_P4].CPU_data[idx] = pl_avr[3];

            idx++;

        su_3 matrix_0,matrix_1,matrix_2;
        su_3 staple;
        coords_4 coords;
    printf("\n");


    for (int i=0; i<lat->NAV; i++){
        printf("\rCPU thermalization [%u]",i);

        for (int dir1 = 0; dir1 < lat->lattice_nd; dir1++) {
            for (int x4 = 0; x4 < lat->lattice_domain_size[3]; x4++)
            for (int x3 = 0; x3 < lat->lattice_domain_size[2]; x3++)
            for (int x2 = 0; x2 < lat->lattice_domain_size[1]; x2++)
            for (int x1 = 0; x1 < lat->lattice_domain_size[0]; x1++) {
                coords.x = x1; coords.y = x2; coords.z = x3; coords.t = x4;
                unsigned int gdi = lattice_coords_to_gid(lat,coords);
                matrix_0 = lattice_table_3(lat,coords,gdi,dir1);
                staple   = lattice_get_staple_cpu(lat,coords,dir1);
                matrix_2 = lattice_heatbath3(lat,staple,matrix_0,beta_effective,lat->prng_pointer,prng_index);
                matrix_1 = lattice_GramSchmidt_3(matrix_2);
                lattice_store_3(lat,matrix_1,gdi,dir1);
            }
        }
    }

    for (int i=1; i<lat->ITER; i++){ // zero measurement - on initial configuration!
        for (int j=0; j<lat->NITER; j++){
        printf("\rCPU working iteration [%u]",i);

        for (int dir1 = 0; dir1 < lat->lattice_nd; dir1++)
            for (int x4 = 0; x4 < lat->lattice_domain_size[3]; x4++)
            for (int x3 = 0; x3 < lat->lattice_domain_size[2]; x3++)
            for (int x2 = 0; x2 < lat->lattice_domain_size[1]; x2++)
            for (int x1 = 0; x1 < lat->lattice_domain_size[0]; x1++) {
                coords.x = x1; coords.y = x2; coords.z = x3; coords.t = x4;
                unsigned int gdi = lattice_coords_to_gid(lat,coords);
                matrix_0 = lattice_table_3(lat,coords,gdi,dir1);
                staple   = lattice_get_staple_cpu(lat,coords,dir1);
                matrix_2 = lattice_heatbath3(lat,staple,matrix_0,beta_effective,lat->prng_pointer,prng_index);
                matrix_1 = lattice_GramSchmidt_3(matrix_2);
                lattice_store_3(lat,matrix_1,gdi,dir1);
            }
        }

        plaq   = lattice_avr_plaquette_cpu(lat);        // Lattice measurement
            lat->Analysis[DM_S_spat].CPU_data[idx] = plaq[0];
            lat->Analysis[DM_S_temp].CPU_data[idx] = plaq[1];
        if (lat->get_plaquettes_avr){
            plaq_plq = lattice_avr_plaquette_plq_cpu(lat);  // Lattice measurement (plaquettes)
                lat->Analysis[DM_Plq_spat].CPU_data[idx] = plaq_plq[0];
                lat->Analysis[DM_Plq_temp].CPU_data[idx] = plaq_plq[1];
        }
        pl_avr = lattice_avr_Polyakov_loop_cpu(lat);    // Lattice Polyakov loop measurement
            lat->Analysis[DM_Polyakov_loop].CPU_data[idx]    = pl_avr[0];
            lat->Analysis[DM_Polyakov_loop_im].CPU_data[idx] = pl_avr[1];
            lat->Analysis[DM_Polyakov_loop_P2].CPU_data[idx] = pl_avr[2];
            lat->Analysis[DM_Polyakov_loop_P4].CPU_data[idx] = pl_avr[3];

        idx++;
    }
    printf("\n\nCPU simulations are done (%f seconds)\n",lat->GPU0->get_timer_CPU(2));

    lat->D_A->lattice_data_analysis_joint_CPU(&lat->Analysis[DM_S_total],  &lat->Analysis[DM_S_spat],  &lat->Analysis[DM_S_temp]);
    lat->D_A->lattice_data_analysis_joint_CPU(&lat->Analysis[DM_Plq_total],&lat->Analysis[DM_Plq_spat],&lat->Analysis[DM_Plq_temp]);

    lattice_analysis_cpu(lat);

    delete[] plaq;
    delete[] pl_avr;

    }
}


}
