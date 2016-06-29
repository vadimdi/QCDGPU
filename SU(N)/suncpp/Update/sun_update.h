/******************************************************************************
 * @file     sun_update.h
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.6
 *
 * @brief    [QCDGPU]
 *           Defines sequential SU(N) update procedures
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

#ifndef sun_update_h
#define sun_update_h

#include "../sunh.h"

template <typename su_n>
su_n staple(modelCPU<su_n> *latCPU, int gid, int dir){
    su_n stap, stap1;
    lattice_zero(&stap);
    
    coords_4 lsize;
    lsize.x = latCPU->lattice_size[0];
    lsize.y = latCPU->lattice_size[1];
    lsize.z = latCPU->lattice_size[2];
    lsize.t = latCPU->lattice_size[3];
    
    su_n m1, m2, m3;
    
    int nd = latCPU->lattice_ndCPU;
    int gid1;
    
    for (int dir1 = 0; dir1 < nd; dir1++)
        if(dir1 != dir){
            gid1 = lattice_neighbours_coords(lsize, gid, dir);
            m1 = latCPU->lattice_tableCPU[gid1 * nd + dir1];
            gid1 = lattice_neighbours_coords(lsize, gid, dir1);
            m2 = Herm(latCPU->lattice_tableCPU[gid1 * nd + dir]);
            m3 = Herm(latCPU->lattice_tableCPU[gid * nd + dir1]);
            stap = stap + (m1 * m2 * m3);
            
            gid1 = lattice_neighbours_coords_backward(lsize, gid, dir1);
            m3 = latCPU->lattice_tableCPU[gid1 * nd + dir1];
            m2 = Herm(latCPU->lattice_tableCPU[gid1 * nd + dir]);
            gid1 = lattice_neighbours_coords(lsize, gid1, dir);
            m1 = Herm(latCPU->lattice_tableCPU[gid1 * nd + dir1]);
            stap = stap + (m1 * m2 * m3);
        }
    
    return stap;
}

template <typename su_n>
void lattice_update_even(modelCPU<su_n> *latCPU, int dir, PRNG_CL::PRNG *prngCPU){
    su_n stap, U;
    int gid;
    
    coords_4 lsize;
    lsize.x = latCPU->lattice_size[0];
    lsize.y = latCPU->lattice_size[1];
    lsize.z = latCPU->lattice_size[2];
    lsize.t = latCPU->lattice_size[3];
    
    for (int i = 0; i < latCPU->lattice_sitesCPU / 2; i++){
        gid = lattice_even_gid(lsize, i);
        stap = staple(latCPU, gid, dir);
        U = latCPU->lattice_tableCPU[gid * latCPU->lattice_ndCPU + dir];
        update_link(&U, stap, (latCPU->beta / latCPU->lattice_group), latCPU->nhit, gid, (latCPU->ints == 2), latCPU->lattice_sitesCPU, 1, prngCPU);
        latCPU->lattice_tableCPU[gid * latCPU->lattice_ndCPU + dir] = U;
    }
}

template <typename su_n>
void lattice_update_odd(modelCPU<su_n> *latCPU, int dir, PRNG_CL::PRNG *prngCPU){
    su_n stap, U;
    int gid;
    
    coords_4 lsize;
    lsize.x = latCPU->lattice_size[0];
    lsize.y = latCPU->lattice_size[1];
    lsize.z = latCPU->lattice_size[2];
    lsize.t = latCPU->lattice_size[3];
    
    for (int i = 0; i < latCPU->lattice_sitesCPU / 2; i++){
        gid = lattice_odd_gid(lsize, i);
        stap = staple(latCPU, gid, dir);
        U = latCPU->lattice_tableCPU[gid * latCPU->lattice_ndCPU + dir];
        update_link(&U, stap, (latCPU->beta / latCPU->lattice_group), latCPU->nhit, gid, (latCPU->ints == 2), latCPU->lattice_sitesCPU, 1, prngCPU);
        latCPU->lattice_tableCPU[gid * latCPU->lattice_ndCPU + dir] = U;
    }
}

#endif
