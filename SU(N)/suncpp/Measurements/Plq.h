/******************************************************************************
 * @file     Plq.h
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.6
 *
 * @brief    [QCDGPU]
 *           Measurement of plaquette, averaged over configuration
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

#ifndef plq_h
#define plq_h

#include "../sunh.h"

#include "../su2/algebra_su2.h"
#include "../su3/algebra_su3.h"

template <typename su_n>
hgpu_complex plqConf(modelCPU<su_n> *latCPU, hgpu_double *pplq){
    hgpu_complex result;
    result.re = 0.0; //spat
    result.im = 0.0;//temp
    
    int nd = latCPU->lattice_ndCPU;
    
    su_n m1, m2, m3, m4;
    int gid1;
    
    coords_4 lsize;
    lsize.x = latCPU->lattice_size[0];
    lsize.y = latCPU->lattice_size[1];
    lsize.z = latCPU->lattice_size[2];
    lsize.t = latCPU->lattice_size[3];
    
    for (int gid = 0; gid < latCPU->lattice_sitesCPU; gid++){
    //--- spat ---------------------------------------------------------------------
    for (int dir = 0; dir < nd - 2; dir++)
        for (int dir1 = dir + 1; dir1 < nd - 1; dir1++){
            m1 = latCPU->lattice_tableCPU[gid * nd + dir];
            gid1 = lattice_neighbours_coords(lsize, gid, dir);
            m2 = latCPU->lattice_tableCPU[gid1 * nd + dir1];
            gid1 = lattice_neighbours_coords(lsize, gid, dir1);
            m3 = Herm(latCPU->lattice_tableCPU[gid1 * nd + dir]);
            m4 = Herm(latCPU->lattice_tableCPU[gid * nd + dir1]);
        
        result.re += ReTr(m1 * m2 * m3 * m4);
    }
    //--- temp ---------------------------------------------------------------------
    for (int dir = 0; dir < nd - 1; dir++){
        m1 = latCPU->lattice_tableCPU[gid * nd + dir];
        gid1 = lattice_neighbours_coords(lsize, gid, dir);
        m2 = latCPU->lattice_tableCPU[gid1 * nd + nd - 1];
        gid1 = lattice_neighbours_coords(lsize, gid, nd - 1);
        m3 = Herm(latCPU->lattice_tableCPU[gid1 * nd + dir]);
        m4 = Herm(latCPU->lattice_tableCPU[gid * nd + nd - 1]);
        
        result.im += ReTr(m1 * m2 * m3 * m4);
    }
    //-------------------------------------------------------------------------------
    }
    
    *pplq = (result.re + result.im) / ((nd - 1) * nd / 2 * latCPU->lattice_sitesCPU);
    result.re /= ((nd - 2) * (nd - 1) / 2 * latCPU->lattice_sitesCPU);
    result.im /= ((nd - 1)* latCPU->lattice_sitesCPU);
    return result;
}
#endif
