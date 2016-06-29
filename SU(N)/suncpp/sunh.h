/******************************************************************************
 * @file     sunh.h
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.6
 *
 * @brief    [QCDGPU]
 *           Defines classes and procedures of lattice initialization and deletion for sequential CPU run of QCDGPU
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

#ifndef sunh_h
#define sunh_h

#include "../suncl/suncl.h"
#include "../kernel/complex.h"
#include "../random/random.h"

#include "coord_work/coord_work.h"

#include "su2/algebra_su2.h"
#include "su2/update_su2.h"
#include "su3/algebra_su3.h"
#include "su3/update_su3.h"

#define X   0
#define Y   1
#define Z   2
#define T   3

#define N_MEAS_QUANTITIES 2

template <typename su_n>
class modelCPU{
public:
    unsigned int lattice_ndCPU;
    int*    lattice_size;
    int     ints;
    int     nhit;
    hgpu_float beta;
    int     lattice_group;
    int nav;
    int iter;
    int niter;
    
    unsigned int lattice_sitesCPU;
    
    su_n *lattice_tableCPU;
    
    void create_latticeCPU(void);
    void delete_latticeCPU(void);
    
    void lattice_initializeCPU(void);
};

class Measurements{
public:
    int iter;
    
    hgpu_complex   *cplq;
    hgpu_double   *tplq;
    
    hgpu_complex   *cs;
    hgpu_double   *ts;
    
    int mask[N_MEAS_QUANTITIES]; //length = number of implemented measurements
    
    Measurements(void){
        for (int k = 0; k < N_MEAS_QUANTITIES; k++) mask[k] = 0;
    };
    ~Measurements(void){};
};
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename su_n>
void modelCPU<su_n>::create_latticeCPU(void){
    lattice_tableCPU = (su_n*)calloc(lattice_ndCPU * lattice_sitesCPU, sizeof(su_n));
};

template <typename su_n>
void modelCPU<su_n>::delete_latticeCPU(void){
    free(lattice_tableCPU);
};

template <typename su_n>
void modelCPU<su_n>::lattice_initializeCPU(void){
    switch (ints){
        case 0: 
            printf("Hot start\n");
            break;
        case 1: 
            {
                for (int gid = 0; gid < lattice_ndCPU * lattice_sitesCPU; gid++)
                    lattice_unity(&lattice_tableCPU[gid]);
            }
            break;
        case 2 :
            {
                for (int gid = 0; gid < lattice_sitesCPU; gid++)
                    for (int dir = 0; dir < lattice_ndCPU; dir++)
                        lattice_matrixGID(&lattice_tableCPU[gid * lattice_ndCPU + dir], gid, dir, lattice_sitesCPU);
            }
            break;
    }
};
#endif
