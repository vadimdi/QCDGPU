/******************************************************************************
 * @file     suncpp.h
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.6
 *
 * @brief    [QCDGPU]
 *           Defines sequential simulation procedure on CPU
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

#ifndef suncpp_h
#define suncpp_h

#include "sunh.h"

#include "../suncl/suncl.h"
#include "../kernel/complex.h"

#include "su2/algebra_su2.h"
#include "su3/algebra_su3.h"

#include "Update/sun_update.h"
#include "Measurements/Plq.h"
#include "Measurements/S.h"
#include "Measurements/analysis_cpp.h"

template <typename su_n>
void lattice_simulateCPU(model_CL::model *lat, su_n *smth);

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++

char*           get_current_datetime(time_t *ltime){
    time_t tim;
    time(&tim);
    char* result = (char*) calloc(256,sizeof(char));
#ifdef _WIN32
    ctime_s(result, 26, &tim);
#else
    sprintf(result,"%s",ctime((const time_t*) &tim));
#endif
    
    *ltime = tim;
    return result;
}

template <typename su_n>
void lattice_simulateCPU(model_CL::model *lat, su_n *smth){
        modelCPU<su_n> *latCPU = new(modelCPU<su_n>);
        
        latCPU->lattice_ndCPU = lat->lattice_nd;
        latCPU->lattice_size = (int*)calloc(latCPU->lattice_ndCPU, sizeof(int));
        for(int i = 0; i < latCPU->lattice_ndCPU; i++)
            latCPU->lattice_size[i] = lat->lattice_full_size[i];
        latCPU->ints = (int)lat->ints;
        latCPU->nhit = (int)lat->NHIT;
        latCPU->beta = (hgpu_float)lat->BETA;
        latCPU->lattice_group = (int)lat->lattice_group;
        latCPU->nav = (int)lat->NAV;
        latCPU->iter = (int)lat->ITER;
        latCPU->niter = (int)lat->NITER;
        
        int nmeas = 0;
        Measurements *meas = new(Measurements);
        meas->iter = (int)lat->ITER;
        if(lat->get_plaquettes_avr){
            meas->cplq = (hgpu_complex*)calloc(latCPU->iter, sizeof(hgpu_complex));
            meas->tplq = (hgpu_double*)calloc(latCPU->iter, sizeof(hgpu_double));
            nmeas += 3; //plq_temp, plq_spat, plq_tot
            meas->mask[0] = 1;
        }
         if(lat->get_actions_avr){
            meas->cs = (hgpu_complex*)calloc(latCPU->iter, sizeof(hgpu_complex));
            meas->ts = (hgpu_double*)calloc(latCPU->iter, sizeof(hgpu_double));
            nmeas += 3; //s_temp, s_spat, s_tot
            meas->mask[1] = 1;
        }
        
        data_analysis_cpp *Analysis = (data_analysis_cpp*)calloc(nmeas, sizeof(data_analysis_cpp));
        
        latCPU->lattice_sitesCPU = latCPU->lattice_size[0];
        for (int i = 1; i < latCPU->lattice_ndCPU; i++)
            latCPU->lattice_sitesCPU *= latCPU->lattice_size[i];
        
        PRNG_CL::PRNG *prngCPU = new(PRNG_CL::PRNG);
        prngCPU->PRNG_generator = lat->PRNG0->PRNG_generator;
        prngCPU->PRNG_randseries = lat->PRNG0->PRNG_randseries;
        prngCPU->PRNG_precision = lat->PRNG0->PRNG_precision;
        prngCPU->initialize_CPU();
        lat->PRNG0->initialize_CPU();
        
        char* header = lat->lattice_make_header();
        printf("%s\n",header);
        
        lat->timestart = get_current_datetime(&lat->ltimestart);
        printf("\nCPU siulations are started\n");
        
        latCPU->create_latticeCPU();
        latCPU->lattice_initializeCPU();

        if(lat->get_plaquettes_avr)
            meas->cplq[0] = plqConf(latCPU, &meas->tplq[0]);
        if(lat->get_actions_avr)
            meas->cs[0] = sConf(latCPU, &meas->ts[0]);
        
        for(int n = 0; n < latCPU->nav; n++){
            lattice_update_odd(latCPU, X, lat->PRNG0);
            lattice_update_odd(latCPU, Y, lat->PRNG0);
            lattice_update_odd(latCPU, Z, lat->PRNG0);
            lattice_update_odd(latCPU, T, lat->PRNG0);
            
            lattice_update_even(latCPU, X, lat->PRNG0);
            lattice_update_even(latCPU, Y, lat->PRNG0);
            lattice_update_even(latCPU, Z, lat->PRNG0);
            lattice_update_even(latCPU, T, lat->PRNG0);
            
            if (n % 10 == 0) printf("\rCPU thermalization [%i]", n);
        }
        
        for(int i = 1; i < latCPU->iter; i++){
            for(int n = 0; n < latCPU->niter; n++){
                lattice_update_odd(latCPU, X, lat->PRNG0);
                lattice_update_odd(latCPU, Y, lat->PRNG0);
                lattice_update_odd(latCPU, Z, lat->PRNG0);
                lattice_update_odd(latCPU, T, lat->PRNG0);
                
                lattice_update_even(latCPU, X, lat->PRNG0);
                lattice_update_even(latCPU, Y, lat->PRNG0);
                lattice_update_even(latCPU, Z, lat->PRNG0);
                lattice_update_even(latCPU, T, lat->PRNG0);
            }
            
            if(lat->get_plaquettes_avr)
                meas->cplq[i] = plqConf(latCPU, &meas->tplq[i]);
            if(lat->get_actions_avr)
                meas->cs[i] = sConf(latCPU, &meas->ts[i]);
            
            if (i % 10 == 0) printf("\rCPU working iteration [%u]",i);
        }
        
        lat->timeend = get_current_datetime(&lat->ltimeend);
        printf("\rCPU simulations are done\n");
        
        lattice_analysis_cpp(meas, Analysis);
        
        int ii = 0;
        if (meas->mask[0]){
            lat->Analysis[DM_Plq_spat].data = (double*)calloc(meas->iter, sizeof(double));
            lat->Analysis[DM_Plq_temp].data = (double*)calloc(meas->iter, sizeof(double));
            lat->Analysis[DM_Plq_total].data = (double*)calloc(meas->iter, sizeof(double));
            
            memcpy(lat->Analysis[DM_Plq_spat].data, Analysis[ii].data, meas->iter * sizeof(double));
            memcpy(lat->Analysis[DM_Plq_temp].data, Analysis[ii + 1].data, meas->iter * sizeof(double));
            memcpy(lat->Analysis[DM_Plq_total].data, Analysis[ii + 2].data, meas->iter * sizeof(double));
            
            lat->Analysis[DM_Plq_spat].data_name = Analysis[ii].data_name;
            lat->Analysis[DM_Plq_spat].mean_value = Analysis[ii].mean_value;
            lat->Analysis[DM_Plq_spat].variance = Analysis[ii].variance;
            
            lat->Analysis[DM_Plq_temp].data_name = Analysis[ii + 1].data_name;
            lat->Analysis[DM_Plq_temp].mean_value = Analysis[ii + 1].mean_value;
            lat->Analysis[DM_Plq_temp].variance = Analysis[ii + 1].variance;
            
            lat->Analysis[DM_Plq_total].data_name = Analysis[ii + 2].data_name;
            lat->Analysis[DM_Plq_total].mean_value = Analysis[ii + 2].mean_value;
            lat->Analysis[DM_Plq_total].variance = Analysis[ii + 2].variance;
            
            ii += 3;
        }
        
         if (meas->mask[1]){
            lat->Analysis[DM_S_spat].data = (double*)calloc(meas->iter, sizeof(double));
            lat->Analysis[DM_S_temp].data = (double*)calloc(meas->iter, sizeof(double));
            lat->Analysis[DM_S_total].data = (double*)calloc(meas->iter, sizeof(double));
            
            memcpy(lat->Analysis[DM_S_spat].data, Analysis[ii].data, meas->iter * sizeof(double));
            memcpy(lat->Analysis[DM_S_temp].data, Analysis[ii + 1].data, meas->iter * sizeof(double));
            memcpy(lat->Analysis[DM_S_total].data, Analysis[ii + 2].data, meas->iter * sizeof(double));
            
            lat->Analysis[DM_S_spat].data_name = Analysis[ii].data_name;
            lat->Analysis[DM_S_spat].mean_value = Analysis[ii].mean_value;
            lat->Analysis[DM_S_spat].variance = Analysis[ii].variance;
            
            lat->Analysis[DM_S_temp].data_name = Analysis[ii + 1].data_name;
            lat->Analysis[DM_S_temp].mean_value = Analysis[ii + 1].mean_value;
            lat->Analysis[DM_S_temp].variance = Analysis[ii + 1].variance;
            
            lat->Analysis[DM_S_total].data_name = Analysis[ii + 2].data_name;
            lat->Analysis[DM_S_total].mean_value = Analysis[ii + 2].mean_value;
            lat->Analysis[DM_S_total].variance = Analysis[ii + 2].variance;
            
            ii += 3;
        }

        delete (prngCPU);
        
        free(Analysis);
        
         if (lat->get_actions_avr){
            free(meas->cs);
            free(meas->ts);
        }
        if (lat->get_plaquettes_avr){
            free(meas->cplq);
            free(meas->tplq);
        }
        delete (meas);
        
        latCPU->delete_latticeCPU();
        free(latCPU->lattice_size);
        delete(latCPU);
}
#endif
