/******************************************************************************
 * @file     analysis_cpp.cpp
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.6
 *
 * @brief    [QCDGPU]
 *           Averaging of measured quantities over run
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

#include "analysis_cpp.h"

void lattice_analysis_cpp(Measurements *meas, data_analysis_cpp *Analysis){
    int i = 0;
    if (meas->mask[0]){ //measurement plq
        Analysis[i].data_name = "Plq_spat";
        Analysis[i + 1].data_name = "Plq_temp";
        Analysis[i + 2].data_name = "Plq_total";
        
        for (int k = 0; k < 3; k++)
            Analysis[i + k].data = (double*)calloc(meas->iter, sizeof(double));
        for (int j = 0; j < meas->iter; j++){
            Analysis[i].data[j] = meas->cplq[j].re;
            Analysis[i + 1].data[j] = meas->cplq[j].im;
            Analysis[i + 2].data[j] = meas->tplq[j];
        }
        
        Analysis[i + 2].data_size = meas->iter;
        lattice_estimates(&Analysis[i + 2]);
        
        i += 3;
    }
    
    if (meas->mask[1]){//measurement action
        Analysis[i].data_name = "S_spat";
        Analysis[i + 1].data_name = "S_temp";
        Analysis[i + 2].data_name = "S_total";
        
        for (int k = 0; k < 3; k++)
            Analysis[i + k].data = (double*)calloc(meas->iter, sizeof(double));
        for (int j = 0; j < meas->iter; j++){
            Analysis[i].data[j] = meas->cs[j].re;
            Analysis[i + 1].data[j] = meas->cs[j].im;
            Analysis[i + 2].data[j] = meas->ts[j];
        }
        
        Analysis[i + 2].data_size = meas->iter;
        lattice_estimates(&Analysis[i + 2]);
        
        i += 3;
    }
}

void lattice_estimates (data_analysis_cpp *Analysis){
    double mean = 0.0;
    double variance = 0.0;
    
    int length = Analysis->data_size - 1;
    
    for (unsigned int i = 1; i < Analysis->data_size; i++)
        mean += Analysis->data[i] / length;
    for (unsigned int i = 1; i < Analysis->data_size; i++)
        variance += pow((Analysis->data[i] - mean), 2) / length;
    Analysis->mean_value = mean;
    Analysis->variance = variance;
}
