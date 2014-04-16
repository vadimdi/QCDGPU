/******************************************************************************
 * @file     data_analysis.cpp
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.4
 *
 * @brief    [QCDGPU]
 *           Data analysis block
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

#include "data_analysis.h"

#define VDELTA_DOUBLE  0.00000000001 // accuracy for results verification (double precision)
#define VDELTA_SINGLE  0.00005       // accuracy for results verification (single precision)
#define CHECKING_PRECISION_SINGLE   (1E-6)  // precision for checking results (single precision)
#define CHECKING_PRECISION_DOUBLE   (1E-12) // precision for checking results (double precision)


#define min(a,b) (((a) < (b)) ? (a) : (b))

namespace analysis_CL{
using analysis_CL::analysis;

           bool analysis::results_verification = true;    // result of CPU-GPU verification

            analysis::data_analysis::data_analysis(void) {
                pointer             = NULL;
                data                = NULL;
                CPU_data            = NULL;
                data_size           = 0;
                pointer_offset      = 0;
                storage_type        = GPU_CL::GPU::GPU_storage_none;
                denominator         = 0.0;
                mean_value          = 0.0;
                CPU_mean_value      = 0.0;
                variance            = 0.0;
                GPU_last_value      = 0.0;
                CPU_last_value      = 0.0;
                CPU_last_variance   = 0.0;
                data_name           = NULL;
                precision_single    = false;
}
            analysis::data_analysis::~data_analysis(void) {
            free(data);
            free(CPU_data);
}

bool        analysis::CPU_GPU_verification_single(double a, double b, const char* err_str){
    bool result = true;
    if ((abs(a)>CHECKING_PRECISION_SINGLE) && (abs(b)>CHECKING_PRECISION_SINGLE))
        if (abs(1.0-b/a)>VDELTA_SINGLE) {
            if (err_str) printf("%s test failed!\n",err_str);
            result = false;
        }
    results_verification &= result;
    return result;
}
bool        analysis::CPU_GPU_verification_double(double a, double b, const char* err_str){
    bool result = true;
    if ((abs(a)>CHECKING_PRECISION_DOUBLE) && (abs(b)>CHECKING_PRECISION_DOUBLE))
        if (abs(1.0-b/a)>VDELTA_DOUBLE) {
            if (err_str) printf("%s test failed!\n",err_str);
            result = false;
        }
    results_verification &= result;
    return result;
}

double inline round(double d){
    return (d>0) ? floor(d + 0.5) : floor(d - 0.5);
}

void        analysis::lattice_data_analysis(data_analysis* data){
    int last_index = (data->data_size - 1);
    if(data->data == NULL) {
        data->data     = (double*) calloc(data->data_size,sizeof(double));
        data->mean_value      = 0.0;
        data->variance        = 0.0;
        double data_value     = 0.0;
        for (unsigned int i=0; i<data->data_size; i++) {
            if (data->storage_type==GPU_CL::GPU::GPU_storage_double2high)
                data_value = GPU_CL::GPU::convert_to_double(data->pointer[4*(i+data->pointer_offset)    ],data->pointer[4*(i+data->pointer_offset) + 1]);
            if (data->storage_type==GPU_CL::GPU::GPU_storage_double2low)
                data_value = GPU_CL::GPU::convert_to_double(data->pointer[4*(i+data->pointer_offset) + 2],data->pointer[4*(i+data->pointer_offset) + 3]);
            if (data->storage_type==GPU_CL::GPU::GPU_storage_double)
                data_value = GPU_CL::GPU::convert_to_double(data->pointer[2*(i+data->pointer_offset)    ],data->pointer[2*(i+data->pointer_offset) + 1]);

            data->data[i]  = data_value / data->denominator;
            if (i>0) data->mean_value += data->data[i];
        }
        if (last_index>0) data->mean_value /= (double) last_index;
        for (unsigned int i=1; i<data->data_size; i++)
            data->variance += pow(data->data[i] - data->mean_value,2);
        if (last_index>0) data->variance   /= (double) last_index;
        data->GPU_last_value = data->data[last_index];
    }
    if (data->precision_single)
        CPU_GPU_verification_single(data->GPU_last_value,data->CPU_last_value,data->data_name);
    else
        CPU_GPU_verification_double(data->GPU_last_value,data->CPU_last_value,data->data_name);
}
void        analysis::lattice_data_analysis_CPU(data_analysis* data){
    int last_index = (data->data_size - 1);
    if(data->CPU_data==NULL) {
        data->CPU_data = (double*) calloc(data->data_size,sizeof(double));
        data->CPU_mean_value  = 0.0;
        for (unsigned int i=1; i<data->data_size; i++)
            data->CPU_mean_value += data->CPU_data[i];
        if (last_index>0) data->CPU_mean_value /= (double) last_index;
        data->CPU_last_value = data->CPU_data[last_index];
    }
}
void        analysis::lattice_data_analysis_joint(data_analysis* data,data_analysis* data1,data_analysis* data2){
    data->data_size = min(data1->data_size,data2->data_size);
    data->storage_type = GPU_CL::GPU::GPU_storage_joint;
    int last_index = (data->data_size - 1);
    if (data->data==NULL) {
        data->data = (double*) calloc(data->data_size,sizeof(double));
        data->mean_value     = 0.0;
        data->CPU_mean_value = 0.0;
        data->variance       = 0.0;
        for (unsigned int i=0; i<data->data_size; i++) {
            data->data[i] = 0.5 * (data1->data[i] + data2->data[i]);
            if (i>0) data->mean_value += data->data[i];
        }
        if (last_index>0) data->mean_value /= (double) last_index;
        for (unsigned int i=1; i<data->data_size; i++) {
            data->variance += pow(data->data[i] - data->mean_value,2);
        }
        if (last_index>0) data->variance   /= (double) last_index;
        data->GPU_last_value = data->data[last_index];
    }
    if (data->precision_single)
        CPU_GPU_verification_single(data->GPU_last_value,data->CPU_last_value,data->data_name);
    else
        CPU_GPU_verification_double(data->GPU_last_value,data->CPU_last_value,data->data_name);
}
void        analysis::lattice_data_analysis_joint_CPU(data_analysis* data,data_analysis* data1,data_analysis* data2){
    data->data_size = min(data1->data_size,data2->data_size);
    data->storage_type = GPU_CL::GPU::GPU_storage_joint;
    int last_index = (data->data_size - 1);
    if(data->CPU_data==NULL) {
        data->CPU_data = (double*) calloc(data->data_size,sizeof(double));
        for (unsigned int i=0; i<data->data_size; i++) {
            data->CPU_data[i] += 0.5*(data1->CPU_data[i] + data2->CPU_data[i]);
            if (i>0) data->CPU_mean_value += data->CPU_data[i];
        }
        if (last_index>0) data->CPU_mean_value /= (double) last_index;
        data->CPU_last_value = data->CPU_data[last_index];
    }
}
void        analysis::lattice_data_analysis_joint3(data_analysis* data,data_analysis* data1,
                                                   data_analysis* data2,data_analysis* data3){
    data->data_size = min(min(data1->data_size,data2->data_size),data3->data_size);
    data->storage_type = GPU_CL::GPU::GPU_storage_joint;
    int last_index = (data->data_size - 1);
    if (data->data==NULL) {
        data->data = (double*) calloc(data->data_size,sizeof(double));
        data->mean_value     = 0.0;
        data->CPU_mean_value = 0.0;
        data->variance       = 0.0;
        for (unsigned int i=0; i<data->data_size; i++) {
            data->data[i] = sqrt (data1->data[i] * data1->data[i] + data2->data[i] * data2->data[i] + data3->data[i] * data3->data[i]);
            if (i>0) data->mean_value += data->data[i];
        }
        if (last_index>0) data->mean_value /= (double) last_index;
        for (unsigned int i=1; i<data->data_size; i++) {
            data->variance += pow(data->data[i] - data->mean_value,2);
        }
        if (last_index>0) data->variance   /= (double) last_index;
        data->GPU_last_value = data->data[last_index];
    }
    if (data->precision_single)
        CPU_GPU_verification_single(data->GPU_last_value,data->CPU_last_value,data->data_name);
    else
        CPU_GPU_verification_double(data->GPU_last_value,data->CPU_last_value,data->data_name);


}

}
