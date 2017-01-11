/******************************************************************************
 * @file     data_analysis.cpp
 * @author   Vadim Demchik <vadimdi@yahoo.com>
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.0
 *
 * @brief    [QCDGPU]
 *           Data analysis module
 *
 * @section  LICENSE
 *
 * Copyright (c) 2013, Vadim Demchik, Natalia Kolomoyets
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

#define VDELTA_DOUBLE  0.0000000001 // accuracy for results verification (double precision)
#define VDELTA_SINGLE  0.00005      // accuracy for results verification (single precision)

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
                skip_checking       = false;
                number_of_series    = 1;
                part_number         = 0;
}
            analysis::data_analysis::~data_analysis(void) {
            free(data);
            free(CPU_data);
}

bool        analysis::CPU_GPU_verification_single(double a, double b, const char* err_str){
    bool result = true;
    if ((a != a)||(b != b)||(abs(1.0-b/a)>VDELTA_SINGLE)) {
        if (err_str) printf("%s test failed!\n",err_str);
        result = false;
    }
    results_verification &= result;
    return result;
}
bool        analysis::CPU_GPU_verification_double(double a, double b, const char* err_str){
    bool result = true;
    if ((a != a)||(b != b)||(abs(1.0-b/a)>VDELTA_DOUBLE)) {
        if (err_str) printf("%s test failed!\n",err_str);
        result = false;
    }
    results_verification &= result;
    return result;
}

bool        analysis::lattice_data_verify(data_analysis* data){
    bool result = true;
    if (!data->skip_checking) {
        if (data->precision_single)
            result = CPU_GPU_verification_single(data->GPU_last_value,data->CPU_last_value,data->data_name);
        else
            result = CPU_GPU_verification_double(data->GPU_last_value,data->CPU_last_value,data->data_name);
    }
    return result;
}

double inline round(double d){
    return (d>0) ? floor(d + 0.5) : floor(d - 0.5);
}

void        analysis::lattice_data_analysis(data_analysis* data){
    int last_index = (data->data_size - 1);
    data->part_number++;
    if(data->data == NULL) data->data     = (double*) calloc((data->data_size * data->number_of_series)+1,sizeof(double));

    double data_value     = 0.0;
    data->mean_value      = 0.0;
    data->variance        = 0.0;

    for (unsigned int i=0; i<data->data_size; i++) {
        if (data->storage_type==GPU_CL::GPU::GPU_storage_double2high)
            data_value = GPU_CL::GPU::convert_to_double(data->pointer[4*(i+data->pointer_offset)    ],data->pointer[4*(i+data->pointer_offset) + 1]);
        if (data->storage_type==GPU_CL::GPU::GPU_storage_double2low)
            data_value = GPU_CL::GPU::convert_to_double(data->pointer[4*(i+data->pointer_offset) + 2],data->pointer[4*(i+data->pointer_offset) + 3]);
        if (data->storage_type==GPU_CL::GPU::GPU_storage_double)
            data_value = GPU_CL::GPU::convert_to_double(data->pointer[2*(i+data->pointer_offset)    ],data->pointer[2*(i+data->pointer_offset) + 1]);

        data->data[i + (data->part_number-1)*data->data_size] = data_value / data->denominator;
    }
    for (unsigned int j=0;j<data->part_number;j++)
        for (unsigned int i=1; i<data->data_size; i++)
            data->mean_value += data->data[i + j*data->data_size];
    if (last_index>0) data->mean_value /= (double) (last_index * data->part_number);

    for (unsigned int j=0;j<data->part_number;j++)
        for (unsigned int i=1; i<data->data_size; i++)
            data->variance += pow(data->data[i + j*data->data_size] - data->mean_value,2);
    if (last_index>0) data->variance /= (double)  (last_index * data->part_number);

    data->GPU_last_value = data_value / data->denominator;
    lattice_data_verify(data);
}
void        analysis::lattice_data_analysis_sub(data_analysis* data,data_analysis* data2){
    int last_index = (data->data_size - 1);
    data->part_number++;
    if(data->data == NULL) data->data     = (double*) calloc((data->data_size * data->number_of_series)+1,sizeof(double));

    double data_value     = 0.0;
    double data2_value    = 0.0;
    data->mean_value      = 0.0;
    data->variance        = 0.0;

    for (unsigned int i=0; i<data->data_size; i++) {
        if (data->storage_type==GPU_CL::GPU::GPU_storage_double2high)
            data_value = GPU_CL::GPU::convert_to_double(data->pointer[4*(i+data->pointer_offset)    ],data->pointer[4*(i+data->pointer_offset) + 1]);
        if (data->storage_type==GPU_CL::GPU::GPU_storage_double2low)
            data_value = GPU_CL::GPU::convert_to_double(data->pointer[4*(i+data->pointer_offset) + 2],data->pointer[4*(i+data->pointer_offset) + 3]);
        if (data->storage_type==GPU_CL::GPU::GPU_storage_double)
            data_value = GPU_CL::GPU::convert_to_double(data->pointer[2*(i+data->pointer_offset)    ],data->pointer[2*(i+data->pointer_offset) + 1]);

        data2_value = data2->data[i+(data->part_number-1)*data2->data_size];
        data->data[i + (data->part_number-1)*data->data_size] = data_value / data->denominator - data2_value * data2_value;
    }
    for (unsigned int j=0;j<data->part_number;j++)
        for (unsigned int i=1; i<data->data_size; i++)
            data->mean_value += data->data[i + j*data->data_size];
    if (last_index>0) data->mean_value /= (double) (last_index * data->part_number);

    for (unsigned int j=0;j<data->part_number;j++)
        for (unsigned int i=1; i<data->data_size; i++)
            data->variance += pow(data->data[i + j*data->data_size] - data->mean_value,2);
    if (last_index>0) data->variance /= (double)  (last_index * data->part_number);

    data->GPU_last_value  = data_value / data->denominator - data2_value * data2_value;
    data->CPU_last_value -= data2_value * data2_value;
    lattice_data_verify(data);
}

void        analysis::lattice_data_analysis_CPU(data_analysis* data){
    int last_index = (data->data_size - 1);
    if(data->CPU_data==NULL) data->CPU_data = (double*) calloc((data->data_size * data->number_of_series)+1,sizeof(double));
    
    data->CPU_mean_value  = 0.0;
    for (unsigned int j=0;j<data->part_number;j++)
        for (unsigned int i=1; i<data->data_size; i++)
            data->CPU_mean_value += data->CPU_data[i + j*data->data_size];
    if (last_index>0) data->CPU_mean_value /= (double) (last_index * data->part_number);
    data->CPU_last_value = data->CPU_data[last_index];
}
void        analysis::lattice_data_analysis_joint(data_analysis* data,data_analysis* data1,data_analysis* data2){
    data->data_size = _MIN(data1->data_size,data2->data_size);
    data->number_of_series = _MAX(data1->number_of_series,data2->number_of_series);
    data->part_number = _MIN(data1->part_number,data2->part_number);
    data->storage_type = GPU_CL::GPU::GPU_storage_joint;
    int last_index = (data->data_size - 1);
    if (data->data==NULL) data->data = (double*) calloc((data->data_size*data->number_of_series)+1,sizeof(double));

    data->mean_value     = 0.0;
    data->CPU_mean_value = 0.0;
    data->variance       = 0.0;
    for (unsigned int j=0;j<data->part_number;j++)
        for (unsigned int i=0; i<data->data_size; i++) {
            data->data[i+j*data->data_size] = 0.5 * (data1->data[i+j*data->data_size] + data2->data[i+j*data->data_size]);
            if (i>0) data->mean_value += data->data[i+j*data->data_size];
        }
    if (last_index>0) data->mean_value /= (double) (last_index * data->part_number);

    for (unsigned int j=0;j<data->part_number;j++)
        for (unsigned int i=1; i<data->data_size; i++)
            data->variance += pow(data->data[i+j*data->data_size] - data->mean_value,2);
    if (last_index>0) data->variance   /= (double) (last_index * data->part_number);
    data->GPU_last_value = data->data[last_index + (data->part_number-1)*data->data_size];
    lattice_data_verify(data);
}
void        analysis::lattice_data_analysis_joint_variance(data_analysis* data,data_analysis* data1,data_analysis* data2){

    data->data_size = _MIN(data1->data_size,data2->data_size);
    data->number_of_series = _MAX(data1->number_of_series,data2->number_of_series);
    data->part_number = _MIN(data1->part_number,data2->part_number);
    data->storage_type = GPU_CL::GPU::GPU_storage_joint;
    int last_index = (data->data_size - 1);
    if (data->data==NULL) data->data = (double*) calloc((data->data_size*data->number_of_series)+1,sizeof(double));

    data->mean_value     = 0.0;
    data->CPU_mean_value = 0.0;
    data->variance       = 0.0;
    for (unsigned int j=0;j<data->part_number;j++)
        for (unsigned int i=0; i<data->data_size; i++) {
            data->data[i+j*data->data_size] = (data2->data[i+j*data->data_size] - data1->data[i+j*data->data_size]*data1->data[i+j*data->data_size]);
        }
    data->mean_value = data2->mean_value - data1->mean_value*data1->mean_value;

    for (unsigned int j=0;j<data->part_number;j++)
        for (unsigned int i=1; i<data->data_size; i++)
            data->variance += pow(data->data[i+j*data->data_size] - data->mean_value,2);
    if (last_index>0) data->variance   /= (double) (last_index * data->part_number);
    data->GPU_last_value = data->data[last_index + (data->part_number-1)*data->data_size];
}

void        analysis::lattice_data_analysis_joint_CPU(data_analysis* data,data_analysis* data1,data_analysis* data2){
    data->data_size = _MIN(data1->data_size,data2->data_size);
    data->number_of_series = _MAX(data1->number_of_series,data2->number_of_series);
    data->part_number = _MIN(data1->part_number,data2->part_number);
    data->storage_type = GPU_CL::GPU::GPU_storage_joint;
    int last_index = (data->data_size - 1);
    if(data->CPU_data==NULL) data->CPU_data = (double*) calloc((data->data_size*data->number_of_series)+1,sizeof(double));

    data->CPU_mean_value = 0.0;
    for (unsigned int j=0;j<data->part_number;j++)
        for (unsigned int i=0; i<data->data_size; i++) {
            data->CPU_data[i+j*data->data_size] = 0.5*(data1->CPU_data[i+j*data->data_size] + data2->CPU_data[i+j*data->data_size]);
            if (i>0) data->CPU_mean_value += data->CPU_data[i+j*data->data_size];
        }
    if (last_index>0) data->CPU_mean_value /= (double) (last_index * data->part_number);
    data->CPU_last_value = data->CPU_data[last_index + (data->part_number-1)*data->data_size];
}
void        analysis::lattice_data_analysis_joint3(data_analysis* data,data_analysis* data1,
                                                   data_analysis* data2,data_analysis* data3){
    data->data_size = _MIN(_MIN(data1->data_size,data2->data_size),data3->data_size);
    data->number_of_series = _MAX(_MAX(data1->number_of_series,data2->number_of_series),data3->number_of_series);
    data->part_number = _MIN(_MIN(data1->part_number,data2->part_number),data3->part_number);
    data->storage_type = GPU_CL::GPU::GPU_storage_joint;
    int last_index = (data->data_size - 1);
    if (data->data==NULL) data->data = (double*) calloc((data->data_size*data->number_of_series)+1,sizeof(double));

    data->mean_value = 0.0;
    data->variance   = 0.0;
    for (unsigned int j=0;j<data->part_number;j++)
        for (unsigned int i=0; i<data->data_size; i++) {
            data->data[i+j*data->data_size] = sqrt (data1->data[i+j*data->data_size] * data1->data[i+j*data->data_size] + 
                                                    data2->data[i+j*data->data_size] * data2->data[i+j*data->data_size] +
                                                    data3->data[i+j*data->data_size] * data3->data[i+j*data->data_size]);
            if (i>0) data->mean_value += data->data[i+j*data->data_size];
        }
    if (last_index>0) data->mean_value /= (double) (last_index * data->part_number);
    for (unsigned int j=0;j<data->part_number;j++)
        for (unsigned int i=1; i<data->data_size; i++)
            data->variance += pow(data->data[i+j*data->data_size] - data->mean_value,2);
    if (last_index>0) data->variance   /= (double) (last_index * data->part_number);
    data->GPU_last_value = data->data[last_index + (data->part_number-1)*data->data_size];
    lattice_data_verify(data);
}

}