/******************************************************************************
 * @file     data_analysis.h
 * @author   Vadim Demchik <vadimdi@yahoo.com>
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.0
 *
 * @brief    [QCDGPU]
 *           Data analysis module (header)
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

#ifndef data_analysis_h
#define data_analysis_h

#include "../clinterface/clinterface.h"

namespace analysis_CL{
class analysis {
        public:
            typedef struct data_analysis {
               unsigned int*   pointer;
               unsigned int    pointer_offset;
               unsigned int    data_size;
GPU_CL::GPU::GPU_storage_type  storage_type;
                       bool    precision_single;
                       bool    skip_checking;
                     double    denominator;
               unsigned int    part_number;
               unsigned int    number_of_series;
                 const char*   data_name;
                     // calculating parameters
                     double*   data;
                     double*   CPU_data;
                     double    mean_value;
                     double    CPU_mean_value;
                     double    CPU_last_variance;
                     double    variance;
                     double    CPU_variance;
                     double    GPU_last_value;
                     double    CPU_last_value;
                      data_analysis(void);
                     ~data_analysis(void);
            } data_analysis;

     static bool    results_verification;

            void    lattice_data_analysis(data_analysis* data);
            void    lattice_data_analysis_sub(data_analysis* data,data_analysis* data2);
            void    lattice_data_analysis_CPU(data_analysis* data);
            void    lattice_data_analysis_joint(data_analysis* data,data_analysis* data1,data_analysis* data2);
            void    lattice_data_analysis_joint_CPU(data_analysis* data,data_analysis* data1,data_analysis* data2);
            void    lattice_data_analysis_joint3(data_analysis* data,data_analysis* data1,data_analysis* data2,data_analysis* data3);
            void    lattice_data_analysis_joint_variance(data_analysis* data,data_analysis* data1,data_analysis* data2);
            bool    lattice_data_verify(data_analysis* data);

        private:
             bool   CPU_GPU_verification_single(double a, double b, const char* err_str);   // CPU-GPU results verification (a?=?b with relative accuracy VDELTA_SINGLE)
             bool   CPU_GPU_verification_double(double a, double b, const char* err_str);   // CPU-GPU results verification (a?=?b with relative accuracy VDELTA_DOUBLE)
    double inline   round(double d);

};
};

#endif
