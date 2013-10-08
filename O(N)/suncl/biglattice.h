/******************************************************************************
 * @file     biglattice.h
 * @author   Vadim Demchik <vadimdi@yahoo.com>
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.0
 *
 * @brief    [QCDGPU]
 *           Dividing lattice into several parts (header)
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

#ifndef biglattice_h
#define biglattice_h

#include "../clinterface/clinterface.h"
#include "../suncl/suncl.h"
#include "../random/random.h"

namespace BIG_LAT{
class BL {

    public:
                    class lattice_data_buffers{
                        public:
                               cl_float4*  plattice_table_float;
                              cl_double4*  plattice_table_double;
                            unsigned int   size_lattice_table;
                            unsigned int   host;
                            unsigned int   platform;
                            unsigned int   device;
                                    bool   device_select;

                                    bool   one_device_one_part; // if TRUE then plattice_table_float=plattice_table_double=NULL
                                                                // -> all are stored data in models[i]
                            unsigned int   models_index;        // index of corresponding models[i] object
                            unsigned int   part_number;         // number of lattice part

                            lattice_data_buffers(void);
                           ~lattice_data_buffers(void);
                    };

                    class lattice_devices{
                        public:
                            unsigned int   host;
                            unsigned int   platform;
                            unsigned int   device;
                                    bool   device_select;

                            // after benchmarking
                                  double   performance;         // % of total performance per device

                            // query
                            unsigned int   lattice_parts;       // number of lattice parts to be simulated on device
                                     int*  lattice_domain_size; // lattice domain for simulations

                           lattice_devices(void);
                          ~lattice_devices(void);
                    };


                      model_CL::model** models;         // array of lattice parts
                          SUN_CPU::SU** SUNcpu;         // array of lattice parts
                 lattice_data_buffers** lattice_data;   // arrays for lattice parts
                      lattice_devices** compute_devices;// array of compute devices in cluster

                      // for big lattices
                      bool     single_device;

                  model_CL::model::run_parameters* global_run;

                       int     big_lattice_parts;       // Number of part for big lattice
                       int     compute_devices_number;  // Number of compute devices for MC simulations

                       int     big_lattice_red_passes;  // Number of red passes (red-green scheme)
                       int     big_lattice_green_passes;// Number of green passes (red-green scheme)

    // functions
                    BL(void);       // constructor
                   ~BL(void);       // destructor
            void  prepare(void);    // 1) prepare lattice
            void  init(void);       // 2) initialization
            void  simulate(void);   // 3) perform lattice simulations
            void  lattice_copy_GPU_to_host(lattice_data_buffers* lat);
            void  lattice_copy_host_to_GPU(lattice_data_buffers* lat);
            void  lattice_get_low_boundary(int* i);
            void  lattice_get_high_boundary(int* i);

            void  lattice_setup_lattice_pointer_initial(int* i);
            void  lattice_setup_lattice_pointer_last(int* i);

            void  lattice_wait_queue_red(int i);
            void  lattice_wait_queue_green(int i);

            void  lattice_wait_kernels_red(int i);
            void  lattice_wait_kernels_green(int i);

            void  lattice_wait_table_read_red(int i);
            void  lattice_wait_table_read_green(int i);

            void  lattice_wait_table_write_red(int i);
            void  lattice_wait_table_write_green(int i);

            void  lattice_copy_boundaries_red(int i);
            void  lattice_copy_boundaries_green(int i);

            void  lattice_measure_red(int i);
            void  lattice_measure_green(int i);

            void  lattice_save_state_red(int i);
            void  lattice_save_state_green(int i);

            void  lattice_update_red(int i);
            void  lattice_update_green(int i);

            void  lattice_orthogonalization_red(int i);
            void  lattice_orthogonalization_green(int i);



    private:
           char*  str_parameter_init(char* str_source);
};
}

#endif
