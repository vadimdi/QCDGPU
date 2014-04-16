/******************************************************************************
 * @file     suncl.h
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.4
 *
 * @brief    [QCDGPU]
 *           Lattice simulation procedures (header)
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
#ifndef suncl_h
#define suncl_h

#include "../clinterface/clinterface.h"
#include "../data_analysis/data_analysis.h"
#include "../random/random.h"
#include "../kernel/complex.h"

namespace model_CL{
class model {

        public:
            typedef enum enum_model_starts{
                model_start_hot,                   // hot start
                model_start_cold,                  // cold start
                model_start_gid                    // fill the lattice with GIDs (for testing purposes)
            } model_starts;

            typedef enum enum_model_precision{
                model_precision_single,            // float (32 bit)
                model_precision_double,            // double (64 bit)
                model_precision_mixed              // mixed precision (32 bit + 32 bit)
            } model_precision;

                       //
                      char*    version;            // version of MC programm
                      char*    path;               // path for output files
                      char*    finishpath;         // path for files start.txt and finish.txt
                      char*    fprefix;
                      char*    fstate;             // file to be loaded

                       // timers section
                    time_t     ltimestart;
                    time_t     ltimeend;
                      char*    timestart;
                      char*    timeend;

                       int     lattice_nd;            // Lattice dimension
                       int*    lattice_full_size;     // Lattice size
                       int*    lattice_domain_size;   // Lattice size

                       int     lattice_group;         // Lattice group
                       int     lattice_full_n1n2;     // n1 * n2
                       int     lattice_full_n2n3;     // n2 * n3
                       int     lattice_full_n1n2n3;   // n1 * n2 * n3
                       int     lattice_full_n2n3n4;   // n2 * n3 * n4
              unsigned int     lattice_full_site;     // Total number of lattice sites
                       int     lattice_domain_n1;     // (n1+B.C.)
                       int     lattice_domain_n1n2;   // (n1+B.C.) * n2
                       int     lattice_domain_n2n3;   // n2 * n3
                       int     lattice_domain_n1n2n3; // (n1+B.C.) * n2 * n3
                       int     lattice_domain_n2n3n4; // n2 * n3 * n4
              unsigned int     lattice_domain_site;   // Total number of lattice sites in simulation domain
                       int     lattice_domain_exact_n1n2n3; // exact n1 * n2 * n3
              unsigned int     lattice_domain_exact_site;   // Total number of lattice sites in simulation domain (exact, without boundary links)


                       // calculational flags
                      bool     big_lattice;        // simulate big lattice (lattice is divided into several parts)

                      bool     get_actions_avr;    // calculate mean actions
                      bool     get_plaquettes_avr; // calculate mean plaquettes
                      bool     get_wilson_loop;    // calculate wilson loop
                      bool     get_Fmunu;          // calculate Fmunu tensor for H field
                      bool     get_F0mu;           // calculate Fmunu tensor for E field

                      bool     get_Fmunu1;         // get Fmunu for lambda1 instead of lambda3
                      bool     get_Fmunu2;         // get Fmunu for lambda2 instead of lambda3
                      bool     get_Fmunu4;         // get Fmunu for lambda4 instead of lambda3
                      bool     get_Fmunu5;         // get Fmunu for lambda5 instead of lambda8
                      bool     get_Fmunu6;         // get Fmunu for lambda6 instead of lambda8
                      bool     get_Fmunu7;         // get Fmunu for lambda7 instead of lambda8

                       int     Fmunu_index1;        // index for first  Fmunu
                       int     Fmunu_index2;        // index for second Fmunu

                      bool     check_prngs;         // check PRNG production

                      bool     turnoff_config_save;// do not write configurations
                      bool     turnoff_updates;    // do not update links
                      bool     turnoff_prns;       // do not produce prns
                      bool     turnoff_gramschmidt;// do not orthogonalize lattice variables

              model_starts     ints;               // type of start
                       int     INIT;               // init : 1=start, 0=continue
                       int     ITER;               // iter : number of samples
                       int     NITER;              // niter: number of bulk simulations between measurements
              unsigned int*    lattice_data;       // Lattice data
                       int     NHIT;               // parameter for multihit
                    double     BETA;               // beta
                       int     NAV;                // number of thermalization cycles
                       int     wilson_R;           // R size for Wilson loop
                       int     wilson_T;           // T size for Wilson loop
           model_precision     precision;          // precision to be used
              unsigned int     PL_level;           // level for calculation Polyakov loops (0 - do not calculate PL, 1 - PL only, 2 - PL, PL^2, PL^4)
                    double     PHI;                // phi angle (lambda_3)
                    double     OMEGA;              // omega angle (lambda_8)

                    double     write_lattice_state_every_secs;  // write lattice state every ... seconds

              // runtime counters
              unsigned int     PRNG_counter;       // counter runs of subroutine PRNG_produce (for load_state purposes)
              unsigned int     NAV_counter;        // number of performed thermalization cycles
              unsigned int     ITER_counter;       // number of performed working cycles
              unsigned int     LOAD_state;         // current load state

              // additional recalculating data
              unsigned int     lattice_table_size;      // Length of lattice table
              unsigned int     lattice_boundary_size;   // Length of lattice boundary
              unsigned int     rowsize;                 // Length of row in lattice
              unsigned int     rowsize4;                // Length of row in lattice (*4)
              unsigned int     halfrowsize;             // Length of half row in lattice
              unsigned int     halfrowsize4;            // Length of half row in lattice (*4)
              unsigned int     lattice_table_row_size;
              unsigned int     lattice_table_exact_row_size;

              unsigned int     size_lattice_table;          // size of buffer lattice_table
              unsigned int     size_lattice_measurement;    // size of buffer lattice_measurement
              unsigned int     size_lattice_energies;       // size of buffer lattice_energies
              unsigned int     size_lattice_wilson_loop;    // size of buffer lattice_wilson_loop
              unsigned int     size_lattice_energies_plq;   // size of buffer lattice_energies_plq
              unsigned int     size_lattice_polyakov_loop;  // size of buffer lattice_polyakov_loop
              unsigned int     size_lattice_boundary;       // size of buffer lattice_boundary
              unsigned int     size_lattice_parameters;     // size of buffer lattice_parameters

              // additional recalculating data
              unsigned int     lattice_boundary_exact_size;
              unsigned int     lattice_energies_size;
              unsigned int     lattice_energies_offset;
              unsigned int     lattice_energies_size_F;
              unsigned int     lattice_polyakov_loop_size;
              unsigned int     lattice_polyakov_loop_offset;
              unsigned int*    lattice_pointer_initial;
              unsigned int*    lattice_pointer_last;
              unsigned int*    lattice_pointer_save;
              unsigned int*    lattice_pointer_measurements;
              unsigned int     lattice_measurement_size;
              unsigned int     lattice_measurement_offset;
              unsigned int     lattice_measurement_size_F;
                 cl_float4*    prng_pointer;
              unsigned int     prngstep;            // step between two threads in prng table

              size_t           local_size_intel;

              unsigned int     desired_platform;
              unsigned int     desired_device;
                      bool     device_select;
               GPU_CL::GPU*    GPU0;                // pointer to GPU instance


#define MODEL_parameter_size    6   // number of parameters for parameters buffer
#define MODEL_energies_size     7   // number of measurements in energy buffer

#define DM_Wilson_loop       0 // index for data measurement for Wilson_loop
#define DM_S_total           1 // index for data measurement for S_total
#define DM_Plq_total         2 // index for data measurement for Plq_total
#define DM_S_spat            3 // index for data measurement for S_spat
#define DM_S_temp            4 // index for data measurement for S_temp
#define DM_Plq_spat          5 // index for data measurement for Plq_spat
#define DM_Plq_temp          6 // index for data measurement for Plq_temp
#define DM_Polyakov_loop     7 // index for data measurement for Polyakov_loop
#define DM_Polyakov_loop_im  8 // index for data measurement for Polyakov_loop_im
#define DM_Polyakov_loop_P2  9 // index for data measurement for Polyakov_loop_P2
#define DM_Polyakov_loop_P4 10 // index for data measurement for Polyakov_loop_P4
#define DM_Fmunu_3          11 // index for data measurement for Fmunu_3's
#define DM_Fmunu_xy_3_re    11 // index for data measurement for Fmunu_xy_3_re
#define DM_Fmunu_xy_3_im    12 // index for data measurement for Fmunu_xy_3_im
#define DM_Fmunu_xz_3_re    13 // index for data measurement for Fmunu_xz_3_re
#define DM_Fmunu_xz_3_im    14 // index for data measurement for Fmunu_xz_3_im
#define DM_Fmunu_yz_3_re    15 // index for data measurement for Fmunu_yz_3_re
#define DM_Fmunu_yz_3_im    16 // index for data measurement for Fmunu_yz_3_im
#define DM_Fmunu_8          17 // index for data measurement for Fmunu_8's
#define DM_Fmunu_xy_8_re    17 // index for data measurement for Fmunu_xy_8_re
#define DM_Fmunu_xy_8_im    18 // index for data measurement for Fmunu_xy_8_im
#define DM_Fmunu_xz_8_re    19 // index for data measurement for Fmunu_xz_8_re
#define DM_Fmunu_xz_8_im    20 // index for data measurement for Fmunu_xz_8_im
#define DM_Fmunu_yz_8_re    21 // index for data measurement for Fmunu_yz_8_re
#define DM_Fmunu_yz_8_im    22 // index for data measurement for Fmunu_yz_8_im
#define DM_Fmunu_abs_3_re   23 // index for data measurement for Fmunu_abs_3_re
#define DM_Fmunu_abs_3_im   24 // index for data measurement for Fmunu_abs_3_im
#define DM_Fmunu_abs_8_re   25 // index for data measurement for Fmunu_abs_8_re
#define DM_Fmunu_abs_8_im   26 // index for data measurement for Fmunu_abs_8_im
#define DM_max              26 // max index for data measurements

             // measurements
analysis_CL::analysis::data_analysis*   Analysis;            // array for all measurements

                      char*    header;
                       int     header_index;

             // PRNG section
            PRNG_CL::PRNG*     PRNG0;                          // pointer to PRNG instance

             // data analysis section
            analysis_CL::analysis*     D_A;                    // pointer to data_analysis instance


            model(void);
           ~model(void);

             int*   lattice_group_elements;

            // identificators for kernels
             int    sun_init_id;
             int    sun_init_X_id;
             int    sun_init_Y_id;
             int    sun_init_Z_id;
             int    sun_init_T_id;
             int    sun_GramSchmidt_id;
             int    sun_update_odd_X_id;
             int    sun_update_odd_Y_id;
             int    sun_update_odd_Z_id;
             int    sun_update_odd_T_id;
             int    sun_update_even_X_id;
             int    sun_update_even_Y_id;
             int    sun_update_even_Z_id;
             int    sun_update_even_T_id;
             int    sun_clear_measurement_id;
             int    sun_get_boundary_low_id;
             int    sun_put_boundary_low_id;
             int    sun_measurement_id;
             int    sun_measurement_reduce_id;
             int    sun_measurement_plq_id;
             int    sun_measurement_plq_reduce_id;
             int    sun_measurement_wilson_id;
             int    sun_wilson_loop_reduce_id;
             int    sun_polyakov_id;
             int    sun_polyakov_reduce_id;
             int    sun_update_indices_id;

             int    argument_wilson_index;
             int    argument_plq_index;
             int    argument_polyakov_index;
             int    argument_measurement_index;
             int    wilson_index;
             int    plq_index;
             int    polyakov_index;
             int    measurement_index;

            // identificators for buffers
    unsigned int    lattice_table;
    unsigned int    lattice_boundary;
    unsigned int    lattice_parameters;
    unsigned int    lattice_measurement;
    unsigned int    lattice_lds;
    unsigned int    lattice_energies;
    unsigned int    lattice_energies_plq;
    unsigned int    lattice_wilson_loop;
    unsigned int    lattice_polyakov_loop;

            // pointers for buffers
    cl_float4*      plattice_table_float;
    cl_double4*     plattice_table_double;
    cl_float4*      plattice_boundary_float;
    cl_double4*     plattice_boundary_double;
    cl_double2*     plattice_measurement;
    cl_double2*     plattice_energies;
    cl_double2*     plattice_energies_plq;
    cl_double*      plattice_wilson_loop;
    cl_double2*     plattice_polyakov_loop;


            // functions
            void    lattice_init(void);
            void    lattice_simulate(void);
            void    lattice_analysis(void);
            void    lattice_print_measurements(void);
            void    lattice_write_results(void);
            void    lattice_get_init_file(char* file);
            void    lattice_save_state(void);
            void    lattice_load_state(void);

            char*   lattice_make_header(void);
            char*   lattice_make_header2(void);
    unsigned int*   lattice_make_bin_header(void);
            bool    lattice_load_bin_header(unsigned int* head);
             int    model_make_header(char* header,int header_size);
            void    lattice_make_programs(void);
            void    lattice_create_buffers(void);

    unsigned int        convert_str_uint(const char* str,unsigned int offset);
    unsigned int        convert_start_to_uint(model::model_starts start);
 model::model_starts    convert_uint_to_start(unsigned int start);
    unsigned int        convert_precision_to_uint(model::model_precision precision);
 model::model_precision convert_uint_to_precision(unsigned int precision);

// PRIVATE STUFF ____________________________________________________________________________________________________
        private:
           static char model_source[FILENAME_MAX];          // path to model-description OpenCL-kernels
           static char complex_source[FILENAME_MAX];        // path to complex OpenCL-kernels
           static char misc_source[FILENAME_MAX];           // path to misc OpenCL-kernels
           static char sun_source[FILENAME_MAX];            // path to SU(N) OpenCL-kernels
           static char sun_measure_source[FILENAME_MAX];    // path to SU(N) measurements OpenCL-kernels
           static char su3_source[FILENAME_MAX];            // path to SU(3) common OpenCL-kernels
           static char su3_update_source[FILENAME_MAX];     // path to SU(3) OpenCL-kernels for update procedures
           static char su3_measure_source[FILENAME_MAX];    // path to SU(3) OpenCL-kernels for measurement procedures

           static char path_suncl[FILENAME_MAX];            // Relative path to SUNCL kernels
           static char path_kernel[FILENAME_MAX];           // Relative path to KERNELS kernels

           int lattice_full_link;
           int lattice_domain_link;
           int lattice_domain_exact_link;
           unsigned int lattice_table_row_size_half;
           unsigned int lattice_table_exact_row_size_half;
           unsigned int lattice_table_group;
           unsigned int lattice_table_exact_group;
           unsigned int lattice_polyakov_size;
           unsigned int lattice_parameters_size;

           void    model_create(void);          // subroutine for tunning particular SU(N) model
           void    model_lattice_init(void);    // subroutine for tunning lattice initialization for SU(N) model

};
};

#endif
