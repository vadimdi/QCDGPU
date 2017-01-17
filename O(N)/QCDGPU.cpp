/******************************************************************************
 * @file     QCDGPU.cpp
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.0
 *
 * @brief    [QCDGPU]
 *           The main program file
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

#include "QCDGPU.h"

char version[] = "1.4";     // version of programm
char fil[]     = "d:\\init.dat";

#if (MODEL_ON==1)
char fprefix[] = "o1-";    // prefix for data-files
#else //if (MODEL_ON>=2)
char fprefix[] = "su3-";    // prefix for data-files
#endif

char fload[]   = "su3-13-Oct-08-00-00-00.qcg"; // file to be load

#ifdef _WIN32
            char path[]   = "d:\\";                         // path for output
#else
            char path[]   = "/home/science/scalar/";        // path for file with last configuration output
#endif

char*           str_parameter_init(char* str_source){
       char* str_destination = (char*) calloc((strlen_s(str_source) + 1),sizeof(char));
       if (str_destination) {
           strcpy_s(str_destination, (strlen_s(str_source) + 1), str_source);
       }
       return str_destination;
}

#ifdef _WIN32
int setenv(const char* name, const char* value, int overwrite)
{
    int errcode = 0;
    if (!overwrite) {
        size_t envsize = 0;
        errcode = getenv_s(&envsize, NULL, 0, name);
        if (errcode || envsize) return errcode;
    }
    return _putenv_s(name, value);
}
#endif



int main(int argc, char ** argv)
{
    using  GPU_CL::GPU;
    using  model_CL::model;
    using  SUN_CPU::SU;
    using  BIG_LAT::BL;

    int result = 0;

    BL* lattice = new(BL);

    setenv("CUDA_CACHE_DISABLE", "1", 1);

    // setup debug_level
    lattice->global_run->GPU_debug->local_run             = true;
    lattice->global_run->GPU_debug->rebuild_binary        = false;
    if (argc>1) lattice->global_run->GPU_debug->local_run = false;
    if (lattice->global_run->GPU_debug->local_run) {
#ifdef _WIN32
        lattice->global_run->GPU_debug->wait_for_keypress = true;
#endif
        lattice->global_run->GPU_debug->profiling         = false;
        lattice->global_run->GPU_debug->brief_report      = true;
        lattice->global_run->GPU_debug->show_stage        = false;
    }

    lattice->global_run->version    = str_parameter_init(version); // setup version for output files
    lattice->global_run->finishpath = str_parameter_init(path);    // setup path to start.txt and finish.txt files
    lattice->global_run->path       = str_parameter_init(path);    // setup path to output files
    lattice->global_run->fprefix    = str_parameter_init(fprefix); // setup prefix for output files
    lattice->global_run->fstate     = str_parameter_init(fload);   // setup filename for state file

    lattice->global_run->PRNG_generator   = PRNG_CL::PRNG::PRNG_generator_RANLUX3; // CONSTANT;
    lattice->global_run->PRNG_randseries  = 7;                     // constant random series

    lattice->global_run->lattice_nd   = 4;
    lattice->global_run->lattice_type = MODEL_ON;   // sites or links
    lattice->global_run->lattice_full_size[0] = 18;  // X
    lattice->global_run->lattice_full_size[1] = 18;  // Y
    lattice->global_run->lattice_full_size[2] = 18;  // Z
    lattice->global_run->lattice_full_size[3] =  4;  // T
#if (MODEL_ON == 1)
    // O(1) params
    lattice->global_run->lattice_group = 1;     // O(1) group
    lattice->global_run->ON_SK_group   = 4;     // Skalozub's O(4) model
    lattice->global_run->ON_z          = 0.1561573990981998;
    lattice->global_run->ON_lambda     = 0.25;
    lattice->global_run->ON_zeta       = 0.25;
    lattice->global_run->ON_eta        = 0.1;
    lattice->global_run->ON_b          = 9.661840260273483;

    lattice->global_run->correlator_X   = 4; // offset for correlator in X direction
    lattice->global_run->correlator_Y   = 4; // offset for correlator in Y direction
    lattice->global_run->correlator_Z   = 4; // offset for correlator in Z direction
    lattice->global_run->correlator_T   = 4; // offset for correlator in T direction

    lattice->global_run->get_correlators       = true;
#else
    // SU(3) params
    lattice->global_run->lattice_group = 3;     // SU(3) group
    lattice->global_run->BETA          = 3.0;
    lattice->global_run->PHI           = 0.0;   // lambda_3
    lattice->global_run->OMEGA         = 0.0;   // lambda_8
    lattice->global_run->PL_level      = 0;     // level of Polyakov loop calculation (0 = do not calculate PL, 1 = calculate PL only, 2 = calculate PL, PL^2, PL^4)
    lattice->global_run->wilson_R      = 1;
    lattice->global_run->wilson_T      = 2;
#endif

    lattice->global_run->precision = model::model_precision_double;
    lattice->global_run->ints      = model::model_start_cold;  // setup start type
    lattice->global_run->INIT      = 1;                        // start simulations (0 - continue, 1 - start)
    lattice->global_run->NAV       = 50;                       // number of thermalization cycles
    lattice->global_run->ITER      = 200;                      // number of working iterations (note: first measurement is performing on the initial configuration!!!)
    lattice->global_run->NITER     = 10;                       // Number of bulk iterations between measurements
    lattice->global_run->NHIT      = 10;

lattice->global_run->write_lattice_state_every_secs = 2.0;     // write lattice configuration every ... seconds (default: 15 minutes)
lattice->global_run->turnoff_state_save             = true;    // do not write lattice states (binary .qcg files)
lattice->global_run->turnoff_config_save            = false;   // do not write configurations (text .cnf files)
lattice->global_run->turnoff_updates                = false;   // turn off lattice updates
lattice->global_run->turnoff_prns                   = false;   // turn off prn production
lattice->global_run->get_plaquettes_avr             = true;    // calculate mean plaquette values
lattice->global_run->get_acceptance_rate            = true;
lattice->global_run->turnoff_boundary_extraction    = true;    // do not extract boundaries
#if (MODEL_ON == 1)
#else
lattice->global_run->PL_level            = 2;                  // level of Polyakov loop calculation (0 = do not calculate PL, 1 = calculate PL only, 2 = calculate PL, PL^2, PL^4)
lattice->global_run->turnoff_gramschmidt = false;                // turn off orthogonalization
lattice->global_run->get_wilson_loop     = false;                // calculate Wilson loop values
lattice->global_run->get_Fmunu           = false;                // calculate Fmunu tensor for H field
lattice->global_run->get_F0mu            = false;                // calculate Fmunu tensor for E field
lattice->global_run->Fmunu_index1        = 1;
lattice->global_run->Fmunu_index2        = 6;
#endif

lattice->global_run->check_prngs         = true;   // check PRNG production

    if (argc>1) {

        for (int i=1;i<argc;i++) {
            unsigned int paramstart = (unsigned int) strcspn(argv[i],"-");
            // get command line parameters!!!
            if (paramstart>0)
                model_CL::model::lattice_get_init_file(argv[i],lattice->global_run);
            else {
                // get parameter
                unsigned int param_name_end = (unsigned int) strcspn(argv[i],"=");
                unsigned int param_len = (unsigned int) strlen_s(argv[i]);
                if (param_name_end < strlen_s(argv[i])){
                    char* param_name = (char*) calloc(param_name_end,sizeof(char));
                    char* param_text = (char*) calloc(param_len-param_name_end+1,sizeof(char));
                    int iVarVal;
                    double fVarVal;
                    for(unsigned int j=1;j<param_name_end;j++) param_name[j-1]=argv[i][j];
                    for(unsigned int j=param_name_end+1;j<param_len;j++) param_text[j-param_name_end-1]=argv[i][j];
                    sscanf_s(param_text,"%d", &iVarVal);
                    sscanf_s(param_text,"%lf", &fVarVal);
                    model_CL::model::parameters_setup(param_name,&iVarVal,&fVarVal,param_text,lattice->global_run);
                    FREE(param_text);
                    FREE(param_name);
                }
            }

        }
    }

    lattice->global_run->device_select    = false;
    lattice->big_lattice_parts      = 1; // 0=autoselection, other=number of sublattices
    lattice->compute_devices_number = 1; // 0=autoselection, other=number of compute devices
    lattice->prepare();      // prepare lattice for simulations


    // available device 0
    lattice->compute_devices[0]->device_select = true;
    lattice->compute_devices[0]->host          = 0;
    lattice->compute_devices[0]->device        = 0;
    lattice->compute_devices[0]->platform      = 0;
    lattice->compute_devices[0]->performance   = 1.0;
    lattice->compute_devices[0]->lattice_parts = 1;             // if compute_devices_number=1, then is copied from big_lattice_parts
    lattice->compute_devices[0]->lattice_domain_size[0] = 18;   // X (other directions are copied from global_run_lattice_full_site)

/*    // available device 1
    lattice->compute_devices[1]->device_select = true;
    lattice->compute_devices[1]->host          = 0;
    lattice->compute_devices[1]->device        = 1;
    lattice->compute_devices[1]->platform      = 0;
    lattice->compute_devices[1]->performance   = 1.0;
    lattice->compute_devices[1]->lattice_parts = 2;   // if compute_devices_number=1, then is copied from big_lattice_parts
    lattice->compute_devices[1]->lattice_domain_size[0] =  8;    // X (other directions are copied from global_run_lattice_full_site)
*/


    lattice->init();      // lattice initialization

    lattice->simulate();  // MC simulations

    // new code____________________
    if (lattice) delete lattice;

    printf("\n[...Program is completed successfully...]\n\n");
    return result;
}
