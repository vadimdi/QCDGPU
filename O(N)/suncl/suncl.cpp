/******************************************************************************
 * @file     suncl.cpp
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.0
 *
 * @brief    [QCDGPU]
 *           Lattice simulation procedures
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

#include "suncl.h"
#include "suncpu.h"

#include "../QCDGPU.h"

namespace model_CL{
using model_CL::model;
using SUN_CPU::SU;
using analysis_CL::analysis;

const char* FXYZ[] = {"x","y","z","t"};   // markers for axis
const char* REIM[] = {"re","im"};         // markers for re / im

#define TIMER_FOR_ELAPSED      0  // index of timer for elapsed time calculation
#define TIMER_FOR_SIMULATIONS  1  // index of timer for  simulation time calculation
#define TIMER_FOR_SAVE         2  // index of timer for saving lattice states during simulation
#define ND_MAX                32  // maximum dimensions
#define BIN_HEADER_SIZE       64  // length of binary header (in dwords)

#if (MODEL_ON==1)
#define SOURCE_UPDATE       "oncl/oncl.cl"
#define SOURCE_BOUNDARY     "suncl/boundary.cl"
#define SOURCE_MEASUREMENTS "oncl/on_measurements_cl.cl"
#define SOURCE_POLYAKOV     "oncl/polyakov.cl"
#define SOURCE_WILSON_LOOP  "oncl/wilson_loop.cl"
#else
#define SOURCE_UPDATE       "suncl/suncl.cl"
#define SOURCE_BOUNDARY     "suncl/boundary.cl"
#define SOURCE_MEASUREMENTS "suncl/sun_measurements_cl.cl"
#define SOURCE_POLYAKOV     "suncl/polyakov.cl"
#define SOURCE_WILSON_LOOP  "suncl/wilson_loop.cl"
#endif

char model::path_oncl[FNAME_MAX_LENGTH]   = "oncl/";
char model::path_suncl[FNAME_MAX_LENGTH]  = "suncl/";
char model::path_kernel[FNAME_MAX_LENGTH] = "kernel/";


// SU(N) section __________________________________________________________________________________________
#define DATA_MEASUREMENTS   48 // number of elements for measurements

// end of SU(N) section -----------------------------------------------------------------------------------

// common section ____________________________________
            model::model(void) {
        PRNG0 = new(PRNG_CL::PRNG);         // PRNG module
        D_A   = new(analysis_CL::analysis); // Data Analysis module
        GPU0  = new(GPU_CL::GPU);           // GPU module
        run   = new(run_parameters);        // run parameters

        run->write_lattice_state_every_secs = 15.0 * 60.0; // write lattice configuration every 15 minutes

        PRNG_counter         = 0;   // counter runs of subroutine PRNG_produce (for load_state purposes)
        NAV_counter          = 0;   // number of performed thermalization cycles
        ITER_counter         = 0;   // number of performed working cycles
        LOAD_state           = 0;
        GramSchmidt_iterator = 0;

        // clear pointers
        lattice_pointer_save         = NULL;
        lattice_pointer_last         = NULL;
        lattice_pointer_initial      = NULL;
        lattice_pointer_measurements = NULL;
        prng_pointer                 = NULL;

        model_create(); // tune particular model

        Analysis = new analysis_CL::analysis::data_analysis[DATA_MEASUREMENTS];
}
            model::~model(void) {
        delete[] Analysis;
        FREE(lattice_group_elements);

        if (GPU0->GPU_debug->profiling) GPU0->print_time_detailed();

        // output elapsed time
        printf("Elapsed time: %f seconds\n",GPU0->get_timer_CPU(TIMER_FOR_ELAPSED));

        GPU0->make_finish_file(run->finishpath);
        GPU0->device_finalize(0);

        delete run;
        delete GPU0;
        delete D_A;
        delete PRNG0;
}
            model::run_parameters::run_parameters(void){
        desired_platform    = 0;
        desired_device      = 0;
        device_select       = false;

        run_PRNG  = new(PRNG_CL::PRNG::PRNG_parameters);
        GPU_debug = new(GPU_CL::GPU::GPU_debug_flags);

        lattice_nd    = 4; // default lattice dimension
        lattice_type  = 2; // default lattice type - links
        lattice_group = 3; // default gauge group

        lattice_full_size   = new int[ND_MAX];
        lattice_domain_size = new int[ND_MAX];

        turnoff_state_save  = false; // do not write lattice states (binary .qcg files)
        turnoff_config_save = false; // do not write lattice configurations (text .cnf files)
        turnoff_prns        = false; // turn off prn production
        turnoff_updates     = false; // turn off lattice updates
        turnoff_gramschmidt = true;  // turn off Gram-Schmidt orthogonalization
        turnoff_boundary_extraction = true; // turn off boundary extraction for external purposes

        get_acceptance_rate = false; // calculate mean acceptance rate
        get_correlators     = false; // calculate correlators
        get_actions_avr     = true;  // calculate mean action values
        check_prngs         = false; // check PRNG production

        get_plaquettes_avr  = false; // calculate mean plaquette values
        get_wilson_loop     = false; // calculate Wilson loop values
        get_Fmunu           = false; // calculate Fmunu tensor for H field
        get_F0mu            = false; // calculate Fmunu tensor for E field

        correlator_X        = 1;     // default offset for correlator in X direction
        correlator_Y        = 1;     // default offset for correlator in Y direction
        correlator_Z        = 1;     // default offset for correlator in Z direction
        correlator_T        = 1;     // default offset for correlator in T direction
        
        Fmunu_index1 = 3;   // get Fmunu for lambda3
        Fmunu_index2 = 8;   // get Fmunu for lambda8

        PL_level     = 0;

        BETA     = 1.0;
        PHI      = 0.0;
        OMEGA    = 0.0;
        wilson_R = 1;
        wilson_T = 1;
}
            model::run_parameters::~run_parameters(void){
        if (run_PRNG)  delete run_PRNG;
        if (GPU_debug) delete GPU_debug;
        if (lattice_domain_size) delete[] lattice_domain_size;
        if (lattice_full_size)   delete[] lattice_full_size;

}
void        model::run_parameters::run_parameters_copy(run_parameters* src,run_parameters* dst){
        GPU_CL::GPU::copy_debug_flags(src->GPU_debug,dst->GPU_debug);

        dst->version    = str_parameter_init(src->version);
        dst->finishpath = str_parameter_init(src->finishpath);
        dst->path       = str_parameter_init(src->path);
        dst->fprefix    = str_parameter_init(src->fprefix);
        dst->fstate     = str_parameter_init(src->fstate);

        dst->run_PRNG->PRNG_generator  = src->PRNG_generator;
        dst->run_PRNG->PRNG_randseries = src->PRNG_randseries;

        dst->lattice_type  = src->lattice_type;
        dst->lattice_nd    = src->lattice_nd;
        dst->lattice_group = src->lattice_group;

        dst->precision  = src->precision;
        dst->ints       = src->ints;
        dst->INIT       = src->INIT;
        dst->NAV        = src->NAV;
        dst->ITER       = src->ITER;
        dst->NITER      = src->NITER;
        dst->NHIT       = src->NHIT;
        dst->BETA       = src->BETA;
        dst->PHI        = src->PHI;
        dst->OMEGA      = src->OMEGA;
        dst->PL_level   = src->PL_level;
        dst->wilson_R   = src->wilson_R;
        dst->wilson_T   = src->wilson_T;

        // for O(N) models
        dst->ON_SK_group = src->ON_SK_group;
        dst->ON_z        = src->ON_z;
        dst->ON_lambda   = src->ON_lambda;
        dst->ON_zeta     = src->ON_zeta;
        dst->ON_eta      = src->ON_eta;
        dst->ON_b        = src->ON_b;

        dst->correlator_X = src->correlator_X;
        dst->correlator_Y = src->correlator_Y;
        dst->correlator_Z = src->correlator_Z;
        dst->correlator_T = src->correlator_T;

        dst->Fmunu_index1        = src->Fmunu_index1;
        dst->Fmunu_index2        = src->Fmunu_index2;

        dst->get_acceptance_rate = src->get_acceptance_rate;
        dst->get_actions_avr     = src->get_actions_avr;
        dst->get_correlators     = src->get_correlators;
        dst->get_plaquettes_avr  = src->get_plaquettes_avr;
        dst->get_wilson_loop     = src->get_wilson_loop;
        dst->get_Fmunu           = src->get_Fmunu;
        dst->get_F0mu            = src->get_F0mu;

        dst->turnoff_state_save  = src->turnoff_state_save;
        dst->turnoff_config_save = src->turnoff_config_save;
        dst->turnoff_updates     = src->turnoff_updates;
        dst->turnoff_prns        = src->turnoff_prns;
        dst->turnoff_gramschmidt = src->turnoff_gramschmidt;
        dst->turnoff_boundary_extraction = src->turnoff_boundary_extraction;

        dst->check_prngs         = src->check_prngs;

        for (int j=0;j<src->lattice_nd;j++){
            dst->lattice_full_size[j]   = src->lattice_full_size[j];
            dst->lattice_domain_size[j] = src->lattice_domain_size[j];
        }

}

void        model::parameters_setup(char* parameter,int* ivalue,double* fvalue,char* text_value,run_parameters* run){
            if (!strcmp(parameter,"PLATFORM"))  run->desired_platform   = (*ivalue);
            if (!strcmp(parameter,"DEVICE"))   {run->desired_device     = (*ivalue); run->device_select = true;}

            if (!strcmp(parameter,"GROUP"))  run->lattice_group   = (*ivalue);
            if (!strcmp(parameter,"TYPE"))   run->lattice_type    = (*ivalue);
            if (!strcmp(parameter,"ND"))     run->lattice_nd      = (*ivalue);
            if (!strcmp(parameter,"N1"))     run->lattice_domain_size[0] = (*ivalue);
            if (!strcmp(parameter,"N2"))     run->lattice_domain_size[1] = (*ivalue);
            if (!strcmp(parameter,"N3"))     run->lattice_domain_size[2] = (*ivalue);
            if (!strcmp(parameter,"N4"))     run->lattice_domain_size[3] = (*ivalue);
            if (!strcmp(parameter,"NS"))    {run->lattice_domain_size[0] = (*ivalue);
                                                                         run->lattice_domain_size[1] = (*ivalue);
                                                                         run->lattice_domain_size[2] = (*ivalue); }
            if (!strcmp(parameter,"NT"))     run->lattice_domain_size[run->lattice_nd - 1] = (*ivalue);
            if (!strcmp(parameter,"L1"))     run->lattice_full_size[0] = (*ivalue);
            if (!strcmp(parameter,"L2"))     run->lattice_full_size[1] = (*ivalue);
            if (!strcmp(parameter,"L3"))     run->lattice_full_size[2] = (*ivalue);
            if (!strcmp(parameter,"L4"))     run->lattice_full_size[3] = (*ivalue);
            if (!strcmp(parameter,"LS"))    {run->lattice_full_size[0] = (*ivalue);
                                                                         run->lattice_full_size[1] = (*ivalue);
                                                                         run->lattice_full_size[2] = (*ivalue); }
            if (!strcmp(parameter,"LT"))     run->lattice_full_size[run->lattice_nd - 1] = (*ivalue);

            if (!strcmp(parameter,"SKGROUP"))run->ON_SK_group     = (*ivalue);
            if (!strcmp(parameter,"ITER"))   run->ITER            = (*ivalue);
            if (!strcmp(parameter,"NITER"))  run->NITER           = (*ivalue);
            if (!strcmp(parameter,"NHIT"))   run->NHIT            = (*ivalue);
            if (!strcmp(parameter,"BETA"))   run->BETA            = (*fvalue);
            if (!strcmp(parameter,"PHI"))    run->PHI             = (*fvalue);
            if (!strcmp(parameter,"OMEGA"))  run->OMEGA           = (*fvalue);
            if (!strcmp(parameter,"NAV"))    run->NAV             = (*ivalue);

            if (!strcmp(parameter,"ON_Z"))     run->ON_z          = (*fvalue);
            if (!strcmp(parameter,"ON_ZETA"))  run->ON_zeta       = (*fvalue);
            if (!strcmp(parameter,"ON_ETA"))   run->ON_eta        = (*fvalue);
            if (!strcmp(parameter,"ON_B"))     run->ON_b          = (*fvalue);
            if (!strcmp(parameter,"ON_LAMBDA"))run->ON_lambda     = (*fvalue);

            if (!strcmp(parameter,"RANDSERIES")) run->PRNG_randseries = (*ivalue);
            if (!strcmp(parameter,"PRNG"))       run->PRNG_generator  = PRNG_CL::PRNG::get_PRNG_by_name(text_value);
            PRNG_CL::PRNG::parameters_setup(parameter,(*ivalue),text_value,run->run_PRNG);

            if (!strcmp(parameter,"INTS"))            run->ints       = model_CL::model::convert_uint_to_start((*ivalue));
            if (!strcmp(parameter,"OUTPUTPATH"))      run->path       = model_CL::model::str_parameter_init(text_value);
            if (!strcmp(parameter,"FINISHPATH"))      run->finishpath = model_CL::model::str_parameter_init(text_value);
            if (!strcmp(parameter,"TURNOFFWAITING"))  run->GPU_debug->local_run = true;
            if (!strcmp(parameter,"TURNOFFKEYPRESS")) run->GPU_debug->wait_for_keypress = false;
            if (!strcmp(parameter,"REBUILDBINARY"))   run->GPU_debug->rebuild_binary = true;
            if (!strcmp(parameter,"GETWILSON"))       run->get_wilson_loop = true;
            if (!strcmp(parameter,"GETRETRACE"))      run->get_plaquettes_avr = true;
            if (!strcmp(parameter,"TURNOFFFMUNU"))  {
                run->get_Fmunu = false;
                run->get_F0mu  = false;
            }
            if (!strcmp(parameter,"F0MU"))  {
                run->get_F0mu  = true;
                run->get_Fmunu = false;
            }
            if (!strcmp(parameter,"FMUNU1"))  run->Fmunu_index1 = 1;
            if (!strcmp(parameter,"FMUNU2"))  run->Fmunu_index1 = 2;
            if (!strcmp(parameter,"FMUNU3"))  run->Fmunu_index1 = 3;
            if (!strcmp(parameter,"FMUNU4"))  run->Fmunu_index1 = 4;
            if (!strcmp(parameter,"FMUNU5"))  run->Fmunu_index2 = 5;
            if (!strcmp(parameter,"FMUNU6"))  run->Fmunu_index2 = 6;
            if (!strcmp(parameter,"FMUNU7"))  run->Fmunu_index2 = 7;
            if (!strcmp(parameter,"FMUNU8"))  run->Fmunu_index2 = 8;
            if (!strcmp(parameter,"WILSONR")) run->wilson_R = (*ivalue);
            if (!strcmp(parameter,"WILSONT")) run->wilson_T = (*ivalue);

}
void        model::lattice_get_init_file(char* file,run_parameters* run){
        int parameters_items = 0;
        GPU_CL::GPU::GPU_init_parameters* parameters = GPU_CL::GPU::get_init_file(file);
        if (parameters==NULL) return;

        bool parameters_flag = false;
        while(!parameters_flag){
            parameters_setup(parameters[parameters_items].Variable,&parameters[parameters_items].iVarVal,&parameters[parameters_items].fVarVal,parameters[parameters_items].txtVarVal,run);
            parameters_flag = parameters[parameters_items].final;
            parameters_items++;
        }
}
char*       model::lattice_make_header(void){
    int header_size = 16384;
    header = (char*) calloc(header_size, sizeof(char));
    int j = 0;

    if (run->lattice_type==1){
        j  += sprintf_s(header+j,header_size-j, " GPU O(%u) simulator (Skalozub action for O(%u)) %s\n\n",run->lattice_group,run->ON_SK_group,run->version);
        j  += sprintf_s(header+j,header_size-j, " Monte Carlo simulation of %uD O(%u) LGT\n\n",run->lattice_nd,run->lattice_group);
    }
    if (run->lattice_type==2){
        j  += sprintf_s(header+j,header_size-j, " GPU SU(%u) simulator %s\n\n",run->lattice_group,run->version);
        j  += sprintf_s(header+j,header_size-j, " Monte Carlo simulation of %uD SU(%u) LGT\n\n",run->lattice_nd,run->lattice_group);
    }
    j  += sprintf_s(header+j,header_size-j, " ***************************************************\n");
    j  += sprintf_s(header+j,header_size-j, " Active OpenCL platform   : %s\n",GPU0->platform_get_name(GPU0->GPU_platform));
    j  += sprintf_s(header+j,header_size-j, " Active OpenCL device     : %s\n",GPU0->GPU_info.device_name);
    j  += sprintf_s(header+j,header_size-j, " ***************************************************\n");
    j  += sprintf_s(header+j,header_size-j, " lattice size                : %3u x %3u x %3u x %3u\n",run->lattice_full_size[0],run->lattice_full_size[1],run->lattice_full_size[2],run->lattice_full_size[3]);
    j  += sprintf_s(header+j,header_size-j, " init                        : %i\n",run->INIT);
    if (run->ints == model_start_hot)
        j  += sprintf_s(header+j,header_size-j," ints (0=hot, 1=cold)        : 0\n");
       else
        j  += sprintf_s(header+j,header_size-j," ints (0=hot, 1=cold)        : 1\n");
    j  += PRNG0->print_generator((header+j),(header_size-j));
    j  += sprintf_s(header+j,header_size-j, " nav                         : %i\n",run->NAV);
    j  += sprintf_s(header+j,header_size-j, " niter                       : %i\n",run->NITER);
    j  += sprintf_s(header+j,header_size-j, " iter (# of samples)         : %i\n",run->ITER);
    j  += sprintf_s(header+j,header_size-j, " nhit                        : %i\n",run->NHIT);
    if (run->precision == model::model_precision_double) j  += sprintf_s(header+j,header_size-j, " precision                   : double\n");
    if (run->precision == model::model_precision_single) j  += sprintf_s(header+j,header_size-j, " precision                   : single\n");
    if (run->precision == model::model_precision_mixed)  j  += sprintf_s(header+j,header_size-j, " precision                   : mixed\n");
    j  += sprintf_s(header+j,header_size-j, " ***************************************************\n");
    j  += model_make_header((header+j),(header_size-j));

    header_index = j;
    return header;
}
void        model::lattice_init(void){
    if (!GPU0->GPU_debug->local_run) GPU0->make_start_file(run->finishpath);

    bool supported_devices = false;
    if (!run->device_select) {
        // auto-select first GPU_vendor_nVidia device
        supported_devices = GPU0->device_auto_select(GPU_CL::GPU::GPU_vendor_any,GPU_CL::GPU::GPU_vendor_any);
    } else {
        // manual selection of platform and device
        supported_devices = GPU0->device_select(run->desired_platform,run->desired_device);
    }
    if(!supported_devices){
        printf("There are no any available OpenCL devices\n");
        exit(0);
    }

    // initialize selected device & show hardwares
    GPU0->device_initialize();
    GPU0->print_available_hardware();
    GPU0->print_stage("device initialized");

    if (run->INIT==0) lattice_load_state();          // load state file if needed

    if (!big_lattice){
        if (!run->lattice_domain_size){
            run->lattice_domain_size = new int[run->lattice_nd];
            for (int i=0;i<run->lattice_nd;i++) run->lattice_domain_size[i] = run->lattice_full_size[i];
        }
    }
    GPU_CL::GPU::copy_debug_flags(run->GPU_debug,GPU0->GPU_debug);

    model_lattice_init();   // model initialization
}
unsigned int model::convert_str_uint(const char* str,unsigned int offset){
    unsigned int result = 0;
    unsigned int pointer = offset;
    for (unsigned int i = 0; i < 4; ++i){
        if (pointer < strlen_s(str)) result += (str[pointer++]<<(i*8));
    }
    return result;
}
unsigned int model::convert_start_to_uint(model::model_starts start){
    if (start == model::model_start_hot)  return 0;
    if (start == model::model_start_cold) return 1;
    if (start == model::model_start_gid)  return 2;
    return 0; // return hot start otherwise
}
unsigned int model::convert_precision_to_uint(model::model_precision precision){
    if (precision == model::model_precision_single) return 1;
    if (precision == model::model_precision_double) return 2;
    if (precision == model::model_precision_mixed)  return 3;
    return 1; // return single precision otherwise
}
model::model_starts    model::convert_uint_to_start(unsigned int start){
    if (start == 0) return model::model_start_hot;
    if (start == 1) return model::model_start_cold;
    if (start == 2) return model::model_start_gid;
    return model::model_start_hot; // return hot start otherwise
}
model::model_precision model::convert_uint_to_precision(unsigned int precision){
    if (precision == 1) return model::model_precision_single;
    if (precision == 2) return model::model_precision_double;
    if (precision == 3) return model::model_precision_mixed;
    return model::model_precision_single; // return single precision otherwise
}
char*        model::str_parameter_init(char* str_source){
       char* str_destination = (char*) calloc((strlen_s(str_source) + 1),sizeof(char));
       if (str_destination) {
           strcpy_s(str_destination, (strlen_s(str_source) + 1), str_source);
       }
       return str_destination;
}

// model-dependent section ___________________________
void        model::model_create(void){
        // parameters for SU(N) model
        lattice_group_elements = (int*) calloc(6,sizeof(int));
        lattice_group_elements[0] =  0;
        lattice_group_elements[1] =  1; //  U(1)
        lattice_group_elements[2] =  4; // SU(2)
        lattice_group_elements[3] = 12; // SU(3)

        big_lattice = false;         // simulate big lattice (lattice is divided into several parts)

        sun_init_id                     = 0;
        sun_init_X_id                   = 0;
        sun_init_Y_id                   = 0;
        sun_init_Z_id                   = 0;
        sun_init_T_id                   = 0;
        sun_GramSchmidt_id              = 0;
        sun_update_odd_X_id             = 0;
        sun_update_odd_Y_id             = 0;
        sun_update_odd_Z_id             = 0;
        sun_update_odd_T_id             = 0;
        sun_update_even_X_id            = 0;
        sun_update_even_Y_id            = 0;
        sun_update_even_Z_id            = 0;
        sun_update_even_T_id            = 0;
        sun_clear_measurement_id        = 0;
        sun_get_boundary_low_id         = 0;
        sun_put_boundary_low_id         = 0;
        sun_measurement_id              = 0;
        sun_measurement_reduce_id       = 0;
        sun_measurement_plq_id          = 0;
        sun_measurement_plq_reduce_id   = 0;
        sun_measurement_corr_id         = 0;
        sun_measurement_corr_reduce_id  = 0;
        sun_measurement_wilson_id       = 0;
        sun_wilson_loop_reduce_id       = 0;
        sun_polyakov_id                 = 0;
        sun_polyakov_reduce_id          = 0;
        sun_update_indices_id           = 0;
        sun_reduce_acceptance_rate_id   = 0;
}
int         model::model_make_header(char* header,int header_size){
    int j = 0;
#if (MODEL_ON == 1)
    j  += sprintf_s(header+j,header_size-j, " eta                         : %16.13e\n",run->ON_eta);
    j  += sprintf_s(header+j,header_size-j, " z                           : %16.13e\n",run->ON_z);
    j  += sprintf_s(header+j,header_size-j, " b                           : %16.13e\n",run->ON_b);
    j  += sprintf_s(header+j,header_size-j, " lambda                      : %16.13e\n",run->ON_lambda);
    j  += sprintf_s(header+j,header_size-j, " zeta                        : %16.13e\n",run->ON_zeta);
    j  += sprintf_s(header+j,header_size-j, " ***************************************************\n");
#else
    j  += sprintf_s(header+j,header_size-j, " Wilson loop R               : %i\n",run->wilson_R);
    j  += sprintf_s(header+j,header_size-j, " Wilson loop T               : %i\n",run->wilson_T);
    if (run->get_Fmunu)
        j  += sprintf_s(header+j,header_size-j," FMUNU(%u, %u)\n",run->Fmunu_index1,run->Fmunu_index2);
    if (run->get_F0mu)
        j  += sprintf_s(header+j,header_size-j," F0MU(%u, %u)\n", run->Fmunu_index1,run->Fmunu_index2);

    j  += sprintf_s(header+j,header_size-j, " BETA                        : %16.13e\n",run->BETA);
    j  += sprintf_s(header+j,header_size-j, " PHI   (lambda_3)            : %16.13e\n",run->PHI);
    j  += sprintf_s(header+j,header_size-j, " OMEGA (lambda_8)            : %16.13e\n",run->OMEGA);
    j  += sprintf_s(header+j,header_size-j, " ***************************************************\n");
#endif

    return j;
}
char*       model::lattice_make_header2(void){
    int header_size = 16384;
    int j = header_index;

    int elapsdays,elapshours,elapsminites,elapsseconds;
    double dif;
    dif=difftime(ltimeend,ltimestart);
    elapsdays=(int)(dif/24/3600);
    elapshours=(((int) dif/3600) % 3600);
    elapsminites=(((int) dif/60) % 60);
    elapsseconds=(((int) dif) % 60);

    j  += sprintf_s(header+j,header_size-j, " Start time               : %s",timestart);
    j  += sprintf_s(header+j,header_size-j, " Finish time              : %s",timeend);
    j  += sprintf_s(header+j,header_size-j, " Elapsed time             : %i:%2.2i:%2.2i:%2.2i\n",elapsdays,elapshours,elapsminites,elapsseconds);
    j  += sprintf_s(header+j,header_size-j, " ***************************************************\n");
    j  += sprintf_s(header+j,header_size-j, " Mean %-20s: % 16.13e\n",    Analysis[DM_S_total].data_name,Analysis[DM_S_total].mean_value);
    j  += sprintf_s(header+j,header_size-j, " Variance %-16s: % 16.13e\n",Analysis[DM_S_total].data_name,Analysis[DM_S_total].variance);
    j  += sprintf_s(header+j,header_size-j, " Mean %-20s: % 16.13e\n",    Analysis[DM_Acc_rate_total].data_name,Analysis[DM_Acc_rate_total].mean_value);
#if (MODEL_ON == 1)
    j  += sprintf_s(header+j,header_size-j, " Mean %-20s: % 16.13e\n",    Analysis[DM_Plq_spat].data_name,Analysis[DM_Plq_spat].mean_value);
    j  += sprintf_s(header+j,header_size-j, " Variance %-16s: % 16.13e\n",Analysis[DM_Plq_total].data_name,Analysis[DM_Plq_total].mean_value);
    j  += sprintf_s(header+j,header_size-j, " Mean %-20s: % 16.13e\n",    Analysis[DM_Correlator1].data_name,Analysis[DM_Correlator1].mean_value);
    j  += sprintf_s(header+j,header_size-j, " Mean %-20s: % 16.13e\n",    Analysis[DM_Correlator2].data_name,Analysis[DM_Correlator2].mean_value);
#else
    j  += sprintf_s(header+j,header_size-j, " Mean %-20s: % 16.13e\n",    Analysis[DM_Plq_total].data_name,Analysis[DM_Plq_total].mean_value);
    j  += sprintf_s(header+j,header_size-j, " Variance %-16s: % 16.13e\n",Analysis[DM_Plq_total].data_name,Analysis[DM_Plq_total].variance);
    j  += sprintf_s(header+j,header_size-j, " Mean %-20s: % 16.13e\n",    Analysis[DM_Polyakov_loop].data_name,Analysis[DM_Polyakov_loop].mean_value);
    j  += sprintf_s(header+j,header_size-j, " Variance %-16s: % 16.13e\n",Analysis[DM_Polyakov_loop].data_name,Analysis[DM_Polyakov_loop].variance);
    j  += sprintf_s(header+j,header_size-j, " Mean %-20s: % 16.13e\n",    Analysis[DM_Polyakov_loop_im].data_name,Analysis[DM_Polyakov_loop_im].mean_value);
    j  += sprintf_s(header+j,header_size-j, " Variance %-16s: % 16.13e\n",Analysis[DM_Polyakov_loop_im].data_name,Analysis[DM_Polyakov_loop_im].variance);
    j  += sprintf_s(header+j,header_size-j, " Mean %-20s: % 16.13e\n",    Analysis[DM_Polyakov_loop_P2].data_name,Analysis[DM_Polyakov_loop_P2].mean_value);
    j  += sprintf_s(header+j,header_size-j, " Variance %-16s: % 16.13e\n",Analysis[DM_Polyakov_loop_P2].data_name,Analysis[DM_Polyakov_loop_P2].variance);
    j  += sprintf_s(header+j,header_size-j, " Mean %-20s: % 16.13e\n",    Analysis[DM_Polyakov_loop_P4].data_name,Analysis[DM_Polyakov_loop_P4].mean_value);
    j  += sprintf_s(header+j,header_size-j, " Variance %-16s: % 16.13e\n",Analysis[DM_Polyakov_loop_P4].data_name,Analysis[DM_Polyakov_loop_P4].variance);
    j  += sprintf_s(header+j,header_size-j, " Mean %-20s: % 16.13e\n",    Analysis[DM_Wilson_loop].data_name,Analysis[DM_Wilson_loop].mean_value);
    j  += sprintf_s(header+j,header_size-j, " Variance %-16s: % 16.13e\n",Analysis[DM_Wilson_loop].data_name,Analysis[DM_Wilson_loop].variance);
#endif
    j  += sprintf_s(header+j,header_size-j, " ***************************************************\n");

#if (MODEL_ON == 1)
#else
    for (int i=0;i<((run->lattice_nd-1)*2+2)*2;i++)
        j  += sprintf_s(header+j,header_size-j, " Mean %-20s: % 16.13e\n",Analysis[DM_Fmunu_3+i].data_name,Analysis[DM_Fmunu_3+i].mean_value);
    for (int i=0;i<((run->lattice_nd-1)*2+2)*2;i++)
        j  += sprintf_s(header+j,header_size-j, " Variance %-16s: % 16.13e\n",Analysis[DM_Fmunu_3+i].data_name,Analysis[DM_Fmunu_3+i].variance);
    j  += sprintf_s(header+j,header_size-j, " ***************************************************\n");
    for (int jj=0;jj<((run->lattice_nd-2)*(run->lattice_nd-1)+2)*2;jj++){
        j  += sprintf_s(header+j,header_size-j, " GPU last %-16s: % 16.13e\n",Analysis[jj+DM_Fmunu_3].data_name,Analysis[jj+DM_Fmunu_3].GPU_last_value);
        j  += sprintf_s(header+j,header_size-j, " CPU last %-16s: % 16.13e\n",Analysis[jj+DM_Fmunu_3].data_name,Analysis[jj+DM_Fmunu_3].CPU_last_value);
    }
    j  += sprintf_s(header+j,header_size-j, " ***************************************************\n");
    for (int jj=0;jj<(run->lattice_nd-2)*(run->lattice_nd-1)*2;jj++)
        j  += sprintf_s(header+j,header_size-j, " CPU l.varnc %-13s: % 16.13e\n",Analysis[jj+DM_Fmunu_3].data_name,Analysis[jj+DM_Fmunu_3].CPU_last_variance);
    j  += sprintf_s(header+j,header_size-j, " ***************************************************\n");
#endif
    j  += sprintf_s(header+j,header_size-j, " CPU last %-16s: % 16.13e\n",Analysis[DM_S_total].data_name,Analysis[DM_S_total].CPU_last_value);
    j  += sprintf_s(header+j,header_size-j, " GPU last %-16s: % 16.13e\n",Analysis[DM_S_total].data_name,Analysis[DM_S_total].GPU_last_value);
#if (MODEL_ON == 1)
    j  += sprintf_s(header+j,header_size-j, " CPU last %-16s: % 16.13e\n",Analysis[DM_Plq_spat].data_name,Analysis[DM_Plq_spat].CPU_last_value);
    j  += sprintf_s(header+j,header_size-j, " GPU last %-16s: % 16.13e\n",Analysis[DM_Plq_spat].data_name,Analysis[DM_Plq_spat].GPU_last_value);
    j  += sprintf_s(header+j,header_size-j, " CPU last %-16s: % 16.13e\n",Analysis[DM_Correlator1].data_name,Analysis[DM_Correlator1].CPU_last_value);
    j  += sprintf_s(header+j,header_size-j, " GPU last %-16s: % 16.13e\n",Analysis[DM_Correlator1].data_name,Analysis[DM_Correlator1].GPU_last_value);
    j  += sprintf_s(header+j,header_size-j, " CPU last %-16s: % 16.13e\n",Analysis[DM_Correlator2].data_name,Analysis[DM_Correlator2].CPU_last_value);
    j  += sprintf_s(header+j,header_size-j, " GPU last %-16s: % 16.13e\n",Analysis[DM_Correlator2].data_name,Analysis[DM_Correlator2].GPU_last_value);
#else
    j  += sprintf_s(header+j,header_size-j, " CPU last %-16s: % 16.13e\n",Analysis[DM_Plq_total].data_name,Analysis[DM_Plq_total].CPU_last_value);
    j  += sprintf_s(header+j,header_size-j, " GPU last %-16s: % 16.13e\n",Analysis[DM_Plq_total].data_name,Analysis[DM_Plq_total].GPU_last_value);
    j  += sprintf_s(header+j,header_size-j, " CPU last %-16s: % 16.13e\n",Analysis[DM_Polyakov_loop].data_name,Analysis[DM_Polyakov_loop].CPU_last_value);
    j  += sprintf_s(header+j,header_size-j, " GPU last %-16s: % 16.13e\n",Analysis[DM_Polyakov_loop].data_name,Analysis[DM_Polyakov_loop].GPU_last_value);
    j  += sprintf_s(header+j,header_size-j, " CPU last %-16s: % 16.13e\n",Analysis[DM_Polyakov_loop_im].data_name,Analysis[DM_Polyakov_loop_im].CPU_last_value);
    j  += sprintf_s(header+j,header_size-j, " GPU last %-16s: % 16.13e\n",Analysis[DM_Polyakov_loop_im].data_name,Analysis[DM_Polyakov_loop_im].GPU_last_value);
    j  += sprintf_s(header+j,header_size-j, " CPU last %-16s: % 16.13e\n",Analysis[DM_Polyakov_loop_P2].data_name,Analysis[DM_Polyakov_loop_P2].CPU_last_value);
    j  += sprintf_s(header+j,header_size-j, " GPU last %-16s: % 16.13e\n",Analysis[DM_Polyakov_loop_P2].data_name,Analysis[DM_Polyakov_loop_P2].GPU_last_value);
    j  += sprintf_s(header+j,header_size-j, " CPU last %-16s: % 16.13e\n",Analysis[DM_Polyakov_loop_P4].data_name,Analysis[DM_Polyakov_loop_P4].CPU_last_value);
    j  += sprintf_s(header+j,header_size-j, " GPU last %-16s: % 16.13e\n",Analysis[DM_Polyakov_loop_P4].data_name,Analysis[DM_Polyakov_loop_P4].GPU_last_value);
    j  += sprintf_s(header+j,header_size-j, " CPU last %-16s: % 16.13e\n",Analysis[DM_Wilson_loop].data_name,Analysis[DM_Wilson_loop].CPU_last_value);
    j  += sprintf_s(header+j,header_size-j, " GPU last %-16s: % 16.13e\n",Analysis[DM_Wilson_loop].data_name,Analysis[DM_Wilson_loop].GPU_last_value);
#endif
    if (analysis_CL::analysis::results_verification)
        j  += sprintf_s(header+j,header_size-j, " *** Verification successfully passed! *************\n");
    else
        j  += sprintf_s(header+j,header_size-j, " --- Verification failed! --------------------------\n");
    j  += sprintf_s(header+j,header_size-j, " Data fields:\n");
    j  += sprintf_s(header+j,header_size-j, "    #, ");
    if (run->get_plaquettes_avr){
                           j += sprintf_s(header+j,header_size-j, "%-22s",Analysis[DM_Plq_spat].data_name);
                           j += sprintf_s(header+j,header_size-j, "%-22s",Analysis[DM_Plq_temp].data_name);
                           j += sprintf_s(header+j,header_size-j, "%-22s",Analysis[DM_Plq_total].data_name);
    }
    if (run->get_wilson_loop)
                           j += sprintf_s(header+j,header_size-j, "%-22s",Analysis[DM_Wilson_loop].data_name);
    if (run->get_actions_avr){
                           j += sprintf_s(header+j,header_size-j, "%-22s",Analysis[DM_S_spat].data_name);
                           j += sprintf_s(header+j,header_size-j, "%-22s",Analysis[DM_S_temp].data_name);
                           j += sprintf_s(header+j,header_size-j, "%-22s",Analysis[DM_S_total].data_name);
    }
    if (run->get_correlators){
                           j += sprintf_s(header+j,header_size-j, "%-22s",Analysis[DM_Correlator1].data_name);
                           j += sprintf_s(header+j,header_size-j, "%-22s",Analysis[DM_Correlator2].data_name);
    }
    if (run->PL_level > 0){
                           j += sprintf_s(header+j,header_size-j, "%-22s",Analysis[DM_Polyakov_loop].data_name);
                           j += sprintf_s(header+j,header_size-j, "%-22s",Analysis[DM_Polyakov_loop_im].data_name);
    }
    if (run->PL_level > 1){
                           j += sprintf_s(header+j,header_size-j, "%-22s",Analysis[DM_Polyakov_loop_P2].data_name);
                           j += sprintf_s(header+j,header_size-j, "%-22s",Analysis[DM_Polyakov_loop_P4].data_name);
    }
    if (run->get_acceptance_rate){
                           j += sprintf_s(header+j,header_size-j, "%-22s",Analysis[DM_Acc_rate_even].data_name);
                           j += sprintf_s(header+j,header_size-j, "%-22s",Analysis[DM_Acc_rate_odd].data_name);
                           j += sprintf_s(header+j,header_size-j, "%-22s",Analysis[DM_Acc_rate_total].data_name);
    }
    if ((run->get_Fmunu)||(run->get_F0mu))
        for (int i=0;i<((run->lattice_nd-1)*2+2)*2;i++)
            j  += sprintf_s(header+j,header_size-j, "%-22s",Analysis[DM_Fmunu_3+i].data_name);
    j  += sprintf_s(header+j,header_size-j, "\n");
    j  += sprintf_s(header+j,header_size-j, " ***************************************************\n");

    header_index = j;

    return header;
}
void        model::lattice_print_measurements(void){
    printf("\n #        ");
    if (run->get_plaquettes_avr){
        printf(" %-11s",Analysis[DM_Plq_spat].data_name);
        printf(" %-11s",Analysis[DM_Plq_temp].data_name);
    }
    if (run->get_wilson_loop)
        printf(" %-11s",Analysis[DM_Wilson_loop].data_name);
    if (run->get_actions_avr){
        printf(" %-11s",Analysis[DM_S_spat].data_name);
        printf(" %-11s",Analysis[DM_S_temp].data_name);
    }
    if (run->PL_level > 0){
        printf(" %-11s",Analysis[DM_Polyakov_loop].data_name);
        printf(" %-11s",Analysis[DM_Polyakov_loop_im].data_name);
    }
    if (run->get_acceptance_rate){
        printf(" %-11s",Analysis[DM_Acc_rate_even].data_name);
        printf(" %-11s",Analysis[DM_Acc_rate_odd].data_name);
    }

    printf("\n");
    unsigned int offs = 0;
    for (int i=0; i<run->ITER; i++) {
        if ((i<5)||(i==(run->ITER - 1))){
            printf("[%5u]: ",(i+offs));
            if (run->get_plaquettes_avr){
                printf(" % 10.8f",Analysis[DM_Plq_spat].data[i]);
                printf(" % 10.8f",Analysis[DM_Plq_temp].data[i]);
            }
            if (run->get_wilson_loop)
                printf(" % 10.8f",Analysis[DM_Wilson_loop].data[i]);
            if (run->get_actions_avr){
                printf(" % 10.8f",Analysis[DM_S_spat].data[i]);
                printf(" % 10.8f",Analysis[DM_S_temp].data[i]);
            }
            if (run->PL_level > 0){
                printf(" % 10.8f",Analysis[DM_Polyakov_loop].data[i]);
                printf(" % 10.8f",Analysis[DM_Polyakov_loop_im].data[i]);
            }
            if (run->get_acceptance_rate){
                printf(" % 10.8f",Analysis[DM_Acc_rate_even].data[i]);
                printf(" % 10.8f",Analysis[DM_Acc_rate_odd].data[i]);
            }

            printf("\n");
        }
    }
}
void        model::lattice_analysis(void){
        unsigned int number_spat = (run->lattice_nd - 1) * (run->lattice_nd - 2) / 2;   // number of spatial plaquettes
#if (MODEL_ON == 1)
#else
        unsigned int number_temp = (run->lattice_nd - 1);                               // number of temporal plquettes
#endif

        for (int i=0; i<=DM_max;i++){
            Analysis[i].data_size       = run->ITER;
            Analysis[i].pointer_offset  = 0;
            if (run->precision==model_precision_double) Analysis[i].precision_single = false;
                else                                    Analysis[i].precision_single = true;
            if ((i&1) == 0) Analysis[i].storage_type = GPU_CL::GPU::GPU_storage_double2low;
            else            Analysis[i].storage_type = GPU_CL::GPU::GPU_storage_double2high;
            Analysis[i].number_of_series = run->number_of_parts;
        }

        Analysis[DM_Wilson_loop].storage_type         = GPU_CL::GPU::GPU_storage_double;

    if (run->get_actions_avr) {
        // S_spat
        Analysis[DM_S_spat].pointer         = GPU0->buffer_map(lattice_energies);
#if (MODEL_ON == 1)
        Analysis[DM_S_spat].denominator     = ((double) (lattice_full_site));
#else
        Analysis[DM_S_spat].denominator     = ((double) (lattice_full_site   * number_spat));
#endif
        Analysis[DM_S_spat].data_name       = "S_spat";
        D_A->lattice_data_analysis(&Analysis[DM_S_spat]);

        // S_temp
        Analysis[DM_S_temp].pointer         = Analysis[DM_S_spat].pointer;
#if (MODEL_ON == 1)
        Analysis[DM_S_temp].denominator     = ((double) (lattice_full_site));
#else
        Analysis[DM_S_temp].denominator     = ((double) (lattice_full_site   * number_temp));
#endif
        Analysis[DM_S_temp].data_name       = "S_temp";
        D_A->lattice_data_analysis(&Analysis[DM_S_temp]);

        // S_total
        Analysis[DM_S_total].data_name      = "S_total";
#if (MODEL_ON == 1)
        D_A->lattice_data_analysis_joint(&Analysis[DM_S_total],&Analysis[DM_S_spat],&Analysis[DM_S_spat]);
#else
        D_A->lattice_data_analysis_joint(&Analysis[DM_S_total],&Analysis[DM_S_spat],&Analysis[DM_S_temp]);
#endif
    }
    if (run->get_plaquettes_avr) {
        // Plq_spat
        Analysis[DM_Plq_spat].pointer         = GPU0->buffer_map(lattice_energies_plq);
#if (MODEL_ON == 1)
        Analysis[DM_Plq_spat].data_name       = "Field";
        Analysis[DM_Plq_spat].denominator     = ((double) (lattice_full_site));
#else
        Analysis[DM_Plq_spat].data_name       = "Plq_spat";
        Analysis[DM_Plq_spat].denominator     = ((double) (lattice_full_site   * number_spat));
#endif
        D_A->lattice_data_analysis(&Analysis[DM_Plq_spat]);

        // Plq_temp
        Analysis[DM_Plq_temp].pointer         = Analysis[DM_Plq_spat].pointer;
#if (MODEL_ON == 1)
        Analysis[DM_Plq_temp].data_name       = "Field^2";
        Analysis[DM_Plq_temp].denominator     = ((double) (lattice_full_site));
#else
        Analysis[DM_Plq_temp].data_name       = "Plq_temp";
        Analysis[DM_Plq_temp].denominator     = ((double) (lattice_full_site   * number_temp));
#endif
        D_A->lattice_data_analysis(&Analysis[DM_Plq_temp]);

        // Plq_total or (Mean_field and Variance)
#if (MODEL_ON == 1)
        Analysis[DM_Plq_total].data_name      = "Field_variance";
        D_A->lattice_data_analysis_joint_variance(&Analysis[DM_Plq_total],&Analysis[DM_Plq_spat],&Analysis[DM_Plq_temp]);
#else
        Analysis[DM_Plq_total].data_name      = "Plq_total";
        D_A->lattice_data_analysis_joint(&Analysis[DM_Plq_total],&Analysis[DM_Plq_spat],&Analysis[DM_Plq_temp]);
#endif
    }
    if (run->get_correlators) {
        // Adjustable correlator
        Analysis[DM_Correlator1].pointer      = GPU0->buffer_map(lattice_correlators);
        Analysis[DM_Correlator1].data_name    = "Correlator";
        Analysis[DM_Correlator1].denominator  = ((double) (lattice_full_site));

        if (run->get_plaquettes_avr) 
            D_A->lattice_data_analysis_sub(&Analysis[DM_Correlator1],&Analysis[DM_Plq_spat]);
        else
            D_A->lattice_data_analysis(&Analysis[DM_Correlator1]);

        // Correlator (+1,+1,+1,+1)
        Analysis[DM_Correlator2].pointer      = Analysis[DM_Correlator1].pointer;
        Analysis[DM_Correlator2].data_name    = "Correlator(+1)";
        Analysis[DM_Correlator2].denominator  = ((double) (lattice_full_site));
        if (run->get_plaquettes_avr) 
            D_A->lattice_data_analysis_sub(&Analysis[DM_Correlator2],&Analysis[DM_Plq_spat]);
        else
            D_A->lattice_data_analysis(&Analysis[DM_Correlator2]);
    }    
    if (run->PL_level > 0) {
        // Polyakov_loop
        Analysis[DM_Polyakov_loop].pointer            = GPU0->buffer_map(lattice_polyakov_loop);
        Analysis[DM_Polyakov_loop].denominator        = ((double) (lattice_full_n1n2n3 * run->lattice_group));
        Analysis[DM_Polyakov_loop].data_name          = "Polyakov_loop";
        D_A->lattice_data_analysis(&Analysis[DM_Polyakov_loop]);

        // Polyakov_loop_im
        Analysis[DM_Polyakov_loop_im].pointer         = Analysis[DM_Polyakov_loop].pointer;
        Analysis[DM_Polyakov_loop_im].denominator     = ((double) (lattice_full_n1n2n3 * run->lattice_group));
        Analysis[DM_Polyakov_loop_im].data_name       = "Polyakov_loop_im";
        D_A->lattice_data_analysis(&Analysis[DM_Polyakov_loop_im]);
    }
    if (run->PL_level > 1){
        // Polyakov_loop_P2
        Analysis[DM_Polyakov_loop_P2].pointer         = Analysis[DM_Polyakov_loop].pointer;
        Analysis[DM_Polyakov_loop_P2].pointer_offset  = lattice_polyakov_loop_size;
        Analysis[DM_Polyakov_loop_P2].denominator     = ((double) (lattice_full_n1n2n3 * run->lattice_group * run->lattice_group));
        Analysis[DM_Polyakov_loop_P2].data_name       = "Polyakov_loop_P2";
        D_A->lattice_data_analysis(&Analysis[DM_Polyakov_loop_P2]);

        // Polyakov_loop_P4
        Analysis[DM_Polyakov_loop_P4].pointer         = Analysis[DM_Polyakov_loop].pointer;
        Analysis[DM_Polyakov_loop_P4].pointer_offset  = lattice_polyakov_loop_size;
        Analysis[DM_Polyakov_loop_P4].denominator     = ((double) (lattice_full_n1n2n3 * run->lattice_group * run->lattice_group * run->lattice_group * run->lattice_group));
        Analysis[DM_Polyakov_loop_P4].data_name       = "Polyakov_loop_P4";
        D_A->lattice_data_analysis(&Analysis[DM_Polyakov_loop_P4]);
    }
    if (run->get_wilson_loop) {
        // Wilson_loop
        Analysis[DM_Wilson_loop].pointer         = GPU0->buffer_map(lattice_wilson_loop);
        Analysis[DM_Wilson_loop].denominator     = ((double) (lattice_full_site * number_spat));
        Analysis[DM_Wilson_loop].data_name       = "Wilson_loop";
        if (big_lattice) Analysis[DM_Wilson_loop].skip_checking = true;
        D_A->lattice_data_analysis(&Analysis[DM_Wilson_loop]);
    }
    if ((run->get_Fmunu)||(run->get_F0mu)) {
        // Fmunu_xy_3_re
        unsigned int* F_pointr;
        if (!Analysis[DM_Plq_spat].pointer)
            F_pointr = GPU0->buffer_map(lattice_energies_plq);
        else
            F_pointr = Analysis[DM_Plq_spat].pointer;

        int index_x = 0;
        int index_y = 1;
        for (int i=0;i<(run->lattice_nd-1)*2;i++){   // loop for re and im
            Analysis[i+DM_Fmunu_3].denominator = ((double) (lattice_full_site));
            Analysis[i+DM_Fmunu_8].denominator = ((double) (lattice_full_site));

            Analysis[i+DM_Fmunu_3].pointer = F_pointr;
            Analysis[i+DM_Fmunu_8].pointer = F_pointr;

            Analysis[i+DM_Fmunu_3].pointer_offset = lattice_energies_offset * (1 + (i >> 1));
            Analysis[i+DM_Fmunu_8].pointer_offset = lattice_energies_offset * (4 + (i >> 1));

            Analysis[i+DM_Fmunu_3].data_name = (char*) calloc(14,sizeof(char));
            Analysis[i+DM_Fmunu_8].data_name = (char*) calloc(14,sizeof(char));

            if (run->get_Fmunu) {
                sprintf_s((char*) Analysis[i+DM_Fmunu_3].data_name,14,"Fmunu_%s%s_%u_%s",FXYZ[index_x],FXYZ[index_y],run->Fmunu_index1,REIM[(i&1)]);
                sprintf_s((char*) Analysis[i+DM_Fmunu_8].data_name,14,"Fmunu_%s%s_%u_%s",FXYZ[index_x],FXYZ[index_y],run->Fmunu_index2,REIM[(i&1)]);
            } else {
                sprintf_s((char*) Analysis[i+DM_Fmunu_3].data_name,14,"Fmunu_%st_%u_%s",FXYZ[((i>>1)&3)],run->Fmunu_index1,REIM[(i&1)]);
                sprintf_s((char*) Analysis[i+DM_Fmunu_8].data_name,14,"Fmunu_%st_%u_%s",FXYZ[((i>>1)&3)],run->Fmunu_index2,REIM[(i&1)]);
            }


            D_A->lattice_data_analysis(&Analysis[i+DM_Fmunu_3]);
            D_A->lattice_data_analysis(&Analysis[i+DM_Fmunu_8]);

            if ((i&1)==1) {
                index_y++;
                if (index_y>(run->lattice_nd-2)) {
                    index_x++;
                    index_y=index_x+1;
                }
            }
        }

        Analysis[DM_Fmunu_abs_3_re].data_name = (char*) calloc(15,sizeof(char));
        Analysis[DM_Fmunu_abs_3_im].data_name = (char*) calloc(15,sizeof(char));
        Analysis[DM_Fmunu_abs_8_re].data_name = (char*) calloc(15,sizeof(char));
        Analysis[DM_Fmunu_abs_8_im].data_name = (char*) calloc(15,sizeof(char));

        sprintf_s((char*) Analysis[DM_Fmunu_abs_3_re].data_name,15,"Fmunu_abs_%1u_re",run->Fmunu_index1);
        sprintf_s((char*) Analysis[DM_Fmunu_abs_3_im].data_name,15,"Fmunu_abs_%1u_im",run->Fmunu_index1);
        sprintf_s((char*) Analysis[DM_Fmunu_abs_8_re].data_name,15,"Fmunu_abs_%1u_re",run->Fmunu_index2);
        sprintf_s((char*) Analysis[DM_Fmunu_abs_8_im].data_name,15,"Fmunu_abs_%1u_im",run->Fmunu_index2);

        // Fmunu_abs_3
        D_A->lattice_data_analysis_joint3(&Analysis[DM_Fmunu_abs_3_re],&Analysis[DM_Fmunu_xy_3_re],&Analysis[DM_Fmunu_xz_3_re],&Analysis[DM_Fmunu_yz_3_re]);
        D_A->lattice_data_analysis_joint3(&Analysis[DM_Fmunu_abs_3_im],&Analysis[DM_Fmunu_xy_3_im],&Analysis[DM_Fmunu_xz_3_im],&Analysis[DM_Fmunu_yz_3_im]);

        // Fmunu_abs_8
        D_A->lattice_data_analysis_joint3(&Analysis[DM_Fmunu_abs_8_re],&Analysis[DM_Fmunu_xy_8_re],&Analysis[DM_Fmunu_xz_8_re],&Analysis[DM_Fmunu_yz_8_re]);
        D_A->lattice_data_analysis_joint3(&Analysis[DM_Fmunu_abs_8_im],&Analysis[DM_Fmunu_xy_8_im],&Analysis[DM_Fmunu_xz_8_im],&Analysis[DM_Fmunu_yz_8_im]);
    }
    if (run->get_acceptance_rate) {
        // Acceptance_rate_even
        Analysis[DM_Acc_rate_even].pointer  = GPU0->buffer_map(lattice_acceptance_rate);
        Analysis[DM_Acc_rate_even].skip_checking   = true;
#if (MODEL_ON == 1)
        Analysis[DM_Acc_rate_even].denominator     = ((double) (lattice_full_site * run->NHIT * run->NITER));
#else
        Analysis[DM_Acc_rate_even].denominator     = ((double) (lattice_full_site * run->NHIT * run->NITER * run->lattice_nd * lattice_group_elements[run->lattice_group]) / 4.0 / 2.0); // 4 -> number of atom updates, 2->even/odd=half of lattice
#endif
        Analysis[DM_Acc_rate_even].data_name       = "AR_even";
        D_A->lattice_data_analysis(&Analysis[DM_Acc_rate_even]);

        // Acceptance_rate_odd
        Analysis[DM_Acc_rate_odd].pointer         = Analysis[DM_Acc_rate_even].pointer;
        Analysis[DM_Acc_rate_odd].skip_checking   = true;
#if (MODEL_ON == 1)
        Analysis[DM_Acc_rate_odd].denominator     = ((double) (lattice_full_site * run->NHIT * run->NITER));
#else
        Analysis[DM_Acc_rate_odd].denominator     = ((double) (lattice_full_site * run->NHIT * run->NITER * run->lattice_nd * lattice_group_elements[run->lattice_group]) / 4.0 / 2.0); // 4 -> number of atom updates, 2->even/odd=half of lattice
#endif
        Analysis[DM_Acc_rate_odd].data_name       = "AR_odd";
        D_A->lattice_data_analysis(&Analysis[DM_Acc_rate_odd]);

        // S_total
        Analysis[DM_Acc_rate_total].data_name     = "AR_total";
        Analysis[DM_Acc_rate_total].skip_checking = true;
        D_A->lattice_data_analysis_joint(&Analysis[DM_Acc_rate_total],&Analysis[DM_Acc_rate_even],&Analysis[DM_Acc_rate_odd]);
    }
    
}
void        model::lattice_write_results(void) {
    FILE *stream;
    char buffer[HGPU_MAX_STRINGLEN];
    int j;

    char* header2 = lattice_make_header2();
    printf("%s\n",header2);

    j  = sprintf_s(buffer  ,sizeof(buffer),  "%s",run->path);
    j += sprintf_s(buffer+j,sizeof(buffer)-j,"%s",run->fprefix);
    j += sprintf_s(buffer+j,sizeof(buffer)-j,"%.2s-%.3s-%.2s-%.2s-%.2s-%.2s",timeend+22,timeend+4,timeend+8,timeend+11,timeend+14,timeend+17);
    j += sprintf_s(buffer+j,sizeof(buffer)-j,".txt");

    fopen_s(&stream,buffer,"w+");
    if(stream)
    {
        fprintf(stream,header);

        // write plaquette data
        for (int i=0; i<run->ITER; i++) {
                fprintf(stream, "%5i",i);
            if (run->get_plaquettes_avr){
                fprintf(stream, " % 16.13e",Analysis[DM_Plq_spat].data[i]);
                fprintf(stream, " % 16.13e",Analysis[DM_Plq_temp].data[i]);
                fprintf(stream, " % 16.13e",Analysis[DM_Plq_total].data[i]);
            }
            if (run->get_wilson_loop)
                fprintf(stream, " % 16.13e",Analysis[DM_Wilson_loop].data[i]);
            if (run->get_actions_avr){
                fprintf(stream, " % 16.13e",Analysis[DM_S_spat].data[i]);
                fprintf(stream, " % 16.13e",Analysis[DM_S_temp].data[i]);
                fprintf(stream, " % 16.13e",Analysis[DM_S_total].data[i]);
            }
            if (run->get_correlators){
                fprintf(stream, " % 16.13e",Analysis[DM_Correlator1].data[i]);
                fprintf(stream, " % 16.13e",Analysis[DM_Correlator2].data[i]);
            }
            if (run->PL_level > 0){
                fprintf(stream, " % 16.13e",Analysis[DM_Polyakov_loop].data[i]);
                fprintf(stream, " % 16.13e",Analysis[DM_Polyakov_loop_im].data[i]);
            }
            if (run->PL_level > 1){
                fprintf(stream, " % 16.13e",Analysis[DM_Polyakov_loop_P2].data[i]);
                fprintf(stream, " % 16.13e",Analysis[DM_Polyakov_loop_P4].data[i]);
            }
            if (run->get_acceptance_rate){
                fprintf(stream, " % 16.13e",Analysis[DM_Acc_rate_even].data[i]);
                fprintf(stream, " % 16.13e",Analysis[DM_Acc_rate_odd].data[i]);
                fprintf(stream, " % 16.13e",Analysis[DM_Acc_rate_total].data[i]);
            }
            if ((run->get_Fmunu)||(run->get_F0mu))
                for (int jj=0;jj<((run->lattice_nd-2)*(run->lattice_nd-1)+2)*2;jj++){
                    fprintf(stream, " % 16.13e",Analysis[jj+DM_Fmunu_3].data[i]);
                }
            fprintf(stream, "\n");


        }

        if ( fclose(stream) ) {printf( "The file was not closed!\n" ); }
    }
}
void        model::lattice_write_configuration(void) {
    FILE *stream;
    char buffer[HGPU_MAX_STRINGLEN];
    int j;

    j  = sprintf_s(buffer  ,sizeof(buffer),  "%s",run->path);
    j += sprintf_s(buffer+j,sizeof(buffer)-j,"%s",run->fprefix);
    j += sprintf_s(buffer+j,sizeof(buffer)-j,"%.2s-%.3s-%.2s-%.2s-%.2s-%.2s",timeend+22,timeend+4,timeend+8,timeend+11,timeend+14,timeend+17);
    j += sprintf_s(buffer+j,sizeof(buffer)-j,".cnf");

    fopen_s(&stream,buffer,"w+");
    if(stream)
    {
        fprintf(stream,header);

#if (MODEL_ON == 1)
        // write configuration data
        fprintf(stream, " Data fields:\n");
        fprintf(stream, "[ T, Z, Y, X]\n");
        fprintf(stream, " ***************************************************\n");
        for (int x1 = 0; x1 < run->lattice_domain_size[0]; x1++)
        for (int x4 = 0; x4 < run->lattice_domain_size[3]; x4++)
        for (int x3 = 0; x3 < run->lattice_domain_size[2]; x3++)
        for (int x2 = 0; x2 < run->lattice_domain_size[1]; x2++)
        {
            unsigned int gid = x2 + x3 * run->lattice_domain_size[1] + x4 * lattice_domain_n2n3 + x1 * lattice_domain_n2n3n4;
            double site = 0.0;
            if (run->precision == model::model_precision_single) {
                site = GPU0->convert_to_double(lattice_pointer_last[gid]);
            } else {
                site = GPU0->convert_to_double(lattice_pointer_last[gid*2],lattice_pointer_last[gid*2 + 1]);
            }
                fprintf(stream, "[%2u,%2u,%2u,%2u] % 17.16e\n",x4,x3,x2,x1,site);
        }
#endif

        if ( fclose(stream) ) {printf( "The file was not closed!\n" ); }
    }
}

void        model::lattice_save_state(void){
    time_t ltimesave;
    time(&ltimesave);
    char* timesave   = GPU0->get_current_datetime();

    lattice_pointer_save       = GPU0->buffer_map(lattice_table);
    unsigned int* lattice_measurement_save   = GPU0->buffer_map(lattice_measurement);
    unsigned int* lattice_energies_save      = GPU0->buffer_map(lattice_energies);
    unsigned int* lattice_energies_plq_save  = NULL;
    unsigned int* lattice_wilson_loop_save   = NULL;
    unsigned int* lattice_polyakov_loop_save = NULL;

    if ((run->get_plaquettes_avr) || (run->get_Fmunu) || (run->get_F0mu))
        lattice_energies_plq_save  = GPU0->buffer_map(lattice_energies_plq);
    if (run->get_wilson_loop)
        lattice_wilson_loop_save   = GPU0->buffer_map(lattice_wilson_loop);
    if (run->PL_level > 0)
        lattice_polyakov_loop_save = GPU0->buffer_map(lattice_polyakov_loop);

    FILE *stream;
    char buffer[HGPU_MAX_STRINGLEN];
    int j = 0;

    unsigned int* head = lattice_make_bin_header();

    j  = sprintf_s(buffer  ,sizeof(buffer),  "%s",run->path);
    j += sprintf_s(buffer+j,sizeof(buffer)-j,"%s",run->fprefix);
    j += sprintf_s(buffer+j,sizeof(buffer)-j,"%.2s-%.3s-%.2s-%.2s-%.2s-%.2s",timesave+22,timesave+4,timesave+8,timesave+11,timesave+14,timesave+17);
    if (big_lattice)
        j += sprintf_s(buffer+j,sizeof(buffer)-j,"_%u",run->part_number);
    j += sprintf_s(buffer+j,sizeof(buffer)-j,".qcg");

    fopen_s(&stream,buffer,"wb");
    if(stream)
    {
        fwrite(head,sizeof(unsigned int),BIN_HEADER_SIZE,stream);                                       // write header
        fwrite(lattice_measurement_save, sizeof(cl_double2), lattice_measurement_size_F, stream);       // write measurements
        fwrite(lattice_energies_save,    sizeof(cl_double2), lattice_energies_size, stream);            // write energies
        if ((run->get_plaquettes_avr) || (run->get_Fmunu) || (run->get_F0mu))
            fwrite(lattice_energies_plq_save,  sizeof(cl_double2), lattice_energies_size_F, stream);    // write energies_plq
        if (run->get_wilson_loop)
            fwrite(lattice_wilson_loop_save,   sizeof(cl_double),  lattice_energies_size, stream);      // write wilson loop
        if (run->PL_level > 0)
            fwrite(lattice_polyakov_loop_save, sizeof(cl_double2), lattice_polyakov_loop_size, stream); // write polyakov loop
        if (run->precision == model_precision_single)                                                        // write configuration
            fwrite(lattice_pointer_save, sizeof(cl_float4), lattice_table_size, stream);
        else
            fwrite(lattice_pointer_save, sizeof(cl_double4), lattice_table_size, stream);

        unsigned int hlen  = BIN_HEADER_SIZE*sizeof(unsigned int);
            if (GPU0->GPU_debug->brief_report) printf("Header: 0x%X-0x%X\n",0,hlen);
        unsigned int hlen2 = lattice_measurement_size_F * sizeof(cl_double2);
            if (GPU0->GPU_debug->brief_report) printf("Lattice mesurements: 0x%X-0x%X\n",hlen,(hlen+hlen2));
            hlen += hlen2;
        hlen2 = lattice_energies_size * sizeof(cl_double2);
            if (GPU0->GPU_debug->brief_report) printf("Lattice energies: 0x%X-0x%X\n",hlen,(hlen+hlen2));
            hlen += hlen2;
        if ((run->get_plaquettes_avr) || (run->get_Fmunu) || (run->get_F0mu)){
            hlen2 = lattice_energies_size_F * sizeof(cl_double2);
            if (GPU0->GPU_debug->brief_report) printf("Lattice plaq av.: 0x%X-0x%X\n",hlen,(hlen+hlen2));
            hlen += hlen2;
        }
        if (run->get_wilson_loop){
            hlen2 = lattice_energies_size * sizeof(cl_double);
            if (GPU0->GPU_debug->brief_report) printf("Lattice Wilson loop: 0x%X-0x%X\n",hlen,(hlen+hlen2));
            hlen += hlen2;
        }
        if (run->PL_level > 0){
            hlen2 = lattice_polyakov_loop_size * sizeof(cl_double2);
            if (GPU0->GPU_debug->brief_report) printf("Lattice Polyakov loop: 0x%X-0x%X\n",hlen,(hlen+hlen2));
            hlen += hlen2;
        }
        if (run->precision == model_precision_single)                                                        // write configuration
            hlen2 = lattice_table_size * sizeof(cl_float4);
        else
            hlen2 = lattice_table_size * sizeof(cl_double4);
        if (GPU0->GPU_debug->brief_report) printf("Lattice data: 0x%X-0x%X\n",hlen,(hlen+hlen2));
        hlen += hlen2;


        if ( fclose(stream) ) printf( "The file was not closed!\n" );
    }
    FREE(head);
}
unsigned int*   model::lattice_make_bin_header(void){
    int k;
    // bin header structure:
    // 0 - 1 - prefix (8 bytes)
    // 2 - 3 - version
    // 4 - ...
    unsigned int* result = (unsigned int*) calloc(BIN_HEADER_SIZE,sizeof(unsigned int));
    const char* bin_prefix = "QCDGPU";
    k =  0; for(int i=0; i<ceil((double) strlen_s(bin_prefix)/4); i++)   result[k++] = convert_str_uint(bin_prefix,i*4);
    k =  2; for(int i=0; i<ceil((double) strlen_s(run->version)/4); i++) result[k++] = convert_str_uint(run->version,i*4);
    k =  4;
    result[k++] = run->INIT;                        // 0x10
    result[k++] = convert_start_to_uint(run->ints); // 0x14
    result[k++] = run->run_PRNG->PRNG_randseries;   // 0x18
    result[k++] = PRNG_CL::PRNG::convert_generator_to_uint(run->run_PRNG->PRNG_generator);  // 0x1C
    result[k++] = run->run_PRNG->RL_nskip;          // 0x20
    result[k++] = PRNG0->PRNG_counter;              // 0x24
    result[k++] = run->NAV;                         // 0x28
    result[k++] = NAV_counter;                      // 0x2C
    result[k++] = run->NITER;                       // 0x30
    result[k++] = run->ITER;                        // 0x34
    result[k++] = ITER_counter;                     // 0x38
    result[k++] = run->NHIT;                        // 0x3C
    result[k++] = run->wilson_R;                    // 0x40
    result[k++] = run->wilson_T;                    // 0x44
    result[k++] = convert_precision_to_uint(run->precision);// 0x48
    result[k++] = GPU0->convert_to_uint_LOW( run->BETA);    // 0x4C
    result[k++] = GPU0->convert_to_uint_HIGH(run->BETA);    // 0x50
    result[k++] = GPU0->convert_to_uint_LOW( run->PHI);     // 0x54
    result[k++] = GPU0->convert_to_uint_HIGH(run->PHI);     // 0x58
    result[k++] = GPU0->convert_to_uint_LOW( run->OMEGA);   // 0x5C
    result[k++] = GPU0->convert_to_uint_HIGH(run->OMEGA);   // 0x60
    result[k++] = run->PL_level;                    // 0x64
    result[k++] = run->Fmunu_index1;                // 0x68
    result[k++] = run->Fmunu_index2;                // 0x6C
    result[k++] = run->get_Fmunu;                   // 0x70
    result[k++] = run->get_F0mu;                    // 0x74
    result[k++] = run->get_actions_avr;             // 0x78
    result[k++] = run->get_plaquettes_avr;          // 0x7C
    result[k++] = run->get_wilson_loop;             // 0x80
    result[k++] = run->get_acceptance_rate;         // 0x84
    result[k++] = lattice_measurement_size_F;       // 0x88 - measurements initialization
    result[k++] = lattice_energies_size;            // 0x8C - energies initialization
    result[k++] = lattice_energies_size_F;          // 0x90 - energies initialization (with Fmunu tensor, if needed)
    result[k++] = lattice_polyakov_loop_size;       // 0x94 - polyakov loop initialization
    result[k++] = run->lattice_type;                // 0x98
    result[k++] = run->lattice_group;               // 0x9C
    result[k++] = run->lattice_nd;                  // 0xA0
    for (int i=0; i<run->lattice_nd; i++) result[k++] = run->lattice_full_size[i];  // 0xA4 0xA8 0xAC 0xB0
    for (int i=0; i<run->lattice_nd; i++) result[k++] = run->lattice_domain_size[i];// 0xB4 0xB8 0xBC 0xC0

    return result;
}
bool        model::lattice_load_bin_header(unsigned int* head){
    bool result = true;
    int k;
    // bin header structure:
    // 0 - 1 - prefix (8 bytes)
    // 2 - 3 - version
    // 4 - ...
    const char* bin_prefix = "QCDGPU";
    k =  0; for(int i=0; i<ceil((double) strlen_s(bin_prefix)/4); i++)   result &= (head[k++] == convert_str_uint(bin_prefix,i*4));
    k =  2; for(int i=0; i<ceil((double) strlen_s(run->version)/4); i++) result &= (head[k++] == convert_str_uint(run->version,i*4));
    if (!result) return result;
    k =  4;
      k++;
    run->ints = convert_uint_to_start(head[k++]);       // 0x14
    run->run_PRNG->PRNG_randseries = head[k++];         // 0x18
    run->run_PRNG->PRNG_generator = PRNG0->convert_uint_to_generator(head[k++]);    // 0x1C
    run->run_PRNG->RL_nskip = head[k++];                // 0x20
    PRNG_counter = head[k++];                           // 0x24
    run->NAV = head[k++];                               // 0x28
    NAV_counter = head[k++];                            // 0x2C
    run->NITER = head[k++];                             // 0x30
    run->ITER = head[k++];                              // 0x34
    ITER_counter = head[k++];                           // 0x38
    run->NHIT = head[k++];                              // 0x3C
    run->wilson_R = head[k++];                          // 0x40
    run->wilson_T = head[k++];                          // 0x44
    run->precision = convert_uint_to_precision(head[k++]);  // 0x48
    unsigned int get_low  = head[k++];                  // 0x4C
    unsigned int get_high = head[k++];                  // 0x50
    run->BETA  = GPU0->convert_to_double(get_low,get_high);
    get_low  = head[k++]; get_high = head[k++];         // 0x54, 0x58
    run->PHI   = GPU0->convert_to_double(get_low,get_high);
    get_low  = head[k++]; get_high = head[k++];         // 0x5C, 0x60
    run->OMEGA = GPU0->convert_to_double(get_low,get_high);
    run->PL_level = head[k++];                          // 0x64
    run->Fmunu_index1 = head[k++];                      // 0x68
    run->Fmunu_index2 = head[k++];                      // 0x6C
    run->get_Fmunu = (head[k++]==0 ? false : true);     // 0x70
    run->get_F0mu  = (head[k++]==0 ? false : true);     // 0x74
    run->get_actions_avr = (head[k++]==0 ? false : true);    // 0x78
    run->get_plaquettes_avr = (head[k++]==0 ? false : true); // 0x7C
    run->get_wilson_loop = (head[k++]==0 ? false : true);    // 0x80
    run->get_acceptance_rate = (head[k++]==0 ? false : true);// 0x84
    lattice_measurement_size_F = head[k++];             // 0x88 - measurements initialization
    lattice_energies_size = head[k++];                  // 0x8C - energies initialization
    lattice_energies_size_F = head[k++];                // 0x90 - energies initialization (with Fmunu tensor, if needed)
    lattice_polyakov_loop_size = head[k++];             // 0x94 - polyakov loop initialization
    run->lattice_type = head[k++];                      // 0x98
    run->lattice_group = head[k++];                     // 0x9C
    run->lattice_nd = head[k++];                        // 0xA0
    for (int i=0; i<run->lattice_nd; i++) run->lattice_full_size[i] = head[k++];   // 0xA4, 0xA8, 0xAC, 0xB0
    for (int i=0; i<run->lattice_nd; i++) run->lattice_domain_size[i] = head[k++]; // 0xB4, 0xB8, 0xBC, 0xC0

    return result;
}
void        model::lattice_load_state(void){
    FILE *stream;
    char buffer[HGPU_MAX_STRINGLEN];
    int j = 0;

    unsigned int* head = (unsigned int*) calloc(BIN_HEADER_SIZE,sizeof(unsigned int));
    bool result = true;

    j  = sprintf_s(buffer  ,sizeof(buffer),  "%s",run->path);
    j += sprintf_s(buffer+j,sizeof(buffer)-j,"%s",run->fstate);

    fopen_s(&stream,buffer,"rb");
    if(stream)
    {
        fread(head,sizeof(unsigned int),BIN_HEADER_SIZE,stream);    // load header
        if (LOAD_state==0) result = lattice_load_bin_header(head);
        if (!result) printf("[ERROR in header!!!]\n");
        if (LOAD_state==1) {
            fread(plattice_measurement, sizeof(cl_double2), lattice_measurement_size_F, stream);       // load measurements
            fread(plattice_energies,    sizeof(cl_double2), lattice_energies_size, stream);            // load energies
            if ((run->get_plaquettes_avr) || (run->get_Fmunu) || (run->get_F0mu))
                fread(plattice_energies_plq,  sizeof(cl_double2), lattice_energies_size_F, stream);    // load energies_plq
            if (run->get_wilson_loop)
                fread(plattice_wilson_loop,   sizeof(cl_double),  lattice_energies_size, stream);      // load wilson loop
            if (run->PL_level > 0)
                fread(plattice_polyakov_loop, sizeof(cl_double2), lattice_polyakov_loop_size, stream); // load polyakov loop
            if (run->precision == model_precision_single)                                              // load configuration
                fread(plattice_table_float,   sizeof(cl_float4),  lattice_table_size, stream);
            else
                fread(plattice_table_double,  sizeof(cl_double4), lattice_table_size, stream);

            if ( fclose(stream) ) printf( "The file was not closed!\n" );
        }
        LOAD_state++;
    }
    FREE(head);
}
 
void        model::model_lattice_init(void){
    if ((run->get_Fmunu)&&(run->get_F0mu)) run->get_Fmunu = false;  // only one field (H or E) may be calculated

    run->ON_sqrt_z_lambda_zeta = (0.25*sqrt(run->ON_z/(run->ON_lambda*run->ON_zeta)));
    // setup U_max for O(N)
    if (run->ON_eta > 0.0){
        double maxU1 = 1.0 - exp((1.0-sqrt(1.0+4.0*sqrt(2.0)*run->ON_eta/sqrt(run->ON_b)))/run->ON_eta);
        double maxU2 = 1.0 - exp((1.0+sqrt(1.0+4.0*sqrt(2.0)*run->ON_eta/sqrt(run->ON_b)))/run->ON_eta);
        if (maxU1<0.0) maxU1=maxU2;
        run->ON_max_U = maxU1;
    } else {
        run->ON_max_U = 1.0 - exp(-1.9999999999/sqrt(run->ON_b));
    }

    local_size_intel = (GPU0->GPU_info.device_vendor == GPU0->GPU::GPU_vendor_Intel) ? 64 : 0;
    if (GPU0->GPU_limit_max_workgroup_size) local_size_intel = GPU0->GPU_limit_max_workgroup_size;
    size_t workgroup_factor = (local_size_intel) ? local_size_intel : 32;

    if (big_lattice) {
        // lattice is divided into domains
        lattice_domain_n1     = (run->lattice_domain_size[0]+2);
    } else {
        // lattice is not divided into domains
        lattice_domain_n1     = run->lattice_domain_size[0];
    }
    lattice_domain_n1n2   = lattice_domain_n1           * run->lattice_domain_size[1];
    lattice_domain_n2n3   = run->lattice_domain_size[1] * run->lattice_domain_size[2];
    lattice_domain_n1n2n3 = lattice_domain_n1           * lattice_domain_n2n3;
    lattice_domain_exact_n1n2n3 = run->lattice_domain_size[0] * run->lattice_domain_size[1] * run->lattice_domain_size[2];
    lattice_domain_n2n3n4 = run->lattice_domain_size[1] * run->lattice_domain_size[2] * run->lattice_domain_size[3];

    lattice_domain_site       = lattice_domain_n1;
    lattice_domain_exact_site = run->lattice_domain_size[0];
    lattice_full_site         = run->lattice_full_size[0];
    for (int i=1;i<run->lattice_nd;i++) {
        lattice_domain_site       *= run->lattice_domain_size[i];
        lattice_domain_exact_site *= run->lattice_domain_size[i];
        lattice_full_site         *= run->lattice_full_size[i];
    }

    lattice_full_n1n2   = run->lattice_full_size[0] * run->lattice_full_size[1];
    lattice_full_n2n3   = run->lattice_full_size[1] * run->lattice_full_size[2];
    lattice_full_n1n2n3 = run->lattice_full_size[0] * run->lattice_full_size[1] * run->lattice_full_size[2];
    lattice_full_n2n3n4 = run->lattice_full_size[1] * run->lattice_full_size[2] * run->lattice_full_size[3];

    lattice_full_link         = run->lattice_nd * lattice_full_site;
    lattice_domain_link       = run->lattice_nd * lattice_domain_site;
    lattice_domain_exact_link = run->lattice_nd * lattice_domain_exact_site;

    lattice_boundary_exact_size = run->lattice_domain_size[1];    // Size of boundary slice (x=1 in depth)
    for (int i=2;i<run->lattice_nd;i++) lattice_boundary_exact_size *= run->lattice_domain_size[i];

    lattice_parameters_size         = GPU0->buffer_size_align(MODEL_parameter_size);
    lattice_energies_size           = GPU0->buffer_size_align(run->ITER);                      // number of working iterations
    lattice_energies_size_F         = lattice_energies_size * MODEL_energies_size;        // number of working iterations (for tensor Fmunu)
    lattice_energies_offset         = lattice_energies_size;

    lattice_boundary_size           = GPU0->buffer_size_align(lattice_boundary_exact_size);
    lattice_table_row_size          = GPU0->buffer_size_align(lattice_domain_site);
    lattice_table_row_size_half     = GPU0->buffer_size_align(lattice_domain_site / 2);
    lattice_table_exact_row_size    = GPU0->buffer_size_align(lattice_domain_exact_site);
    lattice_table_exact_row_size_half=GPU0->buffer_size_align(lattice_domain_exact_site / 2);
    if (run->lattice_type==1) {
        if (run->lattice_group==1)
            lattice_table_size      = GPU0->buffer_size_align(lattice_table_row_size * lattice_group_elements[run->lattice_group]);
        if (run->lattice_group==2)
            lattice_table_size      = GPU0->buffer_size_align(lattice_table_row_size * lattice_group_elements[run->lattice_group] / 2);
        if (run->lattice_group==3)
            lattice_table_size      = GPU0->buffer_size_align(lattice_table_row_size * lattice_group_elements[run->lattice_group] / 3);
        if (run->lattice_group>4)
            lattice_table_size      = GPU0->buffer_size_align(lattice_table_row_size * lattice_group_elements[run->lattice_group] / 4);

        lattice_table_group         = GPU0->buffer_size_align(lattice_table_row_size);
        lattice_table_exact_group   = GPU0->buffer_size_align(lattice_table_exact_row_size);
    } else {
        lattice_table_size          = GPU0->buffer_size_align(lattice_table_row_size * run->lattice_nd * lattice_group_elements[run->lattice_group] / 4);
        lattice_table_group         = GPU0->buffer_size_align(lattice_table_row_size * run->lattice_nd);
        lattice_table_exact_group   = GPU0->buffer_size_align(lattice_table_exact_row_size * run->lattice_nd);
    }

    lattice_measurement_size        = GPU0->buffer_size_align((unsigned int) ceil((double) lattice_table_exact_row_size / workgroup_factor)); // workgroup_factor is the minimum number of workgroup items
    lattice_measurement_size_F      = lattice_measurement_size * MODEL_energies_size;
    lattice_measurement_offset      = lattice_measurement_size;

    lattice_polyakov_size           = GPU0->buffer_size_align((unsigned int) lattice_domain_exact_n1n2n3);
    lattice_polyakov_loop_size      = GPU0->buffer_size_align(run->ITER);      // number of working iterations
    lattice_polyakov_loop_offset    = GPU0->buffer_size_align((unsigned int) ceil((double) lattice_polyakov_size / workgroup_factor));

    if ((!run->get_Fmunu)&&(!run->get_F0mu)){
        lattice_energies_size_F    = lattice_energies_size;
        lattice_measurement_size_F = lattice_measurement_size;
    }

    //_____________________________________________ PRNG preparation
        run->run_PRNG->PRNG_instances   = 0;    // number of instances of generator (or 0 for autoselect)

#if (MODEL_ON==1)
        run->run_PRNG->PRNG_samples     = GPU0->buffer_size_align((unsigned int) ceil(double(lattice_table_row_size_half * 2 * ceil(0.25*double(run->NHIT))))); // 3*NHIT+1 PRNs per link
#else
        run->run_PRNG->PRNG_samples     = GPU0->buffer_size_align((unsigned int) ceil(double(3 * lattice_table_row_size_half * (run->NHIT + 1)))); // 3*NHIT+1 PRNs per link
#endif
        PRNG_CL::PRNG::PRNG_parameters_copy(run->run_PRNG,PRNG0->run_PRNG);
        PRNG0->GPU0 = GPU0;
        prngstep = lattice_table_row_size_half; // lattice_domain_exact_site / 2;
    //-----------------------------------------------------------------

    if (GPU0->GPU_debug->brief_report) {
        printf("Full lattice: -----------------------\n");
        printf("sites                       = %u\n",lattice_full_site);
        printf("n1n2n3                      = %u\n",lattice_full_n1n2n3);
        printf("Domain lattice: ---------------------\n");
        printf("sites                       = %u\n",lattice_domain_site);
        printf("link                        = %u\n",lattice_domain_link);
        printf("n1n2n3                      = %u\n",lattice_domain_n1n2n3);
        printf("sites (exact)               = %u\n",lattice_domain_exact_site);
        printf("links (exact)               = %u\n",lattice_domain_exact_link);
        printf("n1n2n3 (exact)              = %u\n",lattice_domain_exact_n1n2n3);
        printf("lattice_table_exact_row_size= %u\n",lattice_table_exact_row_size);
        printf("..table_exact_row_size_half = %u\n",lattice_table_exact_row_size_half);
        printf("-------------------------------------\n");
        printf("workgroup_size              = %u\n",(unsigned int) (GPU0->GPU_info.max_workgroup_size));
        printf("measurements                = %u\n",lattice_domain_exact_site / (unsigned int) GPU0->GPU_info.max_workgroup_size);
        printf("lattice_table_row_size      = %u\n",lattice_table_row_size);
        printf("lattice_table_row_size_half = %u\n",lattice_table_row_size_half);
        printf("lattice_table_size          = %u\n",lattice_table_size);
        printf("lattice_table_group         = %u\n",lattice_table_group);
        printf("lattice_measurement_size    = %u\n",lattice_measurement_size);
        printf("lattice_measurement_size_F  = %u\n",lattice_measurement_size_F);
        printf("lattice_parameters_size     = %u\n",lattice_parameters_size);
        printf("lattice_energies_size       = %u\n",lattice_energies_size);
        printf("lattice_energies_size_F     = %u\n",lattice_energies_size_F);
        printf("lattice_polyakov_size       = %u\n",lattice_polyakov_size);
        printf("lattice_polyakov_loop_size  = %u\n",lattice_polyakov_loop_size);
        printf("polykov_loop_offset         = %u\n",lattice_polyakov_loop_offset);
        printf("lattice_energies_offset     = %u\n",lattice_energies_offset);
        printf("local_size_intel            = %u\n", (unsigned int) local_size_intel);
        printf("kernels: ----------------------------\n");
        printf("lattice_init                = %u\n",lattice_table_group);
        printf("lattice_GramSchmidt         = %u\n",lattice_table_group);
        printf("lattice_measurement         = %u\n",lattice_table_row_size);
        printf("lattice_measurement_plq     = %u\n",lattice_table_row_size);
        printf("lattice_measurement_wilson  = %u\n",lattice_table_row_size);
        printf("update                      = %u\n",lattice_table_row_size_half);
        printf("lattice_polyakov            = %u\n",lattice_polyakov_size);
        printf("PRNs (in quads): --------------------\n");
        printf("PRNG_instances              = %u\n",PRNG0->run_PRNG->PRNG_instances);
        printf("PRNG_samples                = %u\n",PRNG0->run_PRNG->PRNG_samples);
        printf("PRNs total                  = %u\n",PRNG0->run_PRNG->PRNG_instances * PRNG0->run_PRNG->PRNG_samples);
        printf("PRNG_step                   = %u\n",prngstep);
#if (MODEL_ON == 1)
        printf("PRNs for hot start          = %u\n",lattice_table_row_size);
        printf("PRNs for one update         = %u\n",int(ceil(double(lattice_table_row_size_half * run->NHIT>>1))));
#else
        printf("PRNs for hot start          = %u\n",lattice_table_row_size * 3);
        printf("PRNs for one update         = %u\n",int(ceil(double(3 * lattice_domain_link * (run->NHIT + 1)))));
#endif
        printf("-------------------------------------\n");
        printf(" ROWSIZE                    = %u\n",lattice_table_row_size);
        printf(" NHIT                       = %u\n",run->NHIT);
        printf(" PRNGSTEP                   = %u\n",lattice_table_row_size_half);
        printf(" PRECISION                  = %u\n",run->precision);
        printf(" PL                         = %u\n",run->PL_level);
        printf(" ITER                       = %u\n",run->ITER);
        if (!((run->PHI==0.0)&&(run->OMEGA==0.0))) printf(" TBC\n");     // turn on TBC
        if (run->get_Fmunu) {
            printf(" FMUNU (%u, %u)\n",run->Fmunu_index1,run->Fmunu_index2);   // calculate tensor Fmunu
        }
        if (run->get_F0mu) {
            printf(" F0MU (%u, %u)\n", run->Fmunu_index1,run->Fmunu_index2);   // calculate tensor Fmunu
        }
        printf("-------------------------------------\n");
        printf(" Root path to .CL files     = %s\n",GPU0->cl_root_path);
        printf("-------------------------------------\n");
    }

    //_____________________________________________ PRNG initialization
        PRNG0->initialize();
            GPU0->print_stage("PRNGs initialized");
    //-----------------------------------------------------------------

    char* header = lattice_make_header();
    printf("%s\n",header);


    lattice_create_buffers();
    lattice_make_programs();


    rowsize      = lattice_table_row_size;
    rowsize4     = 4 * lattice_table_row_size;
    halfrowsize  = lattice_domain_link / 8;
    halfrowsize4 = lattice_domain_link / 2;

    GPU0->print_memory_utilized();

    timestart = GPU0->get_current_datetime();
    time(&ltimestart);

    lattice_pointer_initial      = NULL;
    lattice_pointer_measurements = NULL;

    printf("\nrun kernels on GPU (%f seconds)\n",GPU0->get_timer_CPU(TIMER_FOR_ELAPSED));
    GPU0->start_timer_CPU(TIMER_FOR_SIMULATIONS); // start GPU execution timer
    GPU0->start_timer_CPU(TIMER_FOR_SAVE);        // start timer for lattice_state save
}
void        model::lattice_make_programs(void){
    int j;
    int argument_id;
    wilson_index          = 0;
    plq_index             = 0;
    polyakov_index        = 0;
    measurement_index     = 0;
    acceptance_rate_index = 0;

    char options_common[1024];
    int options_length_common  = sprintf_s(options_common,sizeof(options_common),"-Werror");
    if (GPU0->GPU_info.device_vendor == GPU_CL::GPU::GPU_vendor_Intel)
        options_length_common += sprintf_s(options_common + options_length_common,sizeof(options_common)-options_length_common," -D INTEL_ON");
    if (run->lattice_type==1)
        options_length_common += sprintf_s(options_common + options_length_common,sizeof(options_common)-options_length_common," -D ON_MODEL");
    if (run->get_acceptance_rate)
        options_length_common += sprintf_s(options_common + options_length_common,sizeof(options_common)-options_length_common," -D ACC_RATE");
    options_length_common += sprintf_s(options_common + options_length_common,sizeof(options_common)-options_length_common," -D SUN=%u",         run->lattice_group);
    options_length_common += sprintf_s(options_common + options_length_common,sizeof(options_common)-options_length_common," -D ND=%u",          run->lattice_nd);
    options_length_common += sprintf_s(options_common + options_length_common,sizeof(options_common)-options_length_common," -D N1EXACT=%u",     run->lattice_domain_size[0]);
    options_length_common += sprintf_s(options_common + options_length_common,sizeof(options_common)-options_length_common," -D N1=%u",          lattice_domain_n1);
    options_length_common += sprintf_s(options_common + options_length_common,sizeof(options_common)-options_length_common," -D N2=%u",          run->lattice_domain_size[1]);
    options_length_common += sprintf_s(options_common + options_length_common,sizeof(options_common)-options_length_common," -D N3=%u",          run->lattice_domain_size[2]);
    options_length_common += sprintf_s(options_common + options_length_common,sizeof(options_common)-options_length_common," -D N4=%u",          run->lattice_domain_size[3]);
    if (!((run->PHI==0.0)&&(run->OMEGA==0.0))) options_length_common += sprintf_s(options_common + options_length_common,sizeof(options_common)-options_length_common," -D TBC");     // turn on TBC
#if (MODEL_ON==1)
    options_length_common += sprintf_s(options_common + options_length_common,sizeof(options_common)-options_length_common," -I %s%s",           GPU0->cl_root_path,path_oncl);
    options_length_common += sprintf_s(options_common + options_length_common,sizeof(options_common)-options_length_common," -D SKGROUP=%u",     run->ON_SK_group);
#else
    options_length_common += sprintf_s(options_common + options_length_common,sizeof(options_common)-options_length_common," -I %s%s",           GPU0->cl_root_path,path_suncl);
#endif
    options_length_common += sprintf_s(options_common + options_length_common,sizeof(options_common)-options_length_common," -I %s%s",           GPU0->cl_root_path,path_kernel);
    options_length_common += sprintf_s(options_common + options_length_common,sizeof(options_common)-options_length_common," -D ROWSIZE=%u",     lattice_table_row_size);
    options_length_common += sprintf_s(options_common + options_length_common,sizeof(options_common)-options_length_common," -D PRECISION=%u",   run->precision);

    char options[1024];
    int options_length  = sprintf_s(options,sizeof(options),"%s",options_common);
        options_length += sprintf_s(options + options_length,sizeof(options)-options_length," -D NHIT=%u",        run->NHIT);
        options_length += sprintf_s(options + options_length,sizeof(options)-options_length," -D PRNGSTEP=%u",    lattice_table_row_size_half);

    char buffer_update_cl[FNAME_MAX_LENGTH];
        j = sprintf_s(buffer_update_cl  ,FNAME_MAX_LENGTH,  "%s",GPU0->cl_root_path);
        j+= sprintf_s(buffer_update_cl+j,FNAME_MAX_LENGTH-j,"%s",SOURCE_UPDATE);
    char* update_source       = GPU0->source_read(buffer_update_cl);
                                GPU0->program_create(update_source,options);

#if (MODEL_ON==1)
    // O(1)___________________________________________________________________________________
    const size_t init_global_size[]               = {lattice_table_exact_row_size};
    const size_t init_hot_global_size[]           = {lattice_table_exact_row_size};
    const size_t monte_global_size[]              = {lattice_table_exact_row_size_half};

    const size_t measurement_global_size[]        = {lattice_table_exact_row_size};
    const size_t clear_measurement_global_size[]  = {lattice_measurement_size_F};
    const size_t reduce_measurement_global_size[] = {GPU0->GPU_info.max_workgroup_size};
    const size_t reduce_local_size[3]             = {};

    const size_t local_size_lattice_measurement[] = {local_size_intel};

    if (run->ints==model_start_cold) { // cold init
                sun_init_id = GPU0->kernel_init("lattice_init_cold",1,init_global_size,NULL);
                argument_id = GPU0->kernel_init_buffer(sun_init_id,lattice_table);
    } else {                                  // hot init
                sun_init_id = GPU0->kernel_init("lattice_init_hot",1,init_hot_global_size,NULL);
                argument_id = GPU0->kernel_init_buffer(sun_init_id,lattice_table);
                argument_id = GPU0->kernel_init_buffer(sun_init_id,PRNG0->PRNG_randoms_id);
                argument_id = GPU0->kernel_init_buffer(sun_init_id,lattice_parameters);
    }

    int size_reduce_acceptance_rate_odd  = 0;
    int size_reduce_acceptance_rate_even = 0;

    sun_update_odd_X_id = GPU0->kernel_init("update_odd",1,monte_global_size,NULL);
            argument_id = GPU0->kernel_init_buffer(sun_update_odd_X_id,lattice_table);
            argument_id = GPU0->kernel_init_buffer(sun_update_odd_X_id,lattice_parameters);
            argument_id = GPU0->kernel_init_buffer(sun_update_odd_X_id,PRNG0->PRNG_randoms_id);
    if (run->get_acceptance_rate) {
            argument_id = GPU0->kernel_init_buffer(sun_update_odd_X_id,lattice_lds);
            argument_id = GPU0->kernel_init_buffer(sun_update_odd_X_id,lattice_measurement);
            size_reduce_acceptance_rate_odd  = (int) ceil((double) monte_global_size[0] / GPU0->kernel_get_worksize(sun_update_odd_X_id));
    }

    sun_update_even_X_id = GPU0->kernel_init("update_even",1,monte_global_size,NULL);
             argument_id = GPU0->kernel_init_buffer(sun_update_even_X_id,lattice_table);
             argument_id = GPU0->kernel_init_buffer(sun_update_even_X_id,lattice_parameters);
             argument_id = GPU0->kernel_init_buffer(sun_update_even_X_id,PRNG0->PRNG_randoms_id);
    if (run->get_acceptance_rate) {
            argument_id = GPU0->kernel_init_buffer(sun_update_even_X_id,lattice_lds);
            argument_id = GPU0->kernel_init_buffer(sun_update_even_X_id,lattice_measurement);
            size_reduce_acceptance_rate_even = (int) ceil((double) monte_global_size[0] / GPU0->kernel_get_worksize(sun_update_even_X_id));
    }

    if (run->get_acceptance_rate) {
        int lattice_measurement_size_correction = _MAX(size_reduce_acceptance_rate_odd,size_reduce_acceptance_rate_even);

        sun_reduce_acceptance_rate_id = GPU0->kernel_init("reduce_acceptance_rate",1,reduce_measurement_global_size,reduce_local_size);
                    argument_id = GPU0->kernel_init_buffer(sun_reduce_acceptance_rate_id,lattice_measurement);
                    argument_id = GPU0->kernel_init_buffer(sun_reduce_acceptance_rate_id,lattice_acceptance_rate);
                    argument_id = GPU0->kernel_init_buffer(sun_reduce_acceptance_rate_id,lattice_lds);
                    argument_acceptance_rate_index = GPU0->kernel_init_constant(sun_reduce_acceptance_rate_id,&lattice_measurement_size_correction);
    }



        // for all measurements _____________________________________________________________________________________________________________________________________
    char options_measurements[1024];
    int options_measurement_length  = sprintf_s(options_measurements,sizeof(options_measurements),"%s",options_common);

    char buffer_measurements_cl[FNAME_MAX_LENGTH];
        j = sprintf_s(buffer_measurements_cl  ,FNAME_MAX_LENGTH,  "%s",GPU0->cl_root_path);
        j+= sprintf_s(buffer_measurements_cl+j,FNAME_MAX_LENGTH-j,"%s",SOURCE_MEASUREMENTS);
    char* measurements_source = GPU0->source_read(buffer_measurements_cl);
                                GPU0->program_create(measurements_source,options_measurements);

    sun_clear_measurement_id = GPU0->kernel_init("clear_measurement",1,clear_measurement_global_size,NULL);
                 argument_id = GPU0->kernel_init_buffer(sun_clear_measurement_id,lattice_measurement);

    sun_measurement_id = GPU0->kernel_init("lattice_measurement",1,measurement_global_size,local_size_lattice_measurement);
           argument_id = GPU0->kernel_init_buffer(sun_measurement_id,lattice_table);
           argument_id = GPU0->kernel_init_buffer(sun_measurement_id,lattice_measurement);
           argument_id = GPU0->kernel_init_buffer(sun_measurement_id,lattice_parameters);	
           argument_id = GPU0->kernel_init_buffer(sun_measurement_id,lattice_lds);
    int size_reduce_measurement_double2 = (int) ceil((double) lattice_table_exact_row_size / GPU0->kernel_get_worksize(sun_measurement_id));

    sun_measurement_reduce_id = GPU0->kernel_init("reduce_measurement_double2",1,reduce_measurement_global_size,reduce_local_size);
                  argument_id = GPU0->kernel_init_buffer(sun_measurement_reduce_id,lattice_measurement);
                  argument_id = GPU0->kernel_init_buffer(sun_measurement_reduce_id,lattice_energies);
                  argument_id = GPU0->kernel_init_buffer(sun_measurement_reduce_id,lattice_lds);
                  argument_measurement_index = GPU0->kernel_init_constant(sun_measurement_reduce_id,&size_reduce_measurement_double2);
    int size_reduce_measurement_plq_double2   = 0;
    if (run->get_plaquettes_avr) {
        sun_measurement_plq_id = GPU0->kernel_init("lattice_measurement_plq",1,measurement_global_size,local_size_lattice_measurement);
                   argument_id = GPU0->kernel_init_buffer(sun_measurement_plq_id,lattice_table);
                   argument_id = GPU0->kernel_init_buffer(sun_measurement_plq_id,lattice_measurement);
                   argument_id = GPU0->kernel_init_buffer(sun_measurement_plq_id,lattice_parameters);
                   argument_id = GPU0->kernel_init_buffer(sun_measurement_plq_id,lattice_lds);
        size_reduce_measurement_plq_double2 = (int) ceil((double) lattice_table_exact_row_size / GPU0->kernel_get_worksize(sun_measurement_plq_id));

        sun_measurement_plq_reduce_id = GPU0->kernel_init("reduce_measurement_plq_double2",1,reduce_measurement_global_size,reduce_local_size);
                          argument_id = GPU0->kernel_init_buffer(sun_measurement_plq_reduce_id,lattice_measurement);
                          argument_id = GPU0->kernel_init_buffer(sun_measurement_plq_reduce_id,lattice_energies_plq);
                          argument_id = GPU0->kernel_init_buffer(sun_measurement_plq_reduce_id,lattice_lds);
                          argument_plq_index = GPU0->kernel_init_constant(sun_measurement_plq_reduce_id,&size_reduce_measurement_plq_double2);
    }

    int size_reduce_measurement_corr_double2   = 0;
    cl_uint4 correlator_stepz;
        correlator_stepz.s[0] = run->correlator_X;
        correlator_stepz.s[1] = run->correlator_Y;
        correlator_stepz.s[2] = run->correlator_Z;
        correlator_stepz.s[3] = run->correlator_T;
    if (run->get_correlators) {
        sun_measurement_corr_id = GPU0->kernel_init("lattice_measurement_correlator",1,measurement_global_size,local_size_lattice_measurement);
                    argument_id = GPU0->kernel_init_buffer(sun_measurement_corr_id,lattice_table);
                    argument_id = GPU0->kernel_init_buffer(sun_measurement_corr_id,lattice_measurement);
                    argument_id = GPU0->kernel_init_buffer(sun_measurement_corr_id,lattice_parameters);
                    argument_id = GPU0->kernel_init_buffer(sun_measurement_corr_id,lattice_lds);
                    argument_id = GPU0->kernel_init_constant(sun_measurement_corr_id,&correlator_stepz);
        size_reduce_measurement_corr_double2 = (int) ceil((double) lattice_table_exact_row_size / GPU0->kernel_get_worksize(sun_measurement_corr_id));

        sun_measurement_corr_reduce_id = GPU0->kernel_init("reduce_measurement_plq_double2",1,reduce_measurement_global_size,reduce_local_size);
                           argument_id = GPU0->kernel_init_buffer(sun_measurement_corr_reduce_id,lattice_measurement);
                           argument_id = GPU0->kernel_init_buffer(sun_measurement_corr_reduce_id,lattice_correlators);
                           argument_id = GPU0->kernel_init_buffer(sun_measurement_corr_reduce_id,lattice_lds);
                           argument_correlators_index = GPU0->kernel_init_constant(sun_measurement_corr_reduce_id,&size_reduce_measurement_corr_double2);
    }


#else
    // SU(3)__________________________________________________________________________________
    const size_t init_global_size[]               = {lattice_table_exact_row_size};
    const size_t init_hot_global_size[]           = {3*lattice_table_row_size};
    const size_t monte_global_size[]              = {lattice_table_exact_row_size_half};
    const size_t boundary_global_size[]           = {lattice_boundary_size};

    const size_t measurement3_global_size[]       = {lattice_table_exact_row_size};
    const size_t polyakov3_global_size[]          = {lattice_polyakov_size};
    const size_t clear_measurement_global_size[]  = {lattice_measurement_size_F};

    const size_t reduce_measurement_global_size[] = {GPU0->GPU_info.max_workgroup_size};
    const size_t reduce_polyakov_global_size[]    = {GPU0->GPU_info.max_workgroup_size};
    const size_t reduce_local_size[3]             = {};

    const size_t local_size_lattice_measurement[] = {local_size_intel};
    const size_t local_size_lattice_polyakov[]    = {local_size_intel};
    const size_t local_size_lattice_wilson[]      = {local_size_intel};


    if (run->ints==model_start_gid) {         // gid init
                sun_init_id = GPU0->kernel_init("lattice_init_gid",1,init_global_size,NULL);
                argument_id = GPU0->kernel_init_buffer(sun_init_id,lattice_table);
    } else if (run->ints==model_start_cold) { // cold init
                sun_init_id = GPU0->kernel_init("lattice_init_cold",1,init_global_size,NULL);
                argument_id = GPU0->kernel_init_buffer(sun_init_id,lattice_table);
    } else {                                  // hot init
                sun_init_X_id = GPU0->kernel_init("lattice_init_hot_X",1,init_hot_global_size,NULL);
                argument_id = GPU0->kernel_init_buffer(sun_init_X_id,lattice_table);
                argument_id = GPU0->kernel_init_buffer(sun_init_X_id,PRNG0->PRNG_randoms_id);
                sun_init_Y_id = GPU0->kernel_init("lattice_init_hot_Y",1,init_hot_global_size,NULL);
                argument_id = GPU0->kernel_init_buffer(sun_init_Y_id,lattice_table);
                argument_id = GPU0->kernel_init_buffer(sun_init_Y_id,PRNG0->PRNG_randoms_id);
                sun_init_Z_id = GPU0->kernel_init("lattice_init_hot_Z",1,init_hot_global_size,NULL);
                argument_id = GPU0->kernel_init_buffer(sun_init_Z_id,lattice_table);
                argument_id = GPU0->kernel_init_buffer(sun_init_Z_id,PRNG0->PRNG_randoms_id);
                sun_init_T_id = GPU0->kernel_init("lattice_init_hot_T",1,init_hot_global_size,NULL);
                argument_id = GPU0->kernel_init_buffer(sun_init_T_id,lattice_table);
                argument_id = GPU0->kernel_init_buffer(sun_init_T_id,PRNG0->PRNG_randoms_id);
    }

    sun_GramSchmidt_id = GPU0->kernel_init("lattice_GramSchmidt",1,init_global_size,NULL);
           argument_id = GPU0->kernel_init_buffer(sun_GramSchmidt_id,lattice_table);
           argument_id = GPU0->kernel_init_buffer(sun_GramSchmidt_id,lattice_parameters);

    int size_reduce_acceptance_rate_odd_X  = 0;
    int size_reduce_acceptance_rate_odd_Y  = 0;
    int size_reduce_acceptance_rate_odd_Z  = 0;
    int size_reduce_acceptance_rate_odd_T  = 0;

    int size_reduce_acceptance_rate_even_X = 0;
    int size_reduce_acceptance_rate_even_Y = 0;
    int size_reduce_acceptance_rate_even_Z = 0;
    int size_reduce_acceptance_rate_even_T = 0;


    sun_update_odd_X_id = GPU0->kernel_init("update_odd_X",1,monte_global_size,NULL);
            argument_id = GPU0->kernel_init_buffer(sun_update_odd_X_id,lattice_table);
            argument_id = GPU0->kernel_init_buffer(sun_update_odd_X_id,lattice_parameters);
            argument_id = GPU0->kernel_init_buffer(sun_update_odd_X_id,PRNG0->PRNG_randoms_id);
            if (run->get_acceptance_rate) {
                    argument_id = GPU0->kernel_init_buffer(sun_update_odd_X_id,lattice_lds);
                    argument_id = GPU0->kernel_init_buffer(sun_update_odd_X_id,lattice_measurement);
                    size_reduce_acceptance_rate_odd_X  = (int) ceil((double) monte_global_size[0] / GPU0->kernel_get_worksize(sun_update_odd_X_id));
            }

    sun_update_even_X_id = GPU0->kernel_init("update_even_X",1,monte_global_size,NULL);
             argument_id = GPU0->kernel_init_buffer(sun_update_even_X_id,lattice_table);
             argument_id = GPU0->kernel_init_buffer(sun_update_even_X_id,lattice_parameters);
             argument_id = GPU0->kernel_init_buffer(sun_update_even_X_id,PRNG0->PRNG_randoms_id);
            if (run->get_acceptance_rate) {
                    argument_id = GPU0->kernel_init_buffer(sun_update_even_X_id,lattice_lds);
                    argument_id = GPU0->kernel_init_buffer(sun_update_even_X_id,lattice_measurement);
                    size_reduce_acceptance_rate_even_X = (int) ceil((double) monte_global_size[0] / GPU0->kernel_get_worksize(sun_update_even_X_id));
            }

    sun_update_odd_Y_id = GPU0->kernel_init("update_odd_Y",1,monte_global_size,NULL);
            argument_id = GPU0->kernel_init_buffer(sun_update_odd_Y_id,lattice_table);
            argument_id = GPU0->kernel_init_buffer(sun_update_odd_Y_id,lattice_parameters);
            argument_id = GPU0->kernel_init_buffer(sun_update_odd_Y_id,PRNG0->PRNG_randoms_id);
            if (run->get_acceptance_rate) {
                    argument_id = GPU0->kernel_init_buffer(sun_update_odd_Y_id,lattice_lds);
                    argument_id = GPU0->kernel_init_buffer(sun_update_odd_Y_id,lattice_measurement);
                    size_reduce_acceptance_rate_odd_Y  = (int) ceil((double) monte_global_size[0] / GPU0->kernel_get_worksize(sun_update_odd_Y_id));
            }

    sun_update_even_Y_id = GPU0->kernel_init("update_even_Y",1,monte_global_size,NULL);
             argument_id = GPU0->kernel_init_buffer(sun_update_even_Y_id,lattice_table);
             argument_id = GPU0->kernel_init_buffer(sun_update_even_Y_id,lattice_parameters);
             argument_id = GPU0->kernel_init_buffer(sun_update_even_Y_id,PRNG0->PRNG_randoms_id);
            if (run->get_acceptance_rate) {
                    argument_id = GPU0->kernel_init_buffer(sun_update_even_Y_id,lattice_lds);
                    argument_id = GPU0->kernel_init_buffer(sun_update_even_Y_id,lattice_measurement);
                    size_reduce_acceptance_rate_even_Y = (int) ceil((double) monte_global_size[0] / GPU0->kernel_get_worksize(sun_update_even_Y_id));
            }

    sun_update_odd_Z_id = GPU0->kernel_init("update_odd_Z",1,monte_global_size,NULL);
            argument_id = GPU0->kernel_init_buffer(sun_update_odd_Z_id,lattice_table);
            argument_id = GPU0->kernel_init_buffer(sun_update_odd_Z_id,lattice_parameters);
            argument_id = GPU0->kernel_init_buffer(sun_update_odd_Z_id,PRNG0->PRNG_randoms_id);
            if (run->get_acceptance_rate) {
                    argument_id = GPU0->kernel_init_buffer(sun_update_odd_Z_id,lattice_lds);
                    argument_id = GPU0->kernel_init_buffer(sun_update_odd_Z_id,lattice_measurement);
                    size_reduce_acceptance_rate_odd_Z  = (int) ceil((double) monte_global_size[0] / GPU0->kernel_get_worksize(sun_update_odd_Z_id));
            }

    sun_update_even_Z_id = GPU0->kernel_init("update_even_Z",1,monte_global_size,NULL);
             argument_id = GPU0->kernel_init_buffer(sun_update_even_Z_id,lattice_table);
             argument_id = GPU0->kernel_init_buffer(sun_update_even_Z_id,lattice_parameters);
             argument_id = GPU0->kernel_init_buffer(sun_update_even_Z_id,PRNG0->PRNG_randoms_id);
            if (run->get_acceptance_rate) {
                    argument_id = GPU0->kernel_init_buffer(sun_update_even_Z_id,lattice_lds);
                    argument_id = GPU0->kernel_init_buffer(sun_update_even_Z_id,lattice_measurement);
                    size_reduce_acceptance_rate_even_Z = (int) ceil((double) monte_global_size[0] / GPU0->kernel_get_worksize(sun_update_even_Z_id));
            }

    sun_update_odd_T_id = GPU0->kernel_init("update_odd_T",1,monte_global_size,NULL);
            argument_id = GPU0->kernel_init_buffer(sun_update_odd_T_id,lattice_table);
            argument_id = GPU0->kernel_init_buffer(sun_update_odd_T_id,lattice_parameters);
            argument_id = GPU0->kernel_init_buffer(sun_update_odd_T_id,PRNG0->PRNG_randoms_id);
            if (run->get_acceptance_rate) {
                    argument_id = GPU0->kernel_init_buffer(sun_update_odd_T_id,lattice_lds);
                    argument_id = GPU0->kernel_init_buffer(sun_update_odd_T_id,lattice_measurement);
                    size_reduce_acceptance_rate_odd_T  = (int) ceil((double) monte_global_size[0] / GPU0->kernel_get_worksize(sun_update_odd_T_id));
            }

    sun_update_even_T_id = GPU0->kernel_init("update_even_T",1,monte_global_size,NULL);
             argument_id = GPU0->kernel_init_buffer(sun_update_even_T_id,lattice_table);
             argument_id = GPU0->kernel_init_buffer(sun_update_even_T_id,lattice_parameters);
             argument_id = GPU0->kernel_init_buffer(sun_update_even_T_id,PRNG0->PRNG_randoms_id);
            if (run->get_acceptance_rate) {
                    argument_id = GPU0->kernel_init_buffer(sun_update_even_T_id,lattice_lds);
                    argument_id = GPU0->kernel_init_buffer(sun_update_even_T_id,lattice_measurement);
                    size_reduce_acceptance_rate_even_T = (int) ceil((double) monte_global_size[0] / GPU0->kernel_get_worksize(sun_update_even_T_id));
            }

    if (run->get_acceptance_rate) {
        int lattice_measurement_size_correction = _MAX(_MAX(_MAX(_MAX(_MAX(_MAX(_MAX(
                                                    size_reduce_acceptance_rate_odd_X ,size_reduce_acceptance_rate_even_X),
                                                    size_reduce_acceptance_rate_odd_Y),size_reduce_acceptance_rate_even_Y),
                                                    size_reduce_acceptance_rate_odd_Z),size_reduce_acceptance_rate_even_Z),
                                                    size_reduce_acceptance_rate_odd_T),size_reduce_acceptance_rate_even_T);

        sun_reduce_acceptance_rate_id = GPU0->kernel_init("reduce_acceptance_rate",1,reduce_measurement_global_size,reduce_local_size);
                    argument_id = GPU0->kernel_init_buffer(sun_reduce_acceptance_rate_id,lattice_measurement);
                    argument_id = GPU0->kernel_init_buffer(sun_reduce_acceptance_rate_id,lattice_acceptance_rate);
                    argument_id = GPU0->kernel_init_buffer(sun_reduce_acceptance_rate_id,lattice_lds);
                    argument_acceptance_rate_index = GPU0->kernel_init_constant(sun_reduce_acceptance_rate_id,&lattice_measurement_size_correction);
    }



    // boundary exchange for big lattices _______________________________________________________________________________________________________________________
    if (!run->turnoff_boundary_extraction) {
        char buffer_boundary_cl[FNAME_MAX_LENGTH];
            j = sprintf_s(buffer_boundary_cl  ,FNAME_MAX_LENGTH,  "%s",GPU0->cl_root_path);
            j+= sprintf_s(buffer_boundary_cl+j,FNAME_MAX_LENGTH-j,"%s",SOURCE_BOUNDARY);
        char* boundary_source = GPU0->source_read(buffer_boundary_cl);
                                GPU0->program_create(boundary_source,options_common);

        sun_get_boundary_low_id = GPU0->kernel_init("lattice_get_boundary_low",1,boundary_global_size,NULL);
                    argument_id = GPU0->kernel_init_buffer(sun_get_boundary_low_id,lattice_table);
                    argument_id = GPU0->kernel_init_buffer(sun_get_boundary_low_id,lattice_boundary);
0
        sun_put_boundary_low_id = GPU0->kernel_init("lattice_put_boundary_low",1,boundary_global_size,NULL);
                    argument_id = GPU0->kernel_init_buffer(sun_put_boundary_low_id,lattice_table);
                    argument_id = GPU0->kernel_init_buffer(sun_put_boundary_low_id,lattice_boundary);
    }


    // for all measurements _____________________________________________________________________________________________________________________________________
    char options_measurements[1024];
    int options_measurement_length  = sprintf_s(options_measurements,sizeof(options_measurements),"%s",options_common);

    if (run->get_Fmunu) options_measurement_length += sprintf_s(options_measurements + options_measurement_length,sizeof(options_measurements)-options_measurement_length," -D FMUNU");   // calculate tensor Fmunu for H field
    if (run->get_F0mu)  options_measurement_length += sprintf_s(options_measurements + options_measurement_length,sizeof(options_measurements)-options_measurement_length," -D F0MU");   // calculate tensor Fmunu for E field

    if ((run->get_Fmunu)||(run->get_F0mu)) {
        if (run->Fmunu_index1 == 1) options_measurement_length += sprintf_s(options_measurements + options_measurement_length,sizeof(options_measurements)-options_measurement_length," -D FMUNU1");   // calculate tensor Fmunu for lambda1 matrix
        if (run->Fmunu_index1 == 2) options_measurement_length += sprintf_s(options_measurements + options_measurement_length,sizeof(options_measurements)-options_measurement_length," -D FMUNU2");   // calculate tensor Fmunu for lambda2 matrix
        if (run->Fmunu_index1 == 4) options_measurement_length += sprintf_s(options_measurements + options_measurement_length,sizeof(options_measurements)-options_measurement_length," -D FMUNU4");   // calculate tensor Fmunu for lambda4 matrix
        if (run->Fmunu_index2 == 5) options_measurement_length += sprintf_s(options_measurements + options_measurement_length,sizeof(options_measurements)-options_measurement_length," -D FMUNU5");   // calculate tensor Fmunu for lambda5 matrix
        if (run->Fmunu_index2 == 6) options_measurement_length += sprintf_s(options_measurements + options_measurement_length,sizeof(options_measurements)-options_measurement_length," -D FMUNU6");   // calculate tensor Fmunu for lambda6 matrix
        if (run->Fmunu_index2 == 7) options_measurement_length += sprintf_s(options_measurements + options_measurement_length,sizeof(options_measurements)-options_measurement_length," -D FMUNU7");   // calculate tensor Fmunu for lambda7 matrix
    }
    char buffer_measurements_cl[FNAME_MAX_LENGTH];
        j = sprintf_s(buffer_measurements_cl  ,FNAME_MAX_LENGTH,  "%s",GPU0->cl_root_path);
        j+= sprintf_s(buffer_measurements_cl+j,FNAME_MAX_LENGTH-j,"%s",SOURCE_MEASUREMENTS);
    char* measurements_source       = GPU0->source_read(buffer_measurements_cl);
                                GPU0->program_create(measurements_source,options_measurements);

    int size_reduce_measurement_plq_double2   = 0;
    int offset_reduce_measurement_plq_double2 = 0;

    sun_measurement_plq_id        = 0;
    sun_measurement_reduce_id     = 0;
    sun_measurement_plq_reduce_id = 0;

    sun_clear_measurement_id = GPU0->kernel_init("clear_measurement",1,clear_measurement_global_size,NULL);
                 argument_id = GPU0->kernel_init_buffer(sun_clear_measurement_id,lattice_measurement);

    sun_measurement_id = GPU0->kernel_init("lattice_measurement",1,measurement3_global_size,local_size_lattice_measurement);
           argument_id = GPU0->kernel_init_buffer(sun_measurement_id,lattice_table);
           argument_id = GPU0->kernel_init_buffer(sun_measurement_id,lattice_measurement);
           argument_id = GPU0->kernel_init_buffer(sun_measurement_id,lattice_parameters);	
           argument_id = GPU0->kernel_init_buffer(sun_measurement_id,lattice_lds);
    int size_reduce_measurement_double2 = (int) ceil((double) lattice_table_exact_row_size / GPU0->kernel_get_worksize(sun_measurement_id));

    // for mean averaged plaquettes measurements
    cl_uint4 mesurement_plq_param;
    if ((run->get_plaquettes_avr)||(run->get_Fmunu)||(run->get_F0mu)) {
        sun_measurement_plq_id = GPU0->kernel_init("lattice_measurement_plq",1,measurement3_global_size,local_size_lattice_measurement);
                   offset_reduce_measurement_plq_double2 = GPU0->buffer_size_align((unsigned int) ceil((double) lattice_table_exact_row_size / GPU0->kernel_get_worksize(sun_measurement_plq_id)),GPU0->kernel_get_worksize(sun_measurement_plq_id));
                   argument_id = GPU0->kernel_init_buffer(sun_measurement_plq_id,lattice_table);
                   argument_id = GPU0->kernel_init_buffer(sun_measurement_plq_id,lattice_measurement);
                   argument_id = GPU0->kernel_init_buffer(sun_measurement_plq_id,lattice_parameters);
                   argument_id = GPU0->kernel_init_buffer(sun_measurement_plq_id,lattice_lds);
                   argument_id = GPU0->kernel_init_constant(sun_measurement_plq_id,&offset_reduce_measurement_plq_double2);
        size_reduce_measurement_plq_double2 = (int) ceil((double) lattice_table_exact_row_size / GPU0->kernel_get_worksize(sun_measurement_plq_id));

        sun_measurement_plq_reduce_id = GPU0->kernel_init("reduce_measurement_plq_double2",1,reduce_measurement_global_size,reduce_local_size);
                          argument_id = GPU0->kernel_init_buffer(sun_measurement_plq_reduce_id,lattice_measurement);
                          argument_id = GPU0->kernel_init_buffer(sun_measurement_plq_reduce_id,lattice_energies_plq);
                          argument_id = GPU0->kernel_init_buffer(sun_measurement_plq_reduce_id,lattice_lds);
            mesurement_plq_param.s[0] = size_reduce_measurement_plq_double2;
            mesurement_plq_param.s[1] = offset_reduce_measurement_plq_double2;
            mesurement_plq_param.s[2] = lattice_energies_offset;
            mesurement_plq_param.s[3] = 0;
                          argument_plq_index = GPU0->kernel_init_constant(sun_measurement_plq_reduce_id,&mesurement_plq_param);

    int lattice_measurement_size_correction = _MAX(size_reduce_measurement_double2,size_reduce_measurement_plq_double2);
    if (lattice_measurement_size_F<(unsigned int) lattice_measurement_size_correction){

        printf ("buffer plattice_measurement should be resized!!!\n");
        _getch();
    }

    sun_measurement_reduce_id = GPU0->kernel_init("reduce_measurement_double2",1,reduce_measurement_global_size,reduce_local_size);
                  argument_id = GPU0->kernel_init_buffer(sun_measurement_reduce_id,lattice_measurement);
                  argument_id = GPU0->kernel_init_buffer(sun_measurement_reduce_id,lattice_energies);
                  argument_id = GPU0->kernel_init_buffer(sun_measurement_reduce_id,lattice_lds);
                  argument_measurement_index = GPU0->kernel_init_constant(sun_measurement_reduce_id,&size_reduce_measurement_double2);

    // for Wilson loop measurements _____________________________________________________________________________________________________________________________
    char options_wilson[1024];
    sprintf_s(options_wilson,sizeof(options_wilson),"%s",options_common);

    char buffer_wilson_cl[FNAME_MAX_LENGTH];
        j = sprintf_s(buffer_wilson_cl  ,FNAME_MAX_LENGTH,  "%s",GPU0->cl_root_path);
        j+= sprintf_s(buffer_wilson_cl+j,FNAME_MAX_LENGTH-j,"%s",SOURCE_WILSON_LOOP);
    char* wilson_source       = GPU0->source_read(buffer_wilson_cl);
                                GPU0->program_create(wilson_source,options_wilson);

    int size_reduce_wilson_double2 = 0;

    sun_measurement_wilson_id  = 0;
    sun_wilson_loop_reduce_id  = 0;
    if (run->get_wilson_loop) {
        sun_measurement_wilson_id = GPU0->kernel_init("lattice_measurement_wilson",1,measurement3_global_size,local_size_lattice_wilson);
                      argument_id = GPU0->kernel_init_buffer(sun_measurement_wilson_id,lattice_table);
                      argument_id = GPU0->kernel_init_buffer(sun_measurement_wilson_id,lattice_measurement);
                      argument_id = GPU0->kernel_init_buffer(sun_measurement_wilson_id,lattice_parameters);
                      argument_id = GPU0->kernel_init_buffer(sun_measurement_wilson_id,lattice_lds);
        size_reduce_wilson_double2 = (int) ceil((double) lattice_domain_exact_site / GPU0->kernel_get_worksize(sun_measurement_wilson_id));

        sun_wilson_loop_reduce_id = GPU0->kernel_init("reduce_wilson_double2",1,reduce_measurement_global_size,reduce_local_size);
                      argument_id = GPU0->kernel_init_buffer(sun_wilson_loop_reduce_id,lattice_measurement);
                      argument_id = GPU0->kernel_init_buffer(sun_wilson_loop_reduce_id,lattice_wilson_loop);
                      argument_id = GPU0->kernel_init_buffer(sun_wilson_loop_reduce_id,lattice_lds);
                      argument_wilson_index = GPU0->kernel_init_constant(sun_wilson_loop_reduce_id,&size_reduce_wilson_double2);
                      // setup index for wilson loop is before kernel run
    }

    // for Polyakov loop measurements ___________________________________________________________________________________________________________________________
    char options_polyakov[1024];
    int options_length_polyakov  = sprintf_s(options_polyakov,sizeof(options_polyakov),"%s",options_common);
        options_length_polyakov += sprintf_s(options_polyakov + options_length_polyakov,sizeof(options_polyakov)-options_length_polyakov," -D PL=%u",          run->PL_level);

    char buffer_polyakov_cl[FNAME_MAX_LENGTH];
        j = sprintf_s(buffer_polyakov_cl  ,FNAME_MAX_LENGTH,  "%s",GPU0->cl_root_path);
        j+= sprintf_s(buffer_polyakov_cl+j,FNAME_MAX_LENGTH-j,"%s",SOURCE_POLYAKOV);
    char* polyakov_source       = GPU0->source_read(buffer_polyakov_cl);
                                  GPU0->program_create(polyakov_source,options_polyakov);
    int size_reduce_polyakov_double2   = 0;
    int offset_reduce_polyakov_double2 = 0;
    sun_polyakov_id        = 0;
    sun_polyakov_reduce_id = 0;
    cl_uint4 polyakov_param;
    if (run->PL_level > 0) {
        sun_polyakov_id = GPU0->kernel_init("lattice_polyakov",1,polyakov3_global_size,local_size_lattice_polyakov);
        offset_reduce_polyakov_double2 = GPU0->buffer_size_align((unsigned int) ceil((double) lattice_polyakov_size / GPU0->kernel_get_worksize(sun_polyakov_id)),GPU0->kernel_get_worksize(sun_polyakov_id));
            argument_id = GPU0->kernel_init_buffer(sun_polyakov_id,lattice_table);
            argument_id = GPU0->kernel_init_buffer(sun_polyakov_id,lattice_measurement);
            argument_id = GPU0->kernel_init_buffer(sun_polyakov_id,lattice_parameters);
            argument_id = GPU0->kernel_init_buffer(sun_polyakov_id,lattice_lds);
            argument_id = GPU0->kernel_init_constant(sun_polyakov_id,&offset_reduce_polyakov_double2);
        size_reduce_polyakov_double2 = (int) ceil((double) lattice_polyakov_size / GPU0->kernel_get_worksize(sun_polyakov_id));

        sun_polyakov_reduce_id = GPU0->kernel_init("reduce_polyakov_double2",1,reduce_polyakov_global_size,reduce_local_size);
            argument_id = GPU0->kernel_init_buffer(sun_polyakov_reduce_id,lattice_measurement);
            argument_id = GPU0->kernel_init_buffer(sun_polyakov_reduce_id,lattice_polyakov_loop);
            argument_id = GPU0->kernel_init_buffer(sun_polyakov_reduce_id,lattice_lds);
            polyakov_param.s[0] = size_reduce_polyakov_double2;
            polyakov_param.s[1] = offset_reduce_polyakov_double2;
            polyakov_param.s[2] = lattice_polyakov_loop_size;
            polyakov_param.s[3] = 0;
            argument_polyakov_index = GPU0->kernel_init_constant(sun_polyakov_reduce_id,&polyakov_param);
    }
#endif
}
void        model::lattice_create_buffers(void){
    cl_float*    plattice_parameters_float  = NULL;
    cl_double*   plattice_parameters_double = NULL;

    int fc2 = ((run->PL_level>1)&&(lattice_measurement_size_F <= 2 * lattice_polyakov_loop_size)) ? 2 : 1;
    size_lattice_table           = lattice_table_size;
    size_lattice_measurement     = lattice_measurement_size_F;
    size_lattice_energies        = lattice_energies_size;
    size_lattice_wilson_loop     = lattice_energies_size;
    size_lattice_acceptance_rate = lattice_energies_size;
    size_lattice_correlators     = lattice_energies_size;
#if (MODEL_ON==1)
    size_lattice_energies_plq    = lattice_energies_size;
#else
    size_lattice_energies_plq    = lattice_energies_offset * MODEL_energies_size;
#endif
    size_lattice_polyakov_loop   = fc2 * lattice_polyakov_loop_size * run->PL_level;
    size_lattice_boundary        = lattice_boundary_size;
    size_lattice_parameters      = lattice_parameters_size;

    plattice_table_float       = NULL;
    plattice_table_float_1     = NULL;
    plattice_table_double      = NULL;
    plattice_table_double_1    = NULL;
    plattice_boundary_float    = NULL;
    plattice_boundary_float_1  = NULL;
    plattice_boundary_double   = NULL;
    plattice_boundary_double_1 = NULL;

    plattice_measurement     = (cl_double2*) calloc(size_lattice_measurement,    sizeof(cl_double2));

    plattice_energies        = NULL;
        if (run->get_actions_avr) plattice_energies     = (cl_double2*) calloc(size_lattice_energies,       sizeof(cl_double2));
    plattice_energies_plq    = NULL;
        if ((run->get_plaquettes_avr) || (run->get_Fmunu) || (run->get_F0mu)) plattice_energies_plq  = (cl_double2*) calloc(size_lattice_energies_plq,   sizeof(cl_double2));
    plattice_wilson_loop     = NULL;
        if (run->get_wilson_loop) plattice_wilson_loop = (cl_double*)  calloc(size_lattice_wilson_loop,    sizeof(cl_double));
    plattice_polyakov_loop   = NULL;
        if (run->PL_level > 0) plattice_polyakov_loop  = (cl_double2*) calloc(size_lattice_polyakov_loop,sizeof(cl_double2));
    plattice_acceptance_rate = NULL;
        if (run->get_actions_avr) plattice_acceptance_rate = (cl_double2*) calloc(size_lattice_acceptance_rate,sizeof(cl_double2));
    plattice_correlators     = NULL;
        if (run->get_correlators) plattice_correlators = (cl_double2*) calloc(size_lattice_acceptance_rate,sizeof(cl_double2));

    if (run->precision == model_precision_single) {
#if (MODEL_ON==1)
        // O(1)___________________________________________________________________________________
        plattice_table_float_1          = (cl_float*)  calloc(size_lattice_table,     sizeof(cl_float));
        if (!run->turnoff_boundary_extraction) 
            plattice_boundary_float_1   = (cl_float*)  calloc(size_lattice_boundary,  sizeof(cl_float));
#else
        // SU(3)__________________________________________________________________________________
        plattice_table_float            = (cl_float4*)  calloc(size_lattice_table,     sizeof(cl_float4));
        if (!run->turnoff_boundary_extraction) 
            plattice_boundary_float     = (cl_float4*)  calloc(size_lattice_boundary,  sizeof(cl_float4));
#endif
        plattice_parameters_float       = (cl_float*)   calloc(size_lattice_parameters,sizeof(cl_float));
        plattice_parameters_float[0]    = (float) (run->BETA / run->lattice_group);
        plattice_parameters_float[1]    = (float) (run->PHI);
        plattice_parameters_float[2]    = (float) (run->OMEGA);
        plattice_parameters_float[3]    = (float) (run->wilson_R);
        plattice_parameters_float[4]    = (float) (run->wilson_T);

        plattice_parameters_float[16]   = (float) (run->ON_z);
        plattice_parameters_float[17]   = (float) (run->ON_lambda);
        plattice_parameters_float[18]   = (float) (run->ON_zeta);
        plattice_parameters_float[19]   = (float) (run->ON_eta);
        plattice_parameters_float[20]   = (float) (run->ON_b);

        plattice_parameters_float[21]   = (float) (run->ON_sqrt_z_lambda_zeta);  // 1/4*sqrt(z/(lambda*zeta))
        plattice_parameters_float[22]   = (float) (run->ON_max_U);

    } else {
#if (MODEL_ON==1)
        // O(1)___________________________________________________________________________________
        plattice_table_double_1         = (cl_double*) calloc(size_lattice_table,     sizeof(cl_double));
        if (!run->turnoff_boundary_extraction) 
            plattice_boundary_double_1  = (cl_double*) calloc(size_lattice_boundary,  sizeof(cl_double));
#else
        // SU(3)__________________________________________________________________________________
        plattice_table_double           = (cl_double4*) calloc(size_lattice_table,     sizeof(cl_double4));
        if (!run->turnoff_boundary_extraction) 
            plattice_boundary_double    = (cl_double4*) calloc(size_lattice_boundary,  sizeof(cl_double4));
#endif
        plattice_parameters_double      = (cl_double*)  calloc(size_lattice_parameters,sizeof(cl_double));
        plattice_parameters_double[0]   = (double) (run->BETA / run->lattice_group);
        plattice_parameters_double[1]   = (double) (run->PHI);
        plattice_parameters_double[2]   = (double) (run->OMEGA);
        plattice_parameters_double[3]   = (double) (run->wilson_R);
        plattice_parameters_double[4]   = (double) (run->wilson_T);

        plattice_parameters_double[16]  = (double) (run->ON_z);
        plattice_parameters_double[17]  = (double) (run->ON_lambda);
        plattice_parameters_double[18]  = (double) (run->ON_zeta);
        plattice_parameters_double[19]  = (double) (run->ON_eta);
        plattice_parameters_double[20]  = (double) (run->ON_b);

        plattice_parameters_double[21]  = (double) (run->ON_sqrt_z_lambda_zeta);  // 1/4*sqrt(z/(lambda*zeta))
        plattice_parameters_double[22]  = (double) (run->ON_max_U);

    }

    if (run->INIT==0) lattice_load_state();  // load state file if needed

    int lds_size = (int) GPU0->GPU_info.local_memory_size / 4 / 4 / 2;  // 4 bytes per component, 4 components, double2 type <- use all local memory

    lattice_table = 0;
    if (run->precision == model_precision_single) {
#if (MODEL_ON==1)
        // O(1)___________________________________________________________________________________
        lattice_table           = GPU0->buffer_init(GPU0->buffer_type_IO, size_lattice_table,            plattice_table_float_1,     sizeof(cl_float));  // Lattice data
        if (!run->turnoff_boundary_extraction) 
            lattice_boundary    = GPU0->buffer_init(GPU0->buffer_type_IO, size_lattice_boundary,         plattice_boundary_float_1,  sizeof(cl_float));  // Lattice boundary
#else
        // SU(3)__________________________________________________________________________________
        lattice_table           = GPU0->buffer_init(GPU0->buffer_type_IO, size_lattice_table,            plattice_table_float,       sizeof(cl_float4));  // Lattice data
        if (!run->turnoff_boundary_extraction) 
            lattice_boundary    = GPU0->buffer_init(GPU0->buffer_type_IO, size_lattice_boundary,         plattice_boundary_float,    sizeof(cl_float4));  // Lattice boundary
#endif
        lattice_parameters      = GPU0->buffer_init(GPU0->buffer_type_Constant, size_lattice_parameters, plattice_parameters_float,  sizeof(cl_float));   // Lattice counters and indices
    } else {
#if (MODEL_ON==1)
        // O(1)___________________________________________________________________________________
        lattice_table           = GPU0->buffer_init(GPU0->buffer_type_IO, size_lattice_table,            plattice_table_double_1,    sizeof(cl_double));  // Lattice data
        if (!run->turnoff_boundary_extraction) 
            lattice_boundary    = GPU0->buffer_init(GPU0->buffer_type_IO, size_lattice_boundary,         plattice_boundary_double_1, sizeof(cl_double));  // Lattice boundary
#else
        // SU(3)__________________________________________________________________________________
        lattice_table           = GPU0->buffer_init(GPU0->buffer_type_IO, size_lattice_table,            plattice_table_double,      sizeof(cl_double4));  // Lattice data
        if (!run->turnoff_boundary_extraction) 
            lattice_boundary    = GPU0->buffer_init(GPU0->buffer_type_IO, size_lattice_boundary,         plattice_boundary_double,   sizeof(cl_double4));  // Lattice boundary
#endif
        lattice_parameters      = GPU0->buffer_init(GPU0->buffer_type_Constant, size_lattice_parameters, plattice_parameters_double, sizeof(cl_double));  // Lattice counters and indices
    }
    lattice_measurement         = GPU0->buffer_init(GPU0->buffer_type_IO, size_lattice_measurement,      plattice_measurement,       sizeof(cl_double2)); // Lattice measurement
    lattice_lds                 = GPU0->buffer_init(GPU0->buffer_type_LDS,lds_size,                      0,                          sizeof(cl_double2)); // LDS for reduction
    if (run->get_actions_avr) {
        lattice_energies            = GPU0->buffer_init(GPU0->buffer_type_IO, size_lattice_energies,     plattice_energies,          sizeof(cl_double2)); // Lattice energies
        GPU0->buffer_set_name(lattice_energies, (char*) "lattice_energies");
    }
    if ((run->get_plaquettes_avr) || (run->get_Fmunu) || (run->get_F0mu)){
        lattice_energies_plq    = GPU0->buffer_init(GPU0->buffer_type_IO, size_lattice_energies_plq,     plattice_energies_plq,      sizeof(cl_double2)); // Lattice energies (plaquettes)
        GPU0->buffer_set_name(lattice_energies_plq, (char*) "lattice_energies_plq");
    }
    if (run->get_wilson_loop){
        lattice_wilson_loop     = GPU0->buffer_init(GPU0->buffer_type_IO, size_lattice_wilson_loop,      plattice_wilson_loop,       sizeof(cl_double));  // Wilson loop
        GPU0->buffer_set_name(lattice_wilson_loop, (char*) "lattice_wilson_loop");
    }
    if (run->PL_level > 0){
        lattice_polyakov_loop   = GPU0->buffer_init(GPU0->buffer_type_IO, size_lattice_polyakov_loop,    plattice_polyakov_loop,     sizeof(cl_double2)); // Polyakov loops
        GPU0->buffer_set_name(lattice_polyakov_loop, (char*) "lattice_polyakov_loop");
    }
    if (run->get_acceptance_rate){
        lattice_acceptance_rate     = GPU0->buffer_init(GPU0->buffer_type_IO, size_lattice_acceptance_rate, plattice_acceptance_rate,sizeof(cl_double2)); // Lattice acceptance rate
        GPU0->buffer_set_name(lattice_acceptance_rate, (char*) "lattice_acceptance_rate");
    }
    if (run->get_correlators){
        lattice_correlators         = GPU0->buffer_init(GPU0->buffer_type_IO, size_lattice_correlators,  plattice_correlators,       sizeof(cl_double2)); // Lattice correlators
        GPU0->buffer_set_name(lattice_correlators, (char*) "lattice_correlators");
    }
    if (!run->turnoff_boundary_extraction) 
        GPU0->buffer_set_name(lattice_boundary, (char*) "lattice_boundary");

    GPU0->buffer_set_name(lattice_table,      (char*) "lattice_table");
    GPU0->buffer_set_name(lattice_parameters, (char*) "lattice_parameters");
    GPU0->buffer_set_name(lattice_measurement,(char*) "lattice_measurement");
    GPU0->buffer_set_name(lattice_lds,        (char*) "lattice_lds");
}
void        model::lattice_simulate(void){
    // simulations ______________________________________________________________________________________________________________________________________________
    int NAV_start  = 0;
    int ITER_start = 0;

    if (run->INIT!=0) {
        lattice_simulate_start();

        // measurements
        lattice_measure();

        ITER_counter++;
        if (!run->turnoff_state_save)  lattice_save_state();
        if (!run->turnoff_config_save) lattice_write_configuration();
    }
    lattice_pointer_initial = GPU0->buffer_map(lattice_table); // not exactly - will be rewritten by lattice_pointer_last

    GramSchmidt_iterator = 0;

    NAV_start  = NAV_counter;
    ITER_start = ITER_counter;

    if (run->INIT==0) {
        while (PRNG_counter>PRNG0->PRNG_counter) PRNG0->produce();    // adjust PRNG
        if (GPU0->GPU_debug->brief_report) printf("NAV_start=%u, ITER_start=%u\n",NAV_start,ITER_start);
    }


    // perform thermalization
    for (int i=NAV_start; i<run->NAV; i++){
        lattice_update();
        lattice_orthogonalization();

        if (i % 10 == 0) printf("\rGPU thermalization [%i]",i);
        NAV_counter++;

        lattice_periodic_save_state();
    }

    // perform working cycles
    for (int i=ITER_start; i<run->ITER; i++){ // zero measurement - on initial configuration!!!
        for (int j=0; j<run->NITER; j++){
            lattice_update();
            lattice_orthogonalization();
        }
        if (i % 10 == 0) printf("\rGPU working iteration [%u]",i);

        // measurements
        lattice_measure();

        ITER_counter++;

        lattice_periodic_save_state();
    }
    lattice_print_elapsed_time();

    if (!run->turnoff_state_save) {
        lattice_save_state();
        lattice_pointer_last = lattice_pointer_save;
    } else {
        lattice_pointer_last = GPU0->buffer_map(lattice_table);
    }
    if (!run->turnoff_config_save) lattice_write_configuration();
    prng_pointer = GPU0->buffer_map_float4(PRNG0->PRNG_randoms_id);
}

void        model::lattice_simulate_start(void){
      if (run->ints==model_start_hot) {
            PRNG0->produce();
#if (MODEL_ON==1)
            GPU0->kernel_run_async(sun_init_id);          // HOT Lattice initialization
#else
            GPU0->kernel_run_async(sun_init_X_id);        // HOT Lattice initialization
            if (!run->turnoff_prns) PRNG0->produce();
            GPU0->kernel_run_async(sun_init_Y_id);        // HOT Lattice initialization
            if (!run->turnoff_prns) PRNG0->produce();
            GPU0->kernel_run_async(sun_init_Z_id);        // HOT Lattice initialization
            if (!run->turnoff_prns) PRNG0->produce();
            GPU0->kernel_run_async(sun_init_T_id);        // HOT Lattice initialization
            if (!run->turnoff_prns) PRNG0->produce();
#endif
        } else
            GPU0->kernel_run_async(sun_init_id);          // COLD Lattice initialization
        GPU0->print_stage("lattice initialized");
}
void        model::lattice_measure_plq(void){
        if ((run->get_plaquettes_avr)||(run->get_Fmunu)||(run->get_F0mu)) {
            GPU0->kernel_run_async(sun_measurement_plq_id);          // Lattice measurement (plaquettes)
            GPU0->print_stage("measurement done (plaquettes)");
                plq_index = ITER_counter;
                GPU0->kernel_init_constant_reset(sun_measurement_plq_reduce_id,&plq_index,argument_plq_index);
            GPU0->kernel_run_async(sun_measurement_plq_reduce_id);	// Lattice measurement reduction (plaquettes)
            GPU0->print_stage("measurement reduce done (plaquettes)");
        }
}
void        model::lattice_measure_corr(void){
        if (run->get_correlators) {
            GPU0->kernel_run_async(sun_measurement_corr_id);        // Lattice measurement (correlators)
            GPU0->print_stage("measurement done (correlators)");
                correlators_index = ITER_counter;
                GPU0->kernel_init_constant_reset(sun_measurement_corr_reduce_id,&correlators_index,argument_correlators_index);
            GPU0->kernel_run_async(sun_measurement_corr_reduce_id);	// Lattice measurement reduction (correlators)
            GPU0->print_stage("measurement reduce done (correlators)");
        }
}
void        model::lattice_measure_wilson_loop(void){
        if (run->get_wilson_loop) {
            GPU0->kernel_run_async(sun_measurement_wilson_id);       // Lattice Wilson loop measurement
            GPU0->print_stage("Wilson loop measurement done");
                wilson_index = ITER_counter;
                GPU0->kernel_init_constant_reset(sun_wilson_loop_reduce_id,&wilson_index,argument_wilson_index);
            GPU0->kernel_run_async(sun_wilson_loop_reduce_id);       // Lattice Wilson loop measurement reduction
            GPU0->print_stage("Wilson loop measurement reduce done");
        }
}
void        model::lattice_measure_action(void){
        if (run->get_actions_avr) {
            GPU0->kernel_run_async(sun_measurement_id);                  // Lattice measurement
            GPU0->print_stage("measurement done");
                measurement_index = ITER_counter;
                GPU0->kernel_init_constant_reset(sun_measurement_reduce_id,&measurement_index,argument_measurement_index);

            GPU0->kernel_run_async(sun_measurement_reduce_id);           // Lattice measurement reduction
            GPU0->print_stage("measurement reduce done");
        }
}
void        model::lattice_measure_polyakov_loop(void){
        if (run->PL_level > 0) {
            GPU0->kernel_run_async(sun_polyakov_id);                     // Lattice Polyakov loop measurement
            GPU0->print_stage("Polyakov loop measurement done");
                polyakov_index = ITER_counter;
                GPU0->kernel_init_constant_reset(sun_polyakov_reduce_id,&polyakov_index,argument_polyakov_index);
            GPU0->kernel_run_async(sun_polyakov_reduce_id);              // Lattice Polyakov loop measurement reduction
            GPU0->print_stage("Polyakov loop reduce done");
        }
}
void        model::lattice_update(void){
            if (run->get_acceptance_rate) GPU0->kernel_run_async(sun_clear_measurement_id);
            lattice_update_odd();
            lattice_update_even();
}
void        model::lattice_update_odd(void){
                if (!run->turnoff_prns) PRNG0->produce();
            if (!run->turnoff_updates) GPU0->kernel_run_async(sun_update_odd_X_id);     // Lattice measurement staples
#if (MODEL_ON > 1)
               if (!run->turnoff_prns) PRNG0->produce();
            if (!run->turnoff_updates) GPU0->kernel_run_async(sun_update_odd_Y_id);     // Lattice measurement staples
               if (!run->turnoff_prns) PRNG0->produce();
            if (!run->turnoff_updates) GPU0->kernel_run_async(sun_update_odd_Z_id);     // Lattice measurement staples
               if (!run->turnoff_prns) PRNG0->produce();
            if (!run->turnoff_updates) GPU0->kernel_run_async(sun_update_odd_T_id);     // Lattice measurement staples
#endif
}
void        model::lattice_update_even(void){
                if (!run->turnoff_prns) PRNG0->produce();
            if (!run->turnoff_updates) GPU0->kernel_run_async(sun_update_even_X_id);    // Lattice measurement staples
#if (MODEL_ON > 1)
               if (!run->turnoff_prns) PRNG0->produce();
            if (!run->turnoff_updates) GPU0->kernel_run_async(sun_update_even_Y_id);    // Lattice measurement staples
               if (!run->turnoff_prns) PRNG0->produce();
            if (!run->turnoff_updates) GPU0->kernel_run_async(sun_update_even_Z_id);    // Lattice measurement staples
               if (!run->turnoff_prns) PRNG0->produce();
            if (!run->turnoff_updates) GPU0->kernel_run_async(sun_update_even_T_id);    // Lattice measurement staples
#endif
}

void        model::lattice_orthogonalization(void){
        GramSchmidt_iterator++;
        if (!run->turnoff_gramschmidt) GPU0->kernel_run_async(sun_GramSchmidt_id); // Lattice reunitarization
}
void        model::lattice_periodic_save_state(void){
        // write lattice state every [write_lattice_state_every_secs] seconds
        if ((!run->turnoff_state_save)&&(GPU0->timer_in_seconds_CPU(TIMER_FOR_SAVE)>run->write_lattice_state_every_secs)) {
            lattice_save_state();
            printf("\ncurrent configuration saved\n");
            GPU0->start_timer_CPU(TIMER_FOR_SAVE);      // restart timer for lattice_state save
        }
}
void        model::lattice_print_elapsed_time(void){
        printf("\rGPU simulations are done (%f seconds)\n",GPU0->get_timer_CPU(1));
        time(&ltimeend);
        timeend   = GPU0->get_current_datetime();
}
void        model::lattice_measure_clear(void){
    GPU0->kernel_run_async(sun_clear_measurement_id);
}
void        model::lattice_measure(void){
    if (run->get_acceptance_rate) {
        acceptance_rate_index = ITER_counter;
        GPU0->kernel_init_constant_reset(sun_reduce_acceptance_rate_id,&acceptance_rate_index,argument_acceptance_rate_index);
        if (NAV_counter==(unsigned int) run->NAV) GPU0->kernel_run_async(sun_reduce_acceptance_rate_id);
    }
    lattice_measure_action();
    lattice_measure_plq();
    lattice_measure_corr();
    lattice_measure_polyakov_loop();
    if (!big_lattice) lattice_measure_wilson_loop();
}
void*       model::lattice_table_map(void){
        void* ptr = GPU0->buffer_map_void(lattice_table);
        return ptr;
}
void        model::lattice_table_unmap(void* ptr){
        GPU0->buffer_unmap_void(lattice_table,ptr);
}
void*       model::lattice_table_map_async(void){
        void* ptr = GPU0->buffer_map_void(lattice_table);
        return ptr;
}
void        model::lattice_table_unmap_async(void* ptr){
        GPU0->buffer_unmap_void(lattice_table,ptr);
}
void        model::lattice_wait_for_table_write(void){
        GPU0->buffer_wait_for_write(lattice_table);
}
void        model::lattice_wait_for_table_read(void){
        GPU0->buffer_wait_for_read(lattice_table);
}

void        model::lattice_wait_for_queue_finish(void){
        GPU0->wait_for_queue_finish();
}       

}
