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

namespace model_CL{
using model_CL::model;
using SUN_CPU::SU;
using analysis_CL::analysis;

const char* FXYZ[] = {"x","y","z","t"};   // markers for axis
const char* REIM[] = {"re","im"};         // markers for re / im

#ifndef min
#define min(a,b) ((a) < (b)) ? (a) : (b)
#endif
#ifndef max
#define max(a,b) ((a) > (b)) ? (a) : (b)
#endif

#define TIMER_FOR_ELAPSED      0  // index of timer for elapsed time calculation
#define TIMER_FOR_SIMULATIONS  1  // index of timer for  simulation time calculation
#define TIMER_FOR_SAVE         2  // index of timer for saving lattice states during simulation
#define ND_MAX                32  // maximum dimensions
#define BIN_HEADER_SIZE       64  // length of binary header (in dwords)
#define FNAME_MAX_LENGTH     250  // max length of filename with path

#define SOURCE_UPDATE       "suncl/suncl.cl"
#define SOURCE_MEASUREMENTS "suncl/sun_measurements_cl.cl"
#define SOURCE_POLYAKOV     "suncl/polyakov.cl"
#define SOURCE_WILSON_LOOP  "suncl/wilson_loop.cl"

// SU(N) section __________________________________________________________________________________________
#define DATA_MEASUREMENTS   32 // number of elements for measurements

        char model::path_suncl[FILENAME_MAX]         = "suncl/";
        char model::path_kernel[FILENAME_MAX]        = "kernel/";

// end of SU(N) section -----------------------------------------------------------------------------------

// common section ____________________________________
            model::model(void) {
        PRNG0 = new(PRNG_CL::PRNG);         // PRNG module
        D_A   = new(analysis_CL::analysis); // Data Analysis module
        GPU0  = new(GPU_CL::GPU);           // GPU module

        // setup debug_level
        GPU0->GPU_debug.wait_for_keypress   = false;
        GPU0->GPU_debug.profiling           = false;
        GPU0->GPU_debug.rebuild_binary      = false;
        GPU0->GPU_debug.brief_report        = false;
        GPU0->GPU_debug.show_stage          = false;
        GPU0->GPU_debug.local_run           = false;

        desired_platform    = 0;
        desired_device      = 0;
        device_select       = false;

        write_lattice_state_every_secs = 15.0 * 60.0; // write lattice configuration every 15 minutes

        turnoff_config_save = false; // do not write configurations
        turnoff_prns        = false; // turn off prn production
        turnoff_updates     = false; // turn off lattice updates

        get_actions_avr     = true;  // calculate mean action values
        check_prngs         = false; // check PRNG production

        PRNG_counter = 0;   // counter runs of subroutine PRNG_produce (for load_state purposes)
        NAV_counter  = 0;   // number of performed thermalization cycles
        ITER_counter = 0;   // number of performed working cycles
        LOAD_state   = 0;

        lattice_full_size   = new int[ND_MAX];
        lattice_domain_size = new int[ND_MAX];

        // clear pointers
        lattice_pointer_save         = NULL;
        lattice_pointer_last         = NULL;
        lattice_pointer_initial      = NULL;
        lattice_pointer_measurements = NULL;
        prng_pointer                 = NULL;

        model_create(); // tune particular model

        Analysis = (analysis_CL::analysis::data_analysis*) calloc(DATA_MEASUREMENTS,sizeof(analysis_CL::analysis::data_analysis));
}
            model::~model(void) {
        delete[] lattice_domain_size;
        delete[] lattice_full_size;
        delete PRNG0;
        delete D_A;
        delete Analysis;
        free(lattice_group_elements);

        if (GPU0->GPU_debug.profiling) GPU0->print_time_detailed();

        // output elapsed time
        printf("Elapsed time: %f seconds\n",GPU0->get_timer_CPU(TIMER_FOR_ELAPSED));

        GPU0->make_finish_file(finishpath);
        GPU0->device_finalize(0);

        delete GPU0;
}
void        model::lattice_get_init_file(char* file){
        int parameters_items = 0;
        GPU_CL::GPU::GPU_init_parameters* parameters = GPU0->get_init_file(file);
        if (parameters==NULL) return;

        bool parameters_flag = false;
        while(!parameters_flag){
            if (!strcmp(parameters[parameters_items].Variable,"PLATFORM"))  {desired_platform   = parameters[parameters_items].iVarVal;}
            if (!strcmp(parameters[parameters_items].Variable,"DEVICE"))    {desired_device     = parameters[parameters_items].iVarVal; device_select = true;}

            if (!strcmp(parameters[parameters_items].Variable,"GROUP")) {lattice_group   = parameters[parameters_items].iVarVal;}
            if (!strcmp(parameters[parameters_items].Variable,"ND"))    {lattice_nd      = parameters[parameters_items].iVarVal;}
            if (!strcmp(parameters[parameters_items].Variable,"N1"))    {lattice_domain_size[0] = parameters[parameters_items].iVarVal;}
            if (!strcmp(parameters[parameters_items].Variable,"N2"))    {lattice_domain_size[1] = parameters[parameters_items].iVarVal;}
            if (!strcmp(parameters[parameters_items].Variable,"N3"))    {lattice_domain_size[2] = parameters[parameters_items].iVarVal;}
            if (!strcmp(parameters[parameters_items].Variable,"N4"))    {lattice_domain_size[3] = parameters[parameters_items].iVarVal;}
            if (!strcmp(parameters[parameters_items].Variable,"NS"))    {lattice_domain_size[0] = parameters[parameters_items].iVarVal;
                                                                         lattice_domain_size[1] = parameters[parameters_items].iVarVal;
                                                                         lattice_domain_size[2] = parameters[parameters_items].iVarVal; }
            if (!strcmp(parameters[parameters_items].Variable,"NT"))    {lattice_domain_size[lattice_nd - 1] = parameters[parameters_items].iVarVal;}
            if (!strcmp(parameters[parameters_items].Variable,"L1"))    {lattice_full_size[0] = parameters[parameters_items].iVarVal;}
            if (!strcmp(parameters[parameters_items].Variable,"L2"))    {lattice_full_size[1] = parameters[parameters_items].iVarVal;}
            if (!strcmp(parameters[parameters_items].Variable,"L3"))    {lattice_full_size[2] = parameters[parameters_items].iVarVal;}
            if (!strcmp(parameters[parameters_items].Variable,"L4"))    {lattice_full_size[3] = parameters[parameters_items].iVarVal;}
            if (!strcmp(parameters[parameters_items].Variable,"LS"))    {lattice_full_size[0] = parameters[parameters_items].iVarVal;
                                                                         lattice_full_size[1] = parameters[parameters_items].iVarVal;
                                                                         lattice_full_size[2] = parameters[parameters_items].iVarVal; }
            if (!strcmp(parameters[parameters_items].Variable,"LT"))    {lattice_full_size[lattice_nd - 1] = parameters[parameters_items].iVarVal;}
            if (!strcmp(parameters[parameters_items].Variable,"ITER"))  {ITER            = parameters[parameters_items].iVarVal;}
            if (!strcmp(parameters[parameters_items].Variable,"NITER")) {NITER           = parameters[parameters_items].iVarVal;}
            if (!strcmp(parameters[parameters_items].Variable,"NHIT"))  {NHIT            = parameters[parameters_items].iVarVal;}
            if (!strcmp(parameters[parameters_items].Variable,"BETA"))  {BETA            = parameters[parameters_items].fVarVal;}
            if (!strcmp(parameters[parameters_items].Variable,"PHI"))   {PHI             = parameters[parameters_items].fVarVal;}
            if (!strcmp(parameters[parameters_items].Variable,"OMEGA")) {OMEGA           = parameters[parameters_items].fVarVal;}
            if (!strcmp(parameters[parameters_items].Variable,"NAV"))   {NAV             = parameters[parameters_items].iVarVal;}
            if (!strcmp(parameters[parameters_items].Variable,"INTS"))  {ints            = convert_uint_to_start(parameters[parameters_items].iVarVal);}
            PRNG0->parameters_setup(parameters[parameters_items].Variable,parameters[parameters_items].iVarVal,parameters[parameters_items].txtVarVal);
            if (!strcmp(parameters[parameters_items].Variable,"OUTPUTPATH"))  {
                path = (char*) realloc(path, (strlen(parameters[parameters_items].txtVarVal) + 1) * sizeof(char));
                strcpy_s(path,(strlen(parameters[parameters_items].txtVarVal) + 1),parameters[parameters_items].txtVarVal);
            }
            if (!strcmp(parameters[parameters_items].Variable,"FINISHPATH"))  {
                finishpath = (char*) realloc(finishpath, (strlen(parameters[parameters_items].txtVarVal) + 1) * sizeof(char));
                strcpy_s(finishpath,(strlen(parameters[parameters_items].txtVarVal) + 1),parameters[parameters_items].txtVarVal);
            }
            if (!strcmp(parameters[parameters_items].Variable,"TURNOFFWAITING"))  {
                GPU0->GPU_debug.local_run = true;
            }
            if (!strcmp(parameters[parameters_items].Variable,"REBUILDBINARY"))  {
                GPU0->GPU_debug.rebuild_binary = true;
            }
            if (!strcmp(parameters[parameters_items].Variable,"GETWILSON"))  {
                get_wilson_loop = true;
            }
            if (!strcmp(parameters[parameters_items].Variable,"GETRETRACE"))  {
                get_plaquettes_avr = true;
            }
            if (!strcmp(parameters[parameters_items].Variable,"TURNOFFFMUNU"))  {
                get_Fmunu = false;
                get_F0mu  = false;
            }
            if (!strcmp(parameters[parameters_items].Variable,"F0MU"))  {
                get_F0mu  = true;
                get_Fmunu = false;
            }
            if (!strcmp(parameters[parameters_items].Variable,"FMUNU1"))  {
                get_Fmunu1 = true;
                get_Fmunu2 = false;
                get_Fmunu4 = false;
            }
            if (!strcmp(parameters[parameters_items].Variable,"FMUNU2"))  {
                get_Fmunu2 = true;
                get_Fmunu1 = false;
                get_Fmunu4 = false;
            }
            if (!strcmp(parameters[parameters_items].Variable,"FMUNU4"))  {
                get_Fmunu4 = true;
                get_Fmunu1 = false;
                get_Fmunu2 = false;
            }
            if (!strcmp(parameters[parameters_items].Variable,"FMUNU5"))  {
                get_Fmunu5 = true;
                get_Fmunu6 = false;
                get_Fmunu7 = false;
            }
            if (!strcmp(parameters[parameters_items].Variable,"FMUNU6"))  {
                get_Fmunu6 = true;
                get_Fmunu5 = false;
                get_Fmunu7 = false;
            }
            if (!strcmp(parameters[parameters_items].Variable,"FMUNU7"))  {
                get_Fmunu7 = true;
                get_Fmunu5 = false;
                get_Fmunu6 = false;
            }
            if (!strcmp(parameters[parameters_items].Variable,"WILSONR"))   {wilson_R = parameters[parameters_items].iVarVal;}
            if (!strcmp(parameters[parameters_items].Variable,"WILSONT"))   {wilson_T = parameters[parameters_items].iVarVal;}
            parameters_flag = parameters[parameters_items].final;
            parameters_items++;
        }
}

char*       model::lattice_make_header(void){
    int header_size = 16384;
    header = (char*) calloc(header_size, sizeof(char));
    int j = 0;

    j  += sprintf_s(header+j,header_size-j, " GPU SU(%u) simulator %s\n\n",lattice_group,version);
    j  += sprintf_s(header+j,header_size-j, " Monte Carlo simulation of %uD SU(%u) LGT\n\n",lattice_nd,lattice_group);
    j  += sprintf_s(header+j,header_size-j, " ***************************************************\n");
    j  += sprintf_s(header+j,header_size-j, " Active OpenCL platform   : %s\n",GPU0->platform_get_name(GPU0->GPU_platform));
    j  += sprintf_s(header+j,header_size-j, " Active OpenCL device     : %s\n",GPU0->GPU_info.device_name);
    j  += sprintf_s(header+j,header_size-j, " ***************************************************\n");
    j  += sprintf_s(header+j,header_size-j, " lattice size                : %3u x %3u x %3u x %3u\n",lattice_full_size[0],lattice_full_size[1],lattice_full_size[2],lattice_full_size[3]);
    j  += sprintf_s(header+j,header_size-j, " init                        : %i\n",INIT);
    if (ints == model_start_hot)
        j  += sprintf_s(header+j,header_size-j," ints (0=hot, 1=cold)        : 0\n");
       else
        j  += sprintf_s(header+j,header_size-j," ints (0=hot, 1=cold)        : 1\n");
    j  += PRNG0->print_generator((header+j),(header_size-j));
    j  += sprintf_s(header+j,header_size-j, " nav                         : %i\n",NAV);
    j  += sprintf_s(header+j,header_size-j, " niter                       : %i\n",NITER);
    j  += sprintf_s(header+j,header_size-j, " iter (# of samples)         : %i\n",ITER);
    j  += sprintf_s(header+j,header_size-j, " nhit                        : %i\n",NHIT);
    if (precision == model::model_precision_double) j  += sprintf_s(header+j,header_size-j, " precision                   : double\n");
    if (precision == model::model_precision_single) j  += sprintf_s(header+j,header_size-j, " precision                   : single\n");
    if (precision == model::model_precision_mixed)  j  += sprintf_s(header+j,header_size-j, " precision                   : mixed\n");
    j  += sprintf_s(header+j,header_size-j, " ***************************************************\n");
    j  += model_make_header((header+j),(header_size-j));

    header_index = j;
    return header;
}

void        model::lattice_init(void)
{
    if (!GPU0->GPU_debug.local_run) GPU0->make_start_file(finishpath);

    bool supported_devices = false;
    if (!device_select) {
        // auto-select first GPU_vendor_nVidia device
        supported_devices = GPU0->device_auto_select(GPU_CL::GPU::GPU_vendor_any,GPU_CL::GPU::GPU_vendor_any);
    } else {
        // manual selection of platform and device
        supported_devices = GPU0->device_select(desired_platform,desired_device);
    }
    if(!supported_devices){
        printf("There are no any available OpenCL devices\n");
        exit(0);
    }

    // initialize selected device & show hardwares
    GPU0->device_initialize();
    GPU0->print_available_hardware();
    GPU0->print_stage("device initialized");

    if (INIT==0) lattice_load_state();          // load state file if needed

    if (!big_lattice){
        if (!lattice_domain_size){
            lattice_domain_size = new int[lattice_nd];
            for (int i=0;i<lattice_nd;i++) lattice_domain_size[i] = lattice_full_size[i];
        }
    }

    model_lattice_init();   // model initialization
}

unsigned int model::convert_str_uint(const char* str,unsigned int offset){
    unsigned int result = 0;
    unsigned int pointer = offset;
    for (unsigned int i = 0; i<4; i++){
        if (pointer<strlen(str)) result += (str[pointer++]<<(i*8));
    }
    return result;
}

unsigned int model::convert_start_to_uint(model::model_starts start){
    if (start == model::model_start_hot)  return 0;
    if (start == model::model_start_cold) return 1;
    if (start == model::model_start_gid)  return 2;
    return 0; // return hot start otherwise
}

model::model_starts model::convert_uint_to_start(unsigned int start){
    if (start == 0) return model::model_start_hot;
    if (start == 1) return model::model_start_cold;
    if (start == 2) return model::model_start_gid;
    return model::model_start_hot; // return hot start otherwise
}

unsigned int model::convert_precision_to_uint(model::model_precision precision){
    if (precision == model::model_precision_single) return 1;
    if (precision == model::model_precision_double) return 2;
    if (precision == model::model_precision_mixed)  return 3;
    return 1; // return single precision otherwise
}

model::model_precision model::convert_uint_to_precision(unsigned int precision){
    if (precision == 1) return model::model_precision_single;
    if (precision == 2) return model::model_precision_double;
    if (precision == 3) return model::model_precision_mixed;
    return model::model_precision_single; // return single precision otherwise
}

// model-dependent section ___________________________
void        model::model_create(void){
        // parameters for SU(N) model
        lattice_group_elements = (int*) calloc(4,sizeof(int));
        lattice_group_elements[0] =  1;
        lattice_group_elements[1] =  4;
        lattice_group_elements[2] = 12;

        big_lattice = false;         // simulate big lattice (lattice is divided into several parts)

        turnoff_gramschmidt = false; // turn off Gram-Schmidt orthogonalization

        get_plaquettes_avr  = true;  // calculate mean plaquette values
        get_wilson_loop     = true;  // calculate Wilson loop values
        get_Fmunu           = true;  // calculate Fmunu tensor for H field
        get_F0mu            = false; // calculate Fmunu tensor for E field

        get_Fmunu1          = false; // get Fmunu for lambda1 instead of lambda3
        get_Fmunu2          = false; // get Fmunu for lambda2 instead of lambda3
        get_Fmunu4          = false; // get Fmunu for lambda4 instead of lambda3

        get_Fmunu5          = false; // get Fmunu for lambda5 instead of lambda8
        get_Fmunu6          = false; // get Fmunu for lambda6 instead of lambda8
        get_Fmunu7          = false; // get Fmunu for lambda7 instead of lambda8
}

int         model::model_make_header(char* header,int header_size){
    int j = 0;

    j  += sprintf_s(header+j,header_size-j, " Wilson loop R               : %i\n",wilson_R);
    j  += sprintf_s(header+j,header_size-j, " Wilson loop T               : %i\n",wilson_T);
    if (get_Fmunu)
        j  += sprintf_s(header+j,header_size-j," FMUNU(%u, %u)\n",Fmunu_index1,Fmunu_index2);
    if (get_F0mu)
        j  += sprintf_s(header+j,header_size-j," F0MU(%u, %u)\n",Fmunu_index1,Fmunu_index2);

    j  += sprintf_s(header+j,header_size-j, " BETA                        : %16.13e\n",BETA);
    j  += sprintf_s(header+j,header_size-j, " PHI   (lambda_3)            : %16.13e\n",PHI);
    j  += sprintf_s(header+j,header_size-j, " OMEGA (lambda_8)            : %16.13e\n",OMEGA);
    j  += sprintf_s(header+j,header_size-j, " ***************************************************\n");

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
    j  += sprintf_s(header+j,header_size-j, " ***************************************************\n");

    for (int i=0;i<((lattice_nd-1)*2+2)*2;i++)
        j  += sprintf_s(header+j,header_size-j, " Mean %-20s: % 16.13e\n",Analysis[DM_Fmunu_3+i].data_name,Analysis[DM_Fmunu_3+i].mean_value);
    for (int i=0;i<((lattice_nd-1)*2+2)*2;i++)
        j  += sprintf_s(header+j,header_size-j, " Variance %-16s: % 16.13e\n",Analysis[DM_Fmunu_3+i].data_name,Analysis[DM_Fmunu_3+i].variance);
    j  += sprintf_s(header+j,header_size-j, " ***************************************************\n");
    for (int jj=0;jj<((lattice_nd-2)*(lattice_nd-1)+2)*2;jj++){
        j  += sprintf_s(header+j,header_size-j, " GPU last %-16s: % 16.13e\n",Analysis[jj+DM_Fmunu_3].data_name,Analysis[jj+DM_Fmunu_3].GPU_last_value);
        j  += sprintf_s(header+j,header_size-j, " CPU last %-16s: % 16.13e\n",Analysis[jj+DM_Fmunu_3].data_name,Analysis[jj+DM_Fmunu_3].CPU_last_value);
    }
    j  += sprintf_s(header+j,header_size-j, " ***************************************************\n");
    for (int jj=0;jj<(lattice_nd-2)*(lattice_nd-1)*2;jj++)
        j  += sprintf_s(header+j,header_size-j, " CPU l.varnc %-13s: % 16.13e\n",Analysis[jj+DM_Fmunu_3].data_name,Analysis[jj+DM_Fmunu_3].CPU_last_variance);
    j  += sprintf_s(header+j,header_size-j, " ***************************************************\n");
    j  += sprintf_s(header+j,header_size-j, " CPU last %-16s: % 16.13e\n",Analysis[DM_S_total].data_name,Analysis[DM_S_total].CPU_last_value);
    j  += sprintf_s(header+j,header_size-j, " GPU last %-16s: % 16.13e\n",Analysis[DM_S_total].data_name,Analysis[DM_S_total].GPU_last_value);
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
    if (analysis_CL::analysis::results_verification)
        j  += sprintf_s(header+j,header_size-j, " *** Verification successfully passed! *************\n");
    else
        j  += sprintf_s(header+j,header_size-j, " --- Verification failed! --------------------------\n");
    j  += sprintf_s(header+j,header_size-j, " Data fields:\n");
    j  += sprintf_s(header+j,header_size-j, "    #, ");
    if (get_plaquettes_avr){
                           j += sprintf_s(header+j,header_size-j, "%-22s",Analysis[DM_Plq_spat].data_name);
                           j += sprintf_s(header+j,header_size-j, "%-22s",Analysis[DM_Plq_temp].data_name);
                           j += sprintf_s(header+j,header_size-j, "%-22s",Analysis[DM_Plq_total].data_name);
    }
    if (get_wilson_loop)
                           j += sprintf_s(header+j,header_size-j, "%-22s",Analysis[DM_Wilson_loop].data_name);
    if (get_actions_avr){
                           j += sprintf_s(header+j,header_size-j, "%-22s",Analysis[DM_S_spat].data_name);
                           j += sprintf_s(header+j,header_size-j, "%-22s",Analysis[DM_S_temp].data_name);
                           j += sprintf_s(header+j,header_size-j, "%-22s",Analysis[DM_S_total].data_name);
    }
    if (PL_level > 0){
                           j += sprintf_s(header+j,header_size-j, "%-22s",Analysis[DM_Polyakov_loop].data_name);
                           j += sprintf_s(header+j,header_size-j, "%-22s",Analysis[DM_Polyakov_loop_im].data_name);
    }
    if (PL_level > 1){
                           j += sprintf_s(header+j,header_size-j, "%-22s",Analysis[DM_Polyakov_loop_P2].data_name);
                           j += sprintf_s(header+j,header_size-j, "%-22s",Analysis[DM_Polyakov_loop_P4].data_name);
    }
    if ((get_Fmunu)||(get_F0mu))
        for (int i=0;i<((lattice_nd-1)*2+2)*2;i++)
            j  += sprintf_s(header+j,header_size-j, "%-22s",Analysis[DM_Fmunu_3+i].data_name);
    j  += sprintf_s(header+j,header_size-j, "\n");
    j  += sprintf_s(header+j,header_size-j, " ***************************************************\n");

    header_index = j;

    return header;
}

void        model::lattice_print_measurements(void){
    printf("\n #        ");
    if (get_plaquettes_avr){
        printf(" %-11s",Analysis[DM_Plq_spat].data_name);
        printf(" %-11s",Analysis[DM_Plq_temp].data_name);
    }
    if (get_wilson_loop)
        printf(" %-11s",Analysis[DM_Wilson_loop].data_name);
    if (get_actions_avr){
        printf(" %-11s",Analysis[DM_S_spat].data_name);
        printf(" %-11s",Analysis[DM_S_temp].data_name);
    }
    if (PL_level > 0){
        printf(" %-11s",Analysis[DM_Polyakov_loop].data_name);
        printf(" %-11s",Analysis[DM_Polyakov_loop_im].data_name);
    }
    printf("\n");
    unsigned int offs = 0;
    for (int i=0; i<ITER; i++) {
        if ((i<5)||(i==(ITER - 1))){
            printf("[%5u]: ",(i+offs));
            if (get_plaquettes_avr){
                printf(" % 10.8f",Analysis[DM_Plq_spat].data[i]);
                printf(" % 10.8f",Analysis[DM_Plq_temp].data[i]);
            }
            if (get_wilson_loop)
                printf(" % 10.8f",Analysis[DM_Wilson_loop].data[i]);
            if (get_actions_avr){
                printf(" % 10.8f",Analysis[DM_S_spat].data[i]);
                printf(" % 10.8f",Analysis[DM_S_temp].data[i]);
            }
            if (PL_level > 0){
                printf(" % 10.8f",Analysis[DM_Polyakov_loop].data[i]);
                printf(" % 10.8f",Analysis[DM_Polyakov_loop_im].data[i]);
            }
            printf("\n");
        }
    }
}

void        model::lattice_analysis(void){
        for (int i=0; i<=DM_max;i++){
            Analysis[i].data_size       = ITER;
            Analysis[i].pointer_offset  = 0;
            if (precision==model_precision_double) Analysis[i].precision_single = false;
                else                               Analysis[i].precision_single = true;
            if ((i&1) == 0) Analysis[i].storage_type = GPU_CL::GPU::GPU_storage_double2low;
            else            Analysis[i].storage_type = GPU_CL::GPU::GPU_storage_double2high;
        }

        Analysis[DM_Wilson_loop].storage_type         = GPU_CL::GPU::GPU_storage_double;

    if (get_actions_avr) {
        // S_spat
        Analysis[DM_S_spat].pointer         = GPU0->buffer_map(lattice_energies);
        Analysis[DM_S_spat].denominator     = ((double) (lattice_full_site   * 3));
        Analysis[DM_S_spat].data_name       = "S_spat";
        D_A->lattice_data_analysis(&Analysis[DM_S_spat]);

        // S_temp
        Analysis[DM_S_temp].pointer         = Analysis[DM_S_spat].pointer;
        Analysis[DM_S_temp].denominator     = ((double) (lattice_full_site   * 3));
        Analysis[DM_S_temp].data_name       = "S_temp";
        D_A->lattice_data_analysis(&Analysis[DM_S_temp]);

        // S_total
        Analysis[DM_S_total].data_name      = "S_total";
        D_A->lattice_data_analysis_joint(&Analysis[DM_S_total],&Analysis[DM_S_spat],&Analysis[DM_S_temp]);
    }
    if (get_plaquettes_avr) {
        // Plq_spat
        Analysis[DM_Plq_spat].pointer         = GPU0->buffer_map(lattice_energies_plq);
        Analysis[DM_Plq_spat].denominator     = ((double) (lattice_full_site   * 3));
        Analysis[DM_Plq_spat].data_name       = "Plq_spat";
        D_A->lattice_data_analysis(&Analysis[DM_Plq_spat]);

        // Plq_temp
        Analysis[DM_Plq_temp].pointer         = Analysis[DM_Plq_spat].pointer;
        Analysis[DM_Plq_temp].denominator     = ((double) (lattice_full_site   * 3));
        Analysis[DM_Plq_temp].data_name       = "Plq_temp";
        D_A->lattice_data_analysis(&Analysis[DM_Plq_temp]);

        // Plq_total
        Analysis[DM_Plq_total].data_name      = "Plq_total";
        D_A->lattice_data_analysis_joint(&Analysis[DM_Plq_total],&Analysis[DM_Plq_spat],&Analysis[DM_Plq_temp]);
    }    
    if (PL_level > 0) {
        // Polyakov_loop
        Analysis[DM_Polyakov_loop].pointer            = GPU0->buffer_map(lattice_polyakov_loop);
        Analysis[DM_Polyakov_loop].denominator        = ((double) (lattice_full_n1n2n3 * lattice_group));
        Analysis[DM_Polyakov_loop].data_name          = "Polyakov_loop";
        D_A->lattice_data_analysis(&Analysis[DM_Polyakov_loop]);

        // Polyakov_loop_im
        Analysis[DM_Polyakov_loop_im].pointer         = Analysis[DM_Polyakov_loop].pointer;
        Analysis[DM_Polyakov_loop_im].denominator     = ((double) (lattice_full_n1n2n3 * lattice_group));
        Analysis[DM_Polyakov_loop_im].data_name       = "Polyakov_loop_im";
        D_A->lattice_data_analysis(&Analysis[DM_Polyakov_loop_im]);
    }
    if (PL_level > 1){
        // Polyakov_loop_P2
        Analysis[DM_Polyakov_loop_P2].pointer         = Analysis[DM_Polyakov_loop].pointer;
        Analysis[DM_Polyakov_loop_P2].pointer_offset  = lattice_polyakov_loop_size;
        Analysis[DM_Polyakov_loop_P2].denominator     = ((double) (lattice_full_n1n2n3 * lattice_group * lattice_group));
        Analysis[DM_Polyakov_loop_P2].data_name       = "Polyakov_loop_P2";
        D_A->lattice_data_analysis(&Analysis[DM_Polyakov_loop_P2]);

        // Polyakov_loop_P4
        Analysis[DM_Polyakov_loop_P4].pointer         = Analysis[DM_Polyakov_loop].pointer;
        Analysis[DM_Polyakov_loop_P4].pointer_offset  = lattice_polyakov_loop_size;
        Analysis[DM_Polyakov_loop_P4].denominator     = ((double) (lattice_full_n1n2n3 * lattice_group * lattice_group * lattice_group * lattice_group));
        Analysis[DM_Polyakov_loop_P4].data_name       = "Polyakov_loop_P4";
        D_A->lattice_data_analysis(&Analysis[DM_Polyakov_loop_P4]);
    }
    if (get_wilson_loop) {
        // Wilson_loop
        Analysis[DM_Wilson_loop].pointer         = GPU0->buffer_map(lattice_wilson_loop);
        Analysis[DM_Wilson_loop].denominator     = ((double) (lattice_full_site * 3));
        Analysis[DM_Wilson_loop].data_name       = "Wilson_loop";
        D_A->lattice_data_analysis(&Analysis[DM_Wilson_loop]);
    }
    if ((get_Fmunu)||(get_F0mu)) {
        // Fmunu_xy_3_re
        unsigned int* F_pointr;
        if (!Analysis[DM_Plq_spat].pointer)
            F_pointr = GPU0->buffer_map(lattice_energies_plq);
        else
            F_pointr = Analysis[DM_Plq_spat].pointer;

        int index_x = 0;
        int index_y = 1;
        for (int i=0;i<(lattice_nd-1)*2;i++){   // loop for re and im
            Analysis[i+DM_Fmunu_3].denominator = ((double) (lattice_full_site));
            Analysis[i+DM_Fmunu_8].denominator = ((double) (lattice_full_site));

            Analysis[i+DM_Fmunu_3].pointer = F_pointr;
            Analysis[i+DM_Fmunu_8].pointer = F_pointr;

            Analysis[i+DM_Fmunu_3].pointer_offset = lattice_energies_offset * (1 + (i >> 1));
            Analysis[i+DM_Fmunu_8].pointer_offset = lattice_energies_offset * (4 + (i >> 1));

            Analysis[i+DM_Fmunu_3].data_name = (char*) calloc(14,sizeof(char));
            Analysis[i+DM_Fmunu_8].data_name = (char*) calloc(14,sizeof(char));

            if (get_Fmunu) {
                sprintf_s((char*) Analysis[i+DM_Fmunu_3].data_name,14,"Fmunu_%s%s_%u_%s",FXYZ[index_x],FXYZ[index_y],Fmunu_index1,REIM[(i&1)]);
                sprintf_s((char*) Analysis[i+DM_Fmunu_8].data_name,14,"Fmunu_%s%s_%u_%s",FXYZ[index_x],FXYZ[index_y],Fmunu_index2,REIM[(i&1)]);
            } else {
                sprintf_s((char*) Analysis[i+DM_Fmunu_3].data_name,14,"Fmunu_%st_%u_%s",FXYZ[((i>>1)&3)],Fmunu_index1,REIM[(i&1)]);
                sprintf_s((char*) Analysis[i+DM_Fmunu_8].data_name,14,"Fmunu_%st_%u_%s",FXYZ[((i>>1)&3)],Fmunu_index2,REIM[(i&1)]);
            }


            D_A->lattice_data_analysis(&Analysis[i+DM_Fmunu_3]);
            D_A->lattice_data_analysis(&Analysis[i+DM_Fmunu_8]);

            if ((i&1)==1) {
                index_y++;
                if (index_y>(lattice_nd-2)) {
                    index_x++;
                    index_y=index_x+1;
                }
            }
        }

        Analysis[DM_Fmunu_abs_3_re].data_name = (char*) calloc(15,sizeof(char));
        Analysis[DM_Fmunu_abs_3_im].data_name = (char*) calloc(15,sizeof(char));
        Analysis[DM_Fmunu_abs_8_re].data_name = (char*) calloc(15,sizeof(char));
        Analysis[DM_Fmunu_abs_8_im].data_name = (char*) calloc(15,sizeof(char));

        sprintf_s((char*) Analysis[DM_Fmunu_abs_3_re].data_name,15,"Fmunu_abs_%1u_re",Fmunu_index1);
        sprintf_s((char*) Analysis[DM_Fmunu_abs_3_im].data_name,15,"Fmunu_abs_%1u_im",Fmunu_index1);
        sprintf_s((char*) Analysis[DM_Fmunu_abs_8_re].data_name,15,"Fmunu_abs_%1u_re",Fmunu_index2);
        sprintf_s((char*) Analysis[DM_Fmunu_abs_8_im].data_name,15,"Fmunu_abs_%1u_im",Fmunu_index2);

        // Fmunu_abs_3
        D_A->lattice_data_analysis_joint3(&Analysis[DM_Fmunu_abs_3_re],&Analysis[DM_Fmunu_xy_3_re],&Analysis[DM_Fmunu_xz_3_re],&Analysis[DM_Fmunu_yz_3_re]);
        D_A->lattice_data_analysis_joint3(&Analysis[DM_Fmunu_abs_3_im],&Analysis[DM_Fmunu_xy_3_im],&Analysis[DM_Fmunu_xz_3_im],&Analysis[DM_Fmunu_yz_3_im]);

        // Fmunu_abs_8
        D_A->lattice_data_analysis_joint3(&Analysis[DM_Fmunu_abs_8_re],&Analysis[DM_Fmunu_xy_8_re],&Analysis[DM_Fmunu_xz_8_re],&Analysis[DM_Fmunu_yz_8_re]);
        D_A->lattice_data_analysis_joint3(&Analysis[DM_Fmunu_abs_8_im],&Analysis[DM_Fmunu_xy_8_im],&Analysis[DM_Fmunu_xz_8_im],&Analysis[DM_Fmunu_yz_8_im]);
    }
    
}

void        model::lattice_write_results(void) {
    FILE *stream;
    char buffer[250];
    int j;

    char* header2 = lattice_make_header2();
    printf("%s\n",header2);

    j  = sprintf_s(buffer  ,sizeof(buffer),  "%s",path);
    j += sprintf_s(buffer+j,sizeof(buffer)-j,"%s",fprefix);
    j += sprintf_s(buffer+j,sizeof(buffer)-j,"%.2s-%.3s-%.2s-%.2s-%.2s-%.2s.txt",timeend+22,timeend+4,timeend+8,timeend+11,timeend+14,timeend+17);

    fopen_s(&stream,buffer,"w+");
    if(stream)
    {
        fprintf(stream,header);

        // write plaquette data
        for (int i=0; i<ITER; i++) {
                                         fprintf(stream, "%5i",i);
            if (get_plaquettes_avr){
                fprintf(stream, " % 16.13e",Analysis[DM_Plq_spat].data[i]);
                fprintf(stream, " % 16.13e",Analysis[DM_Plq_temp].data[i]);
                fprintf(stream, " % 16.13e",Analysis[DM_Plq_total].data[i]);
            }
            if (get_wilson_loop)
                fprintf(stream, " % 16.13e",Analysis[DM_Wilson_loop].data[i]);
            if (get_actions_avr){
                fprintf(stream, " % 16.13e",Analysis[DM_S_spat].data[i]);
                fprintf(stream, " % 16.13e",Analysis[DM_S_temp].data[i]);
                fprintf(stream, " % 16.13e",Analysis[DM_S_total].data[i]);
            }
            if (PL_level > 0){
                fprintf(stream, " % 16.13e",Analysis[DM_Polyakov_loop].data[i]);
                fprintf(stream, " % 16.13e",Analysis[DM_Polyakov_loop_im].data[i]);
            }
            if (PL_level > 1){
                fprintf(stream, " % 16.13e",Analysis[DM_Polyakov_loop_P2].data[i]);
                fprintf(stream, " % 16.13e",Analysis[DM_Polyakov_loop_P4].data[i]);
            }
            if ((get_Fmunu)||(get_F0mu))
                for (int jj=0;jj<((lattice_nd-2)*(lattice_nd-1)+2)*2;jj++){
                    fprintf(stream, " % 16.13e",Analysis[jj+DM_Fmunu_3].data[i]);
                }
            fprintf(stream, "\n");
        }

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

    if ((get_plaquettes_avr) || (get_Fmunu) || (get_F0mu))
        lattice_energies_plq_save  = GPU0->buffer_map(lattice_energies_plq);
    if (get_wilson_loop)
        lattice_wilson_loop_save   = GPU0->buffer_map(lattice_wilson_loop);
    if (PL_level > 0)
        lattice_polyakov_loop_save = GPU0->buffer_map(lattice_polyakov_loop);

    FILE *stream;
    char buffer[250];
    int j = 0;

    unsigned int* head = lattice_make_bin_header();

    j  = sprintf_s(buffer  ,sizeof(buffer),  "%s",path);
    j += sprintf_s(buffer+j,sizeof(buffer)-j,"%s",fprefix);
    j += sprintf_s(buffer+j,sizeof(buffer)-j,"%.2s-%.3s-%.2s-%.2s-%.2s-%.2s.qcg",timesave+22,timesave+4,timesave+8,timesave+11,timesave+14,timesave+17);

    fopen_s(&stream,buffer,"wb");
    if(stream)
    {
        fwrite(head,sizeof(unsigned int),BIN_HEADER_SIZE,stream);                                       // write header
        fwrite(lattice_measurement_save, sizeof(cl_double2), lattice_measurement_size_F, stream);       // write measurements
        fwrite(lattice_energies_save,    sizeof(cl_double2), lattice_energies_size, stream);            // write energies
        if ((get_plaquettes_avr) || (get_Fmunu) || (get_F0mu))
            fwrite(lattice_energies_plq_save,  sizeof(cl_double2), lattice_energies_size_F, stream);    // write energies_plq
        if (get_wilson_loop)
            fwrite(lattice_wilson_loop_save,   sizeof(cl_double),  lattice_energies_size, stream);      // write wilson loop
        if (PL_level > 0)
            fwrite(lattice_polyakov_loop_save, sizeof(cl_double2), lattice_polyakov_loop_size, stream); // write polyakov loop
        if (precision == model_precision_single)                                                        // write configuration
            fwrite(lattice_pointer_save, sizeof(cl_float4), lattice_table_size, stream);
        else
            fwrite(lattice_pointer_save, sizeof(cl_double4), lattice_table_size, stream);

        unsigned int hlen  = BIN_HEADER_SIZE*sizeof(unsigned int);
            if (GPU0->GPU_debug.brief_report) printf("Header: 0x%X-0x%X\n",0,hlen);
        unsigned int hlen2 = lattice_measurement_size_F * sizeof(cl_double2);
            if (GPU0->GPU_debug.brief_report) printf("Lattice mesurements: 0x%X-0x%X\n",hlen,(hlen+hlen2));
            hlen += hlen2;
        hlen2 = lattice_energies_size * sizeof(cl_double2);
            if (GPU0->GPU_debug.brief_report) printf("Lattice energies: 0x%X-0x%X\n",hlen,(hlen+hlen2));
            hlen += hlen2;
        if ((get_plaquettes_avr) || (get_Fmunu) || (get_F0mu)){
            hlen2 = lattice_energies_size_F * sizeof(cl_double2);
            if (GPU0->GPU_debug.brief_report) printf("Lattice plaq av.: 0x%X-0x%X\n",hlen,(hlen+hlen2));
            hlen += hlen2;
        }
        if (get_wilson_loop){
            hlen2 = lattice_energies_size * sizeof(cl_double);
            if (GPU0->GPU_debug.brief_report) printf("Lattice Wilson loop: 0x%X-0x%X\n",hlen,(hlen+hlen2));
            hlen += hlen2;
        }
        if (PL_level > 0){
            hlen2 = lattice_polyakov_loop_size * sizeof(cl_double2);
            if (GPU0->GPU_debug.brief_report) printf("Lattice Polyakov loop: 0x%X-0x%X\n",hlen,(hlen+hlen2));
            hlen += hlen2;
        }
        if (precision == model_precision_single)                                                        // write configuration
            hlen2 = lattice_table_size * sizeof(cl_float4);
        else
            hlen2 = lattice_table_size * sizeof(cl_double4);
        if (GPU0->GPU_debug.brief_report) printf("Lattice data: 0x%X-0x%X\n",hlen,(hlen+hlen2));
        hlen += hlen2;


        if ( fclose(stream) ) printf( "The file was not closed!\n" );
    }
    free(head);
}

unsigned int*   model::lattice_make_bin_header(void){
    int k;
    // bin header structure:
    // 0 - 1 - prefix (8 bytes)
    // 2 - 3 - version
    // 4 - ...
    unsigned int* result = (unsigned int*) calloc(BIN_HEADER_SIZE,sizeof(unsigned int));
    const char* bin_prefix = "QCDGPU";
    k =  0; for(int i=0; i<ceil((double) strlen(bin_prefix)/4); i++) result[k++] = convert_str_uint(bin_prefix,i*4);
    k =  2; for(int i=0; i<ceil((double) strlen(version)/4); i++)    result[k++] = convert_str_uint(version,i*4);
    k =  4;
    result[k++] = INIT;
    result[k++] = convert_start_to_uint(ints);
    result[k++] = PRNG0->PRNG_randseries;
    result[k++] = PRNG0->convert_generator_to_uint(PRNG0->PRNG_generator);
    result[k++] = PRNG0->RL_nskip;
    result[k++] = PRNG0->PRNG_counter;
    result[k++] = NAV;
    result[k++] = NAV_counter;
    result[k++] = NITER;
    result[k++] = ITER;
    result[k++] = ITER_counter;
    result[k++] = NHIT;
    result[k++] = wilson_R;
    result[k++] = wilson_T;
    result[k++] = convert_precision_to_uint(precision);
    result[k++] = GPU0->convert_to_uint_LOW( BETA);
    result[k++] = GPU0->convert_to_uint_HIGH(BETA);
    result[k++] = GPU0->convert_to_uint_LOW( PHI);
    result[k++] = GPU0->convert_to_uint_HIGH(PHI);
    result[k++] = GPU0->convert_to_uint_LOW( OMEGA);
    result[k++] = GPU0->convert_to_uint_HIGH(OMEGA);
    result[k++] = PL_level;
    result[k++] = Fmunu_index1;
    result[k++] = Fmunu_index2;
    result[k++] = get_Fmunu;
    result[k++] = get_actions_avr;
    result[k++] = get_plaquettes_avr;
    result[k++] = get_wilson_loop;
    result[k++] = lattice_measurement_size_F; // measurements initialization
    result[k++] = lattice_energies_size;      // energies initialization
    result[k++] = lattice_energies_size_F;    // energies initialization (with Fmunu tensor, if needed)
    result[k++] = lattice_polyakov_loop_size; // polyakov loop initialization
    result[k++] = lattice_group;
    result[k++] = lattice_nd;
    for (int i=0; i<lattice_nd; i++) result[k++] = lattice_full_size[i];
    for (int i=0; i<lattice_nd; i++) result[k++] = lattice_domain_size[i];

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
    k =  0; for(int i=0; i<ceil((double) strlen(bin_prefix)/4); i++) result &= (head[k++] == convert_str_uint(bin_prefix,i*4));
    k =  2; for(int i=0; i<ceil((double) strlen(version)/4); i++)    result &= (head[k++] == convert_str_uint(version,i*4));
    if (!result) return result;
    k =  4;
      k++;
    ints = convert_uint_to_start(head[k++]);            // 0x14
    PRNG0->PRNG_randseries = head[k++];                 // 0x18
    PRNG0->PRNG_generator = PRNG0->convert_uint_to_generator(head[k++]);    // 0x1C
    PRNG0->RL_nskip = head[k++];                        // 0x20
    PRNG_counter = head[k++];                           // 0x24
    NAV = head[k++];                                    // 0x28
    NAV_counter = head[k++];                            // 0x2C
    NITER = head[k++];                                  // 0x30
    ITER = head[k++];                                   // 0x34
    ITER_counter = head[k++];                           // 0x38
    NHIT = head[k++];                                   // 0x3C
    wilson_R = head[k++];                               // 0x40
    wilson_T = head[k++];                               // 0x44
    precision = convert_uint_to_precision(head[k++]);
    unsigned int get_low  = head[k++];                  // 0x48
    unsigned int get_high = head[k++];                  // 0x4C
    BETA  = GPU0->convert_to_double(get_low,get_high);
    get_low  = head[k++]; get_high = head[k++];         // 0x50, 0x54
    PHI   = GPU0->convert_to_double(get_low,get_high);
    get_low  = head[k++]; get_high = head[k++];         // 0x58, 0x5C
    OMEGA = GPU0->convert_to_double(get_low,get_high);
    PL_level = head[k++];                               // 0x60
    Fmunu_index1 = head[k++];                           // 0x64
    Fmunu_index2 = head[k++];                           // 0x68
    get_Fmunu = (head[k++]==0 ? false : true);          // 0x6C
    get_actions_avr = (head[k++]==0 ? false : true);    // 0x70
    get_plaquettes_avr = (head[k++]==0 ? false : true); // 0x74
    get_wilson_loop = (head[k++]==0 ? false : true);    // 0x78
    lattice_measurement_size_F = head[k++];             // 0x7C - measurements initialization
    lattice_energies_size = head[k++];                  // 0x80 - energies initialization
    lattice_energies_size_F = head[k++];                // 0x84 - energies initialization (with Fmunu tensor, if needed)
    lattice_polyakov_loop_size = head[k++];             // 0x88 - polyakov loop initialization
    lattice_group = head[k++];                          // 0x8C
    lattice_nd = head[k++];                             // 0x90
    for (int i=0; i<lattice_nd; i++) lattice_full_size[i] = head[k++];   // 0x98, 0x9C, 0xA0, 0xA4
    for (int i=0; i<lattice_nd; i++) lattice_domain_size[i] = head[k++]; // 0xA8, 0xAC, 0xB0, 0xB4

    return result;
}

void        model::lattice_load_state(void){
    FILE *stream;
    char buffer[250];
    int j = 0;

    unsigned int* head = (unsigned int*) calloc(BIN_HEADER_SIZE,sizeof(unsigned int));
    bool result = true;

    j  = sprintf_s(buffer  ,sizeof(buffer),  "%s",path);
    j += sprintf_s(buffer+j,sizeof(buffer)-j,"%s",fstate);

    fopen_s(&stream,buffer,"rb");
    if(stream)
    {
        fread(head,sizeof(unsigned int),BIN_HEADER_SIZE,stream);    // load header
        if (LOAD_state==0) result = lattice_load_bin_header(head);
        if (!result) printf("[ERROR in header!!!]\n");
        if (LOAD_state==1) {
            fread(plattice_measurement, sizeof(cl_double2), lattice_measurement_size_F, stream);       // load measurements
            fread(plattice_energies,    sizeof(cl_double2), lattice_energies_size, stream);            // load energies
            if ((get_plaquettes_avr) || (get_Fmunu) || (get_F0mu))
                fread(plattice_energies_plq,  sizeof(cl_double2), lattice_energies_size_F, stream);    // load energies_plq
            if (get_wilson_loop)
                fread(plattice_wilson_loop,   sizeof(cl_double),  lattice_energies_size, stream);      // load wilson loop
            if (PL_level > 0)
                fread(plattice_polyakov_loop, sizeof(cl_double2), lattice_polyakov_loop_size, stream); // load polyakov loop
            if (precision == model_precision_single)                                                   // load configuration
                fread(plattice_table_float,   sizeof(cl_float4),  lattice_table_size, stream);
            else
                fread(plattice_table_double,  sizeof(cl_double4), lattice_table_size, stream);

            if ( fclose(stream) ) printf( "The file was not closed!\n" );
        }
        LOAD_state++;
    }
    free(head);
}
 
void        model::model_lattice_init(void){
    if(lattice_group == 2)
     {
        Fmunu_index1 = 3;
        if (get_Fmunu1) Fmunu_index1 = 1;
           Fmunu_index2 = 2;
     }
    if(lattice_group == 3)
     {
        Fmunu_index1 = 3;
            if (get_Fmunu1) Fmunu_index1 = 1;
            if (get_Fmunu2) Fmunu_index1 = 2;
            if (get_Fmunu4) Fmunu_index1 = 4;
        Fmunu_index2 = 8;
            if (get_Fmunu5) Fmunu_index2 = 5;
            if (get_Fmunu6) Fmunu_index2 = 6;
            if (get_Fmunu7) Fmunu_index2 = 7;
     }
    
    if ((get_Fmunu)&&(get_F0mu)) get_F0mu = 0;  // only one field (H or E) may be calculated

    local_size_intel = (GPU0->GPU_info.device_vendor == GPU0->GPU::GPU_vendor_Intel) ? 64 : 0;
    if (GPU0->GPU_limit_max_workgroup_size) local_size_intel = GPU0->GPU_limit_max_workgroup_size;
    size_t workgroup_factor = (local_size_intel) ? local_size_intel : 32;

    if (big_lattice) {
        // lattice is divided into domains
        lattice_domain_n1     = (lattice_domain_size[0]+2);
    } else {
        // lattice is not divided into domains
        lattice_domain_n1     = lattice_domain_size[0];
    }
    lattice_domain_n1n2   = lattice_domain_n1      * lattice_domain_size[1];
    lattice_domain_n2n3   = lattice_domain_size[1] * lattice_domain_size[2];
    lattice_domain_n1n2n3 = lattice_domain_n1      * lattice_domain_n2n3;
    lattice_domain_exact_n1n2n3 = lattice_domain_size[0] * lattice_domain_size[1] * lattice_domain_size[2];
    lattice_domain_n2n3n4 =  lattice_domain_size[1]    * lattice_domain_size[2] * lattice_domain_size[3];

    lattice_domain_site       = lattice_domain_n1;
    lattice_domain_exact_site = lattice_domain_size[0];
    lattice_full_site         = lattice_full_size[0];
    for (int i=1;i<lattice_nd;i++) {
        lattice_domain_site       *= lattice_domain_size[i];
        lattice_domain_exact_site *= lattice_domain_size[i];
        lattice_full_site         *= lattice_full_size[i];
    }

    lattice_full_n1n2   = lattice_full_size[0] * lattice_full_size[1];
    lattice_full_n2n3   = lattice_full_size[1] * lattice_full_size[2];
    lattice_full_n1n2n3 = lattice_full_size[0] * lattice_full_size[1] * lattice_full_size[2];
    lattice_full_n2n3n4 = lattice_full_size[1] * lattice_full_size[2] * lattice_full_size[3];

    lattice_full_link         = lattice_nd * lattice_full_site;
    lattice_domain_link       = lattice_nd * lattice_domain_site;
    lattice_domain_exact_link = lattice_nd * lattice_domain_exact_site;

    lattice_boundary_exact_size = lattice_domain_size[1];    // Size of boundary slice (x=1 in depth)
    for (int i=2;i<lattice_nd;i++) lattice_boundary_exact_size *= lattice_domain_size[i];

    lattice_parameters_size         = GPU0->buffer_size_align(MODEL_parameter_size);
    lattice_energies_size           = GPU0->buffer_size_align(ITER);                      // number of working iterations
    lattice_energies_size_F         = lattice_energies_size * MODEL_energies_size;        // number of working iterations (for tensor Fmunu)
    lattice_energies_offset         = lattice_energies_size;

    lattice_boundary_size           = GPU0->buffer_size_align(lattice_boundary_exact_size);
    lattice_table_row_size          = GPU0->buffer_size_align(lattice_domain_site);
    lattice_table_row_size_half     = GPU0->buffer_size_align(lattice_domain_site / 2);
    lattice_table_exact_row_size    = GPU0->buffer_size_align(lattice_domain_exact_site);
    lattice_table_exact_row_size_half=GPU0->buffer_size_align(lattice_domain_exact_site / 2);
    lattice_table_size              = GPU0->buffer_size_align(lattice_table_row_size * lattice_nd * lattice_group_elements[lattice_group-1] / 4);
    lattice_table_group             = GPU0->buffer_size_align(lattice_table_row_size * lattice_nd);
    lattice_table_exact_group       = GPU0->buffer_size_align(lattice_table_exact_row_size * lattice_nd);

    lattice_measurement_size        = GPU0->buffer_size_align((unsigned int) ceil((double) lattice_table_exact_row_size / workgroup_factor)); // workgroup_factor is the minimum number of workgroup items
    lattice_measurement_size_F      = lattice_measurement_size * MODEL_energies_size;
    lattice_measurement_offset      = lattice_measurement_size;

    lattice_polyakov_size           = GPU0->buffer_size_align((unsigned int) lattice_domain_exact_n1n2n3);
    lattice_polyakov_loop_size      = GPU0->buffer_size_align(ITER);      // number of working iterations
    lattice_polyakov_loop_offset    = GPU0->buffer_size_align((unsigned int) ceil((double) lattice_polyakov_size / workgroup_factor));

    if ((!get_Fmunu)&&(!get_F0mu)){
        lattice_energies_size_F    = lattice_energies_size;
        lattice_measurement_size_F = lattice_measurement_size;
    }

    //_____________________________________________ PRNG preparation
        PRNG0->PRNG_instances   = 0;    // number of instances of generator (or 0 for autoselect)
        // number of samples produced by each generator (quads)
        PRNG0->PRNG_samples     = GPU0->buffer_size_align((unsigned int) ceil(double(3 * lattice_table_row_size_half * (NHIT + 1)))); // 3*NHIT+1 PRNs per link
        PRNG0->GPU0 = GPU0;
        prngstep = lattice_table_row_size_half;
    //-----------------------------------------------------------------

    if (GPU0->GPU_debug.brief_report) {
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
        printf("local_size_intel            = %u\n",local_size_intel);
        printf("kernels: ----------------------------\n");
        printf("lattice_init                = %u\n",lattice_table_group);
        printf("lattice_GramSchmidt         = %u\n",lattice_table_group);
        printf("lattice_measurement         = %u\n",lattice_table_row_size);
        printf("lattice_measurement_plq     = %u\n",lattice_table_row_size);
        printf("lattice_measurement_wilson  = %u\n",lattice_table_row_size);
        printf("update                      = %u\n",lattice_table_row_size_half);
        printf("lattice_polyakov            = %u\n",lattice_polyakov_size);
        printf("PRNs (in quads): --------------------\n");
        printf("PRNG_instances              = %u\n",PRNG0->PRNG_instances);
        printf("PRNG_samples                = %u\n",PRNG0->PRNG_samples);
        printf("PRNs total                  = %u\n",PRNG0->PRNG_instances * PRNG0->PRNG_samples);
        printf("PRNG_step                   = %u\n",prngstep);
        printf("PRNs for hot start          = %u\n",lattice_table_row_size * 3);
        printf("PRNs for one update         = %u\n",int(ceil(double(3 * lattice_domain_link * (NHIT + 1)))));
        printf("-------------------------------------\n");
        printf(" ROWSIZE                    = %u\n",lattice_table_row_size);
        printf(" NHIT                       = %u\n",NHIT);
        printf(" PRNGSTEP                   = %u\n",lattice_table_row_size_half);
        printf(" PRECISION                  = %u\n",precision);
        printf(" PL                         = %u\n",PL_level);
        printf(" ITER                       = %u\n",ITER);
        if (!((PHI==0.0)&&(OMEGA==0.0))) printf(" TBC\n");     // turn on TBC
        if (get_Fmunu) {
            printf(" FMUNU (%u, %u)\n",Fmunu_index1,Fmunu_index2);   // calculate tensor Fmunu
        }
        if (get_F0mu) {
            printf(" F0MU (%u, %u)\n",Fmunu_index1,Fmunu_index2);   // calculate tensor Fmunu
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
}

void        model::lattice_make_programs(void)
{
    int j;
    int argument_id;
    wilson_index      = 0;
    plq_index         = 0;
    polyakov_index    = 0;
    measurement_index = 0;

    char options_common[1024];
    int options_length_common  = sprintf_s(options_common,sizeof(options_common),"-Werror");
    if (GPU0->GPU_info.device_vendor == GPU_CL::GPU::GPU_vendor_Intel)
        options_length_common += sprintf_s(options_common + options_length_common,sizeof(options_common)-options_length_common," -D INTEL_ON");
    options_length_common += sprintf_s(options_common + options_length_common,sizeof(options_common)-options_length_common," -D SUN=%u",         lattice_group);
    options_length_common += sprintf_s(options_common + options_length_common,sizeof(options_common)-options_length_common," -D ND=%u",          lattice_nd);
    options_length_common += sprintf_s(options_common + options_length_common,sizeof(options_common)-options_length_common," -D N1=%u",          lattice_domain_n1);
    options_length_common += sprintf_s(options_common + options_length_common,sizeof(options_common)-options_length_common," -D N2=%u",          lattice_domain_size[1]);
    options_length_common += sprintf_s(options_common + options_length_common,sizeof(options_common)-options_length_common," -D N3=%u",          lattice_domain_size[2]);
    options_length_common += sprintf_s(options_common + options_length_common,sizeof(options_common)-options_length_common," -D N4=%u",          lattice_domain_size[3]);
    if (!((PHI==0.0)&&(OMEGA==0.0))) options_length_common += sprintf_s(options_common + options_length_common,sizeof(options_common)-options_length_common," -D TBC");     // turn on TBC
    options_length_common += sprintf_s(options_common + options_length_common,sizeof(options_common)-options_length_common," -I %s%s",           GPU0->cl_root_path,path_suncl);
    options_length_common += sprintf_s(options_common + options_length_common,sizeof(options_common)-options_length_common," -I %s%s",           GPU0->cl_root_path,path_kernel);
    options_length_common += sprintf_s(options_common + options_length_common,sizeof(options_common)-options_length_common," -D ROWSIZE=%u",     lattice_table_row_size);
    options_length_common += sprintf_s(options_common + options_length_common,sizeof(options_common)-options_length_common," -D PRECISION=%u",   precision);

    char options[1024];
    int options_length  = sprintf_s(options,sizeof(options),"%s",options_common);
        options_length += sprintf_s(options + options_length,sizeof(options)-options_length," -D NHIT=%u",        NHIT);
        options_length += sprintf_s(options + options_length,sizeof(options)-options_length," -D PRNGSTEP=%u",    lattice_table_row_size_half);

    char buffer_update_cl[FNAME_MAX_LENGTH];
        j = sprintf_s(buffer_update_cl  ,FNAME_MAX_LENGTH,  "%s",GPU0->cl_root_path);
        j+= sprintf_s(buffer_update_cl+j,FNAME_MAX_LENGTH-j,"%s",SOURCE_UPDATE);
    char* update_source       = GPU0->source_read(buffer_update_cl);
                                GPU0->program_create(update_source,options);


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


    if (ints==model_start_gid) {            // gid init
                sun_init_id = GPU0->kernel_init("lattice_init_gid",1,init_global_size,NULL);
                argument_id = GPU0->kernel_init_buffer(sun_init_id,lattice_table);
    } else if (ints==model_start_cold) {    // cold init
                sun_init_id = GPU0->kernel_init("lattice_init_cold",1,init_global_size,NULL);
                argument_id = GPU0->kernel_init_buffer(sun_init_id,lattice_table);
    } else {                                // hot init
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

    sun_update_odd_X_id = GPU0->kernel_init("update_odd_X",1,monte_global_size,NULL);
            argument_id = GPU0->kernel_init_buffer(sun_update_odd_X_id,lattice_table);
            argument_id = GPU0->kernel_init_buffer(sun_update_odd_X_id,lattice_parameters);
            argument_id = GPU0->kernel_init_buffer(sun_update_odd_X_id,PRNG0->PRNG_randoms_id);

    sun_update_even_X_id = GPU0->kernel_init("update_even_X",1,monte_global_size,NULL);
             argument_id = GPU0->kernel_init_buffer(sun_update_even_X_id,lattice_table);
             argument_id = GPU0->kernel_init_buffer(sun_update_even_X_id,lattice_parameters);
             argument_id = GPU0->kernel_init_buffer(sun_update_even_X_id,PRNG0->PRNG_randoms_id);

    sun_update_odd_Y_id = GPU0->kernel_init("update_odd_Y",1,monte_global_size,NULL);
            argument_id = GPU0->kernel_init_buffer(sun_update_odd_Y_id,lattice_table);
            argument_id = GPU0->kernel_init_buffer(sun_update_odd_Y_id,lattice_parameters);
            argument_id = GPU0->kernel_init_buffer(sun_update_odd_Y_id,PRNG0->PRNG_randoms_id);

    sun_update_even_Y_id = GPU0->kernel_init("update_even_Y",1,monte_global_size,NULL);
             argument_id = GPU0->kernel_init_buffer(sun_update_even_Y_id,lattice_table);
             argument_id = GPU0->kernel_init_buffer(sun_update_even_Y_id,lattice_parameters);
             argument_id = GPU0->kernel_init_buffer(sun_update_even_Y_id,PRNG0->PRNG_randoms_id);

    sun_update_odd_Z_id = GPU0->kernel_init("update_odd_Z",1,monte_global_size,NULL);
            argument_id = GPU0->kernel_init_buffer(sun_update_odd_Z_id,lattice_table);
            argument_id = GPU0->kernel_init_buffer(sun_update_odd_Z_id,lattice_parameters);
            argument_id = GPU0->kernel_init_buffer(sun_update_odd_Z_id,PRNG0->PRNG_randoms_id);

    sun_update_even_Z_id = GPU0->kernel_init("update_even_Z",1,monte_global_size,NULL);
             argument_id = GPU0->kernel_init_buffer(sun_update_even_Z_id,lattice_table);
             argument_id = GPU0->kernel_init_buffer(sun_update_even_Z_id,lattice_parameters);
             argument_id = GPU0->kernel_init_buffer(sun_update_even_Z_id,PRNG0->PRNG_randoms_id);

    sun_update_odd_T_id = GPU0->kernel_init("update_odd_T",1,monte_global_size,NULL);
            argument_id = GPU0->kernel_init_buffer(sun_update_odd_T_id,lattice_table);
            argument_id = GPU0->kernel_init_buffer(sun_update_odd_T_id,lattice_parameters);
            argument_id = GPU0->kernel_init_buffer(sun_update_odd_T_id,PRNG0->PRNG_randoms_id);

    sun_update_even_T_id = GPU0->kernel_init("update_even_T",1,monte_global_size,NULL);
             argument_id = GPU0->kernel_init_buffer(sun_update_even_T_id,lattice_table);
             argument_id = GPU0->kernel_init_buffer(sun_update_even_T_id,lattice_parameters);
             argument_id = GPU0->kernel_init_buffer(sun_update_even_T_id,PRNG0->PRNG_randoms_id);

    // for all measurements _____________________________________________________________________________________________________________________________________
    char options_measurements[1024];
    int options_measurement_length  = sprintf_s(options_measurements,sizeof(options_measurements),"%s",options_common);

    if (get_Fmunu) options_measurement_length += sprintf_s(options_measurements + options_measurement_length,sizeof(options_measurements)-options_measurement_length," -D FMUNU");   // calculate tensor Fmunu for H field
    if (get_F0mu)  options_measurement_length += sprintf_s(options_measurements + options_measurement_length,sizeof(options_measurements)-options_measurement_length," -D F0MU");   // calculate tensor Fmunu for E field

    if ((get_Fmunu)||(get_F0mu)) {
        if (get_Fmunu1) options_measurement_length += sprintf_s(options_measurements + options_measurement_length,sizeof(options_measurements)-options_measurement_length," -D FMUNU1");   // calculate tensor Fmunu for lambda1 matrix
        if (get_Fmunu2) options_measurement_length += sprintf_s(options_measurements + options_measurement_length,sizeof(options_measurements)-options_measurement_length," -D FMUNU2");   // calculate tensor Fmunu for lambda2 matrix
        if (get_Fmunu4) options_measurement_length += sprintf_s(options_measurements + options_measurement_length,sizeof(options_measurements)-options_measurement_length," -D FMUNU4");   // calculate tensor Fmunu for lambda4 matrix
        if (get_Fmunu5) options_measurement_length += sprintf_s(options_measurements + options_measurement_length,sizeof(options_measurements)-options_measurement_length," -D FMUNU5");   // calculate tensor Fmunu for lambda5 matrix
        if (get_Fmunu6) options_measurement_length += sprintf_s(options_measurements + options_measurement_length,sizeof(options_measurements)-options_measurement_length," -D FMUNU6");   // calculate tensor Fmunu for lambda6 matrix
        if (get_Fmunu7) options_measurement_length += sprintf_s(options_measurements + options_measurement_length,sizeof(options_measurements)-options_measurement_length," -D FMUNU7");   // calculate tensor Fmunu for lambda7 matrix
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
    if ((get_plaquettes_avr)||(get_Fmunu)||(get_F0mu)) {
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
    }

    int lattice_measurement_size_correction = max(size_reduce_measurement_double2,size_reduce_measurement_plq_double2);

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
    if (get_wilson_loop) {
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
        options_length_polyakov += sprintf_s(options_polyakov + options_length_polyakov,sizeof(options_polyakov)-options_length_polyakov," -D PL=%u",          PL_level);

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
    if (PL_level > 0) {
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
}

void        model::lattice_create_buffers(void)
{
    cl_float*    plattice_parameters_float  = NULL;
    cl_double*   plattice_parameters_double = NULL;

    int fc = 2;
    int fc2 = ((PL_level>1)&&(lattice_measurement_size_F <= 2 * lattice_polyakov_loop_size)) ? 2 : 1;
    size_lattice_table         = fc * lattice_table_size;
    size_lattice_measurement   = fc * lattice_measurement_size_F;
    size_lattice_energies      = fc * lattice_energies_size;
    size_lattice_wilson_loop   = fc * lattice_energies_size;
    size_lattice_energies_plq  = fc * lattice_energies_offset * MODEL_energies_size;
    size_lattice_polyakov_loop = fc * fc2 * lattice_polyakov_loop_size * PL_level;
    size_lattice_boundary      = fc * lattice_boundary_size;
    size_lattice_parameters    = fc * lattice_parameters_size;

    plattice_table_float    = NULL;
    plattice_table_double   = NULL;
    plattice_boundary_float = NULL;
    plattice_boundary_double= NULL;

    plattice_measurement    = (cl_double2*) calloc(size_lattice_measurement,  sizeof(cl_double2));
    plattice_energies       = (cl_double2*) calloc(size_lattice_energies,     sizeof(cl_double2));
    plattice_energies_plq   = (cl_double2*) calloc(size_lattice_energies_plq, sizeof(cl_double2));
    plattice_wilson_loop    = (cl_double*)  calloc(size_lattice_wilson_loop,  sizeof(cl_double));
    plattice_polyakov_loop  = NULL;
        if (PL_level > 0) plattice_polyakov_loop = (cl_double2*) calloc(size_lattice_polyakov_loop,sizeof(cl_double2));

    if (precision == model_precision_single) {
        plattice_table_float            = (cl_float4*)  calloc(size_lattice_table,     sizeof(cl_float4));
        plattice_boundary_float         = (cl_float4*)  calloc(size_lattice_boundary,  sizeof(cl_float4));
        plattice_parameters_float       = (cl_float*)   calloc(size_lattice_parameters,sizeof(cl_float));
        plattice_parameters_float[0]    = (float) (BETA / lattice_group);
        plattice_parameters_float[1]    = (float) (PHI);
        plattice_parameters_float[2]    = (float) (OMEGA);
        plattice_parameters_float[3]    = (float) (wilson_R);
        plattice_parameters_float[4]    = (float) (wilson_T);
    } else {
        plattice_table_double           = (cl_double4*) calloc(size_lattice_table,     sizeof(cl_double4));
        plattice_boundary_double        = (cl_double4*) calloc(size_lattice_boundary,  sizeof(cl_double4));
        plattice_parameters_double      = (cl_double*)  calloc(size_lattice_parameters,sizeof(cl_double));
        plattice_parameters_double[0]   = (double) (BETA / lattice_group);
        plattice_parameters_double[1]   = (double) (PHI);
        plattice_parameters_double[2]   = (double) (OMEGA);
        plattice_parameters_double[3]   = (double) (wilson_R);
        plattice_parameters_double[4]   = (double) (wilson_T);
    }

    if (INIT==0) lattice_load_state();  // load state file if needed

    int lds_size = (int) GPU0->GPU_info.local_memory_size / 4 / 4 / 2;  // 4 bytes per component, 4 components, double2 type <- use all local memory

    lattice_table = 0;
    if (precision == model_precision_single) {
        lattice_table           = GPU0->buffer_init(GPU0->buffer_type_IO, size_lattice_table,            plattice_table_float,       sizeof(cl_float4));  // Lattice data
        lattice_boundary        = GPU0->buffer_init(GPU0->buffer_type_IO, size_lattice_boundary,         plattice_boundary_float,    sizeof(cl_float4));  // Lattice boundary
        lattice_parameters      = GPU0->buffer_init(GPU0->buffer_type_Constant, size_lattice_parameters, plattice_parameters_float,  sizeof(cl_float));   // Lattice counters and indices
    } else {
        lattice_table           = GPU0->buffer_init(GPU0->buffer_type_IO, size_lattice_table,            plattice_table_double,      sizeof(cl_double4)); // Lattice data
        lattice_boundary        = GPU0->buffer_init(GPU0->buffer_type_IO, size_lattice_boundary,         plattice_boundary_double,   sizeof(cl_double4)); // Lattice boundary
        lattice_parameters      = GPU0->buffer_init(GPU0->buffer_type_Constant, size_lattice_parameters, plattice_parameters_double, sizeof(cl_double));  // Lattice counters and indices
    }
    lattice_measurement         = GPU0->buffer_init(GPU0->buffer_type_IO, size_lattice_measurement,      plattice_measurement,       sizeof(cl_double2)); // Lattice measurement
    lattice_lds                 = GPU0->buffer_init(GPU0->buffer_type_LDS,lds_size,                      0,                          sizeof(cl_double2)); // LDS for reduction
    lattice_energies            = GPU0->buffer_init(GPU0->buffer_type_IO, size_lattice_energies,         plattice_energies,          sizeof(cl_double2)); // Lattice energies
    if ((get_plaquettes_avr) || (get_Fmunu) || (get_F0mu))
        lattice_energies_plq    = GPU0->buffer_init(GPU0->buffer_type_IO, size_lattice_energies_plq,     plattice_energies_plq,      sizeof(cl_double2)); // Lattice energies (plaquettes)
    if (get_wilson_loop)
        lattice_wilson_loop     = GPU0->buffer_init(GPU0->buffer_type_IO, size_lattice_wilson_loop,      plattice_wilson_loop,       sizeof(cl_double));  // Wilson loop
    if (PL_level > 0)
        lattice_polyakov_loop   = GPU0->buffer_init(GPU0->buffer_type_IO, size_lattice_polyakov_loop,    plattice_polyakov_loop, sizeof(cl_double2)); // Polyakov loops

    plattice_measurement    = (cl_double2*) calloc(size_lattice_measurement,  sizeof(cl_double2));
    plattice_energies       = (cl_double2*) calloc(size_lattice_energies,     sizeof(cl_double2));
    plattice_energies_plq   = (cl_double2*) calloc(size_lattice_energies_plq, sizeof(cl_double2));
    plattice_wilson_loop    = (cl_double*)  calloc(size_lattice_wilson_loop,  sizeof(cl_double));
    plattice_polyakov_loop  = NULL;
        if (PL_level > 0) plattice_polyakov_loop = (cl_double2*) calloc(size_lattice_polyakov_loop,sizeof(cl_double2));
}

void        model::lattice_simulate(void)
{
    lattice_create_buffers();
    lattice_make_programs();


    rowsize      = lattice_table_row_size;
    rowsize4     = 4 * lattice_table_row_size;
    halfrowsize  = lattice_domain_link / 8;
    halfrowsize4 = lattice_domain_link / 2;


    // simulations ______________________________________________________________________________________________________________________________________________
    GPU0->print_memory_utilized();

    timestart = GPU0->get_current_datetime();
    time(&ltimestart);

    lattice_pointer_initial      = NULL;
    lattice_pointer_measurements = NULL;

    printf("\nrun kernels on GPU (%f seconds)\n",GPU0->get_timer_CPU(TIMER_FOR_ELAPSED));
    GPU0->start_timer_CPU(TIMER_FOR_SIMULATIONS); // start GPU execution timer
    GPU0->start_timer_CPU(TIMER_FOR_SAVE);        // start timer for lattice_state save
    if (ints==model_start_hot) PRNG0->produce();

    int NAV_start  = 0;
    int ITER_start = 0;

    if (INIT!=0) {

        if (ints==model_start_hot) {
            GPU0->kernel_run(sun_init_X_id);        // Lattice initialization
            if (!turnoff_prns) PRNG0->produce();
            GPU0->kernel_run(sun_init_Y_id);        // Lattice initialization
            if (!turnoff_prns) PRNG0->produce();
            GPU0->kernel_run(sun_init_Z_id);        // Lattice initialization
            if (!turnoff_prns) PRNG0->produce();
            GPU0->kernel_run(sun_init_T_id);        // Lattice initialization
            if (!turnoff_prns) PRNG0->produce();
        } else
            GPU0->kernel_run(sun_init_id);          // Lattice initialization
        GPU0->print_stage("lattice initialized");

        if ((get_plaquettes_avr)||(get_Fmunu)||(get_F0mu)) {
            GPU0->kernel_run(sun_measurement_plq_id);          // Lattice measurement (plaquettes)
            GPU0->print_stage("measurement done (plaquettes)");
                plq_index = ITER_counter;
                GPU0->kernel_init_constant_reset(sun_measurement_plq_reduce_id,&plq_index,argument_plq_index);
            GPU0->kernel_run(sun_measurement_plq_reduce_id);    // Lattice measurement reduction (plaquettes)
            GPU0->print_stage("measurement reduce done (plaquettes)");
        }
    
        if (get_wilson_loop) {
            GPU0->kernel_run(sun_measurement_wilson_id);       // Lattice Wilson loop measurement
            GPU0->print_stage("Wilson loop measurement done");
                wilson_index = ITER_counter;
                GPU0->kernel_init_constant_reset(sun_wilson_loop_reduce_id,&wilson_index,argument_wilson_index);
            GPU0->kernel_run(sun_wilson_loop_reduce_id);       // Lattice Wilson loop measurement reduction
            GPU0->print_stage("Wilson loop measurement reduce done");
        }

        if (get_actions_avr) {
            GPU0->kernel_run(sun_measurement_id);                  // Lattice measurement
            GPU0->print_stage("measurement done");
                measurement_index = ITER_counter;
                GPU0->kernel_init_constant_reset(sun_measurement_reduce_id,&measurement_index,argument_measurement_index);
            GPU0->kernel_run(sun_measurement_reduce_id);           // Lattice measurement reduction
            GPU0->print_stage("measurement reduce done");
        }

        if (PL_level > 0) {
            GPU0->kernel_run(sun_polyakov_id);                     // Lattice Polyakov loop measurement
            GPU0->print_stage("Polyakov loop measurement done");
                polyakov_index = ITER_counter;
                GPU0->kernel_init_constant_reset(sun_polyakov_reduce_id,&polyakov_index,argument_polyakov_index);
            GPU0->kernel_run(sun_polyakov_reduce_id);              // Lattice Polyakov loop measurement reduction
            GPU0->print_stage("Polyakov loop reduce done");
        }

        ITER_counter++;

        if (!turnoff_config_save) lattice_save_state();
    }
    lattice_pointer_initial = GPU0->buffer_map(lattice_table);

    NAV_start  = NAV_counter;
    ITER_start = ITER_counter;

    if (INIT==0) {
        while (PRNG_counter>PRNG0->PRNG_counter) PRNG0->produce();    // adjust PRNG
        if (GPU0->GPU_debug.brief_report) printf("NAV_start=%u, ITER_start=%u\n",NAV_start,ITER_start);

        wilson_index = ITER_start + 1;
    }


    // perform thermalization
    for (int i=NAV_start; i<NAV; i++){
           if (!turnoff_prns) PRNG0->produce();
        if (!turnoff_updates) GPU0->kernel_run(sun_update_odd_X_id);     // Update odd X links
           if (!turnoff_prns) PRNG0->produce();
        if (!turnoff_updates) GPU0->kernel_run(sun_update_odd_Y_id);     // Update odd Y links
           if (!turnoff_prns) PRNG0->produce();
        if (!turnoff_updates) GPU0->kernel_run(sun_update_odd_Z_id);     // Update odd Z links
           if (!turnoff_prns) PRNG0->produce();
        if (!turnoff_updates) GPU0->kernel_run(sun_update_odd_T_id);     // Update odd T links

           if (!turnoff_prns) PRNG0->produce();
        if (!turnoff_updates) GPU0->kernel_run(sun_update_even_X_id);    // Update even X links
           if (!turnoff_prns) PRNG0->produce();
        if (!turnoff_updates) GPU0->kernel_run(sun_update_even_Y_id);    // Update even Y links
           if (!turnoff_prns) PRNG0->produce();
        if (!turnoff_updates) GPU0->kernel_run(sun_update_even_Z_id);    // Update even Z links
           if (!turnoff_prns) PRNG0->produce();
        if (!turnoff_updates) GPU0->kernel_run(sun_update_even_T_id);    // Update even T links

            if (!turnoff_gramschmidt)
                GPU0->kernel_run(sun_GramSchmidt_id);          // Lattice reunitarization

        if (i % 10 == 0) printf("\rGPU thermalization [%i]",i);
        NAV_counter++;

        // write lattice state every [write_lattice_state_every_secs] seconds
        if ((!turnoff_config_save)&&(GPU0->timer_in_seconds_CPU(TIMER_FOR_SAVE)>write_lattice_state_every_secs)) {
            lattice_save_state();
            printf("\ncurrent configuration saved\n");
            GPU0->start_timer_CPU(TIMER_FOR_SAVE);      // restart timer for lattice_state save
        } 
    }

    // perform working cycles
    for (int i=ITER_start; i<ITER; i++){ // zero measurement - on initial configuration!
        for (int j=0; j<NITER; j++){
               if (!turnoff_prns) PRNG0->produce();
            if (!turnoff_updates) GPU0->kernel_run(sun_update_odd_X_id);     // Lattice measurement staples
               if (!turnoff_prns) PRNG0->produce();
            if (!turnoff_updates) GPU0->kernel_run(sun_update_odd_Y_id);     // Lattice measurement staples
               if (!turnoff_prns) PRNG0->produce();
            if (!turnoff_updates) GPU0->kernel_run(sun_update_odd_Z_id);     // Lattice measurement staples
               if (!turnoff_prns) PRNG0->produce();
            if (!turnoff_updates) GPU0->kernel_run(sun_update_odd_T_id);     // Lattice measurement staples
               if (!turnoff_prns) PRNG0->produce();

            if (!turnoff_updates) GPU0->kernel_run(sun_update_even_X_id);    // Lattice measurement staples
               if (!turnoff_prns) PRNG0->produce();
            if (!turnoff_updates) GPU0->kernel_run(sun_update_even_Y_id);    // Lattice measurement staples
               if (!turnoff_prns) PRNG0->produce();
            if (!turnoff_updates) GPU0->kernel_run(sun_update_even_Z_id);    // Lattice measurement staples
               if (!turnoff_prns) PRNG0->produce();
            if (!turnoff_updates) GPU0->kernel_run(sun_update_even_T_id);    // Lattice measurement staples

            if (!turnoff_gramschmidt)
                GPU0->kernel_run(sun_GramSchmidt_id);           // Lattice reunitarization

            if (i % 10 == 0) printf("\rGPU working iteration [%u]",i);
        }

        if (get_wilson_loop) {
            GPU0->kernel_run(sun_measurement_wilson_id);        // Lattice Wilson loop measurement
                wilson_index = ITER_counter;
                GPU0->kernel_init_constant_reset(sun_wilson_loop_reduce_id,&wilson_index,argument_wilson_index);
            GPU0->kernel_run(sun_wilson_loop_reduce_id);        // Lattice Wilson loop measurement reduction
        }

        // measurements
        if (get_actions_avr) {
            GPU0->kernel_run(sun_measurement_id);                 // Lattice measurement
                measurement_index = ITER_counter;
                GPU0->kernel_init_constant_reset(sun_measurement_reduce_id,&measurement_index,argument_measurement_index);
            GPU0->kernel_run(sun_measurement_reduce_id);          // Lattice measurement reduction
        }

GPU0->kernel_run(sun_clear_measurement_id);
        if (PL_level > 0) {
            GPU0->kernel_run(sun_polyakov_id);                    // Lattice Polyakov loop measurement
                polyakov_index = ITER_counter;
                GPU0->kernel_init_constant_reset(sun_polyakov_reduce_id,&polyakov_index,argument_polyakov_index);
            GPU0->kernel_run(sun_polyakov_reduce_id);             // Lattice Polyakov loop measurement reduction
        }

GPU0->kernel_run(sun_clear_measurement_id);
        // measurements (plaquettes)
        if ((get_plaquettes_avr)||(get_Fmunu)||(get_F0mu)) {
            GPU0->kernel_run(sun_measurement_plq_id);           // Lattice measurement (plaquettes)
                plq_index = ITER_counter;
                GPU0->kernel_init_constant_reset(sun_measurement_plq_reduce_id,&plq_index,argument_plq_index);
            GPU0->kernel_run(sun_measurement_plq_reduce_id);    // Lattice measurement reduction (plaquettes)
        }

        ITER_counter++;

        // write lattice state every [write_lattice_state_every_secs] seconds
        if ((!turnoff_config_save)&&(GPU0->timer_in_seconds_CPU(TIMER_FOR_SAVE)>write_lattice_state_every_secs)) {
            lattice_save_state();
            printf("\ncurrent configuration saved\n");
            GPU0->start_timer_CPU(TIMER_FOR_SAVE);      // restart timer for lattice_state save
        } 
    }
    printf("\rGPU simulations are done (%f seconds)\n",GPU0->get_timer_CPU(1));
    time(&ltimeend);
    timeend   = GPU0->get_current_datetime();

    if (!turnoff_config_save) {
        lattice_save_state();
        lattice_pointer_last = lattice_pointer_save;
    } else {
        lattice_pointer_last = GPU0->buffer_map(lattice_table);
    }
    prng_pointer = GPU0->buffer_map_float4(PRNG0->PRNG_randoms_id);
}

}
