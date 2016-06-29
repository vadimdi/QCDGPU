/******************************************************************************
 * @file     QCDGPU.cpp
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.6
 *
 * @brief    [QCDGPU]
 *           The main program file
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
#ifndef QCDGPU_CPP
#define QCDGPU_CPP
 
#include "QCDGPU.h"

char version[] = "1.5";     /**< version of programm */
char fil[]     = "c:\\init.dat"; /**< init.dat file loading by default */
char fprefix[] = "su3-";     /**< prefix for data-files */

char fload[] = "su3-14-Aug-20-00-00-00.qcg"; /**< file with configuration to be load by default*/

#ifdef _WIN32
            char path[]   = "C:/DATA/TEST/"; /**< path for output file *//*< path for output file (corresponding to Windows/Linux OS) **/
#else
            char path[]   = "/home/science/scalar/";        
#endif

#define PI 3.1415926535897932384626433832795

/** !!!!!!!!!!!!!!!!!!!!
 * @param argc an integer argument
 * @param argv an array of character pointers
 * @return 0
 */
int main(int argc, char ** argv)
{
    int result = 0;
#ifndef CPU_RUN
    using  GPU_CL::GPU;
#endif
    using  model_CL::model;
#ifdef BIGLAT
    using  model_CL::SubLattice;
#endif

    model* model0 = new(model);
#ifndef CPU_RUN
model0->GPU0->GPU_limit_max_workgroup_size = 64; // please, remark it for Intel platform
//#ifndef WIN32
//setenv("CUDA_CACHE_DISABLE", "1", 1); // disable nVidia cach, if needed
//#endif

    // setup debug_level
    model0->GPU0->GPU_debug.local_run             = true;
    if (argc>1) model0->GPU0->GPU_debug.local_run = false;
    int local_run = model0->GPU0->GPU_debug.local_run;

    if (model0->GPU0->GPU_debug.local_run) {
#ifdef _WIN32
        model0->GPU0->GPU_debug.wait_for_keypress   = true;
#endif
        model0->GPU0->GPU_debug.profiling           = false;
        model0->GPU0->GPU_debug.brief_report        = true;
    }

    
        model0->desired_platform = 0;
        model0->desired_device   = 0;
        model0->device_select    = true;
        model0->GPU0->GPU_debug.rebuild_binary = false;
#endif
        model0->PRNG0->PRNG_generator   = PRNG_CL::PRNG::PRNG_generator_RANLUX3;
		model0->PRNG0->PRNG_randseries = 54;//1432996491; // constant random series; 0 <--> system time
        model0->PRNG0->PRNG_precision   = PRNG_CL::PRNG::PRNG_precision_single;
        model0->lattice_nd      = 4;
        model0->lattice_group   = 3;       // SU(N) group
		model0->lattice_full_size[0] = 32;  // Nx
		model0->lattice_full_size[1] = 32;// 38;  // Ny
		model0->lattice_full_size[2] = 32;// 38;  // Nz
		model0->lattice_full_size[3] = 32;// 38;  // Nt

        //model0->big_lattice         = false;
                    model0->lattice_domain_size[0] = model0->lattice_full_size[0];
                    model0->lattice_domain_size[1] = model0->lattice_full_size[1];
                    model0->lattice_domain_size[2] = model0->lattice_full_size[2];
                    model0->lattice_domain_size[3] = model0->lattice_full_size[3];
#ifdef BIGLAT
                    model0->lattice_Nparts = NPARTS;

                    model0->SubLat = (SubLattice*)calloc(model0->lattice_Nparts, sizeof(SubLattice));
			model0->SubLat[0].Nx = 18;
			model0->SubLat[1].Nx = 14;
			//model0->SubLat[2].Nx = 8;
			//model0->SubLat[3].Nx = 8;

                    for(int k = 0; k < model0->lattice_Nparts; k++)
                    {
                        model0->SubLat[k].Ny = model0->lattice_full_size[1];
                        model0->SubLat[k].Nz = model0->lattice_full_size[2];
                        model0->SubLat[k].Nt = model0->lattice_full_size[3];
                    }

                    model0->SubLat[0].desired_platform = 0;
                    model0->SubLat[0].desired_device   = 0;
					//model0->SubLat[0].device_no = 0;
                        model0->SubLat[1].desired_platform = 0;
                        model0->SubLat[1].desired_device   = 0;
						//model0->SubLat[1].device_no = 0;
                    /*model0->SubLat[2].desired_platform = 0;
                    model0->SubLat[2].desired_device   = 0;
					//model0->SubLat[0].device_no = 1;
                        model0->SubLat[3].desired_platform = 0;
                        model0->SubLat[3].desired_device   = 0;
						//model0->SubLat[0].device_no = 1;*/

						//model0->Ndevices = 2;
                        model0->Ndevices = NDEV;
						model0->devParts = (unsigned int*)calloc(model0->Ndevices, sizeof(unsigned int));
                                                for(int k = 0; k < model0->Ndevices; k++)
						//for(int k = 0; k < model0->lattice_Nparts; k++)
						  model0->devParts[k] = 1;
						//model0->devParts[0] = 1;
						//model0->devParts[1] = 1;
						//model0->devParts[2] = 1;
						//model0->devParts[3] = 1;

                        for(int k = 0; k < model0->lattice_Nparts; k++)
                        {
                            model0->SubLat[k].device_select    = true;
                            model0->SubLat[k].GPU0 = new(GPU_CL::GPU);

                            model0->SubLat[k].GPU0->GPU_debug.rebuild_binary = false;
                            model0->SubLat[k].GPU0->GPU_debug.local_run      = local_run;

                            model0->SubLat[k].GPU0->GPU_limit_max_workgroup_size = model0->GPU0->GPU_limit_max_workgroup_size; 

                            if (model0->SubLat[k].GPU0->GPU_debug.local_run) {
#ifdef _WIN32
                                model0->SubLat[k].GPU0->GPU_debug.wait_for_keypress   = true;
#endif
                            }

                            model0->SubLat[k].PRNG0 = new(PRNG_CL::PRNG);
                            model0->SubLat[k].PRNG0->PRNG_generator   = PRNG_CL::PRNG::PRNG_generator_RANLUX3;
	                        model0->SubLat[k].PRNG0->PRNG_randseries  = model0->PRNG0->PRNG_randseries;
                            model0->SubLat[k].PRNG0->PRNG_precision   = PRNG_CL::PRNG::PRNG_precision_single;
                        }
						//model0->SubLat[2].GPU0->GPU_limit_max_workgroup_size = 64;
						//model0->SubLat[3].GPU0->GPU_limit_max_workgroup_size = 64;
#endif

        //model0->precision     = model::model_precision_single;
						model0->precision = model::model_precision_double;
        model0->INIT          = 1;                  // start simulations (0 - continue, 1 - start)
        model0->NAV           = 0;                  // number of thermalization cycles
        model0->ITER          = 2;                // number of working iterations (note: first measurement is performing on the initial configuration!!!)
        model0->NITER         = 1;                  // Number of bulk iterations between measurements
        model0->NHIT          = 10;
        model0->NHITPar       = 1;
        model0->BETA          = 6.5;
        model0->PHI           = 0.0;//2.0E-6;    // lambda_3 (for SU(2) and SU(3) groups)
        model0->OMEGA         = 0.0;    // lambda_8 (for SU(3) group)
        model0->ints          = model::model_start_cold;  // setup start type (model_start_cold or model_start_hot)
#ifndef CPU_RUN
        model0->PL_level      = 3;                  // level of Polyakov loop calculation
                                                    // PL_level = 0 - do not calculate PL
                                                    // PL_level = 1 - calculate PL only
                                                    // PL_level = 2 - calculate PL, PL^2, PL^4
                                                    // PL_level = 3 - calculate PL_diff
#endif
        model0->wilson_R      = 7;
        model0->wilson_T      = 2;
        model0->version       = (char*) calloc((strlen(version) + 1),sizeof(char));
        strcpy_s(model0->version,(strlen(version) + 1),version);

model0->turnoff_updates     = false;                // turn off lattice updates
model0->turnoff_config_save = true;                 // do not write configurations
model0->turnoff_prns        = false;                // turn off prn production
model0->turnoff_gramschmidt = false;                // turn off orthogonalization
model0->get_actions_avr     = true;                 // calculate mean action values
model0->get_plaquettes_avr  = true;                 // calculate mean plaquette values
#ifndef CPU_RUN
model0->get_wilson_loop     = true;                // calculate Wilson loop values
model0->get_Fmunu           = true;                 // calculate Fmunu tensor for H field (3-rd and 8-ht by default)
model0->get_F0mu            = false;                // calculate Fmunu tensor for E field (3-rd and 8-ht by default)
model0->get_Fmunu1          = false;                // flag to measure chomomagnetic field strength corresponding to the 1-st SU(3) generator
model0->get_Fmunu5          = false;                // --"-- for 5-th group generator

model0->write_lattice_state_every_secs = 2.0;       // write lattice configuration every ... seconds (default: 15 minutes)
#endif

model0->get_actions_diff    = false;

model0->check_prngs         = true;   // check PRNG production

        model0->finishpath = (char*) calloc((strlen(path) + 1),sizeof(char));  // setup path to start.txt and finish.txt files
        strcpy_s(model0->finishpath,(strlen(path) + 1),path);

        model0->path = (char*) calloc((strlen(path) + 1),sizeof(char));        // setup path to output files
        strcpy_s(model0->path,(strlen(path) + 1),path);

        model0->fprefix = (char*) calloc((strlen(fprefix) + 1),sizeof(char));  // setup prefix for output files
		if(model0->lattice_group == 3)
                        strcpy_s(model0->fprefix,(strlen(fprefix) + 1),fprefix);
		if(model0->lattice_group == 2)
			strcpy_s(model0->fprefix,(strlen(fprefix) + 1),"su2-");

        model0->fstate = (char*) calloc((strlen(fload) + 1),sizeof(char));  // setup filename for state file
        strcpy_s(model0->fstate,(strlen(fload) + 1),fload);

        if (argc>1) {
           // get command line parameters
           model0->lattice_get_init_file(argv[1]);
        }

	time_t timer1, timer2;
	time(&timer1);
        
#ifdef CPU_RUN
        if(model0->lattice_group == 2)
            lattice_simulateCPU<su_2>(model0, NULL);
        if(model0->lattice_group == 3)
            lattice_simulateCPU<su_3>(model0, NULL);
#else
    model0->lattice_init();
    model0->lattice_simulate();           // MC simulations on GPU
#ifndef BIGLAT
    LatticeCheckCPU(model0);
#endif
    model0->lattice_analysis();
#endif
    model0->lattice_write_results();
    model0->lattice_print_measurements();

	time(&timer2);
	printf("\nEXECUTION TIME: %f sec\n", difftime(timer2, timer1));

    delete model0;

    printf("\n[...Program is completed successfully...]\n\n");

    return result;
}

#endif
