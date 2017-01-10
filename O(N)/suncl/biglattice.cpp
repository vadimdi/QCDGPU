/******************************************************************************
 * @file     biglattice.cpp
 * @author   Vadim Demchik <vadimdi@yahoo.com>
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.0
 *
 * @brief    [QCDGPU]
 *           Dividing lattice into several parts
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
#include "biglattice.h"

namespace BIG_LAT{
using BIG_LAT::BL;
using model_CL::model;
using SUN_CPU::SU;

#define ND_MAX                32  // maximum dimensions


                BL::BL(void){
        big_lattice_parts      = 1; // default number of sublattices
        compute_devices_number = 1; // default number of compute devices
        single_device     = true;   // simulation on single device
        global_run        = new model_CL::model::run_parameters;
}
                BL::~BL(void){

    printf("compute_devices_number = %u\n",compute_devices_number);
    for (int i=0;i<compute_devices_number;i++){
        if (compute_devices[i]) delete compute_devices[i];
        if (models[i]) delete models[i];
        if (SUNcpu[i]) delete SUNcpu[i];
    }

    for (int i=0;i<big_lattice_parts;i++){
        if (lattice_data[i]) delete lattice_data[i];
    }
    if (compute_devices) delete[] compute_devices;
    if (models) delete[] models;
    if (SUNcpu) delete[] SUNcpu;
    if (lattice_data) delete[] lattice_data;
    if (global_run) delete global_run;
}
                BL::lattice_data_buffers::lattice_data_buffers(void){
    plattice_table_float  = NULL;
    plattice_table_double = NULL;
}
                BL::lattice_data_buffers::~lattice_data_buffers(void){

                  
}
                BL::lattice_devices::lattice_devices(void){
        lattice_domain_size = new int[ND_MAX];
        for (int i=0;i<ND_MAX;i++) lattice_domain_size[i] = 0;
}
                BL::lattice_devices::~lattice_devices(void){
        if (lattice_domain_size) delete lattice_domain_size;
}

char*           BL::str_parameter_init(char* str_source){
       char* str_destination = (char*) calloc((strlen(str_source) + 1),sizeof(char));
       strcpy_s(str_destination,(strlen(str_source) + 1),str_source);
       return str_destination;
}
void            BL::prepare(void){
    if (!big_lattice_parts) {
        // autoselection of number of sublattices

        big_lattice_parts = 1;
    }

    if (!compute_devices_number) {
        // autoselection of number of compute devices

        compute_devices_number = 1;
    }

    lattice_data    = new lattice_data_buffers*[big_lattice_parts];
    compute_devices = new lattice_devices*[compute_devices_number];
    models          = new model_CL::model*[compute_devices_number];
    SUNcpu          = new SUN_CPU::SU*[compute_devices_number];

    for (int i=0;i<compute_devices_number;i++){
        compute_devices[i] = new BL::lattice_devices;
    }
}
void            BL::init(void){
    int idx=0;  // index in lattice_data array

    if (compute_devices_number==1) {    // if only one computing device in the system
        compute_devices[0]->lattice_parts = big_lattice_parts;
        compute_devices[0]->performance = 1.0;
        if (compute_devices[0]->lattice_parts==1)
            compute_devices[0]->lattice_domain_size[0] = global_run->lattice_full_size[0];
    }
        
    for (int i=0;i<compute_devices_number;i++){
        for (unsigned int j=0;j<compute_devices[i]->lattice_parts;j++){
            lattice_data[idx] = new lattice_data_buffers;
            lattice_data[idx]->device_select = compute_devices[i]->device_select;
            lattice_data[idx]->host          = compute_devices[i]->host;
            lattice_data[idx]->platform      = compute_devices[i]->platform;
            lattice_data[idx]->device        = compute_devices[i]->device;
            lattice_data[idx]->one_device_one_part = (compute_devices[i]->lattice_parts>1) ? false : true;
            lattice_data[idx]->models_index = i;
            lattice_data[idx]->part_number = idx;
            idx++;
        }
    }
    for (int i=0;i<compute_devices_number;i++){
        models[i] = new model_CL::model;

        GPU_CL::GPU::copy_debug_flags(global_run->GPU_debug,models[i]->GPU0->GPU_debug);

        SUNcpu[i] = new SUN_CPU::SU;   // for CPU check

        model_CL::model::run_parameters::run_parameters_copy(global_run,models[i]->run);

        models[i]->big_lattice = (big_lattice_parts>1) ? true : false;

        models[i]->run->number_of_parts = compute_devices[0]->lattice_parts;

        models[i]->run->device_select    = compute_devices[i]->device_select;
        models[i]->run->desired_platform = compute_devices[i]->platform;
        models[i]->run->desired_device   = compute_devices[i]->device;

        if (compute_devices[i]->performance>0){

        }

        if (compute_devices[i]->lattice_domain_size[0]==0)
            compute_devices[i]->lattice_domain_size[0] = global_run->lattice_full_size[0] / compute_devices[i]->lattice_parts;

        for (int j=0;j<global_run->lattice_nd;j++){
            models[i]->run->lattice_full_size[j]   = global_run->lattice_full_size[j];
            if (compute_devices[i]->lattice_domain_size[j]==0) 
                compute_devices[i]->lattice_domain_size[j] = global_run->lattice_full_size[j];
            models[i]->run->lattice_domain_size[j] = compute_devices[i]->lattice_domain_size[j];
            if (big_lattice_parts==1) models[i]->run->lattice_domain_size[j] = global_run->lattice_full_size[j];
        }
        models[i]->lattice_init();
    }
}

void            BL::simulate(void){
    // prepare shadow buffers for lattice_data[i]
    for (int i=0;i<big_lattice_parts;i++){
        int idx = lattice_data[i]->models_index;
        if (!lattice_data[i]->one_device_one_part) lattice_data[i]->size_lattice_table = models[idx]->size_lattice_table;
    }

    big_lattice_green_passes = (big_lattice_parts/compute_devices_number)>>1;
    big_lattice_red_passes   = (big_lattice_parts/compute_devices_number) - big_lattice_green_passes;

    // lattice initialization
    for (int i=0;i<big_lattice_red_passes;  i++) {
        for (int dv=0;dv<compute_devices_number;dv++){
            int i_part = 2*(dv + i * compute_devices_number);
            int idx = lattice_data[i_part]->models_index;
            if (models[idx]->run->INIT!=0) {
                models[idx]->run->part_number = lattice_data[i_part]->part_number;
                models[idx]->lattice_simulate_start();
            } else {
              
            }
        }
        lattice_wait_queue_red(i);
        lattice_wait_table_read_red(i);
    }
    for (int i=0;i<big_lattice_green_passes;i++) {
        for (int dv=0;dv<compute_devices_number;dv++){
            int i_part = 2*(dv + i * compute_devices_number) + 1;
            int idx = lattice_data[i_part]->models_index;
            if (models[idx]->run->INIT!=0) {
                models[idx]->run->part_number = lattice_data[i_part]->part_number;
                models[idx]->lattice_simulate_start();
            } else {
              
            }
        }
        lattice_wait_queue_green(i);
        lattice_wait_table_read_green(i);
    }

    // initial measurements
    for (int i=0;i<big_lattice_red_passes;  i++) {
        lattice_copy_boundaries_red(i);
        lattice_wait_table_write_red(i);
        lattice_measure_red(i);
        lattice_wait_queue_red(i);
        lattice_save_state_red(i);
        for (int dv=0;dv<compute_devices_number;dv++){
            int i_part = 2*(dv + i * compute_devices_number);
            lattice_setup_lattice_pointer_initial(&i_part);
        }
        lattice_wait_table_read_red(i);
    }
    for (int i=0;i<big_lattice_green_passes;i++) {
        lattice_copy_boundaries_green(i);
        lattice_wait_table_write_green(i);
        lattice_measure_green(i);
        lattice_wait_queue_green(i);
        lattice_save_state_green(i);


        for (int dv=0;dv<compute_devices_number;dv++){
            int i_part = 2*(dv + i * compute_devices_number)+1;
            lattice_setup_lattice_pointer_initial(&i_part);
        }
        lattice_wait_table_read_green(i);
    }

    // update ITER_counter over all sublattices
    for (int i=0;i<big_lattice_parts;i++){
        int idx = lattice_data[i]->models_index;
        if ((models[idx]->run->INIT!=0)&&(models[idx]->ITER_counter==0)) models[idx]->ITER_counter = 1;
    }

    // perform thermalization
    for (unsigned int j=0; j<(unsigned int) global_run->NAV; j++){
        for (int i=0;i<big_lattice_red_passes;  i++) {
            int i_part = i;
            int idx = lattice_data[i_part]->models_index;
            if (models[idx]->NAV_counter==j) {
                lattice_copy_boundaries_red(i);
                lattice_wait_table_write_red(i);
                lattice_update_red(i);
                lattice_wait_kernels_red(i);
                lattice_orthogonalization_red(i);
                lattice_wait_queue_red(i);
                lattice_wait_table_read_red(i);
            }
        }
        for (int i=0;i<big_lattice_green_passes;i++) {
            int i_part = i;
            int idx = lattice_data[i_part]->models_index;
            if (models[idx]->NAV_counter==j) {
                lattice_copy_boundaries_green(i);
                lattice_wait_table_write_green(i);
                lattice_update_green(i);
                lattice_wait_kernels_green(i);
                lattice_orthogonalization_green(i);
                lattice_wait_queue_green(i);
                lattice_wait_table_read_green(i);
            }
        }

        // update NAV_counters
        for (int i=0;i<big_lattice_parts;i++){
            int idx = lattice_data[i]->models_index;
            if (models[idx]->NAV_counter==j) models[idx]->NAV_counter = j+1;
        }
        if (j % 10 == 0) printf("\rGPU thermalization [%i]",j);

        // save lattice state
            
    }

    // perform working cycles
    for (unsigned int t=1; t<(unsigned int) global_run->ITER; t++){ // zero measurement - on initial configuration!!!
        for (int j=0; j<global_run->NITER; j++){
            for (int i=0;i<big_lattice_red_passes;  i++) {
                int i_part = i;
                int idx = lattice_data[i_part]->models_index;
                if (models[idx]->ITER_counter==t) {
                    lattice_copy_boundaries_red(i);
                    lattice_wait_table_write_red(i);
                    lattice_update_red(i);
                    lattice_wait_kernels_red(i);
                    lattice_orthogonalization_red(i);
                    lattice_wait_queue_red(i);
                    lattice_wait_table_read_red(i);
                }
            }
            for (int i=0;i<big_lattice_green_passes;i++) {
                int i_part = i;
                int idx = lattice_data[i_part]->models_index;
                if (models[idx]->ITER_counter==t) {
                    lattice_copy_boundaries_green(i);
                    lattice_wait_table_write_green(i);
                    lattice_update_green(i);
                    lattice_wait_kernels_green(i);
                    lattice_orthogonalization_green(i);
                    lattice_wait_queue_green(i);
                    lattice_wait_table_read_green(i);
                }
            }
        }
        for (int i=0;i<big_lattice_red_passes;  i++) {
            int idx = lattice_data[i]->models_index;
            if (models[idx]->ITER_counter==t) {
                lattice_copy_boundaries_red(i);
                lattice_wait_table_write_red(i);
                lattice_measure_red(i);
                lattice_wait_queue_red(i);
                if (((unsigned int) global_run->ITER - 1)==t) {
                    for (int dv=0;dv<compute_devices_number;dv++){
                        int i_part = 2*(dv + i * compute_devices_number);
                        lattice_setup_lattice_pointer_last(&i_part);
                    }
                }
                lattice_wait_table_read_red(i);
                if (((unsigned int) global_run->ITER - 1)==t)
                    SUNcpu[lattice_data[i]->models_index]->lattice_check_cpu(models[lattice_data[i]->models_index]);
            }
        }
        for (int i=0;i<big_lattice_green_passes;i++) {
            int idx = lattice_data[i]->models_index;
            if (models[idx]->ITER_counter==t) {
                lattice_copy_boundaries_green(i);
                lattice_wait_table_write_green(i);
                lattice_measure_green(i);
                lattice_wait_queue_green(i);
                if (((unsigned int) global_run->ITER - 1)==t) {
                    for (int dv=0;dv<compute_devices_number;dv++){
                        int i_part = 2*(dv + i * compute_devices_number)+1;
                        lattice_setup_lattice_pointer_last(&i_part);
                    }
                }
                lattice_wait_table_read_green(i);
                if (((unsigned int) global_run->ITER - 1)==t)
                    SUNcpu[lattice_data[i]->models_index]->lattice_check_cpu(models[lattice_data[i]->models_index]);
            }
        }

        // update ITER_counters
        for (int i=0;i<big_lattice_parts;i++){
            int idx = lattice_data[i]->models_index;
            if (models[idx]->ITER_counter==t) models[idx]->ITER_counter = t+1;
        }

        if (t % 10 == 0) printf("\rGPU working iteration [%u]",t);

        // save lattice state

    }

    for (int i=0;i<compute_devices_number;i++){
        models[i]->lattice_print_elapsed_time();
        models[i]->lattice_analysis();
        models[i]->lattice_write_results();
        if (!global_run->turnoff_config_save) models[i]->lattice_write_configuration();
        models[i]->lattice_print_measurements();
    }

}

void            BL::lattice_copy_host_to_GPU(lattice_data_buffers* lat){
        // host->GPU
        int idx = lat->models_index;
        if (!lat->one_device_one_part) {
            // shadow copy active part: lattice_data[i]->models[idx]
            if (models[idx]->run->precision==model::model_precision_double){
                models[idx]->lattice_table_unmap_async(lat->plattice_table_double);
            } else {
                models[idx]->lattice_table_unmap_async(lat->plattice_table_float);
            }
        }
        models[idx]->run->part_number = lat->part_number;
}
void            BL::lattice_copy_GPU_to_host(lattice_data_buffers* lat){
        // GPU->host
        int idx = lat->models_index;
        if (!lat->one_device_one_part) {             
            // shadow copy active part: models[idx]->lattice_data[i]
            if (models[idx]->run->precision==model::model_precision_double){
                lat->plattice_table_double = (cl_double4*) models[idx]->lattice_table_map_async();
            } else {
                lat->plattice_table_float = (cl_float4*) models[idx]->lattice_table_map_async();
            }
            
        }
}

void            BL::lattice_get_low_boundary(int* i){
        // copy slice (with width x=1) from LOW->(N1-1)*N2N3N4 to I->(N1+1)*N2N3N4
        int idx = lattice_data[(*i)]->models_index;
        int i_low=(*i)-1;
        if (i_low<0) i_low=big_lattice_parts-1;
        int idx_low = lattice_data[i_low]->models_index;

        int lattice_group = models[idx]->run->lattice_group;
        int lattice_nd = models[idx]->run->lattice_nd;
        size_t boundary_size = models[idx]->lattice_domain_n2n3n4;
        unsigned int offset_src = abs(models[idx_low]->run->lattice_domain_size[0]-1) * models[idx_low]->lattice_domain_n2n3n4;
        unsigned int offset_dst = (models[idx]->run->lattice_domain_size[0]+1) * models[idx]->lattice_domain_n2n3n4;

        if (models[idx]->run->precision==model::model_precision_double){
            boundary_size *= sizeof(cl_double4);
            for (int j=0;j<models[idx]->lattice_group_elements[lattice_group]/4;j++){
                for (int t=0;t<lattice_nd;t++){
                    void* plattice_table_low_boundary = (void*) (lattice_data[i_low]->plattice_table_double + (offset_src + models[idx_low]->lattice_table_row_size * (t + j*lattice_nd)));
                    void* plattice_table_boundary     = (void*) (lattice_data[(*i) ]->plattice_table_double + (offset_dst + models[idx_low]->lattice_table_row_size * (t + j*lattice_nd)));
                    memcpy_s(plattice_table_boundary,boundary_size,plattice_table_low_boundary,boundary_size);
                }
            }
        } else {
            boundary_size *= sizeof(cl_float4);
            for (int j=0;j<models[idx]->lattice_group_elements[lattice_group]/4;j++){
                for (int t=0;t<lattice_nd;t++){
                    void* plattice_table_low_boundary = (void*) (lattice_data[i_low]->plattice_table_float + (offset_src + models[idx_low]->lattice_table_row_size * (t + j*lattice_nd)));
                    void* plattice_table_boundary     = (void*) (lattice_data[(*i) ]->plattice_table_float + (offset_dst + models[idx_low]->lattice_table_row_size * (t + j*lattice_nd)));
                    memcpy_s(plattice_table_boundary,boundary_size,plattice_table_low_boundary,boundary_size);
                }
            }
        }
}
void            BL::lattice_get_high_boundary(int* i){
        // copy slice (with width x=1) from HIGH->0 to I->N1*N2N3N4
        int idx = lattice_data[(*i)]->models_index;
        int i_high=(*i)+1;
        if (i_high>=big_lattice_parts) i_high=0;
        int idx_high = lattice_data[i_high]->models_index;
        int lattice_group = models[idx]->run->lattice_group;
        int lattice_nd = models[idx]->run->lattice_nd;
        size_t boundary_size = models[idx]->lattice_domain_n2n3n4;
        unsigned int offset = models[idx]->run->lattice_domain_size[0] * models[idx]->lattice_domain_n2n3n4;

        if (models[idx]->run->precision==model::model_precision_double){
            boundary_size *= sizeof(cl_double4);
            for (int j=0;j<models[idx]->lattice_group_elements[lattice_group]/4;j++){
                for (int t=0;t<lattice_nd;t++){
                    void* plattice_table_high_boundary = (void*) (lattice_data[(*i)]->plattice_table_double + (offset + models[idx]->lattice_table_row_size * (t + j*lattice_nd)));
                    void* plattice_table_boundary = (void*) (lattice_data[i_high]->plattice_table_double + (models[idx_high]->lattice_table_row_size * (t + j*lattice_nd)));
                    memcpy_s(plattice_table_high_boundary,boundary_size,plattice_table_boundary,boundary_size);
                }
            }
        } else {
            boundary_size *= sizeof(cl_float4);
            for (int j=0;j<models[idx]->lattice_group_elements[lattice_group]/4;j++){
                for (int t=0;t<lattice_nd;t++){
                    void* plattice_table_high_boundary = (void*) (lattice_data[(*i)]->plattice_table_float + (offset + models[idx]->lattice_table_row_size * (t + j*lattice_nd)));
                    void* plattice_table_boundary = (void*) (lattice_data[i_high]->plattice_table_float + (models[idx_high]->lattice_table_row_size * (t + j*lattice_nd)));
                    memcpy_s(plattice_table_high_boundary,boundary_size,plattice_table_boundary,boundary_size);
                }
            }
        }
}

void            BL::lattice_setup_lattice_pointer_initial(int* i){
        models[lattice_data[(*i)]->models_index]->lattice_pointer_initial = (unsigned int*) models[lattice_data[(*i)]->models_index]->lattice_table_map_async();
}
void            BL::lattice_setup_lattice_pointer_last(int* i){
        models[lattice_data[(*i)]->models_index]->lattice_pointer_last = (unsigned int*) models[lattice_data[(*i)]->models_index]->lattice_table_map_async();
}

void            BL::lattice_wait_queue_red(int i){
    for (int dv=0;dv<compute_devices_number;dv++){
        int i_part = 2*(dv + i * compute_devices_number);
        int idx = lattice_data[i_part]->models_index;
        models[idx]->lattice_wait_for_queue_finish();
        if (!lattice_data[i_part]->one_device_one_part) lattice_copy_GPU_to_host(lattice_data[i_part]);  // GPU->host
    }
}
void            BL::lattice_wait_queue_green(int i){
    for (int dv=0;dv<compute_devices_number;dv++){
        int i_part = 2*(dv + i * compute_devices_number) + 1;
        int idx = lattice_data[i_part]->models_index;
        models[idx]->lattice_wait_for_queue_finish();
        if (!lattice_data[i_part]->one_device_one_part) lattice_copy_GPU_to_host(lattice_data[i_part]);  // GPU->host
    }
}

void            BL::lattice_wait_kernels_red(int i){
    for (int dv=0;dv<compute_devices_number;dv++){
        int i_part = 2*(dv + i * compute_devices_number);
        int idx = lattice_data[i_part]->models_index;
        models[idx]->lattice_wait_for_queue_finish();
    }
}
void            BL::lattice_wait_kernels_green(int i){
    for (int dv=0;dv<compute_devices_number;dv++){
        int i_part = 2*(dv + i * compute_devices_number) + 1;
        int idx = lattice_data[i_part]->models_index;
        models[idx]->lattice_wait_for_queue_finish();
    }
}

void            BL::lattice_wait_table_read_red(int i){
        for (int dv=0;dv<compute_devices_number;dv++){
            int i_part = 2*(dv + i * compute_devices_number);
            int idx = lattice_data[i_part]->models_index;
            if (!lattice_data[i_part]->one_device_one_part) models[idx]->lattice_wait_for_table_read();
        }
}
void            BL::lattice_wait_table_read_green(int i){
        for (int dv=0;dv<compute_devices_number;dv++){
            int i_part = 2*(dv + i * compute_devices_number)+1;
            int idx = lattice_data[i_part]->models_index;
            if (!lattice_data[i_part]->one_device_one_part) models[idx]->lattice_wait_for_table_read();
        }
}

void            BL::lattice_wait_table_write_red(int i){
        for (int dv=0;dv<compute_devices_number;dv++){
            int i_part = 2*(dv + i * compute_devices_number);
            int idx = lattice_data[i_part]->models_index;
            if (!lattice_data[i_part]->one_device_one_part) models[idx]->lattice_wait_for_table_write();
        }
}
void            BL::lattice_wait_table_write_green(int i){
        for (int dv=0;dv<compute_devices_number;dv++){
            int i_part = 2*(dv + i * compute_devices_number)+1;
            int idx = lattice_data[i_part]->models_index;
            if (!lattice_data[i_part]->one_device_one_part) models[idx]->lattice_wait_for_table_write();
        }
}

void            BL::lattice_copy_boundaries_red(int i){
        for (int dv=0;dv<compute_devices_number;dv++){
            int i_part = 2*(dv + i * compute_devices_number);
            int idx = lattice_data[i_part]->models_index;
            if (models[idx]->run->INIT!=0) {
                if (!lattice_data[i_part]->one_device_one_part) {
                    // 1) update low boundary links in lat->plattice_table_xxx with lattice_data[i_part-1]
                    lattice_get_low_boundary(&i_part);
                    // 2) update high boundary links in lat->plattice_table_xxx with lattice_data[i_part+1]
                    lattice_get_high_boundary(&i_part);
                    // 3) copy lattice part from [host] to [GPU]
                    lattice_copy_host_to_GPU(lattice_data[i_part]);  // host->GPU
                }
            }
        }
}
void            BL::lattice_copy_boundaries_green(int i){
        for (int dv=0;dv<compute_devices_number;dv++){
            int i_part = 2*(dv + i * compute_devices_number)+1;
            int idx = lattice_data[i_part]->models_index;
            if (models[idx]->run->INIT!=0) {
                if (!lattice_data[i_part]->one_device_one_part) {
                    // 1) update low boundary links in lat->plattice_table_xxx with lattice_data[i_part-1]
                    lattice_get_low_boundary(&i_part);
                    // 2) update high boundary links in lat->plattice_table_xxx with lattice_data[i_part+1]
                    lattice_get_high_boundary(&i_part);
                    // 3) copy lattice part from [host] to [GPU]
                    lattice_copy_host_to_GPU(lattice_data[i_part]);  // host->GPU
                }
            }
        }
}

void            BL::lattice_measure_red(int i){
        for (int dv=0;dv<compute_devices_number;dv++){
            int i_part = 2*(dv + i * compute_devices_number);
            int idx = lattice_data[i_part]->models_index;
            // 4) measure the part
            models[idx]->lattice_measure();
        }
}
void            BL::lattice_measure_green(int i){
        for (int dv=0;dv<compute_devices_number;dv++){
            int i_part = 2*(dv + i * compute_devices_number)+1;
            int idx = lattice_data[i_part]->models_index;
            // 4) measure the part
            models[idx]->lattice_measure();
        }
}

void            BL::lattice_save_state_red(int i){
        for (int dv=0;dv<compute_devices_number;dv++){
            int i_part = 2*(dv + i * compute_devices_number);
            int idx = lattice_data[i_part]->models_index;
            models[idx]->run->part_number = lattice_data[i_part]->part_number;
            if (!models[idx]->run->turnoff_state_save) models[idx]->lattice_save_state();
        }
}
void            BL::lattice_save_state_green(int i){
        for (int dv=0;dv<compute_devices_number;dv++){
            int i_part = 2*(dv + i * compute_devices_number)+1;
            int idx = lattice_data[i_part]->models_index;
            models[idx]->run->part_number = lattice_data[i_part]->part_number;
            if (!models[idx]->run->turnoff_state_save) models[idx]->lattice_save_state();
        }
}


void            BL::lattice_update_red(int i){
        for (int dv=0;dv<compute_devices_number;dv++){
            int i_part = 2*(dv + i * compute_devices_number);
            int idx = lattice_data[i_part]->models_index;
            // 4) measure the part
            models[idx]->lattice_update();
        }
}
void            BL::lattice_update_green(int i){
        for (int dv=0;dv<compute_devices_number;dv++){
            int i_part = 2*(dv + i * compute_devices_number)+1;
            int idx = lattice_data[i_part]->models_index;
            // 4) measure the part
            models[idx]->lattice_update();
        }
}

void            BL::lattice_orthogonalization_red(int i){
        for (int dv=0;dv<compute_devices_number;dv++){
            int i_part = 2*(dv + i * compute_devices_number);
            int idx = lattice_data[i_part]->models_index;
            // 4) measure the part
            models[idx]->lattice_orthogonalization();
        }
}
void            BL::lattice_orthogonalization_green(int i){
        for (int dv=0;dv<compute_devices_number;dv++){
            int i_part = 2*(dv + i * compute_devices_number)+1;
            int idx = lattice_data[i_part]->models_index;
            // 4) measure the part
            models[idx]->lattice_orthogonalization();
        }
}

}
