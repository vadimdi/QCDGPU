/******************************************************************************
 * @file     clinterface.h
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.4
 *
 * @brief    [QCDGPU] 
 *           Interface for OpenCL AMD APP & nVidia SDK environment
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
#ifndef clinterface_h
#define clinterface_h

#include <CL/cl.h>
#include "platform.h"


namespace GPU_CL{
class GPU {
    public:
            typedef enum enum_GPU_error_codes{
                GPU_error_SUCCESS,                          // OK
                GPU_error_no_platform,                      // no OpenCL platform
                GPU_error_no_device,                        // no OpenCL device
                GPU_error_device_initialization_failed,     // error of device initialization
                GPU_error_no_buffer,                        // no buffer
            } GPU_error_codes;

            typedef enum enum_GPU_vendors{
                GPU_vendor_any,             // any vendor
                GPU_vendor_AMD,             // Advanced Micro Devices, Inc.
                GPU_vendor_Apple,           // Apple
                GPU_vendor_nVidia,          // NVIDIA Corporation
                GPU_vendor_Intel,           // Intel(R) Corporation
                GPU_vendor_None             // 
            } GPU_vendors;

            typedef enum enum_GPU_buffer_types{
                buffer_type_None,           // NULL buffer type
                buffer_type_Global,         // Global buffer type
                buffer_type_Input,          // Input buffer type
                buffer_type_IO,             // Like Global buffer, but with initialization
                buffer_type_Output,         // Output buffer type
                buffer_type_Constant,       // Constant buffer type
                buffer_type_LDS,            // Local buffer type
                buffer_type_UAV             // UAV buffer type
            } GPU_buffer_types;


            typedef enum enum_GPU_storage_type{
                GPU_storage_none,           // default unknown storage type
                GPU_storage_joint,          // joint data based on two other data arrays
                GPU_storage_double2high,    // high dword in double2
                GPU_storage_double2low,     // low  dword in double2
                GPU_storage_double          // double
            } GPU_storage_type;

            typedef struct _GPU_time_deviation{
                double mean;                // mean value
                double deviation;           // deviation
                int    number;              // number of elements
            } GPU_time_deviation;

            /**
            * A set of debug options
            */
            typedef struct _GPU_debug_flags{
                bool wait_for_keypress      : 1; /**< Prevents exit from program up to key press*/
                bool profiling              : 1; /**< Turn on/off kernels profiling*/
                bool brief_report           : 1; /**< Turn on/off brief report about quantitative characteristics of run (worgroup size, sizes of buffers, number pseudorandom numbers to produce, etc.)*/
                bool show_stage             : 1; /**< Enable/disable the output of run stages*/
                bool local_run              : 1; /**< Enable/disable local run option*/
                bool rebuild_binary         : 1; /**< Enable/disable rebuilding of binary files for each run*/
            } GPU_debug_flags;
            /**
            * Defines initial conditions for simulation, as physical as hardware
            */
            typedef struct _GPU_init_parameters{
                char   Variable[250];            /**< The name of a member of the structure*/
                int    iVarVal;                  /**< The integer value of a member of the structure*/
                double fVarVal;                  /**< The float value of a member of the structure*/
                char   txtVarVal[250];           /**< The text value of a member of the structure*/
                bool   final;                    /**< The flag defining that a member of the strucure is the last in a set of such members*/
            } GPU_init_parameters;

            /**
            * Information about the desired compute device
            */
            typedef struct _GPU_device_info{
                char*       device_name;                    /**< The name of desired compute device*/   // CL_DEVICE_NAME                          
                cl_ulong    local_memory_size;              /**< The local memory size of desired compute device (in bytes)*/   // CL_DEVICE_LOCAL_MEM_SIZE                
                cl_ulong    max_constant_size;              /**< The maximal size of the constant memory buffer of the desired device*/   // CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE
                cl_ulong    max_memory_size;                /**< The maximal size of memory can be allocated on the desired device*/   // CL_DEVICE_MAX_MEM_ALLOC_SIZE
                size_t      max_memory_width;               /**< The number of PRNG instances*/   // CL_DEVICE_IMAGE2D_MAX_WIDTH
                size_t      max_workgroup_size;             /**< The greatest of the maximal numbers of OpenCL work items in work group for each compute dimension */   // max(CL_DEVICE_MAX_WORK_ITEM_SIZES)
                size_t      memory_align_factor;            /**< The factor for buffers aligment*/   // memory align factor for buffers
             GPU_vendors    platform_vendor;                /**< The vendor of the platform containing the desired device*/   // active platform vendor
             GPU_vendors    device_vendor;                  /**< The vendor of desired compute device*/   // active device vendor
            } GPU_device_info;

// ___________________________________________________________ static variables
            static GPU_debug_flags GPU_debug;              // structure for debuging 
            static char current_path[FILENAME_MAX];        // FILENAME_MAX is definned by <stdio.h>
                  char* cl_root_path;                      // FILENAME_MAX is definned by <stdio.h>

// end of static variables ___________________________________

            GPU_device_info GPU_info;                       // GPU device info

            unsigned int GPU_limit_max_workgroup_size;      // manually limit max workgroup size (0 - do not limit)

            int GPU_current_kernel;                         // current kernel counter
            int GPU_current_buffer;                         // current buffer counter
            int GPU_current_program;                        // current program counter
            int GPU_active_program;                         // current program counter

            int GPU_inf_max_n;                              // max number in .inf-filenames

            unsigned int GPU_platform_id;                  // GPU platform to be used (0 = first, 1 = second, etc)
            unsigned int GPU_device_id;                    // GPU device to be used   (0 = first, 1 = second, etc)

            cl_uint GPU_platforms_number;                  // number of available platforms
            cl_uint GPU_devices_number;                    // number of available devices
            cl_uint GPU_total_devices_number;              // total number of available devices

            cl_uint GPU_max_work_item_dimensions;          // CL_MAX_WORK_ITEM_DIMENSIONS
            size_t* GPU_max_work_item_sizes;               // CL_MAX_WORK_ITEMS_SIZES

            cl_platform_id   GPU_platform;                 // utilized platform
            cl_device_id     GPU_device;                   // utilized platform
            cl_context       GPU_context;                  // utilized context
            cl_command_queue GPU_queue;                    // utilized command queue

            cl_int GPU_error;

            int     CPU_timers;                            // number of reserved CPU timers
            int*    CPU_timer;                             // CPU timer to measure performance
            int     CPU_current_timer_id;                  // current timer id

            GPU(void);
           ~GPU(void);


            // ___________________________________________ debug section
                    void OpenCL_Check_Error(cl_int CL_Error_code,const char * CL_Error_description);
            unsigned int convert_to_uint(float value);
            unsigned int convert_to_uint_HIGH(double x);
            unsigned int convert_to_uint_LOW(double x);
                   float convert_to_float(unsigned int value);
                   float convert_to_float(double value);
                  double convert_to_double(float value);
           static double convert_to_double(unsigned int value);
           static double convert_to_double(unsigned int value_LOW, unsigned int value_HIGH);

            // ___________________________________________ public functions
            int     device_initialize(void);
            int     device_finalize(int error_code);
            bool    device_auto_select(int platform_vendor,int vendor);
            bool    device_select(unsigned int platform_id,unsigned int device_id);
            char*   device_get_name(cl_device_id device);
            char*   platform_get_name(cl_platform_id platform);

            char*   source_read(const char* file_name);
            char*   source_add(char* source, const char* file_name);

            int     program_create(const char* source);
            int     program_create(const char* source,const char* options);
            int     program_set_active(int program_id);
            int     program_get_active(void);

            int     kernel_init(const char* kernel_name, unsigned int work_dimensions, const size_t* global_size, const size_t* local_size);
            int     kernel_init_buffer(int kernel_id,int buffer_id);
            int     kernel_init_constant(int kernel_id,int* host_ptr);
            int     kernel_init_constant_reset(int kernel_id,int* host_ptr,int argument_id);
            int     kernel_init_constant(int kernel_id,cl_uint4* host_ptr);
            int     kernel_init_constant(int kernel_id,float* host_ptr);
            int     kernel_init_constant(int kernel_id,double* host_ptr);
            int     kernel_run(int kernel_id);
            int     kernel_get_worksize(int kernel_id);
 GPU_time_deviation kernel_get_execution_time(int kernel_id);

            int     buffer_init(int buffer_type, int size, void* host_ptr, int size_of);
            int     buffer_write(int buffer_id);
   unsigned int*    buffer_map(int buffer_id);
      cl_float4*    buffer_map_float4(int buffer_id);
      cl_float4*    buffer_read_float4(int buffer_id);
            int     buffer_kill(int buffer_id);
 GPU_time_deviation buffer_write_get_time(int buffer_id);
 GPU_time_deviation buffer_read_get_time(int buffer_id);
   unsigned int     buffer_size_align(unsigned int size);
   unsigned int     buffer_size_align(unsigned int size,unsigned int step);
 GPU_init_parameters* get_init_file(char finitf[]);
            void    make_start_file(char* path);
            void    make_finish_file(char* path);

            bool    is_file_exist(char* path);

            void    print_available_hardware(void);
            void    print_stage(const char * stage);
            void    print_memory_utilized(void);
            int     print_mapped_buffer_uint   (int buffer_id,unsigned int number_of_elements);
            int     print_mapped_buffer_uint2  (int buffer_id,unsigned int number_of_elements);
            int     print_mapped_buffer_uint3  (int buffer_id,unsigned int number_of_elements);
            int     print_mapped_buffer_uint4  (int buffer_id,unsigned int number_of_elements);
            int     print_mapped_buffer_float  (int buffer_id,unsigned int number_of_elements);
            int     print_mapped_buffer_float2 (int buffer_id,unsigned int number_of_elements);
            int     print_mapped_buffer_float3 (int buffer_id,unsigned int number_of_elements);
            int     print_mapped_buffer_float4 (int buffer_id,unsigned int number_of_elements);
            int     print_mapped_buffer_float4 (int buffer_id,unsigned int number_of_elements,unsigned int offset);
            int     print_mapped_buffer_double (int buffer_id,unsigned int number_of_elements);
            int     print_mapped_buffer_double2(int buffer_id,unsigned int number_of_elements);
            int     print_time_detailed(void);

            int     start_timer_CPU(void);
            int     start_timer_CPU(int timer);
            int     stop_timer_CPU(int timer);
          double    timer_in_seconds_CPU(int timer);
          double    get_timer_CPU(int timer);
           char*    get_current_datetime(void);

           char*    trim(char* str);
             int    str_char_replace(char* str, char search, char replace);

    private:

      class kernels_hash {
         public:
                int          argument_id;
                cl_kernel    kernel;                        // kernel
                const char*  kernel_name;                   // ponter to kernel's name declaration
                unsigned int work_dimensions;               // work dimensions
                size_t*      global_size;                   // global size
                size_t*      local_size;                    // local size
                cl_program   program;                       // program
                int          program_id;                    // program_id for GPU_programs array
                // profiling data __________________________
                cl_ulong     kernel_start;                  // kernel last start time
                cl_ulong     kernel_finish;                 // kernel last finish time
                double       kernel_elapsed_time;           // total kernel execution time (in nanoseconds)
                double       kernel_elapsed_time_squared;   // total kernel execution time squared (in nanoseconds) - for deviation calculation
                int          kernel_number_of_starts;       // total number of kernel starts - for deviation calculation
                size_t       kernel_preferred_workgroup_size_multiple; // kernel preferred work group size multiple
                cl_ulong     kernel_local_mem_size;         // kernel local memory size

                kernels_hash(void);
               ~kernels_hash(void);
      };

      
      class buffers_hash{
         public:
                cl_mem       buffer;
                int          buffer_type;
                void*        host_ptr;
                int          size;
                size_t       size_in_bytes;                     // buffer size in bytes ( size*sizeof(...) )
              unsigned int*  mapped_ptr;                        // ptr to corresponding host memory after mapping
                // profiling data __________________________
                cl_ulong     buffer_write_start;                // buffer start write time
                cl_ulong     buffer_write_finish;               // buffer finish write time
                double       buffer_write_elapsed_time;         // total buffer write time (in nanoseconds)
                double       buffer_write_elapsed_time_squared; // total buffer write time squared (in nanoseconds) - for deviation calculation
                int          buffer_write_number_of;            // total number of writes - for deviation calculation
                cl_ulong     buffer_read_start;                 // buffer start read time
                cl_ulong     buffer_read_finish;                // buffer finish read time
                double       buffer_read_elapsed_time;          // total buffer read time (in nanoseconds)
                double       buffer_read_elapsed_time_squared;  // total buffer read time squared (in nanoseconds) - for deviation calculation
                int          buffer_read_number_of;             // total number of reads - for deviation calculation
                
                buffers_hash(void);
               ~buffers_hash(void);
      };
      
      class programs_hash{
         public:
                cl_program   program;
                const char*  source_ptr;             // ponter to source's declaration
                const char*  options;                // pointer to compiling options (-D... in clBuildProgram)
                const char*  md5;                    // md5 for source code
                      char*  datetime;               // datetime of kernel compiling
                      char*  build_log;              // program build log
                const char*  device;                 // device
                const char*  platform;               // platform
                int          source_length;

                programs_hash(void);
                ~programs_hash(void);
      };
    
            typedef union _Uint_and_Float{          // Uint <---> Float converter
                unsigned int uint_value[1];
                float        float_value;
            } Uint_and_Float;

            typedef union _Int_and_Double{          // Int <---> Double converter
                unsigned int int_value[2];
                double       double_value;
            } Int_to_Double;

            typedef union _Float_and_Double{        // Float <---> Double converter
                float        intval[2];
                double       dbVal;
            } Float_and_Double;

            kernels_hash*    GPU_kernels;           // Hash for kernels
            buffers_hash*    GPU_buffers;           // Hash for buffers pointers
            programs_hash*   GPU_programs;          // Hash for programs pointers 

            // ____________________________________ MD5 section
            #define MD5_blocksize   64
            unsigned char MD5_buffer[MD5_blocksize];
            unsigned int  MD5_T[MD5_blocksize];
            unsigned int  MD5_s[MD5_blocksize];
            unsigned int  MD5_state[4];
            unsigned int  MD5_bytes[2];
            unsigned char MD5_result[16];
            bool          MD5_finalized;

            static inline unsigned int MD5_rol(unsigned int x,int s);
            void MD5_init(void);
            void MD5_finalize(void);
            void MD5_step(const unsigned char* MD5_block);
            void MD5_getword(const unsigned char* buffer, unsigned int* x, unsigned int len);
            void MD5_setword(const unsigned int* x, unsigned char* buffer, unsigned int len);
            void MD5_update(const unsigned char* input, unsigned int len);
            void MD5_update(const char*          input, unsigned int len);
           char* MD5_getresult(void);
           char* MD5(const char* str);

            // private functions
            cl_uint     clGetDeviceInfoUint(cl_device_id device,cl_device_info inf);
            cl_ulong    clGetDeviceInfoUlong(cl_device_id device,cl_device_info device_info);
   GPU_time_deviation   time_get_deviation(double elapsed_time, double elapsed_time_squared, int number);
            double      time_gigabytes_per_second(double elapsed_time,double size);
            double      time_megabytes_per_second(double elapsed_time,double size);
            int         inf_file_delete(int index);
            int         inf_file_rename(int index_old,int index_new);
            int         inf_file_get_max_n(void);

};
};

#endif
