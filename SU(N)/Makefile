all: amd intel nvidia

CC = g++
CFLAGS = -Wextra -W -O2

SRCS =  QCDGPU.cpp \
	clinterface/clinterface.cpp \
	random/random.cpp \
	suncl/suncl.cpp \
	suncl/suncpu.cpp \
	suncl/su2cpu.cpp \
	suncl/su3cpu.cpp \
	data_analysis/data_analysis.cpp

HDRS =  QCDGPU.h \
	clinterface/platform.h \
	clinterface/clinterface.h \
	random/random.h \
	kernel/complex.h \
	suncl/suncl.h \
	suncl/suncpu.h \
	suncl/su2cpu.h \
	suncl/su3cpu.h \
	data_analysis/data_analysis.h

is_64=$(shell s=`uname -m`; if (echo $$s | grep x86_64 > /dev/null); then echo 1; fi)

ifndef AMDAPPSDKROOT
AMDAPPSDKROOT=/opt/AMDAPP/
endif
		
# amd
SDK_INC_AMD=$(AMDAPPSDKROOT)/include
ifeq ($(is_64), 1)
SDK_LIB_AMD=$(AMDAPPSDKROOT)/lib/x86_64/
else
SDK_LIB_AMD=$(AMDAPPSDKROOT)/lib/x86/
endif

# intel
SDK_INC_INTEL=/usr/include/
ifeq ($(is_64), 1)
SDK_LIB_INTEL=/usr/lib64/
endif


# nvidia
SDK_INC_NVIDIA=/usr/local/cuda/include/
SDK_LIB_NVIDIA=/usr/lib/

OBJS_AMD=$(patsubst %.c,%.amd.o,$(SRCS))
OBJS_INTEL=$(patsubst %.c,%.intel.o,$(SRCS))
OBJS_NVIDIA=$(patsubst %.c,%.nvidia.o,$(SRCS))

TARGET=$(patsubst %.cpp,%,$(wildcard *.cpp))

amd:	$(join $(TARGET),.amd)
intel:	$(join $(TARGET),.intel)
nvidia:	$(join $(TARGET),.nvidia)

%.amd:	$(SRCS) $(HDRS)
	$(CC) $(CFLAGS) $(SRCS) -I$(SDK_INC_AMD) -L$(SDK_LIB_AMD) -Wl,-rpath,$(SDK_LIB_AMD)  -o $@ -lOpenCL

%.intel:$(SRCS) $(HDRS)
	$(CC) $(CFLAGS) $(SRCS) -I$(SDK_INC_INTEL) -L$(SDK_LIB_INTEL) -Wl,-rpath,$(PATH_TO_INTEL_LIB) -o $@ -lOpenCL

%.nvidia:$(SRCS) $(HDRS)
	$(CC) $(CFLAGS) $(SRCS) -I$(SDK_INC_NVIDIA) -L$(SDK_LIB_NVIDIA) -o $@ -lOpenCL

clobber:
	rm -rf $(TARGET).amd $(TARGET).intel $(TARGET).nvidia $(OBJS_AMD) $(OBJS_INTEL) $(OBJS_NVIDIA)

clean:
	rm -f $(TARGET).amd $(TARGET).intel $(TARGET).nvidia
