CC = g++
#CC = gcc

#####################################################
AMD = 0
INTEL = 1
NVIDIA = 2
#PLAT = $(AMD)
PLAT = $(NVIDIA)

# sequential run on CPU:
CPU_RUN = 0

BIGLAT = 0
USE_OPENMP = 0
CHB2 = 0

# If defined BIGLAT:
NPARTS = 2
NDEV = 1
#####################################################

ifndef AMDAPPSDKROOT
AMDAPPSDKROOT=/opt/AMD/CodeXL_1.4.5728/examples/Teapot/AMDTTeaPotLib/AMDTTeaPotLib
endif

NVIDIASDKROOT=/usr/local/cuda

ifeq ($(PLAT), $(AMD))
CFLAGS = -I$(AMDAPPSDKROOT)/include -g -O2
LDFLAGS = -L$(AMDAPPSDKROOT)/lib/x86_64 -lOpenCL
#CFLAGS = -I$(AMDAPPSDKROOT)/inc -g -O2
#LDFLAGS = -lOpenCL
endif

ifeq ($(PLAT), $(NVIDIA))
CFLAGS = -I$(NVIDIASDKROOT)/include -g -O2
LDFLAGS = -L/usr/local/cuda/lib64 -lOpenCL
endif

ifeq ($(CPU_RUN), 1)
CFLAGS += -D CPU_RUN
endif


ifeq ($(BIGLAT), 1)
CFLAGS += -D BIGLAT
CFLAGS += -D NPARTS=$(NPARTS) -D NDEV=$(NDEV)
endif

ifeq ($(USE_OPENMP), 1)
CFLAGS += -D USE_OPENMP
LDFLAGS += -fopenmp
endif

ifeq ($(CHB2), 1)
CFLAGS += -D CHB2
endif

# -Wall

# project name
TARGET = QCDGPU

SRCS =  QCDGPU.cpp \
	random/random.cpp \
	suncl/suncl.cpp 

ifeq ($(CPU_RUN), 0)
SRCS += clinterface/clinterface.cpp \
	suncl/suncpu.cpp \
	suncl/su2cpu.cpp \
	suncl/su3cpu.cpp \
	data_analysis/data_analysis.cpp 
endif
ifeq ($(CPU_RUN), 1)
SRCS += suncpp/su2/algebra_su2.cpp \
	suncpp/su2/update_su2.cpp \
	suncpp/su3/algebra_su3.cpp \
	suncpp/su3/update_su3.cpp \
	suncpp/Measurements/analysis_cpp.cpp \
	suncpp/coord_work/coord_work.cpp \
	suncpp/IO/io.cpp 
endif

HDRS =  QCDGPU.h \
	clinterface/clinterface.h \
	random/random.h \
	kernel/complex.h \
	suncl/suncl.h \
	data_analysis/data_analysis.h
	
ifeq ($(CPU_RUN), 0)
HDRS += clinterface/platform.h \
	suncl/suncpu.h \
	suncl/su2cpu.h \
	suncl/su3cpu.h
endif
ifeq ($(CPU_RUN), 1)
HDRS += suncpp/sunh.h \
	suncpp/suncpp.h \
	suncpp/su2/algebra_su2.h \
	suncpp/su2/update_su2.h \
	suncpp/su3/algebra_su3.h \
	suncpp/su3/update_su3.h \
	suncpp/Measurements/Plq.h \
	suncpp/Measurements/S.h \
	suncpp/Measurements/analysis_cpp.h \
	suncpp/coord_work/coord_work.h \
	suncpp/Update/sun_update.h \
	suncpp/IO/io.h 
endif

OBJS = $(SRCS:.cpp=.o)

$(TARGET): $(SRCS) $(HDRS)

ifndef AMDAPPSDKROOT
	@echo >&2
	@echo "AMD APP SDK not installed" >&2
	@exit 1
else
	$(CC) $(CFLAGS) $(SRCS) -o $(TARGET) $(LDFLAGS)
endif

clobber:
	rm -rf $(TARGET) $(OBJS)

clean:
	rm -f $(TARGET)
