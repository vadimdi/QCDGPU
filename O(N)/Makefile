CC = g++
CFLAGS = -I$(AMDAPPSDKROOT)/include -g
# -Wall
LDFLAGS = -L$(AMDAPPSDKROOT)/lib/x86_64 -lOpenCL

# project name
TARGET = QCDGPU

SRCS =  QCDGPU.cpp \
	clinterface/clinterface.cpp \
	random/random.cpp \
	suncl/suncl.cpp \
	suncl/suncpu.cpp \
	data_analysis/data_analysis.cpp \
	suncl/biglattice.cpp

HDRS =  QCDGPU.h \
	clinterface/platform.h \
	clinterface/clinterface.h \
	random/random.h \
	kernel/complex.h \
	suncl/suncl.h \
	suncl/suncpu.h \
	data_analysis/data_analysis.h \
	suncl/biglattice.h

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
