CXX ?= g++

NETCDF_CFLAGS := $(shell pkg-config --cflags netcdf 2>/dev/null)
NETCDF_LIBS := $(shell pkg-config --libs netcdf 2>/dev/null)

ifeq ($(NETCDF_CFLAGS),)
	ifneq ($(wildcard /opt/homebrew/include/netcdf.h),)
		NETCDF_CFLAGS = -I/opt/homebrew/include
		NETCDF_LIBS = -L/opt/homebrew/lib -lnetcdf
	else ifneq ($(wildcard /usr/local/include/netcdf.h),)
		NETCDF_CFLAGS = -I/usr/local/include
		NETCDF_LIBS = -L/usr/local/lib -lnetcdf
	else
		NETCDF_LIBS = -lnetcdf
	endif
endif

CXXFLAGS ?= -O2 -std=c++11 -Wall -Wextra -pedantic $(NETCDF_CFLAGS)
LDFLAGS ?=
LDLIBS ?= $(NETCDF_LIBS)

TARGET = trace_fieldlines
SRCS = main.cpp field_data.cpp rk4.cpp trace.cpp
OBJS = $(SRCS:.cpp=.o)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
