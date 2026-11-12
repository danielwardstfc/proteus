# Builds generate_box_mesh against MOOSE's MFEM.

MFEM_DIR ?= ../../../../../moose/framework/contrib/mfem/build-opt
MFEM_LIB_SUFFIX ?= -opt

CONFIG_MK = $(MFEM_DIR)/config/config.mk
ifeq ($(wildcard $(CONFIG_MK)),)
  $(error No MFEM config at '$(CONFIG_MK)'. Set MFEM_DIR to your MFEM build, \
          e.g. export MFEM_DIR=/path/to/mfem/build-opt)
endif
include $(CONFIG_MK)

TARGET = mfem_mesh_generator

CXX      = $(MFEM_CXX)
CXXFLAGS = $(MFEM_FLAGS)
LDLIBS   = $(subst -lmfem,-lmfem$(MFEM_LIB_SUFFIX),$(MFEM_LIBS))

all: $(TARGET)

$(TARGET): mfem_mesh_generator.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDLIBS)

clean:
	rm -f $(TARGET)

.PHONY: all clean
