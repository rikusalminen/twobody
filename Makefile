CFLAGS+=-std=gnu99
CFLAGS+=-MMD
CFLAGS+=-W -Wall -Wextra

ifneq ($(DEBUG), 1)
CFLAGS+=-O3 -ffast-math -march=native
CFLAGS+=-DNDEBUG
else
CFLAGS+=-O0 -g -ggdb
CFLAGS+=-DDEBUG
endif

CFLAGS+=-I$(SRC_DIR)/include

CFLAGS+=-Wno-psabi # GCC warnings about AVX ABI (simd)
CFLAGS+=-Wno-unknown-warning-option

LDLIBS+=-lm
LDFLAGS+=

SRCS= \
	src/twobody/conic.c \
	src/twobody/anomaly.c \
	src/twobody/true_anomaly.c \
	src/twobody/eccentric_anomaly.c \
	src/twobody/orientation.c \
	src/twobody/orbit.c \
	test/twobody/conic_test.c \
	test/twobody/anomaly_test.c \
	test/twobody/true_anomaly_test.c \
	test/twobody/eccentric_anomaly_test.c \
	test/twobody/orientation_test.c \
	test/twobody/orbit_test.c \
	test/twobody/twobody_test.c \
	test/numtest.c \
	src/twobody/twobody.c

TARGETS= \
	test/twobody/twobody_test \
	libtwobody.a

libtwobody.a: \
	src/twobody/conic.o \
	src/twobody/anomaly.o \
	src/twobody/true_anomaly.o \
	src/twobody/eccentric_anomaly.o \
	src/twobody/orientation.o \
	src/twobody/orbit.o \
	src/twobody/twobody.o

test/twobody/twobody_test: \
	test/twobody/conic_test.o \
	test/twobody/anomaly_test.o \
	test/twobody/true_anomaly_test.o \
	test/twobody/eccentric_anomaly_test.o \
	test/twobody/orientation_test.o \
	test/twobody/orbit_test.o \
	test/twobody/twobody_test.o \
	test/numtest.o \
	libtwobody.a

.DEFAULT_GOAL=all
.PHONY: all
all: $(TARGETS)

SRC_DIR ?= $(patsubst %/,%, $(dir $(abspath $(firstword $(MAKEFILE_LIST)))))

.PHONY: clean
.SILENT: clean
clean:
	$(RM) $(TARGETS)
	$(RM) $(OBJS)
	$(RM) $(DEPS)
ifneq ($(SRC_DIR), $(CURDIR))
	-@rmdir --ignore-fail-on-non-empty -p $(OBJDIRS)
endif
	$(RM) cscope.out cscope.out.in cscope.out.po
	$(RM) tags

OBJS=$(SRCS:.c=.o)
DEPS=$(OBJS:.o=.d)

# Object file subdirectories
ifneq ($(SRC_DIR), $(CURDIR))
vpath %.c $(SRC_DIR)

OBJDIRS=$(filter-out ./, $(sort $(dir $(OBJS))))

$(OBJDIRS): ; @mkdir -p $@
$(DEPS): | $(OBJDIRS)
$(OBJS): | $(OBJDIRS)
endif

-include $(DEPS)

# implicit rules for building archives not parallel safe (e.g. make -j 3)
%.a: ; $(AR) rcs $@ $^

# cscope.out
cscope.out: $(SRCS)
	cscope -f $@ -I$(SRC_DIR)/include -bq $^

# ctags
tags: $(SRCS)
	ctags -f $@ -R $(SRC_DIR)/include $^
