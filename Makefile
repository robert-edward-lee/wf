PROJECT_NAME = wf

WORK_DIRS = .

INCLUDE_DIRS =
DEFINES =
OPT_LEVEL = 3

SHARED_LIB = lib$(PROJECT_NAME).dll
IMP_LIB = lib$(PROJECT_NAME).dll.a

CC = gcc
LD = $(CC)
LDFLAGS = -shared -Wl,--out-implib,$(IMP_LIB)
STRIP = strip

INC_FLAGS = $(addprefix -I,$(INCLUDE_DIRS))
DEF_FLAGS = $(addprefix -D,$(DEFINES))
OPT_FLAGS = $(addprefix -O,$(OPT_LEVEL))
STDC_FLAGS = -ansi
WARN_FLAGS = -Wall -Wextra -pedantic
DEPEND_FLAGS = -MMD -MP

CFLAGS = \
	$(INC_FLAGS) \
	$(DEF_FLAGS) \
	$(OPT_FLAGS) \
	$(STDC_FLAGS) \
	$(WARN_FLAGS) \
	$(DEPEND_FLAGS) \
	$(EXTRA_FLAGS)

all: test

shared: $(SHARED_LIB)

$(SHARED_LIB): $(PROJECT_NAME).o
	@echo '  LD      ' $@
	@$(LD) $(LDFLAGS) -o $@ $^
	@echo '  STRIP   ' $@
	@$(STRIP) $@

$(IMP_LIB): $(SHARED_LIB)

$(PROJECT_NAME).o: $(PROJECT_NAME).c
	@echo '  CC      ' $@
	@$(CC) -c $(CFLAGS) -o $@ $<

clean:
	@$(RM) \
		$(foreach dir,$(WORK_DIRS),$(addsuffix /*.a,$(dir))) \
		$(foreach dir,$(WORK_DIRS),$(addsuffix /*.d,$(dir))) \
		$(foreach dir,$(WORK_DIRS),$(addsuffix /*.dll,$(dir))) \
		$(foreach dir,$(WORK_DIRS),$(addsuffix /*.exe,$(dir))) \
		$(foreach dir,$(WORK_DIRS),$(addsuffix /*.o,$(dir))) \
		$(foreach dir,$(WORK_DIRS),$(addsuffix /*.obj,$(dir))) \
		$(foreach dir,$(WORK_DIRS),$(addsuffix /*.tds,$(dir))) \
		$(foreach dir,$(WORK_DIRS),$(addsuffix /*.d,$(dir)))

test: shared
	@python pytest

-include $(PROJECT_NAME).d
