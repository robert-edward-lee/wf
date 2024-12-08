################################################################################
#                              НАСТРОЙКА ПРОЕКТА                               #
################################################################################
PROJECT_NAME = wf

WORK_DIRS = .
BUILD_DIR = .
SOURCE_DIRS = .
SOURCES = $(wildcard $(SOURCE_DIRS)/*.c)
OBJECTS = $(addprefix $(BUILD_DIR)/,$(notdir $(patsubst %.c,%.o,$(SOURCES))))
DEPENDS = $(patsubst %.o,%.d,$(OBJECTS))
-include $(DEPENDS)

SHARED_LIB = $(BUILD_DIR)/lib$(PROJECT_NAME).dll
IMP_LIB = $(BUILD_DIR)/lib$(PROJECT_NAME).dll.a

CC = gcc
LD = $(CC)
LDFLAGS = -shared -Wl,--out-implib,$(IMP_LIB)
STRIP = $(TOOLCHAIN_PREFIX)strip

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
	$(WARN_FLAGS) \
	$(STDC_FLAGS) \
	$(EXTRA_FLAGS) \
	$(DEPEND_FLAGS)

all: test

$(SHARED_LIB): $(OBJECTS)
	@echo '  LD      ' $@
	@$(LD) $(LDFLAGS) -o $@ $^
	@echo '  STRIP   ' $@
	@$(STRIP) $@

$(IMP_LIB): $(SHARED_LIB)

$(BUILD_DIR)/%.o: %.c
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

test: $(SHARED_LIB)
	@python pytest
