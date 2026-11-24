ROOT_CONFIG ?= root-config
ifeq ($(origin CXX),default)
CXX := $(shell $(ROOT_CONFIG) --cxx)
endif

CXXFLAGS ?= -O2 -Wall -Wextra -Wpedantic
CPPFLAGS += -Iinclude -Ibuild -DTCTSCAN_PROJECT_ROOT='"$(CURDIR)"'
ROOT_CFLAGS := $(shell $(ROOT_CONFIG) --cflags)
ROOT_LIBS := $(shell $(ROOT_CONFIG) --libs)
LDLIBS += $(ROOT_LIBS) -lMinuit2 -lHistPainter

SOURCES := $(wildcard src/*.cpp)
OBJECTS := $(patsubst src/%.cpp,build/%.o,$(SOURCES))
APPS := trc_to_root average_waveforms scan_amplitude scan_xy_maps peak_timing propagation_xy_maps
TARGETS := $(addprefix bin/,$(APPS))
TEST_SOURCES := $(wildcard tests/test_*.cpp)
TESTS := $(patsubst tests/%.cpp,build/%,$(TEST_SOURCES))

.PHONY: all test test-propagation test-template check-template-fit regression clean FORCE
all: $(TARGETS)

build/source_revision.h: FORCE
	python3 tools/source_revision.py $@

build/result_writer.o: build/source_revision.h

build/%.o: src/%.cpp
	@mkdir -p build
	$(CXX) $(CPPFLAGS) $(ROOT_CFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

bin/%: app/%.cpp $(OBJECTS)
	@mkdir -p bin
	$(CXX) $(CPPFLAGS) $(ROOT_CFLAGS) $(CXXFLAGS) -MMD -MP -MF build/app_$*.d -MT $@ $< $(OBJECTS) -o $@ $(LDLIBS)

build/test_%: tests/test_%.cpp $(OBJECTS)
	$(CXX) $(CPPFLAGS) $(ROOT_CFLAGS) $(CXXFLAGS) -MMD -MP -MF build/test_$*.d -MT $@ $< $(OBJECTS) -o $@ $(LDLIBS)

test: $(TESTS)
	@set -e; for test in $(TESTS); do ./$$test; done

regression: all
	python3 tests/regression.py

test-template: all
	python3 tests/regression_template.py

# Optional inspection utility, deliberately not part of APPS/all.
check-template-fit: build/check_initial_fit

build/check_initial_fit: checks/template_fit/check_initial_fit.cpp $(OBJECTS) build/source_revision.h
	$(CXX) $(CPPFLAGS) $(ROOT_CFLAGS) $(CXXFLAGS) -MMD -MP -MF build/check_initial_fit.d -MT $@ $< $(OBJECTS) -o $@ $(LDLIBS)

test-propagation: all
	python3 tests/regression_propagation.py

clean:
	rm -f build/*.o build/*.d build/test_* $(TARGETS)

-include $(wildcard build/*.d)
