CXX      := g++
CXXFLAGS := -std=c++20 -O3 -march=native -Wall -Wextra -pedantic

TARGET_BIN := insan_calc
TARGET_LIB := libinsan_calc.so
SRC        := insan_calc_core.cpp

all: $(TARGET_BIN) $(TARGET_LIB)

$(TARGET_BIN): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET_BIN)

$(TARGET_LIB): $(SRC)
	$(CXX) $(CXXFLAGS) -shared -fPIC $(SRC) -o $(TARGET_LIB)

clean:
	rm -f $(TARGET_BIN) $(TARGET_LIB)

.PHONY: all clean
