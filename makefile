# # 编译器及选项
# CXX = g++
# CXXFLAGS = -std=c++17 -O2 -fopenmp -Iinclude \
#            -I/opt/intel/oneapi/mkl/latest/include  # MKL 头文件路径

# # 链接器选项
# LDFLAGS = -fopenmp
# LDLIBS = -L/opt/intel/oneapi/mkl/latest/lib/intel64 \
#          -lmkl_rt -lpthread -lm -ldl  # MKL 动态库
# RPATH = -Wl,-rpath=/opt/intel/oneapi/mkl/latest/lib/intel64  # 运行时库路径

# # 可执行文件
# TARGET = build/bem_solver

# # 所有源文件（可以扩展）
# SRCDIR = src
# SRCFILES = $(wildcard $(SRCDIR)/*.cpp)
# OBJFILES = $(patsubst $(SRCDIR)/%.cpp, build/%.o, $(SRCFILES))

# # 默认目标
# all: $(TARGET)

# # 编译可执行文件
# $(TARGET): $(OBJFILES)
# 	@echo "Linking executable with MKL..."
# 	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $@ $(LDLIBS) $(RPATH)

# # 对每个 .cpp 编译成 .o
# build/%.o: $(SRCDIR)/%.cpp
# 	@mkdir -p build
# 	@echo "Compiling $<..."
# 	$(CXX) $(CXXFLAGS) -c $< -o $@

# # 清理命令
# clean:
# 	@echo "Cleaning build files..."
# 	rm -rf build/*.o $(TARGET)

# # 验证 MKL 链接的测试目标
# test-mkl:
# 	@echo "#include <mkl.h>\n#include <stdio.h>\nint main() { printf(\"MKL Version: %d.%d.%d\\n\", __INTEL_MKL__, __INTEL_MKL_MINOR__, __INTEL_MKL_UPDATE__); return 0; }" > test_mkl.c
# 	$(CXX) test_mkl.c -o test_mkl $(LDLIBS) $(RPATH)
# 	./test_mkl
# 	@rm test_mkl.c test_mkl

# # 运行 Python 脚本（保留原始功能）
# .PHONY: print
# print:
# 	@echo "=== Running all Python scripts in py/ ==="
# 	@for f in py/*.py; do \
# 	  echo ">>> $$f"; \
# 	  python3 "$$f"; \
# 	done

# .PHONY: all clean test-mkl print

# 编译器及选项
CXX       = g++
# CXXFLAGS  = -std=c++17 -O2 -fopenmp -Iinclude -DEIGEN_USE_LAPACKE
CXXFLAGS  = -std=c++17 -O2 -fopenmp -Iinclude -DEIGEN_USE_LAPACKE \
			-I$(CURDIR)/include/eigen-3.4.0 \
			-I$(CURDIR)/include/eigen-3.4.0/spectra
LDFLAGS   = -llapacke -llapack -lblas

# 可执行文件
TARGET    = build/bem_solver

SRCDIR    = src
SRCFILES  = $(wildcard $(SRCDIR)/*.cpp)
OBJFILES  = $(patsubst $(SRCDIR)/%.cpp, build/%.o, $(SRCFILES))

# 测试
TEST_SRC    = src/test.cpp
TEST_OBJ    = build/test.o
TEST_TARGET = build/test_solver

all: $(TARGET)

# 链接主程序，注意加上 $(LDFLAGS)
$(TARGET): $(OBJFILES)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# 编译每个 .cpp 为 .o
build/%.o: $(SRCDIR)/%.cpp
	@mkdir -p build
	$(CXX) $(CXXFLAGS) -c $< -o $@

# 测试程序的编译和链接
$(TEST_OBJ): $(TEST_SRC)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TEST_TARGET): $(TEST_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

test: $(TEST_TARGET)

clean:
	rm -rf build/*.o $(TARGET) $(TEST_TARGET)
	rm -rf data/*

.PHONY: print
print:
	@echo "=== Running all Python scripts in py/ ==="
	@for f in py/*.py; do \
	  echo ">>> $$f"; \
	  python3 "$$f"; \
	done