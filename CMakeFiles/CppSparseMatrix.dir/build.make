# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/quentin/PFE/CppCode/cppmpi

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/quentin/PFE/CppCode/cppmpi

# Include any dependencies generated for this target.
include CMakeFiles/CppSparseMatrix.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/CppSparseMatrix.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/CppSparseMatrix.dir/flags.make

CMakeFiles/CppSparseMatrix.dir/matrix.cpp.o: CMakeFiles/CppSparseMatrix.dir/flags.make
CMakeFiles/CppSparseMatrix.dir/matrix.cpp.o: matrix.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/quentin/PFE/CppCode/cppmpi/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/CppSparseMatrix.dir/matrix.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/CppSparseMatrix.dir/matrix.cpp.o -c /home/quentin/PFE/CppCode/cppmpi/matrix.cpp

CMakeFiles/CppSparseMatrix.dir/matrix.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CppSparseMatrix.dir/matrix.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/quentin/PFE/CppCode/cppmpi/matrix.cpp > CMakeFiles/CppSparseMatrix.dir/matrix.cpp.i

CMakeFiles/CppSparseMatrix.dir/matrix.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CppSparseMatrix.dir/matrix.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/quentin/PFE/CppCode/cppmpi/matrix.cpp -o CMakeFiles/CppSparseMatrix.dir/matrix.cpp.s

CMakeFiles/CppSparseMatrix.dir/matrix.cpp.o.requires:

.PHONY : CMakeFiles/CppSparseMatrix.dir/matrix.cpp.o.requires

CMakeFiles/CppSparseMatrix.dir/matrix.cpp.o.provides: CMakeFiles/CppSparseMatrix.dir/matrix.cpp.o.requires
	$(MAKE) -f CMakeFiles/CppSparseMatrix.dir/build.make CMakeFiles/CppSparseMatrix.dir/matrix.cpp.o.provides.build
.PHONY : CMakeFiles/CppSparseMatrix.dir/matrix.cpp.o.provides

CMakeFiles/CppSparseMatrix.dir/matrix.cpp.o.provides.build: CMakeFiles/CppSparseMatrix.dir/matrix.cpp.o


# Object files for target CppSparseMatrix
CppSparseMatrix_OBJECTS = \
"CMakeFiles/CppSparseMatrix.dir/matrix.cpp.o"

# External object files for target CppSparseMatrix
CppSparseMatrix_EXTERNAL_OBJECTS =

CppSparseMatrix: CMakeFiles/CppSparseMatrix.dir/matrix.cpp.o
CppSparseMatrix: CMakeFiles/CppSparseMatrix.dir/build.make
CppSparseMatrix: CMakeFiles/CppSparseMatrix.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/quentin/PFE/CppCode/cppmpi/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable CppSparseMatrix"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/CppSparseMatrix.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/CppSparseMatrix.dir/build: CppSparseMatrix

.PHONY : CMakeFiles/CppSparseMatrix.dir/build

CMakeFiles/CppSparseMatrix.dir/requires: CMakeFiles/CppSparseMatrix.dir/matrix.cpp.o.requires

.PHONY : CMakeFiles/CppSparseMatrix.dir/requires

CMakeFiles/CppSparseMatrix.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/CppSparseMatrix.dir/cmake_clean.cmake
.PHONY : CMakeFiles/CppSparseMatrix.dir/clean

CMakeFiles/CppSparseMatrix.dir/depend:
	cd /home/quentin/PFE/CppCode/cppmpi && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/quentin/PFE/CppCode/cppmpi /home/quentin/PFE/CppCode/cppmpi /home/quentin/PFE/CppCode/cppmpi /home/quentin/PFE/CppCode/cppmpi /home/quentin/PFE/CppCode/cppmpi/CMakeFiles/CppSparseMatrix.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/CppSparseMatrix.dir/depend

