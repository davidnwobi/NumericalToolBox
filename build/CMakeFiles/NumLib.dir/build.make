# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.27

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\CMake\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\CMake\bin\cmake.exe" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = C:\Users\user\OneDrive\Desktop\object_oriented_numerical_anaylsis\NumericalToolbox

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = C:\Users\user\OneDrive\Desktop\object_oriented_numerical_anaylsis\NumericalToolbox\build

# Include any dependencies generated for this target.
include CMakeFiles/NumLib.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/NumLib.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/NumLib.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/NumLib.dir/flags.make

# Object files for target NumLib
NumLib_OBJECTS =

# External object files for target NumLib
NumLib_EXTERNAL_OBJECTS =

libNumLib.a: CMakeFiles/NumLib.dir/build.make
libNumLib.a: CMakeFiles/NumLib.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=C:\Users\user\OneDrive\Desktop\object_oriented_numerical_anaylsis\NumericalToolbox\build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Linking CXX static library libNumLib.a"
	$(CMAKE_COMMAND) -P CMakeFiles\NumLib.dir\cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\NumLib.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/NumLib.dir/build: libNumLib.a
.PHONY : CMakeFiles/NumLib.dir/build

CMakeFiles/NumLib.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\NumLib.dir\cmake_clean.cmake
.PHONY : CMakeFiles/NumLib.dir/clean

CMakeFiles/NumLib.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" C:\Users\user\OneDrive\Desktop\object_oriented_numerical_anaylsis\NumericalToolbox C:\Users\user\OneDrive\Desktop\object_oriented_numerical_anaylsis\NumericalToolbox C:\Users\user\OneDrive\Desktop\object_oriented_numerical_anaylsis\NumericalToolbox\build C:\Users\user\OneDrive\Desktop\object_oriented_numerical_anaylsis\NumericalToolbox\build C:\Users\user\OneDrive\Desktop\object_oriented_numerical_anaylsis\NumericalToolbox\build\CMakeFiles\NumLib.dir\DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/NumLib.dir/depend

