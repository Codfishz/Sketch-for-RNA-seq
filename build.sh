#!/bin/bash

# Set the project root directory to the current working directory.
PROJECT_ROOT=$(pwd)

# Set the source directory (where .cpp files are located) and the include directory (where header files are located).
SRC_DIR="$PROJECT_ROOT/src"
INCLUDE_DIR="$PROJECT_ROOT/include"

# Set the output directory where the build (executable and object files) will be placed.
OUTPUT_DIR="$PROJECT_ROOT/build"

# Create the build directory if it does not exist.
mkdir -p $OUTPUT_DIR

# Set the C++ compiler and the compilation options.
# - std=c++17: use the C++17 standard.
# - -I$INCLUDE_DIR: add the project include directory.
# - -I/usr/local/include: include additional include directories if needed (e.g., for external libraries).
# - -g: generate debug information.
# - -pg: enable profiling support (for gprof).
CXX=g++
CXXFLAGS="-std=c++17 -I$INCLUDE_DIR -I/usr/local/include -g -pg"

# Find all the source (.cpp) files in the SRC_DIR directory (and its subdirectories).
SOURCES=$(find $SRC_DIR -name '*.cpp')

# Set the name (and path) for the output executable.
OUTPUT="$OUTPUT_DIR/test"

# Compile all source files and link with the nthash library.
# -L/usr/local/lib: specify the directory where the nthash library is located.
# -lnthash: link against the nthash library.
$CXX $CXXFLAGS $SOURCES -L/usr/local/lib -lnthash -o $OUTPUT

# Check whether the compilation was successful.
if [ $? -eq 0 ]; then
    echo "Build successful. The executable is located at $OUTPUT_DIR/test"
    # echo "Running the program..."
    
    # # Run the program in "index" mode:
    # # Mode: -o index, followed by reference genome FASTA file and the desired output index file.
    # $OUTPUT -o index "$PROJECT_ROOT/Test_Data/gencode.v45.transcripts.fa" "$PROJECT_ROOT/Test_Data/gencode_v45.index"
    
    # # Run the program in "quant" mode:
    # # Mode: -o quant, followed by the index file, the FASTQ reads file, and the desired output CSV file.
    # $OUTPUT -o quant "$PROJECT_ROOT/Test_Data/gencode_v45.index" "$PROJECT_ROOT/Test_Data/sd_02_099.fastq" "$PROJECT_ROOT/output_0406_02_31_0.9frac.csv"
    
    # # Uncomment the following line if you want to run another quant mode with different parameters:
    # # $OUTPUT -o quant "$PROJECT_ROOT/Test_Data/gencode.v45.index" "$PROJECT_ROOT/Test_Data/sd_02_099.fastq" "$PROJECT_ROOT/output_0324_02_67_0.75frac.csv"
    
    # Check if the profiling output file (gmon.out) exists.
    if [ -f gmon.out ]; then
        echo "Profiling data generated. Generating analysis report..."
        # Generate a profiling analysis report using gprof.
        gprof $OUTPUT gmon.out > $OUTPUT_DIR/analysis.txt
        echo "Analysis report generated at $OUTPUT_DIR/analysis.txt"
    else
        echo "gmon.out file not found. Profiling might not have been successful."
    fi
else
    echo "Build failed."
fi
