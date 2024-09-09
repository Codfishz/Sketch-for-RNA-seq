#!/bin/bash

# 设置项目根目录
PROJECT_ROOT=$(pwd)

# 设置源文件目录和头文件目录
SRC_DIR="$PROJECT_ROOT/src"
INCLUDE_DIR="$PROJECT_ROOT/include"
OUTPUT_DIR="$PROJECT_ROOT/build"

# 创建build目录如果它不存在
mkdir -p $OUTPUT_DIR

# 设置编译器和编译选项
CXX=g++
CXXFLAGS="-std=c++11 -I$INCLUDE_DIR -g -pg"

# 查找所有的源文件
SOURCES=$(find $SRC_DIR -name '*.cpp')

# 设置输出可执行文件的名称
OUTPUT="$OUTPUT_DIR/test"

# 编译所有源文件并生成可执行文件
$CXX $CXXFLAGS $SOURCES -o $OUTPUT

# 检查编译是否成功
if [ $? -eq 0 ]; then
    echo "Build successful. The executable is located at $OUTPUT_DIR/test"
    echo "Running the program..."
    # 运行生成的可执行文件
    $OUTPUT "$PROJECT_ROOT/Test_Data/gencode.v45.transcripts.fa" "$PROJECT_ROOT/Test_Data/sd_0249.fastq"  "$PROJECT_ROOT/Test_Data/gencode.v45.chr_patch_hapl_scaff.annotation.gtf" "$PROJECT_ROOT/output.csv"
    #$OUTPUT "$PROJECT_ROOT/Test_Data/reference_genome.fasta" "$PROJECT_ROOT/Test_Data/reads.fastq" "$PROJECT_ROOT/Test_Data/Test.gtf"
    
    # 检查是否生成了 gmon.out 文件
    if [ -f gmon.out ]; then
        echo "Profiling data generated. Generating analysis report..."
        gprof $OUTPUT gmon.out > $OUTPUT_DIR/analysis.txt  # 生成分析报告
        echo "Analysis report generated at $OUTPUT_DIR/analysis.txt"
    else
        echo "gmon.out file not found. Profiling might not have been successful."
    fi
else
    echo "Build failed."
fi

