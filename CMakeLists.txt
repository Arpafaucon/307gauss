cmake_minimum_required (VERSION 2.6)
project (gauss)



add_executable(bench_pc
main_pc.cpp
bench_pc.cpp
bench_common.cpp
io.cpp
gauss_var.cpp
gauss_fixed.cpp
)

add_executable(bench_hls
bench_hls.cpp
bench_pc.cpp
bench_common.cpp
io.cpp
gauss_var.cpp
gauss_fixed.cpp
gauss_hls_stream.cpp
gauss_fixed_hls.cpp
)

target_include_directories(bench_hls
PUBLIC
"/home/arpad/install/vivado/Vivado/2018.2/include"
"/usr/include/c++/7"
)

