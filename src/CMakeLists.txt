add_library(nbody_lib particle.cpp)
target_compile_features(nbody_lib PUBLIC cxx_std_17)
target_include_directories(nbody_lib PUBLIC ../include)

find_package(Eigen3 3.4 REQUIRED)
find_package(OpenMP REQUIRED)

target_link_libraries(nbody_lib PUBLIC Eigen3::Eigen OpenMP::OpenMP_CXX)