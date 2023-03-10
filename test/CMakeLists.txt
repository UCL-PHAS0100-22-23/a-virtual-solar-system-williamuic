add_executable(tests test.cpp)
find_package(Catch2 3 REQUIRED)
target_include_directories(tests PUBLIC ../include)
target_link_libraries(tests PUBLIC Catch2::Catch2WithMain nbody_lib)

include(Catch)
catch_discover_tests(tests)