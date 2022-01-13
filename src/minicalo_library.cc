//
// Created by Shah Rukh Qasim on 07.01.22.
//
#include <pybind11/pybind11.h>
#include "Calorimeter.hh"

namespace py = pybind11;

float cpp_test_function(float a, float b) {
    return a*2 + b;
}


py::tuple cpp_test_function_2(float a, float b) {
    return py::make_tuple(a*2+1, a*3+b);
}



PYBIND11_MODULE(minicalo, handle){
    handle.doc() = "Test c++ function";
    handle.def("cpp_test_function", &cpp_test_function);
    handle.def("cpp_test_function_2", &cpp_test_function_2);
    handle.def("dict_check", &dict_check);

    handle.def("simulate_pu", &simulate_pu);
    handle.def("simulate_particle", &simulate_particle);
    handle.def("get_sensor_data", &get_sensor_data);
    handle.def("initialize", &initialize);
    handle.def("initialize_test", &initialize_test);
    handle.def("wrap_up", &wrap_up);
}