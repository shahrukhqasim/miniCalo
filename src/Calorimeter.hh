//
// Created by Shah Rukh Qasim on 07.01.22.
//

#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>


#ifndef B4A_CALORIMETER_HH
#define B4A_CALORIMETER_HH

void initialize(long rseed, std::string detector_specs);
void wrap_up();
pybind11::dict simulate_pu();

pybind11::dict simulate_particle(double position_x, double position_y, double position_z, double direction_x,
                           double direction_y,
                           double direction_z, int pdgid, double energy);

pybind11::dict dict_check(int number_of_events, long rseed, std::string detector_specs);

#endif //B4A_CALORIMETER_HH
