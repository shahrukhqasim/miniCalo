cmake_minimum_required(VERSION 3.13)
project(miniCalo)

set(CMAKE_CXX_STANDARD 14)

include_directories(include)

add_executable(miniCalo
        include/B4aActionInitialization.hh
        include/B4aEventAction.hh
        include/B4Analysis.hh
        include/B4aSteppingAction.hh
        include/B4DetectorConstruction.hh
        include/B4PrimaryGeneratorAction.hh
        include/B4RunAction.hh
        include/sensorContainer.h
        src/B4aActionInitialization.cc
        src/B4aEventAction.cc
        src/B4aSteppingAction.cc
        src/B4DetectorConstruction.cc
        src/B4PrimaryGeneratorAction.cc
        src/B4RunAction.cc
        src/sensorContainer.cc
        exampleB4a.cc
        plotHisto.C
        plotNtuple.C)
