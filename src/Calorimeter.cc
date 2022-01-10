//
// Created by Shah Rukh Qasim on 07.01.22.
//

#include "Calorimeter.hh"

#include "B4DetectorConstruction.hh"
#include "B4aActionInitialization.hh"
#include <string>

#include "defines.h"

#ifdef G4MULTITHREADED
#undef G4MULTITHREADED
#endif

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else

#include "G4RunManager.hh"

#endif

#include<string>

#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "FTFP_BERT.hh"
#include "QGSP_FTFP_BERT.hh"

#include "Randomize.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4RandomTools.hh"
#include "B4DetectorConstructionFP.hh"

//limiter stuff

#include "G4PhysicsListHelper.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4UserSpecialCuts.hh"
#include "G4StepLimiter.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "B4PartGeneratorBase.hh"
#include "B4aEventAction.hh"
#include "PrimariesGenerator.hh"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/embed.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace py::literals;



G4RunManager* runManager;
B4aEventAction* eventAction;
G4UImanager *ui_manager;
PrimariesGenerator *primariesGenerator;

void wrap_up() {
    delete runManager;
}

void initialize(long rseed, std::string detector_specs) {
    rseed++;

    B4PartGeneratorBase::seedsoffset_ = 800 * rseed;
    G4cout << "random seed " << rseed << G4endl;

    runManager = new G4RunManager;

    Json::Value root;
    std::ifstream ifs;
    ifs.open(detector_specs);

    Json::CharReaderBuilder builder;
    builder["collectComments"] = true;
    JSONCPP_STRING errs;
    if (!parseFromStream(builder, ifs, &root, &errs)) {
        std::cout << errs << std::endl;
        throw std::runtime_error("Error in parsing detector_specs");
    }
    auto detConstruction = new B4DetectorConstruction(root);
    runManager->SetUserInitialization(detConstruction);

    auto physicsList = new FTFP_BERT;
    physicsList->RegisterPhysics(new G4StepLimiterPhysics());
    runManager->SetUserInitialization(physicsList);

    eventAction = new B4aEventAction(std::string(), false);
    primariesGenerator = new PrimariesGenerator;
    auto actionInitialization = new B4aActionInitialization(detConstruction, eventAction, primariesGenerator);
    runManager->SetUserInitialization(actionInitialization);


    // Get the pointer to the User Interface manager
    ui_manager = G4UImanager::GetUIpointer();

    ui_manager->ApplyCommand(std::string("/run/initialize"));
    ui_manager->ApplyCommand(std::string("/run/printProgress 100"));
}

py::dict simulate_pu() {
    if (ui_manager == nullptr) {
        G4cout<<"Forgot to call initialize?\n";
        throw std::runtime_error("Forgot to call initialize?");
    }

    primariesGenerator->SetNextToPU();
    ui_manager->ApplyCommand(std::string("/run/beamOn ") + std::to_string(1));

    py::array rechitEnergy = py::cast(eventAction->getRechitEnergy());
    py::array rechitX = py::cast(eventAction->getRechitX());
    py::array rechitY = py::cast(eventAction->getRechitY());
    py::array rechitZ = py::cast(eventAction->getRechitZ());

    py::array hitsParticles = py::cast(eventAction->getHitsParticles());
    py::array hitsDeposits = py::cast(eventAction->getHitsParticles());
    py::array hitsSensorIndex = py::cast(eventAction->getHitsParticlesSensorIdx());


//    delete ev;
    py::dict d = py::dict("rechit_energy"_a =rechitEnergy,
                          "rechit_x"_a = rechitX,
                          "rechit_y"_a = rechitY,
                          "rechit_z"_a = rechitZ,
                          "rechit_particle_id"_a = hitsParticles,
                          "rechit_particle_deposit"_a = hitsDeposits,
                          "rechit_particle_sensor_idx"_a = hitsSensorIndex);
    return d;
}

py::dict simulate_particle(double position_x, double position_y, double position_z, double direction_x,
                           double direction_y,
                           double direction_z, int pdgid, double energy) {
    if (ui_manager == nullptr) {
        G4cout<<"Forgot to call initialize?\n";
        throw std::runtime_error("Forgot to call initialize?");
    }

    primariesGenerator->SetNextToParticle(position_x,position_y,position_z,direction_x,direction_y,direction_z,pdgid,energy);
    ui_manager->ApplyCommand(std::string("/run/beamOn ") + std::to_string(1));
    std::vector<float> rechitEnergy = eventAction->getRechitEnergy();
    std::vector<float> rechitX = eventAction->getRechitX();
    std::vector<float> rechitY = eventAction->getRechitY();
    std::vector<float> rechitZ = eventAction->getRechitZ();

//    delete ev;
    py::dict d = py::dict("rechit_energy"_a = rechitEnergy,
                          "rechit_x"_a = rechitX,
                          "rechit_y"_a = rechitY,
                          "rechit_z"_a = rechitZ);
    return d;
}


pybind11::dict dict_check(int number_of_events, long rseed, std::string detector_specs) {
//    py::dict d;
    std::vector<float> rechitEnergy;
    rechitEnergy.push_back(1.0);
    rechitEnergy.push_back(2.0);
    rechitEnergy.push_back(3.0);

    py::dict d = py::dict("rechit_energy"_a = rechitEnergy);

    return d;
}
