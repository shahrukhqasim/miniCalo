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
B4DetectorConstruction * detector;

void wrap_up() {
    delete runManager;
}
void initialize_test(long rseed, std::string detector_specs) {
    rseed++;

    B4PartGeneratorBase::seedsoffset_ = 800 * rseed;
    G4cout << "random seed " << rseed << G4endl;

    runManager = new G4RunManager;

    Json::Value root(detector_specs);

    auto layers_upper = root["layers"];

}

void initialize(long rseed, std::string detector_specs) {
    rseed++;

    B4PartGeneratorBase::seedsoffset_ = 800 * rseed;
    G4cout << "random seed " << rseed << G4endl;

    runManager = new G4RunManager;

    Json::Value root(detector_specs);
    Json::Reader reader;
    reader.parse(detector_specs, root);

    detector = new B4DetectorConstruction(root);
    runManager->SetUserInitialization(detector);

    auto physicsList = new FTFP_BERT;
    physicsList->RegisterPhysics(new G4StepLimiterPhysics());
    runManager->SetUserInitialization(physicsList);

    eventAction = new B4aEventAction(std::string(), false);
    primariesGenerator = new PrimariesGenerator;
    auto actionInitialization = new B4aActionInitialization(detector, eventAction, primariesGenerator);
    runManager->SetUserInitialization(actionInitialization);


    // Get the pointer to the User Interface manager
    ui_manager = G4UImanager::GetUIpointer();

    ui_manager->ApplyCommand(std::string("/run/initialize"));
    ui_manager->ApplyCommand(std::string("/run/printProgress 100"));
}

py::dict collect() {
    const std::vector<double>& rechit_energy = eventAction->getRechitEnergy();
    const std::vector<double>& rechit_x= eventAction->getRechitX();
    const std::vector<double>& rechit_y= eventAction->getRechitY();
    const std::vector<double>& rechit_z= eventAction->getRechitZ();
    const std::vector<int>& rechit_stupid_pid= eventAction->getRechitStupidPid();

    std::vector<int> rechit_idx;
    std::vector<double> rechit_energy_2;
    std::vector<double> rechit_x_2;
    std::vector<double> rechit_y_2;
    std::vector<double> rechit_z_2;
    std::vector<int> rechit_stupid_pid_2;

    for (int i=0;i<rechit_energy.size();i++) {
        if (rechit_energy[i]>0) {
            rechit_idx.push_back(i);
            rechit_energy_2.push_back(rechit_energy[i]);
            rechit_x_2.push_back(rechit_x[i]);
            rechit_y_2.push_back(rechit_y[i]);
            rechit_z_2.push_back(rechit_z[i]);
            rechit_stupid_pid_2.push_back(rechit_stupid_pid[i]);
        }
    }

    py::array rechitEnergy = py::cast(rechit_energy_2);
    py::array rechitX = py::cast(rechit_x_2);
    py::array rechitY = py::cast(rechit_y_2);
    py::array rechitZ = py::cast(rechit_z_2);
    py::array rechitIdx = py::cast(rechit_idx);
    py::array rechitStupidPid = py::cast(rechit_stupid_pid_2);

    auto&hits_particles = eventAction->getHitsParticles();
    auto&hits_deposits = eventAction->getHitsDeposits();
    auto&hits_sensor_idx = eventAction->getHitsParticlesSensorIdx();

    std::vector<int> hits_particles_2;
    std::vector<double> hits_deposits_2;
    std::vector<int> hits_sensor_idx_2;

    std::unordered_map<std::string, int> ps_to_idx;
    int max_k = 0;
    for(int i=0;i < hits_particles.size();i++) {
        if(hits_deposits[i]==0) {
            continue;
        }
        std::string k = std::to_string(hits_particles[i]) + std::string("_") + std::to_string(hits_sensor_idx[i]);
        if (ps_to_idx.find(k) != ps_to_idx.end()) {
            hits_deposits_2[ps_to_idx[k]] += hits_deposits[i];
        }
        else {
            hits_particles_2.push_back(hits_particles[i]);
            hits_deposits_2.push_back(hits_deposits[i]);
            hits_sensor_idx_2.push_back(hits_sensor_idx[i]);
            ps_to_idx[k] = max_k;
            max_k += 1;
        }
    }

    py::array np_hits_particles_2  = py::cast(hits_particles_2);
    py::array np_hits_deposits_2  = py::cast(hits_deposits_2);
    py::array np_hits_sensor_idx_2  = py::cast(hits_sensor_idx_2);

    py::array hitsParticles = py::cast(eventAction->getHitsParticles());
    py::array hitsDeposits = py::cast(eventAction->getHitsDeposits());
    py::array hitsSensorIndex = py::cast(eventAction->getHitsParticlesSensorIdx());

    std::vector<int> a;
    a.push_back(1);
    a.push_back(2);
    a.push_back(3);

    py::array particles_vertex_position_x = py::cast(eventAction->particles_caught.particles_vertex_position_x);
    py::array particles_vertex_position_y = py::cast(eventAction->particles_caught.particles_vertex_position_y);
    py::array particles_vertex_position_z = py::cast(eventAction->particles_caught.particles_vertex_position_z);
    py::array particles_momentum_direction_x = py::cast(eventAction->particles_caught.particles_momentum_direction_x);
    py::array particles_momentum_direction_y = py::cast(eventAction->particles_caught.particles_momentum_direction_y);
    py::array particles_momentum_direction_z = py::cast(eventAction->particles_caught.particles_momentum_direction_z);
    py::array particles_kinetic_energy = py::cast(eventAction->particles_caught.particles_kinetic_energy);
    py::array particles_pdgid = py::cast(eventAction->particles_caught.particles_pdgid);
    py::array particles_total_energy_deposited_active = py::cast(eventAction->particles_caught.particles_total_energy_deposited_active);
    py::array particles_total_energy_deposited_all = py::cast(eventAction->particles_caught.particles_total_energy_deposited_all);
    py::array particles_parent_id = py::cast(eventAction->particles_caught.particles_parent_idx);

//    delete ev;
    py::dict d = py::dict("rechit_energy"_a =rechitEnergy,
                          "rechit_x"_a = rechitX,
                          "rechit_y"_a = rechitY,
                          "rechit_z"_a = rechitZ,
                          "rechit_idx"_a = rechitIdx,
                          "rechit_stupid_pid"_a = rechitStupidPid,

//                          "rechit_particle_id"_a = hitsParticles,
//                          "rechit_particle_deposit"_a = hitsDeposits,
//                          "rechit_particle_sensor_idx"_a = hitsSensorIndex,

                          "hit_particle_id"_a = np_hits_particles_2,
                          "hit_particle_deposit"_a = np_hits_deposits_2,
                          "hit_particle_sensor_idx"_a = np_hits_sensor_idx_2,

//                          "particle_positions_x"_a = (py::array) py::cast(eventAction->positions_tracked_x),
//                          "particle_positions_y"_a = (py::array) py::cast(eventAction->positions_tracked_y),
//                          "particle_positions_z"_a = (py::array) py::cast(eventAction->positions_tracked_z),
//                          "particle_positions_particle_id"_a = (py::array) py::cast(eventAction->positions_tracked_particle_index),

                          "particles_vertex_position_x"_a = particles_vertex_position_x,
                          "particles_vertex_position_y"_a = particles_vertex_position_y,
                          "particles_vertex_position_z"_a = particles_vertex_position_z,
                          "particles_momentum_direction_x"_a = particles_momentum_direction_x,
                          "particles_momentum_direction_y"_a = particles_momentum_direction_y,
                          "particles_momentum_direction_z"_a = particles_momentum_direction_z,
                          "particles_kinetic_energy"_a = particles_kinetic_energy,
                          "particles_pdgid"_a = particles_pdgid,
                          "particles_total_energy_deposited_active"_a = particles_total_energy_deposited_active,
                          "particles_total_energy_deposited_all"_a = particles_total_energy_deposited_all,
                          "particles_parent_idx"_a = particles_parent_id
    );

    for (auto& it: eventAction->result_arrays_int) {
        py::array np = py::cast(*it.second);
        d[it.first.c_str()] = np;
    }
    for (auto& it: eventAction->result_arrays_double) {
        py::array np = py::cast(*it.second);
        d[it.first.c_str()] = np;
    }

    return d;
}

py::dict simulate_pu() {
    if (ui_manager == nullptr) {
        G4cout<<"Forgot to call initialize?\n";
        throw std::runtime_error("Forgot to call initialize?");
    }

    primariesGenerator->SetNextToPU();
    ui_manager->ApplyCommand(std::string("/run/beamOn ") + std::to_string(1));

    return collect();
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
    return collect();
}


py::dict get_sensor_data() {
    if (ui_manager == nullptr) {
        G4cout<<"Forgot to call initialize?\n";
        throw std::runtime_error("Forgot to call initialize?");
    }

    const auto &activesensors = detector->getActiveSensors();

    std::vector<double> rechit_x;
    std::vector<double> rechit_y;
    std::vector<double> rechit_z;

    for(int i =0;i<activesensors->size();i++) {
        rechit_x.push_back(activesensors->at(i)->getPosx());
        rechit_y.push_back(activesensors->at(i)->getPosy());
        rechit_z.push_back(activesensors->at(i)->getPosz());
    }

    py::dict d = py::dict("rechit_x"_a = (py::array)py::cast(rechit_x),
                          "rechit_y"_a = (py::array) py::cast(rechit_y),
                          "rechit_z"_a = (py::array) py::cast(rechit_z));

    return d;

}

pybind11::dict dict_check(int number_of_events, long rseed, std::string detector_specs) {
//    py::dict d;
    std::vector<float> rechitEnergy;
    rechitEnergy.push_back(1.0);
    rechitEnergy.push_back(2.0);
    rechitEnergy.push_back(3.0);


    std::vector<float> rechitDummy;
    rechitDummy.push_back(1.0);
    rechitDummy.push_back(2.2);
    rechitDummy.push_back(3.0);

    std::string name("dummy_dummy");

    py::dict d = py::dict("rechit_energy"_a = rechitEnergy,
                          py::arg(name.c_str()) = rechitDummy
                          );

    d["damn"] = rechitDummy;

    return d;
}
