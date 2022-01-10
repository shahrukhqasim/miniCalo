//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: B4aEventAction.hh 75215 2013-10-29 16:07:06Z gcosmo $
// 
/// \file B4aEventAction.hh
/// \brief Definition of the B4aEventAction class

#ifndef B4aEventAction_h
#define B4aEventAction_h 1

#include "defines.h"

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "B4PartGeneratorBase.hh"
#include "G4StepPoint.hh"
#include "sensorContainer.h"
#include "B4DetectorConstruction.hh"
#include "G4Step.hh"
#include "B4RunAction.hh"
#include "B4PrimaryGeneratorAction.hh"

/// Event action class
///
/// It defines data members to hold the energy deposit and track lengths
/// of charged particles in Absober and Gap layers:
/// - fEnergyAbs, fEnergyGap, fTrackLAbs, fTrackLGap
/// which are collected step by step via the functions
/// - AddAbs(), AddGap()


class SOAParticles {
public:
    std::vector<double> particles_vertex_position_x;
    std::vector<double> particles_vertex_position_y;
    std::vector<double> particles_vertex_position_z;

    std::vector<double> particles_momentum_direction_x;
    std::vector<double> particles_momentum_direction_y;
    std::vector<double> particles_momentum_direction_z;

    std::vector<double> particles_kinetic_energy;
    std::vector<int> particles_pdgid;

    std::vector<double> particles_total_energy_deposited_active;
    std::vector<double> particles_total_energy_deposited_all;

    std::unordered_map<int, int> trackid_to_idx;
};


class G4VPhysicalVolume;
class B4aEventAction : public G4UserEventAction
{
	friend B4RunAction;
  public:
    B4aEventAction(std::string output_folder_name, bool do_root);
    virtual ~B4aEventAction();

    virtual void  BeginOfEventAction(const G4Event* event);
    virtual void    EndOfEventAction(const G4Event* event);
    
    void AddEnergy(G4double de, G4double dl);


    template<typename T>
    void writeToBinaryFile(std::string file, std::vector<T> &x, int id, int num_elements);

    void accumulateStepData(G4VPhysicalVolume *volume, const G4Step *step);

    void clear(){
    	rechit_energy_.clear();
    	allvolumes_.clear();
    	rechit_absorber_energy_.clear();
        rechit_x_.clear();
        rechit_y_.clear();
        rechit_z_.clear();
        rechit_eta_.clear();
        rechit_phi_.clear();
        rechit_vxy_.clear();
        rechit_layer_.clear();
        rechit_detid_.clear();
        nsteps_=0;
        totalen_=0;
    }


    void setGenerator(B4PartGeneratorBase  * generator){
    	generator_=generator;
    	navail_parts = generator_->generateAvailableParticles().size();
    }
    void setDetector(B4DetectorConstructionBase *detector){
    	detector_=detector;
    }

  private:
    G4double  fEnergyAbs;
    G4double totalen_;
public:
    std::vector<float>  rechit_energy_,rechit_absorber_energy_;
    std::vector<float>  rechit_x_;
    std::vector<float>  rechit_y_;
    std::vector<float>  rechit_z_;
    std::vector<float>  rechit_layer_;
    std::vector<float>  rechit_eta_;
    std::vector<float>  rechit_phi_;
    std::vector<float>  rechit_vxy_;
public:
    const std::vector<float> &getRechitEnergy() const;

    const std::vector<float> &getRechitX() const;

    const std::vector<float> &getRechitY() const;

    const std::vector<float> &getRechitZ() const;

private:
    std::vector<int>       rechit_detid_;
    std::vector<const G4VPhysicalVolume * > allvolumes_;


//    std::vector<G4ThreeVector> particles_vertex_position;
//    std::vector<G4ThreeVector> particles_position;
//    std::vector<double> particles_kinetic_energy;
//    std::vector<G4ThreeVector> particles_momentum_direction;
//    std::vector<int> particles_pdgid;

//    std::vector<double> particles_total_energy_deposited;
//    std::vector<double> particles_total_energy_deposited_2;
//    std::unordered_map<int, int> particles_buckets;
//    std::unordered_map<int, int> tracks_buckets;

    SOAParticles particlesEntering;

    std::vector<int> hits_particles_id;
    std::vector<float> hits_particles_deposits;
public:
    const std::vector<float> &getHitsParticlesSensorIdx() const;

private:
    std::vector<float> hits_particles_sensor_idx;
public:
    const std::vector<int> &getHitsParticles() const;

    const std::vector<float> &getHitsDeposits() const;

public:
    int particle_index_not_found = 0;
    double total_deposit = 0;
    double total_deposit_2 = 0;
private:

    std::vector<int>  rechit_idx__;

    G4double  fEnergyGap;
    G4double  fTrackLAbs; 
    G4double  fTrackLGap;


    B4PartGeneratorBase * generator_;
    B4DetectorConstructionBase * detector_;

    size_t nsteps_;
    size_t navail_parts;

    std::string output_bin_folder;

    bool do_root;

};

// inline functions



inline void B4aEventAction::AddEnergy(G4double de, G4double dl) {
  fEnergyGap += de; 
  fTrackLGap += dl;
}
                     
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
