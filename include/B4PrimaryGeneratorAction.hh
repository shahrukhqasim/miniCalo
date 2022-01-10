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
// $Id: B4PrimaryGeneratorAction.hh 94808 2015-12-10 08:22:26Z gcosmo $
// 
/// \file B4PrimaryGeneratorAction.hh
/// \brief Definition of the B4PrimaryGeneratorAction class

#ifndef B4PrimaryGeneratorAction_h
#define B4PrimaryGeneratorAction_h 1


#include "defines.h"

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include <vector>
#include "json/json.h"
#include "B4PartGeneratorBase.hh"

class G4ParticleGun;

class G4Event;

/// The primary generator action class with particle gum.
///
/// It defines a single particle which hits the calorimeter 
/// perpendicular to the input face. The type of the particle
/// can be changed via the G4 build-in commands of G4ParticleGun class 
/// (see the macros provided with this example).



class B4PrimaryGeneratorAction : public B4PartGeneratorBase {
public:
    B4PrimaryGeneratorAction(std::string jsonfile);

    virtual ~B4PrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event *event);

    virtual void GeneratePrimariesx(G4Event *event);

    // set methods
    void SetRandomFlag(G4bool value);

    G4ParticleGun *getGun() { return fParticleGun; }


    virtual bool isJetGenerator() { return false; }


    std::vector<G4String> generateAvailableParticles() const;

    particles getParticle() const { return particleid_; }

    int isParticle(int i) const {
        return i == particleid_;
    }

    G4String getParticleName(enum particles) const;

private:
    G4ParticleGun *fParticleGun; // G4 particle gun

    G4String setParticleID(enum particles);
    G4String setParticleID(std::string particles);

    particles particleid_;
    bool from_beamspot_;
    Json::Value m_particles;


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
