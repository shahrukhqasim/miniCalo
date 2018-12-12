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
// $Id: B4PrimaryGeneratorAction.cc 100946 2016-11-03 11:28:08Z gcosmo $
// 
/// \file B4PrimaryGeneratorAction.cc
/// \brief Implementation of the B4PrimaryGeneratorAction class

#include "B4PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4INCLRandom.hh"
#include <G4INCLGeant4Random.hh>
#include <G4INCLRandomSeedVector.hh>
#include<ctime>
#include<sys/types.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

 B4PrimaryGeneratorAction * B4PrimaryGeneratorAction::globalgen=0;

B4PrimaryGeneratorAction::B4PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(nullptr)
{
  G4int nofParticles = 1;
  fParticleGun = new G4ParticleGun(nofParticles);

  // default particle kinematic
  //
  auto particleDefinition 
    = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
 // G4ParticleTable::GetParticleTable()->DumpTable();
  fParticleGun->SetParticleDefinition(particleDefinition);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fParticleGun->SetParticleEnergy(100.*GeV);


  G4INCL::Random::setGenerator( new G4INCL::Geant4RandomGenerator());

  globalgen=this;

  xorig_=0;
  yorig_=0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4PrimaryGeneratorAction::~B4PrimaryGeneratorAction()
{
  delete fParticleGun;
}
std::vector<G4String> B4PrimaryGeneratorAction::generateAvailableParticles(){
	std::vector<G4String> out;
	auto oldid=particleid_;
	for(int i=0;i<particles_size;i++)
		out.push_back(setParticleID((B4PrimaryGeneratorAction::particles)i));
	setParticleID(oldid);
	return out;
}

G4String B4PrimaryGeneratorAction::setParticleID(enum particles p){
	particleid_=p;
	G4ParticleDefinition * particleDefinition =0;

	if(p==elec){
		particleDefinition  = G4ParticleTable::GetParticleTable()->FindParticle("e-");
		fParticleGun->SetParticleDefinition(particleDefinition);
		return "isElectron";
	}
	if(p==muon){
		particleDefinition  = G4ParticleTable::GetParticleTable()->FindParticle("mu-");
		fParticleGun->SetParticleDefinition(particleDefinition);
		return "isMuon";
	}
	if(p==pioncharged){
		particleDefinition  = G4ParticleTable::GetParticleTable()->FindParticle("pi+");
		fParticleGun->SetParticleDefinition(particleDefinition);
		return "isPionCharged";
	}
	if(p==pionneutral){
		particleDefinition  = G4ParticleTable::GetParticleTable()->FindParticle("pi0");
		fParticleGun->SetParticleDefinition(particleDefinition);
		return "isPionNeutral";
	}
	if(p==klong){
		particleDefinition  = G4ParticleTable::GetParticleTable()->FindParticle("kaon0L");
		fParticleGun->SetParticleDefinition(particleDefinition);
		return "isK0Long";
	}
	if(p==kshort){
		particleDefinition  = G4ParticleTable::GetParticleTable()->FindParticle("kaon0S");
		fParticleGun->SetParticleDefinition(particleDefinition);
		return "isK0Short";
	}
	if(p==gamma){
		particleDefinition  = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
		fParticleGun->SetParticleDefinition(particleDefinition);
		return "isGamma";
	}

	return "isInvalid";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // This function is called at the begining of event

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get world volume 
  // from G4LogicalVolumeStore
  //
  G4double worldZHalfLength = 0.;
  auto worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");

  // Check that the world volume has box shape
  G4Box* worldBox = nullptr;
  if (  worldLV ) {
    worldBox = dynamic_cast<G4Box*>(worldLV->GetSolid());
  }

  if ( worldBox ) {
    worldZHalfLength = worldBox->GetZHalfLength();  
  }
  else  {
    G4ExceptionDescription msg;
    msg << "World volume of box shape not found." << G4endl;
    msg << "Perhaps you have changed geometry." << G4endl;
    msg << "The gun will be place in the center.";
    G4Exception("B4PrimaryGeneratorAction::GeneratePrimaries()",
      "MyCode0002", JustWarning, msg);
  } 
  
  // Set gun position

  //generate a few of them

  particles mainid=particleid_;
  int nshots=1;
  G4double sign=1;

  G4double zposition = -1*cm;
  double energy_max=100;

  for(int i=0;i<nshots;i++){

	  if(particleid_==gamma)
		  particleid_=pionneutral;
	  else
		  particleid_=gamma;

	  // particleid_=pioncharged;
	  setParticleID(particleid_);

	  energy_=1101;
	  while(energy_>energy_max){//somehow sometimes the random gen shoots >1??
		  G4double rand =  G4INCL::Random::shoot();
		  energy_=(energy_max-1)*rand+1;
	  }
	  //G4cout << "shooting particle at " ;
	  double xpos=10;
	  while(fabs(xpos)>5){
		  xpos=10*G4INCL::Random::shoot() -5;
		  //G4cout << xpos <<  G4endl;
	  }
	  double ypos=10;
	  while(fabs(ypos)>5){
		  ypos=10*G4INCL::Random::shoot() -5;
	  }

	  G4ThreeVector position(xpos*cm, ypos*cm, zposition);
	  xorig_=xpos;
	  yorig_=ypos;

	  //G4cout << position <<  G4endl;

	  if(i==0){
		  mainid=particleid_;
	  }
	  else{
		  G4double ringsize=10*cm;
		  //make a ring
		  auto xycoords=G4INCL::Random::correlatedUniform(-1);
		  G4double magnitude=xycoords.first*xycoords.first+xycoords.second*xycoords.second;
		  magnitude=sqrt(magnitude);
		  xycoords.first*=ringsize/magnitude*sign;
		  xycoords.second*=ringsize/magnitude*sign;
		  position=G4ThreeVector(xycoords.first, xycoords.second, zposition);
		  //G4cout << "adding particle at " << position << " " << id <<  G4endl;
		  sign*=-1.;
	  }

	  fParticleGun->SetParticleEnergy(energy_ * GeV);
	  fParticleGun->SetParticlePosition(position);
	  fParticleGun->GeneratePrimaryVertex(anEvent);

  }
  particleid_=mainid;

}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

