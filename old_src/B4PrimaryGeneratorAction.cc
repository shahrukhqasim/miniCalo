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

template<class T>
static G4String createString(const T &i) {
    std::stringstream ss;
    ss << i;
    std::string number = ss.str();
    return number;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


B4PrimaryGeneratorAction::B4PrimaryGeneratorAction(std::string jsonfile)
        : B4PartGeneratorBase(),
          fParticleGun(nullptr) {
    G4int nofParticles = 1;
    fParticleGun = new G4ParticleGun(nofParticles);

    // default particle kinematic
    //
    auto particleDefinition
            = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
    // G4ParticleTable::GetParticleTable()->DumpTable();
    fParticleGun->SetParticleDefinition(particleDefinition);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
    fParticleGun->SetParticleEnergy(100. * GeV);


    G4INCL::Random::setGenerator(new G4INCL::Geant4RandomGenerator());

    for (int i = 0; i < seedsoffset_; i++)
        G4double rand = G4INCL::Random::shoot();

    xorig_ = 0;
    yorig_ = 0;
    setParticleID(gamma);

    from_beamspot_ = false;

#ifdef FROMBEAMSPOT
    from_beamspot_=true;
#endif
    std::ifstream ifs;
    ifs.open(jsonfile);

    Json::CharReaderBuilder builder;
    builder["collectComments"] = true;
    JSONCPP_STRING errs;
    if (!parseFromStream(builder, ifs, &m_particles, &errs)) {
        std::cout << errs << std::endl;
        std::cout<<"Error in reading particles json file "<<jsonfile<<std::endl;
        exit(0);
    }



}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4PrimaryGeneratorAction::~B4PrimaryGeneratorAction() {
    delete fParticleGun;
}

std::vector<G4String> B4PrimaryGeneratorAction::generateAvailableParticles() const {
    std::vector<G4String> out;
    for (int i = 0; i < particles_size; i++)
        out.push_back(getParticleName((B4PrimaryGeneratorAction::particles) i));

    return out;
}


G4String B4PrimaryGeneratorAction::getParticleName(enum particles p) const {
    if (p == elec) {
        return "isElectron";
    }
    if (p == positron) {
        return "isPositron";
    }
    if (p == muon) {
        return "isMuon";
    }
    if (p == pioncharged) {
        return "isPionCharged";
    }
    if (p == pionneutral) {
        return "isPionNeutral";
    }
    if (p == klong) {
        return "isK0Long";
    }
    if (p == kshort) {
        return "isK0Short";
    }
    if (p == gamma) {
        return "isGamma";
    }

    return "isInvalid" + createString((int) p);
}

G4String B4PrimaryGeneratorAction::setParticleID(enum particles p) {
    particleid_ = p;
    G4ParticleDefinition *particleDefinition = 0;

    if (p == elec) {
        particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("e-");
        fParticleGun->SetParticleDefinition(particleDefinition);
    }
    if (p == positron) {
        particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("e+");
        fParticleGun->SetParticleDefinition(particleDefinition);
    }
    if (p == muon) {
        particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("mu-");
        fParticleGun->SetParticleDefinition(particleDefinition);
    }
    if (p == pioncharged) {
        particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("pi+");
        fParticleGun->SetParticleDefinition(particleDefinition);
    }
    if (p == pionneutral) {
        particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("pi0");
        fParticleGun->SetParticleDefinition(particleDefinition);
    }
    if (p == klong) {
        particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("kaon0L");
        fParticleGun->SetParticleDefinition(particleDefinition);
    }
    if (p == kshort) {
        particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("kaon0S");
        fParticleGun->SetParticleDefinition(particleDefinition);
    }
    if (p == gamma) {
        particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
        fParticleGun->SetParticleDefinition(particleDefinition);
    }

    return getParticleName(p);
    return "isInvalid" + createString((int) p);
}

G4String B4PrimaryGeneratorAction::setParticleID(std::string p) {
//    particleid_ = p;
    G4ParticleDefinition *particleDefinition = 0;

    if (p == "electron") {
        particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("e-");
        fParticleGun->SetParticleDefinition(particleDefinition);
    }
    if (p == "positron") {
        particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("e+");
        fParticleGun->SetParticleDefinition(particleDefinition);
    }
    if (p == "muon") {
        particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("mu-");
        fParticleGun->SetParticleDefinition(particleDefinition);
    }
    if (p == "pion_charged") {
        particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("pi+");
        fParticleGun->SetParticleDefinition(particleDefinition);
    }
    if (p == "pion_neutral") {
        particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("pi0");
        fParticleGun->SetParticleDefinition(particleDefinition);
    }
    if (p == "klong") {
        particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("kaon0L");
        fParticleGun->SetParticleDefinition(particleDefinition);
    }
    if (p == "kshort") {
        particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("kaon0S");
        fParticleGun->SetParticleDefinition(particleDefinition);
    }
    if (p == "gamma") {
        particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
        fParticleGun->SetParticleDefinition(particleDefinition);
    }

    return p;
//    return "isInvalid" + createString((int) p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4ThreeVector generatePosition(G4double dRmin, G4double dRmax, G4double z) {

    G4double x, y;
    while (true) {
        x = dRmax - 2 * dRmax * G4INCL::Random::shoot();
        y = dRmax - 2 * dRmax * G4INCL::Random::shoot();
        G4ThreeVector dir(x, y, 0);
        G4double dr = dir.r();
        if (dr > dRmin && dr <= dRmax)
            break;
    }
    return G4ThreeVector(x, y, z);
}

double gen_etaToR(const G4double &eta, const G4double &z) {
    return z * exp(-eta);
}

bool particleTraversesFully(G4double etamin, G4double etamax, const G4ThreeVector &position,
                            const G4ThreeVector &direction) {

    G4double calo_z = 30 * cm + 320 * cm;//should hit all layers
    G4double dz = calo_z - position.z();
    G4ThreeVector dir = direction;
    dir = dir.unit();
    dir *= dz / dir.z();

    G4ThreeVector newpos = position + dir;

    G4double eta = newpos.eta();
    if (eta <= etamax && eta >= etamin)
        return true;
    return false;
}

//just make sure it hits the first calo layer
G4ThreeVector generateDirection(G4double etamin, G4double etamax, const G4ThreeVector &position) {

    while (true) {
        G4double x_dir = 1. - 2. * G4INCL::Random::shoot();
        G4double y_dir = 1. - 2. * G4INCL::Random::shoot();
        G4double z_dir = G4INCL::Random::shoot() + 0.01;

        G4ThreeVector dir(x_dir, y_dir, z_dir);
        if (particleTraversesFully(etamin, etamax, position, dir))
            return dir.unit();
    }
}

//is allowed to change position in x and y
G4ThreeVector generateDirectionAndPosition(G4double etamin, G4double etamax, G4ThreeVector &position,
                                           G4double angle, G4double pos_drmin, G4double pos_drmax, G4double pos_z) {

    position = generatePosition(pos_drmin, pos_drmax, pos_z);

    auto rotaxis = position.orthogonal().unit();
    auto poscp = position;

    G4ThreeVector direction;

    int try_counter = 0;

    while (true) {
        G4double rotation = 2. * M_PI * G4INCL::Random::shoot();
        //G4cout << "rotaxis pre" << rotaxis << G4endl;
        rotaxis = rotaxis.rotate(position, rotation).unit(); //defines the rotation for direction
        //G4cout << "rotaxis post" << rotaxis << G4endl;

        poscp = position;
        direction = poscp.rotate(rotaxis, angle).unit();

        if (particleTraversesFully(etamin, etamax, position, direction))
            return direction;

        if (try_counter > 20) {
            position = generatePosition(pos_drmin, pos_drmax, pos_z);
            poscp = position;
        }
        if (try_counter > 2000) {
            throw std::runtime_error("no direction with angle found");
        }
        try_counter++;
    }

}

void B4PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent) {


    // This function is called at the begining of event

    // In order to avoid dependence of PrimaryGeneratorAction
    // on DetectorConstruction class we get world volume
    // from G4LogicalVolumeStore
    //
    G4double worldZHalfLength = 0.;
    auto worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");

    // Check that the world volume has box shape
    G4Box *worldBox = nullptr;
    if (worldLV) {
        worldBox = dynamic_cast<G4Box *>(worldLV->GetSolid());
    }

    if (worldBox) {
        worldZHalfLength = worldBox->GetZHalfLength();
    } else {
        G4ExceptionDescription msg;
        msg << "World volume of box shape not found." << G4endl;
        msg << "Perhaps you have changed geometry." << G4endl;
        msg << "The gun will be place in the center.";
        G4Exception("B4PrimaryGeneratorAction::GeneratePrimaries()",
                    "MyCode0002", JustWarning, msg);
    }
    // Set gun position
    int id = anEvent->GetEventID();

    double position_x = m_particles[id]["position"][0].asDouble();
    double position_y = m_particles[id]["position"][1].asDouble();
    double position_z = m_particles[id]["position"][2].asDouble();


    double direction_x = m_particles[id]["direction"][0].asDouble();
    double direction_y = m_particles[id]["direction"][1].asDouble();
    double direction_z = m_particles[id]["direction"][2].asDouble();


    double energy = m_particles[id]["energy"].asDouble();
    energy_ = energy;

    std::string pid = m_particles[id]["id"].asString();

    std::string particleName = setParticleID(pid);
    G4cout << "shooting " << particleName << " with event id " << anEvent->GetEventID() << G4endl;

    G4ThreeVector position(position_x * m, position_y * m, position_z * m);
    G4ThreeVector direction(direction_x * m, direction_y * m, direction_z * m);

    fParticleGun->SetParticleMomentumDirection(direction);
    fParticleGun->SetParticleEnergy(energy * GeV);
    fParticleGun->SetParticlePosition(position);
    fParticleGun->GeneratePrimaryVertex(anEvent);


    G4cout << "Energy: " << energy_ << ", position: " << position << " direction " << direction << G4endl;


}

void B4PrimaryGeneratorAction::GeneratePrimariesx(G4Event *anEvent) {
    // This function is called at the begining of event

    // In order to avoid dependence of PrimaryGeneratorAction
    // on DetectorConstruction class we get world volume
    // from G4LogicalVolumeStore
    //
    G4double worldZHalfLength = 0.;
    auto worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");

    // Check that the world volume has box shape
    G4Box *worldBox = nullptr;
    if (worldLV) {
        worldBox = dynamic_cast<G4Box *>(worldLV->GetSolid());
    }

    if (worldBox) {
        worldZHalfLength = worldBox->GetZHalfLength();
    } else {
        G4ExceptionDescription msg;
        msg << "World volume of box shape not found." << G4endl;
        msg << "Perhaps you have changed geometry." << G4endl;
        msg << "The gun will be place in the center.";
        G4Exception("B4PrimaryGeneratorAction::GeneratePrimaries()",
                    "MyCode0002", JustWarning, msg);
    }
    // Set gun position


    double energy_max = 100;
    double energy_min = 1;
    //energy_=15;

    //iterate
//  if(particleid_==gamma)
//      G4cout << setParticleID(elec) <<G4endl;
//  else if(particleid_==elec)
//      G4cout << setParticleID(gamma) <<G4endl;
    // else if(particleid_==gamma)
    //    G4cout << setParticleID(pioncharged) <<G4endl;

    if (particleid_ + 1 < positron)
        setParticleID((particles)(particleid_ + 1));
    else
        setParticleID(elec);
    //positron
    G4cout << "shooting " << getParticleName(particleid_) << " with event id " << anEvent->GetEventID() << G4endl;


// good metrics are: angle w.r.t. projective and energy



    //particleid_=pioncharged;

    energy_ = 10001;
    while (energy_ > energy_max) {//somehow sometimes the random gen shoots >1??
        G4double rand = G4INCL::Random::shoot();
        energy_ = (energy_max - energy_min) * rand + energy_min;
    }


    G4cout << "get position " << G4endl;


    angle_ = M_PI / 60. * G4INCL::Random::shoot0(); //almost pointing, 3 degrees difference allowed

    G4ThreeVector position;//15.*cm,66*cm,xorig_,yorig_);
    G4ThreeVector direction = generateDirectionAndPosition(1.55, 2.95, position, angle_, 20. * cm, 60 * cm, 329 * cm);

    std::cout << "position " << position << ", direction " << direction << ", angle " << angle_ << std::endl;

    xorig_ = position.x();
    yorig_ = position.y();


    if (from_beamspot_) {
        xorig_ = 0;
        yorig_ = 0;
        position = G4ThreeVector(xorig_, yorig_, 0);
    }

    dirx_ = direction.x();
    diry_ = direction.y();
    dirz_ = direction.z();


    G4cout << "projective direction " << position.unit() << G4endl;
    G4cout << "particle direction " << direction.unit() << G4endl;


    diff_proj_phi_ = direction.deltaPhi(position);
    diff_proj_theta_ = direction.theta(position);
    //same as theta for some reason..
    angle_ = direction.angle(position.unit());

    //G4cout << "Dphi, DTheta " <<  diff_proj_phi_ << ", " << diff_proj_theta_ << " angle " << angle_ << G4endl;

    if (from_beamspot_)
        G4cout << "Direction eta,phi " << direction.eta() << ", " << direction.phi() << G4endl;

    G4cout << "new direction " << direction << G4endl;

    //get eta

//0.2*G4INCL::Random::shoot(),0.2*G4INCL::Random::shoot(),0);


//  fParticleGun->setParticleId(xyz);
    fParticleGun->SetParticleMomentumDirection(direction);
    fParticleGun->SetParticleEnergy(energy_ * GeV);
    fParticleGun->SetParticlePosition(position);
    fParticleGun->GeneratePrimaryVertex(anEvent);


    G4cout << "Energy: " << energy_ << ", position: " << position << " direction " << direction << G4endl;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

