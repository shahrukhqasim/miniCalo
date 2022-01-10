//
// Created by Shah Rukh Qasim on 1/6/22.
//

#include "PrimariesGenerator.hh"

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
#include <G4INCLRanecu.hh>
#include<ctime>
#include<sys/types.h>
#include <cstdlib>


PrimariesGenerator::PrimariesGenerator() :
        pythia_("/Users/shahrukhqasim/Workspace/NextCal/miniCalo/pythia8-data"),
        jetDef_(new fastjet::JetDefinition(fastjet::antikt_algorithm, 0.4)) {

    G4INCL::Random::setGenerator(new G4INCL::Ranecu());

    pythia_.readString("Beams:eCM = 14000.");
    pythia_.readString("Init:showChangedParticleData = off");
    pythia_.readString("Init:showChangedSettings = on");
    pythia_.readString("Next:numberShowLHA = 0");
    pythia_.readString("Next:numberShowInfo = 0");
    pythia_.readString("Next:numberShowProcess = 0");
    pythia_.readString("Next:numberShowEvent = 0");

    pythia_.readString("Tune:pp 14");
    pythia_.readString("Tune:ee 7");
    pythia_.readString("MultipartonInteractions:ecmPow=0.03344");
    pythia_.readString("MultipartonInteractions:bProfile=2");
    pythia_.readString("MultipartonInteractions:pT0Ref=1.41");
    pythia_.readString("MultipartonInteractions:coreRadius=0.7634");
    pythia_.readString("MultipartonInteractions:coreFraction=0.63");
    pythia_.readString("ColourReconnection:range=5.176");
    pythia_.readString("SigmaTotal:zeroAXB=off");
    pythia_.readString("SpaceShower:alphaSorder=2");
    pythia_.readString("SpaceShower:alphaSvalue=0.118");
    pythia_.readString("SigmaProcess:alphaSvalue=0.118");
    pythia_.readString("SigmaProcess:alphaSorder=2");
    pythia_.readString("MultipartonInteractions:alphaSvalue=0.118");
    pythia_.readString("MultipartonInteractions:alphaSorder=2");
    pythia_.readString("TimeShower:alphaSorder=2");
    pythia_.readString("TimeShower:alphaSvalue=0.118");
    pythia_.readString("SigmaTotal:mode = 0");
    pythia_.readString("SigmaTotal:sigmaEl = 21.89");
    pythia_.readString("SigmaTotal:sigmaTot = 100.309");
    //pythia_.readString("PDF:pSet=LHAPDF6:NNPDF31_nnlo_as_0118");

    //  pythia_.readString("SoftQCD:nonDiffractive = on");
    //  pythia_.readString("SoftQCD:singleDiffractive = on");
    //  pythia_.readString("SoftQCD:doubleDiffractive = on");
    pythia_.readString("SoftQCD:inelastic = on");

    pythia_.readString("Tune:preferLHAPDF = 2");
    pythia_.readString("Main:timesAllowErrors = 10000");
    pythia_.readString("Check:epTolErr = 0.01");
    pythia_.readString("Beams:setProductionScalesFromLHEF = off");
    pythia_.readString("SLHA:minMassSM = 1000.");
    pythia_.readString("ParticleDecays:limitTau0 = on");
    pythia_.readString("ParticleDecays:tau0Max = 10");
    pythia_.readString("ParticleDecays:allowPhotonRadiation = on");

    pythia_.init();
    std::cout << "Pythia initialized..." << std::endl;

    pythia_.rndm.init(0);

    fParticleGun = new G4ParticleGun();

}


void PrimariesGenerator::generate() {
}


void PrimariesGenerator::GenerateSingleVertex(G4PrimaryVertex *vertex) {
    double const etaTargetMin = 1.5;
    double const etaTargetMax = 3.0;
    double const eMin = 10.;
    double const eMax = 4000.;

    // jets from this position at eta ~ 3.6 will hit the center of the detector
    xorig_ = 0;
    yorig_ = 0;


    // // make a dummy primary proton
    // G4PrimaryParticle* primaryParticle = new G4PrimaryParticle(2212);

    G4double zsign = 1.;
    std::vector<int> primaries;
    primaries.reserve(1024 * 16);

    G4ThreeVector vertex_position(xorig_, yorig_, 0);
    G4double vertex_time(0.);


    while (true) {
        // loop until an event passing a gen filter is generated

        unsigned iTry(0);
        while (true) {
            // loop until pythia generates a consistent event (should succeed)
            if (!pythia_.next()) {
                if (++iTry < 10)
                    continue;
//                pythia_.stat();
            }
            break;
        }

        if (iTry == 10) {
            // failed event generation - make a dummy vertex (or just crash?)
            G4cerr << "EVENT GENERATION FAILED!!!!" << G4endl;


            std::cout << "Error in event gen PrimariesGenerator.cc" << std::endl;
            exit(1);

        }


        primaries.clear();
//            fjinputs_.clear();
        double totalen = 0;

        // cluster ak4 jets from final state particles

        for (int i{0}; i < pythia_.event.size(); ++i) {
            auto &part(pythia_.event[i]);

            if (part.isFinal()) {
                if (part.eta() < 4 && part.eta() > 1) {
                        std::cout<<"B eta "<<part.eta()<<" and e "<<part.e()<<std::endl;

                    primaries.push_back(i);
                    totalen += part.e();
                }
            }
        }
        std::cout << "N nfilt E," << pythia_.event.size() << "," << primaries.size() << "," << totalen << std::endl;

        break;
    }


    for (int ipart: primaries) {
        auto &pj(pythia_.event[ipart]);

        int pdgId = pj.id();
        auto *partDefinition = G4ParticleTable::GetParticleTable()->FindParticle(pdgId);
        if (partDefinition == nullptr)
            continue; //throw std::runtime_error(std::string("Unknown particle ") + std::to_string(pdgId));

        G4PrimaryParticle *particle = new G4PrimaryParticle(pdgId);
        particle->SetMass(pj.m() * GeV);
        particle->SetMomentum(pj.px() * GeV, pj.py() * GeV, pj.pz() * GeV * zsign);
        particle->SetCharge(partDefinition->GetPDGCharge());
        vertex->SetPrimary(particle);


//        std::cout<<"Shooting particle " << particle->GetPx() << " " << particle->GetPy()
//        << " " << particle->GetPz() << " "<< particle->GetKineticEnergy() <<std::endl;

        //      primaryParticle->SetDaughter(particle);
    }
}

bool PrimariesGenerator::isJetGenerator() {
    return false;
}

std::vector<G4String> PrimariesGenerator::generateAvailableParticles() const {
    return {"isMinbias"};
}

B4PartGeneratorBase::particles PrimariesGenerator::getParticle() const {
    return minbias;
}

int PrimariesGenerator::isParticle(int i) const {
    return 0;
}

void PrimariesGenerator::GenerateParticle(G4Event *anEvent) {
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
    G4ThreeVector vertex_position(xorig_, yorig_, 0);
    G4double vertex_time(0.);
    G4PrimaryVertex *vertex = new G4PrimaryVertex(vertex_position, vertex_time);

    G4ThreeVector position(particle_position_x * m, particle_position_y * m, particle_position_z * m);
    G4ThreeVector direction(particle_direction_x * m, particle_direction_y * m, particle_direction_z * m);

    auto *partDefinition = G4ParticleTable::GetParticleTable()->FindParticle(particle_pdgid);
    if (partDefinition == nullptr)
        throw std::runtime_error(std::string("Unknown particle ") + std::to_string(particle_pdgid));


    fParticleGun->SetParticleMomentumDirection(direction);
    fParticleGun->SetParticleEnergy(particle_energy * GeV);
    fParticleGun->SetParticlePosition(position);
    fParticleGun->GeneratePrimaryVertex(anEvent);
}
void PrimariesGenerator::GeneratePU(G4Event *anEvent) {
    G4ThreeVector vertex_position(xorig_, yorig_, 0);
    G4double vertex_time(0.);

    G4PrimaryVertex *vertex = new G4PrimaryVertex(vertex_position, vertex_time);

    GenerateSingleVertex(vertex);

    anEvent->AddPrimaryVertex(vertex);
}
void PrimariesGenerator::GeneratePrimaries(G4Event *anEvent) {
    if (generatePU) {
        GeneratePU(anEvent);
    } else {
        GenerateParticle(anEvent);
    }
}


void PrimariesGenerator::SetNextToPU() {
    generatePU=true;
}

void PrimariesGenerator::SetNextToParticle(double position_x, double position_y, double position_z, double direction_x,
                                           double direction_y,
                                           double direction_z, int pdgid, double energy) {
    generatePU=false;
    particle_position_x = position_x;
    particle_position_y = position_y;
    particle_position_z = position_z;

    particle_direction_x = direction_x;
    particle_direction_y = direction_y;
    particle_direction_z = direction_z;

    particle_pdgid = pdgid;
    particle_energy = energy;
}