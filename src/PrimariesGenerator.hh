//
// Created by Shah Rukh Qasim on 1/6/22.
//

#ifndef B4A_PRIMARIESGENERATOR_HH
#define B4A_PRIMARIESGENERATOR_HH


#include "B4PrimaryGeneratorAction.hh"
#include "G4PrimaryVertex.hh"
#include "B4PartGeneratorBase.hh"
#include "globals.hh"
#include <vector>
#include "Pythia8/Pythia.h"

#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"



class PrimariesGenerator : public B4PartGeneratorBase {
public:
    PrimariesGenerator();
    void generate();
    void GenerateSingleVertex(G4PrimaryVertex* vertex);

    bool isJetGenerator() override;

    std::vector<G4String> generateAvailableParticles() const override;

    particles getParticle() const override;

    int isParticle(int i) const override;

    void GeneratePU(G4Event *anEvent);
    void GenerateParticle(G4Event *anEvent);
    void GeneratePrimaries(G4Event *anEvent) override;

    void SetNextToPU();
    void
    SetNextToParticle(double position_x, double position_y, double position_z, double direction_x, double direction_y,
                      double direction_z, int pdgid, double energy);


protected:
    Pythia8::Pythia pythia_;
    fastjet::JetDefinition* jetDef_;

    double xorig_ = 0;
    double yorig_ = 0;

    bool generatePU=true;

    double particle_position_x = 0;
    double particle_position_y = 0;
    double particle_position_z = 0;

    double particle_direction_x = 0;
    double particle_direction_y = 0;
    double particle_direction_z = 0;
    double particle_energy = 0;
    int particle_pdgid = 0;

    G4ParticleGun* fParticleGun;


};



#endif //B4A_PRIMARIESGENERATOR_HH
