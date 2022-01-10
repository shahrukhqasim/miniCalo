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
// $Id: B4aActionInitialization.hh 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file B4aActionInitialization.hh
/// \brief Definition of the B4aActionInitialization class

#ifndef B4aActionInitialization_h
#define B4aActionInitialization_h 1

#include "defines.h"
#include "G4VUserActionInitialization.hh"
#include "G4String.hh"
#include "B4aEventAction.hh"

class B4DetectorConstruction;
class B4DetectorConstructionBase;

/// Action initialization class.
///

class B4aActionInitialization : public G4VUserActionInitialization
{
  public:
    B4aActionInitialization(B4DetectorConstructionBase *detConstruction,
                            B4aEventAction *event_action,
                            B4PartGeneratorBase *primaries_generator);
    virtual ~B4aActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;

    void setFilename(G4String fname){
    	fname_=fname;
    }

  private:
    B4DetectorConstructionBase* fDetConstruction;
    G4String fname_;
    B4aEventAction* eventAction;
    B4PartGeneratorBase *primaries_generator;
};

#endif

    
