////
//// ********************************************************************
//// * License and Disclaimer                                           *
//// *                                                                  *
//// * The  Geant4 software  is  copyright of the Copyright Holders  of *
//// * the Geant4 Collaboration.  It is provided  under  the terms  and *
//// * conditions of the Geant4 Software License,  included in the file *
//// * LICENSE and available at  http://cern.ch/geant4/license .  These *
//// * include a list of copyright holders.                             *
//// *                                                                  *
//// * Neither the authors of this software system, nor their employing *
//// * institutes,nor the agencies providing financial support for this *
//// * work  make  any representation or  warranty, express or implied, *
//// * regarding  this  software system or assume any liability for its *
//// * use.  Please see the license in the file  LICENSE  and URL above *
//// * for the full disclaimer and the limitation of liability.         *
//// *                                                                  *
//// * This  code  implementation is the result of  the  scientific and *
//// * technical work of the GEANT4 collaboration.                      *
//// * By using,  copying,  modifying or  distributing the software (or *
//// * any work based  on the software)  you  agree  to acknowledge its *
//// * use  in  resulting  scientific  publications,  and indicate your *
//// * acceptance of all terms of the Geant4 Software license.          *
//// ********************************************************************
////
//// $Id: B4DetectorConstruction.cc 101905 2016-12-07 11:34:39Z gunter $
////
///// \file B4DetectorConstruction.cc
///// \brief Implementation of the B4DetectorConstruction class
//
//#include "B4DetectorConstructionFP.hh"
//
//#include "G4Material.hh"
//#include "G4NistManager.hh"
//
//#include "G4RotationMatrix.hh"
//#include "G4Box.hh"
//#include "G4LogicalVolume.hh"
//#include "G4PVPlacement.hh"
//#include "G4PVReplica.hh"
//#include "G4GlobalMagFieldMessenger.hh"
//#include "G4AutoDelete.hh"
//#include "G4UserLimits.hh"
//#include "G4PVDivision.hh"
//#include "B4DetectorConstruction.hh"
//
//#include "G4GeometryManager.hh"
//#include "G4PhysicalVolumeStore.hh"
//#include "G4LogicalVolumeStore.hh"
//#include "G4SolidStore.hh"
//
//#include "G4VisAttributes.hh"
//#include "G4Colour.hh"
//
//#include "G4PhysicalConstants.hh"
//#include "G4SystemOfUnits.hh"
//#include "G4FieldManager.hh"
//#include "G4TransportationManager.hh"
//#include "G4PropagatorInField.hh"
//#include "json/json.h"
//
//#include "G4Cons.hh"
//#include <cmath>
//
//#include "sensorContainer.h"
//
//#include <cstdlib>
//#include <Geant4/G4PVParameterised.hh>
////#include "Math/Vector3D.h"
//
////#define USEDIVISION
//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//G4ThreadLocal
//G4GlobalMagFieldMessenger* B4DetectorConstructionFP::fMagFieldMessenger = nullptr;
//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//template<class T>
//static G4String createString(const T& i){
//	std::stringstream ss;
//	ss << i;
//	std::string number=ss.str();
//	return number;
//}
//
//
//B4DetectorConstructionFP::B4DetectorConstructionFP(Json::Value& value)
//: B4DetectorConstructionBase(),
//  fCheckOverlaps(true),
//  m_vacuum(0),
//  m_pb(0),
//  m_pbtungsten(0),
//  m_silicon(0),
//  m_cu(0),
//  m_detector_specs(value)
//{
//
//
//    limit_in_calo_time_max_=500*ns;//this could be more low energy stuff
////    limit_in_calo_energy_max_=500*keV;
//	limit_in_calo_energy_max_=50*MeV;
//
//    limit_world_time_max_=500*ns; //either got there or not (30ns should be easily sufficient
//    limit_world_energy_max_=1000*eV;
//
//}
//
//G4VPhysicalVolume* B4DetectorConstructionFP::Construct()
//{
//	//DefineGeometry(homogenous_ecal_only);
//	//DefineGeometry(hcal_only_irregular);
//	DefineGeometry(fullendcap);
//	// Define materials
//	DefineMaterials();
//
//	// Define volumes
//	return DefineVolumes();
//}
//
//void  B4DetectorConstructionFP::DefineGeometry(geometry geo){
//
//    layerAbsorberFractions.clear();
//    layerThicknesses.clear();
//    layerGranularity.clear();
//	if(geo == fullendcap){
//		calorThickness=2200*mm;
//
//        nofEELayers = 14;
//        Ncalowedges=96;
//        nofHB=18;
//        int etasegments=24;
//
//        Calo_start_eta=1.5;
//        Calo_end_eta=3.0;
//        m_calo_start_z=320*cm;
//
//        for(int i=0;i<nofEELayers+nofHB;i++){
//            if(i<nofEELayers){
//                layerThicknesses.push_back(10.382*mm + 0.300*mm);////absorber plus silicon, this makes 1.85 X0 per layer
//                layerAbsorberFractions.push_back(0.9719153717);
//            }
//            else{//copper absorber: 15.32 cm == 1 l_0 , ECal part already has about 1. Add 9 more// 1/2 per layer
//                layerThicknesses.push_back(15.32*cm / 2. + 0.300*mm);
//                layerAbsorberFractions.push_back(0.9960988296);
//            }
//
//		    layerGranularity.push_back(etasegments);
//		}
//
//		layerThicknessEE=26*cm / (float)nofEELayers;
//		layerThicknessHB=(calorThickness-nofEELayers*layerThicknessEE)/(float)nofHB; //100*mm;
//
//	}
//	else{
//	    G4ExceptionDescription msg;
//	            msg << "Geometry not supported in this branch";
//	            G4Exception("B4DetectorConstructionFP::DefineGeometry()",
//	                    "MyCode0001", FatalException, msg);
//	}
//}
//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//B4DetectorConstructionFP::~B4DetectorConstructionFP()
//{
//}
//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//
//
//
//G4Material* B4DetectorConstructionFP::getMaterial(std::string name) {
//	if (name ==std::string("G4_Cu")) {
//		return m_cu;
//	}
//	else if (name ==std::string("G4_Si")) {
//		return m_silicon;
//	}
//	if (name ==std::string("G4_PbWO4")) {
//		return m_pbtungsten;
//	}
//	if (name ==std::string("Galactic")) {
//		return m_vacuum;
//	}
//	if (name ==std::string("Air")) {
//		return m_air;
//
//	}
//	if (name ==std::string("StainlessSteel")) {
//		return m_stainlesssteel;
//
//	}
//	if (name ==std::string("G4_Pb")) {
//		return m_pb;
//
//	}
//
//	return nullptr;
//}
//
//void B4DetectorConstructionFP::createCalo(G4LogicalVolume * worldLV,
//										G4ThreeVector position,G4String name){
//
//
//
//
//	m_worldLV = worldLV;
//	createDataVectors();
//
//
//	G4double  calo_start_z = m_detector_specs["calo_start_z"].asDouble() * m;
//	G4double  calo_end_z = m_detector_specs["calo_end_z"].asDouble() * m;
//	G4double midz = (calo_end_z + calo_start_z) /2;
//	G4double dz = ((calo_end_z - calo_start_z)/2);
//	G4double calo_start_eta = m_detector_specs["calo_start_eta"].asDouble();
//	G4double calo_end_eta = m_detector_specs["calo_end_eta"].asDouble();
//
//	G4ThreeVector positionCalo = G4ThreeVector(0,0,midz);
//
//	G4double pRmin1 = etaToR(calo_end_eta, calo_start_z);
//	G4double pRmax1 = etaToR(calo_start_eta, calo_start_z);
//	G4double pRmin2 = etaToR(calo_end_eta, calo_end_z);
//	G4double pRmax2 = etaToR(calo_start_eta, calo_end_z);
//
//
//	std::cout<<"Calo start z "<<calo_start_z<<std::endl;
//	std::cout<<"Calo end z "<<calo_end_z<<std::endl;
//	std::cout<<"Calo start eta "<<calo_start_eta<<std::endl;
//	std::cout<<"Calo end eta "<<calo_end_eta<<std::endl;
//	std::cout<<"R min 1 "<<pRmin1<<std::endl;
//	std::cout<<"R max 1 "<<pRmax1<<std::endl;
//	std::cout<<"R min 2 "<<pRmin2<<std::endl;
//	std::cout<<"R max 2 "<<pRmax2<<std::endl;
//
////	G4Cons* caloS
////			= new G4Cons("Calorimeter",pRmin1, pRmax1,
////			pRmin2, pRmax2,
////			dz,
////			0, 2*pi);
////	G4LogicalVolume* caloLV
////			= new G4LogicalVolume(caloS, m_vacuum, "Calorimeter",0,0,0);
////	new G4PVPlacement(0,               // no rotation
////					  positionCalo, // at (x,y,z)
////					  caloLV,       // its logical volume
////					  "Calorimeter",       // its name
////					  worldLV,         // its mother  volume
////					  false,           // no boolean operations
////					  0,               // copy number
////					  fCheckOverlaps); // checking overlaps
//
//
////	for (int i = 0; i< m_parameterized_data.m_rmin1.size();i++) {
////		printParameterizedData(i);
//////	exit(0);
////		G4Cons *chamberS
////				= new G4Cons("x", m_parameterized_data.m_rmin1[i], m_parameterized_data.m_rmax1[i],
////							 m_parameterized_data.m_rmin2[i], m_parameterized_data.m_rmax2[i],
////							 m_parameterized_data.m_dz[i], m_parameterized_data.m_sphi[i], m_parameterized_data.m_dphi[i]);
////		auto fLogicChamber
////				= new G4LogicalVolume(chamberS, m_cu, "y", 0, 0, 0);
////
////		G4ThreeVector positionChamber = G4ThreeVector(0, 0, m_parameterized_data.m_midz[i]);
////		new G4PVPlacement(0,               // no rotation
////						  positionChamber, // at (x,y,z)
////						  fLogicChamber,       // its logical volume
////						  "Calorimeter",       // its name
////						  worldLV,         // its mother  volume
////						  false,           // no boolean operations
////						  0,               // copy number
////						  fCheckOverlaps); // checking overlaps
////
////		G4VisAttributes *boxVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
////		worldLV->SetVisAttributes(boxVisAtt);
//////	caloLV ->SetVisAttributes(boxVisAtt);
////		fLogicChamber->SetVisAttributes(boxVisAtt);
////	}
//
//
//	G4VPVParameterisation* chamberParam =
//			new B4FParameterisation(&m_parameterized_data);
//
//	// dummy value : kZAxis -- modified by parameterised volume
//
//	int i=0;
//	G4Cons *chamberS
//			= new G4Cons("x", m_parameterized_data.m_rmin1[i], m_parameterized_data.m_rmax1[i],
//						 m_parameterized_data.m_rmin2[i], m_parameterized_data.m_rmax2[i],
//						 m_parameterized_data.m_dz[i], m_parameterized_data.m_sphi[i], m_parameterized_data.m_dphi[i]);
//	auto fLogicChamber
//			= new G4LogicalVolume(chamberS, m_cu, "y", 0, 0, 0);
//
//	std::cout<<"m_parameterized_data.m_rmin1.size() "<<m_parameterized_data.m_rmin1.size()<<std::endl;
//
//	new G4PVParameterised("z",       // their name
//						  fLogicChamber,   // their logical volume
//						  worldLV,       // Mother logical volume
//						  kZAxis,          // Are placed along this axis
//						  m_parameterized_data.m_rmin1.size(),    // Number of chambers
//						  chamberParam,    // The parametrisation
//						  fCheckOverlaps); // checking overlaps
//
//	int layernum = 0;
//	auto layers_upper = m_detector_specs["layers"];
//	for (auto &layer_upper: layers_upper) {
//		for (auto &layer: layer_upper) {
//			layernum += 1;
//		}
//	}
//}
//
//void B4DetectorConstructionFP::DefineMaterials()
//{
//	// Lead material defined using NIST Manager
//	auto nistManager = G4NistManager::Instance();
//	nistManager->FindOrBuildMaterial("G4_Pb");
//    nistManager->FindOrBuildMaterial("G4_Cu");
//	nistManager->FindOrBuildMaterial("G4_PbWO4");
//    nistManager->FindOrBuildMaterial("G4_Si");
//   // nistManager->FindOrBuildMaterial("G4_Air");
//
//	//nistManager->ListMaterials("all");
//
//	// Liquid argon material
//	G4double a;  // mass of a mole;
//	G4double z;  // z=mean number of protons;
//	G4double density;
//	new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
//	// The argon by NIST Manager is a gas with a different density
//
//	// Vacuum
//	new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
//			kStateGas, 2.73*kelvin, 3.e-18*pascal);
//
//	// Print materials
//	//G4cout << *(G4Material::GetMaterialTable()) << G4endl;
//
//
//	// Get materials
//	m_vacuum = G4Material::GetMaterial("Galactic"); // "Galactic"
//	m_pb = G4Material::GetMaterial("G4_Pb");
//	m_pbtungsten = G4Material::GetMaterial("G4_PbWO4");
//	m_silicon = G4Material::GetMaterial("G4_Si");
//	m_cu = G4Material::GetMaterial("G4_Cu");
////	m_cuw = G4Material::GetMaterial("G4_CuW");
//
//
////	m_stainlesssteel = G4Material::GetMaterial("G4_STAINLESS-STEEL");
//
//	G4Element* Mn  = new G4Element("Manganese", "Mn", z=25, a=    54.94*g/mole);
//	G4Element* Cr  = new G4Element("Chromium",  "Cr", z=24, a=    52.00*g/mole);
//	G4Element* Ni  = new G4Element("Nickel",    "Ni", z=28, a=    58.70*g/mole);
//
//	a=16.*g/mole;
//	auto elO=new G4Element("Oxygen","O2",8.,a);
//	a=14.01*g/mole;
//	auto elN=new G4Element("Nitrogen","N2",7.,a);
//
//	density = 1.290*mg/cm3;
//
//
//	m_air = new G4Material("Air",density, 2);
//	m_air->AddElement(elN, 0.7);
//	m_air->AddElement(elO, 0.3);
//
//	m_stainlesssteel = new G4Material("StainlessSteel",8.02*g/cm3, 5);
//	m_stainlesssteel->AddElement(Mn, 0.02);
//	m_stainlesssteel->AddMaterial(nistManager->FindOrBuildMaterial("G4_Si"), 0.01);
//	m_stainlesssteel->AddElement(Cr, 0.19);
//	m_stainlesssteel->AddElement(Ni, 0.10);
//	m_stainlesssteel->AddMaterial(nistManager->FindOrBuildMaterial("G4_Fe"), 0.68);
//
//
//	if ( ! m_vacuum || ! m_pb || ! m_pbtungsten || !m_silicon ) {
//		G4ExceptionDescription msg;
//		msg << "Cannot retrieve materials already defined.";
//		G4Exception("B4DetectorConstructionFP::DefineVolumes()",
//				"MyCode0001", FatalException, msg);
//	}
//
//}
//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//G4VPhysicalVolume* B4DetectorConstructionFP::DefineVolumes()
//{
//	// Geometry parameters
//	//calorSizeXY  = 200*cm;
//
//
//
//
//
//
//    //auto calorThickness = nofEELayers * layerThicknessEE + nofHB*layerThicknessHB;
//    G4double worldSizeXY = 4.9 * m;
//	G4double worldSizeZ  = 15. * m;
//
//
//
//	//
//	// World
//	//
//	auto worldS
//	= new G4Box("World",           // its name
//			worldSizeXY/2, worldSizeXY/2, worldSizeZ); // its size
//
//	auto worldLV
//	= new G4LogicalVolume(
//			worldS,           // its solid
//			m_vacuum,  // its material
//			"World");         // its name
//
//	worldLV->SetUserLimits(new G4UserLimits(
//	        worldSizeZ/10., //max step length
//            worldSizeZ*50., //max track length
//            limit_world_time_max_, //max track time
//            limit_world_energy_max_)); //min track energy
//
//
//	auto worldPV
//	= new G4PVPlacement(
//			0,                // no rotation
//			G4ThreeVector(),  // at (0,0,0)
//			worldLV,          // its logical volume
//			"World",          // its name
//			0,                // its mother  volume
//			false,            // no boolean operation
//			0,                // copy number
//			fCheckOverlaps);  // checking overlaps
//
//	//
//	// Calorimeter
//	//
//	createCalo(worldLV,G4ThreeVector(0,0,m_calo_start_z),"");
//
//
//
//	G4cout << "created in total "<< activecells_.size()<<" sensors" <<G4endl;
//
//	//
//	// Visualization attributes
//	//
//
//
////	worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());
////
////	double g=0;
////	for(auto& v: activecells_){
////	    auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(.5,g,.0));
////	    simpleBoxVisAtt->SetVisibility(true);
////	    simpleBoxVisAtt->SetForceSolid(true);
////		v.getVol()->GetLogicalVolume()->SetVisAttributes(simpleBoxVisAtt);
////		g+= 1./(double)activecells_.size();
////	}
//	//
//	// Always return the physical World
//	//
//	return worldPV;
//}
//
//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//void B4DetectorConstructionFP::ConstructSDandField()
//{
//	// Create global magnetic field messenger.
//	// Uniform magnetic field is then created automatically if
//	// the field value is not zero.
//	G4ThreeVector fieldValue(0,0,0);// no need for field here 1.*tesla);
//	fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
//	fMagFieldMessenger->SetVerboseLevel(2);
//
//	// Register the field messenger for deleting
//	G4AutoDelete::Register(fMagFieldMessenger);
//
//
//	auto* fieldprop = G4TransportationManager::GetTransportationManager()->GetPropagatorInField();
//	fieldprop->SetMaxLoopCount(100) ;//check 100, 10 is bad, less is bad, default is 1000, maybe a bit less works, too
//
//
//}
//
//void B4DetectorConstructionFP::createDataVectors() {
//	int numCons = 0;
//
//	auto layers_upper = m_detector_specs["layers"];
//	for (auto &layer_upper: layers_upper) {
//		for (auto &layer: layer_upper) {
//
//			double start_z = layer["start_z"].asDouble();
//			double end_z = layer["end_z"].asDouble();
//			double z_length = end_z - start_z;
//			G4Material* material = getMaterial(layer["material"].asString());
//			for (auto &eta_segment : layer["eta_segments"]) {
//				int numPhi = eta_segment["phi_segments"].asInt();
//				numCons += numPhi;
//				double start_phi = 0;
//				double delta_phi = 2 * pi / (float)numPhi;
//				double start_eta = eta_segment["start_eta"].asDouble();
//				double end_eta = eta_segment["end_eta"].asDouble();
//				double eta_width = end_eta - start_eta;
//
//
//				for (int p=0; p < numPhi; p++) {
//					m_parameterized_data.m_midz.push_back(((end_z + start_z)/2)*m);
//					m_parameterized_data.m_rmin1.push_back(etaToR(start_eta+eta_width ,start_z)*m);
//					m_parameterized_data.m_rmax1.push_back(etaToR(start_eta ,start_z)*m);
//					m_parameterized_data.m_rmin2.push_back(etaToR(start_eta+eta_width ,start_z+z_length)*m);
//					m_parameterized_data.m_rmax2.push_back(etaToR(start_eta ,start_z+z_length)*m);
//					m_parameterized_data.m_dz.push_back((z_length/2.)*m);
//					m_parameterized_data.m_sphi.push_back(start_phi);
//					m_parameterized_data.m_dphi.push_back(delta_phi);
//					m_parameterized_data.m_material.push_back(material);
//					start_phi += delta_phi;
////					break;
//				}
//			}
////			break;
//		}
//	}
//}
//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//
//B4FParameterisation::B4FParameterisation(ParameterizedData *parameterized_data) {
//	m_parameterized_data = parameterized_data;
//}
//
//B4FParameterisation::~B4FParameterisation() {
//
//}
//
//void B4FParameterisation::ComputeTransformation(const G4int copyNo, G4VPhysicalVolume *physVol) const {
//
//	G4double z_position = (m_parameterized_data->m_midz[copyNo]);
//	G4ThreeVector origin(0,0,z_position);
//	physVol->SetTranslation(origin);
//	physVol->SetRotation(0);
//
//}
//
//void B4FParameterisation::ComputeDimensions(G4Cons &box, const G4int copyNo,
//											const G4VPhysicalVolume *physVol) const {
//
//	box.SetInnerRadiusMinusZ (m_parameterized_data->m_rmin1[copyNo]);
//	box.SetOuterRadiusMinusZ (m_parameterized_data->m_rmax1[copyNo]);
//	box.SetInnerRadiusPlusZ  (m_parameterized_data->m_rmin2[copyNo]);
//	box.SetOuterRadiusPlusZ  (m_parameterized_data->m_rmax2[copyNo]);
//	box.SetZHalfLength       (m_parameterized_data->m_dz[copyNo]);
////	box.SetStartPhiAngle(m_parameterized_data->m_sphi[copyNo]);
////	box.SetDeltaPhiAngle(m_parameterized_data->m_dphi[copyNo]);
//	box.SetStartPhiAngle(0);
//	box.SetDeltaPhiAngle(2*pi);
//
//}
//
////G4VSolid* B4FParameterisation::ComputeSolid(const G4int copyNo, G4VPhysicalVolume *physVol) {
////	return new G4Cons(
////			std::string("copy_") + std::to_string(copyNo),
////			m_parameterized_data->m_rmin1[copyNo],
////			m_parameterized_data->m_rmax1[copyNo],
////			m_parameterized_data->m_rmin2[copyNo],
////			m_parameterized_data->m_rmax2[copyNo],
////			m_parameterized_data->m_dz[copyNo],
////			m_parameterized_data->m_sphi[copyNo],
////			m_parameterized_data->m_dphi[copyNo]);
////}
//
//
//G4Material* B4FParameterisation::ComputeMaterial(const G4int copyNo, G4VPhysicalVolume *physVol,
//												 const G4VTouchable *parentTouch) {
//	return m_parameterized_data->m_material[copyNo];
//
//}
//
//
//void B4DetectorConstructionFP::printParameterizedData(const size_t &index) {
//	std::cout<<"At index "<<index<<std::endl;
//	std::cout<<"Mid Z: "<<m_parameterized_data.m_midz[index]<<std::endl;
//	std::cout<<"R min1: "<<m_parameterized_data.m_rmin1[index]<<std::endl;
//	std::cout<<"R max1: "<<m_parameterized_data.m_rmax1[index]<<std::endl;
//	std::cout<<"R min2: "<<m_parameterized_data.m_rmin2[index]<<std::endl;
//	std::cout<<"R max2: "<<m_parameterized_data.m_rmax2[index]<<std::endl;
//	std::cout<<"dz/2: "<<m_parameterized_data.m_dz[index]<<std::endl;
//	std::cout<<"Start phi: "<<m_parameterized_data.m_sphi[index]<<std::endl;
//	std::cout<<"End phi: "<<m_parameterized_data.m_dphi[index]<<std::endl;
////	std::cout<<"Material: "<<m_parameterized_data.m_material[index]<<std::endl;
//}