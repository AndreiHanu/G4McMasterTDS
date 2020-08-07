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
// ********************************************************************
// G4SiDetector.cc
//
// Description: Definition of the Eljen Technologies M550-20x8-1 (P/N: 4467-2233)
// plastic scintillator detector (EJ-204) that is used by McMaster University to
// perform dosimetry measurements for the lens of eye.
//
// ********************************************************************

#include "DetectorConstruction.hh"
#include <cmath>

// Units and constants
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

// Manager classes
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4GeometryManager.hh"
#include "G4SDManager.hh"

// Store classes
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

// Geometry classes
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4VisAttributes.hh"

// Primitive geometry types
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Polycone.hh"

// Boolean operations on volumes
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"

// Regions
#include "G4Region.hh"

// Messenger classes
#include "G4GenericMessenger.hh"

// Scoring Components
#include "G4MultiFunctionalDetector.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSTrackLength.hh"
#include "G4PSPassageTrackLength.hh"
#include "G4PSSphereSurfaceCurrent.hh"
#include "G4PSSphereSurfaceCurrent3D.hh"
#include "G4PSSphereSurfaceFlux.hh"
#include "G4PSSphereSurfaceFlux3D.hh"
#include "G4PSIncidentKineticEnergy.hh"
#include "G4SDParticleFilter.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction(): G4VUserDetectorConstruction(), fCheckOverlaps(true),
WorldPhysical(0)
{	
	// Source Radius
	sourceRadius = 25.*cm;
			 
	// Define Materials
	DefineMaterials();
	
	// Define commands to control the geometry
   	DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
    // NIST Manager
	G4NistManager* nistManager = G4NistManager::Instance();
	nistManager->SetVerbose(0);

	G4double a;  // atomic mass
  	G4double z;  // atomic number
  	G4double density,ncomponents,fractionmass;
  
  	G4Element* La = new G4Element("Lanthanum", "La", z=57.,a=138.9055*g/mole);
  	G4Element* Br = new G4Element("Bromium", "Br", z=35., a=79.904*g/mole);
  	G4Element* Ce = new G4Element("Cerium", "Ce", z=58., a=140.116*g/mole);

	G4Material* LaBr3 = new G4Material("LaBr3", density = 5.29*g/cm3, ncomponents=2);
  	LaBr3->AddElement(La, fractionmass=25.*perCent);
  	LaBr3->AddElement(Br, fractionmass=75.*perCent);

	G4Material* LaBr3_Ce = new G4Material("LaBr3_Ce", density = 5.29*g/cm3, ncomponents=2);
	LaBr3_Ce->AddMaterial(LaBr3, fractionmass=99.5*perCent);
	LaBr3_Ce->AddElement(Ce, fractionmass=0.5*perCent);
	  
  	// Set the materials for the Geometry
	fMatWorld = nistManager->FindOrBuildMaterial("G4_Galactic");
	//fMatWorld = nistManager->FindOrBuildMaterial("G4_AIR");
	fMatDMB = nistManager->FindOrBuildMaterial("G4_NYLON-8062");

	// LaBr3 Detector Materials
	fMatLaBr3Housing = nistManager->FindOrBuildMaterial("G4_Al");
	fMatLaBr3Interior = nistManager->FindOrBuildMaterial("G4_AIR");
	fMatLaBr3Reflector = nistManager->FindOrBuildMaterial("G4_TEFLON");
	fMatLaBr3Crystal = LaBr3_Ce;
	fMatLaBr3LightGuide = nistManager->FindOrBuildMaterial("G4_PLEXIGLASS");
	fMatLaBr3PMT = nistManager->FindOrBuildMaterial("G4_Pyrex_Glass");
	fMatLaBr3PMTInterior = nistManager->FindOrBuildMaterial("G4_Galactic");

	// Plastic Detector Materials
	fMatPlasticHousing = nistManager->FindOrBuildMaterial("G4_Al");
	fMatPlasticInterior = nistManager->FindOrBuildMaterial("G4_AIR");
	fMatPlasticEntranceWindow = nistManager->FindOrBuildMaterial("G4_MYLAR");
	fMatPlasticCrystal = nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
	fMatPlasticLightGuide = nistManager->FindOrBuildMaterial("G4_PLEXIGLASS");
	fMatPlasticPMT = nistManager->FindOrBuildMaterial("G4_Pyrex_Glass");
	fMatPlasticPMTInterior = nistManager->FindOrBuildMaterial("G4_Galactic");

	// PIPS Detector Materials
    fMatPIPSHousing = nistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL");
	fMatPIPSInterior = nistManager->FindOrBuildMaterial("G4_AIR");
    fMatPIPSDetector = nistManager->FindOrBuildMaterial("G4_Si");
    fMatPIPSElastomer = nistManager->FindOrBuildMaterial("G4_POLYETHYLENE");
    fMatPIPSLBI = nistManager->FindOrBuildMaterial("G4_POLYETHYLENE");
    fMatPIPSRC = nistManager->FindOrBuildMaterial("G4_BRASS");
    fMatPIPSRI = nistManager->FindOrBuildMaterial("G4_POLYETHYLENE");
	fMatPIPSBNCIns = nistManager->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
  	
  	// Print materials
	//G4cout << *(G4Material::GetMaterialTable()) << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{ 	
	// Cleanup old geometry
  	G4GeometryManager::GetInstance()->OpenGeometry();
  	G4PhysicalVolumeStore::GetInstance()->Clean();
  	G4LogicalVolumeStore::GetInstance()->Clean();
  	G4SolidStore::GetInstance()->Clean();	
  	
	////////////////////////////////////////////////////////////////////////
	// Construct The World Volume

	G4double world_X = 2*(sourceRadius + 1.*cm);
	G4double world_Y = world_X;
	G4double world_Z = world_X;
	
	G4Box* WorldSolid = new G4Box("World", world_X/2, world_Y/2, world_Z/2);
  
	WorldLogical = 
		new G4LogicalVolume(WorldSolid,						// The Solid
							fMatWorld,						// Material
							"World");						// Name
  
	WorldPhysical = 
		new G4PVPlacement(	0,								// Rotation
							G4ThreeVector(),				// Translation vector
							WorldLogical,					// Logical volume
							"World",						// Name
							0,								// Mother volume
							false,							// Unused boolean parameter
							0,								// Copy number
							fCheckOverlaps);				// Overlap Check
							

	////////////////////////////////////////////////////////////////////////
	// Construct the Source sphere
	// Note: The actual radius of the Source solid will be slightly smaller (1 mm) than
	// specified in the macro files in order to allow tracking the incident kinetic energy
	// of particles.
	// NOTE: Make sure the outer radius is smaller than the source radius otherwise you'll
	// get "track stuck" warnings during navigation.

	G4Sphere* SourceSolid = new G4Sphere("SourceSolid", sourceRadius - 1.5*mm, sourceRadius - 0.5*mm, 0., 360.0*degree, 0., 180.0*degree);
	SetSourceInnerRadius(SourceSolid->GetInnerRadius());

	SourceLogical = 
		new G4LogicalVolume(SourceSolid,					// The Solid
							fMatWorld,		    			// Material
							"SourceLogical");	     		// Name

	SourcePhysical = 
		new G4PVPlacement(	0,								// No Rotation
							G4ThreeVector(0,0,0),
							SourceLogical,					// Logical volume
							"SourcePhysical",				// Name
							WorldLogical,					// Mother volume
							false,							// Unused boolean parameter
							0,								// Copy number
							fCheckOverlaps);				// Overlap Check

	////////////////////////////////////////////////////////////////////////
	// Construct the Detector Mounting Bracket (DMB)
	//
	// Material: Delrin-150

	// Create the base solid
	G4VSolid* DMB_Solid = new G4UnionSolid("DMB Solid", 
											new G4Box("DMB Base Solid", 0.5*25.4/2*mm, 10*25.4/2*mm, 6*25.4/2*mm), 
											new G4Box("DMB Side Solid", (2.5+3)*25.4/2*mm, 10*25.4/2*mm, 0.75*25.4/2*mm), 
											0, 
											G4ThreeVector((-2.5-3-0.5)*25.4/2*mm,0,(-6+0.75)*25.4/2*mm));

	// Create the cutout for the Plastic Detector
	DMB_Solid = new G4SubtractionSolid("DMB Solid", 
										DMB_Solid, 
										new G4Tubs("Plastic Detector Cutout Solid", 0., 2.66*25.4/2*mm, 1*25.4/2*mm, 0., 360.*deg), 
										0, 
										G4ThreeVector(-(3-0.5/2)*25.4*mm, 2*25.4*mm, (-6+0.75)*25.4/2*mm));
	
	// Create the cutout for the LaBr3 detector
	DMB_Solid = new G4SubtractionSolid("DMB Solid", 
										DMB_Solid, 
										new G4Tubs("LaBr3 Detector Cutout Solid", 0., 3.35*25.4/2*mm, 1*25.4/2*mm, 0., 360.*deg), 
										0, 
										G4ThreeVector(-(3-0.5/2)*25.4*mm, -1.5*25.4*mm, (-6+0.75)*25.4/2*mm));

	// Create the cutout for the PIPS Detector
	DMB_Solid = new G4SubtractionSolid("DMB Solid", 
										DMB_Solid, 
										new G4Tubs("PIPS Detector Cutout Solid", 0., 32.0/2*mm, 1*25.4/2*mm, 0., 360.*deg), 
										0, 
										G4ThreeVector(-(3+3.5/2)*25.4*mm, 0.25*25.4*mm, (-6+0.75)*25.4/2*mm));	

	G4LogicalVolume* DMB_Base_Logical = new G4LogicalVolume(DMB_Solid, fMatDMB, "DMB Logical");

	new G4PVPlacement(G4Transform3D(G4RotationMatrix(0,0,0), 
									G4ThreeVector((2.5+0.5/2)*25.4*mm,0,6*25.4/2*mm)),
					  DMB_Base_Logical,					// Logical volume
					  "DMB Physical",					// Name
					  WorldLogical,						// Mother volume
					  false,							// Unused boolean parameter
					  0,								// Copy number
					  fCheckOverlaps);					// Overlap Check
	

	////////////////////////////////////////////////////////////////////////
	// LaBr3 - Detector Housing
	//
	// Material: Aluminium

	G4double startPhi_LaBr3_Housing = 0.0*deg;
   	G4double endPhi_LaBr3_Housing = 360.0*deg;
   	G4int nrRZ_LaBr3_Housing = 6;
   	G4double zPlane_LaBr3_Housing[]={0, 1.53*25.4*mm, 2.29*25.4*mm, 4.11*25.4*mm, 5.21*25.4*mm, 8.82*25.4*mm};
   	G4double rInner_LaBr3_Housing[]={0, 0, 0, 0, 0, 0};
	G4double rOuter_LaBr3_Housing[]={2.2*25.4/2*mm, 2.2*25.4/2*mm, 3.15*25.4/2*mm, 3.15*25.4/2*mm, 2.31*25.4/2*mm, 2.31*25.4/2*mm};

	G4VSolid* LaBr3_Housing_Solid = new G4Polycone("LaBr3_Housing_Solid",
													startPhi_LaBr3_Housing,
													endPhi_LaBr3_Housing,
													nrRZ_LaBr3_Housing,
													zPlane_LaBr3_Housing,
													rInner_LaBr3_Housing,
													rOuter_LaBr3_Housing);

	G4LogicalVolume* LaBr3_Housing_Logical = new G4LogicalVolume(LaBr3_Housing_Solid, fMatLaBr3Housing, "LaBr3_Housing_Logical");

	new G4PVPlacement(G4Transform3D(G4RotationMatrix(0,0,0), 
									G4ThreeVector(0,-1.5*25.4*mm,-2.5*25.4*mm)),
					  LaBr3_Housing_Logical,			// Logical volume
					  "LaBr3_Housing_Physical",			// Name
					  WorldLogical,						// Mother volume
					  false,							// Unused boolean parameter
					  0,								// Copy number
					  fCheckOverlaps);					// Overlap Check

	////////////////////////////////////////////////////////////////////////
	// LaBr3 - Detector Interior
	//
	// Material: Air

	G4double startPhi_LaBr3_Interior = 0.0*deg;
   	G4double endPhi_LaBr3_Interior = 360.0*deg;
   	G4int nrRZ_LaBr3_Interior = 6;
   	G4double zPlane_LaBr3_Interior[]={0.5*mm, 1.53*25.4*mm, 2.29*25.4*mm, 4.11*25.4*mm, 5.21*25.4*mm, 8.82*25.4*mm - 0.5*mm};
   	G4double rInner_LaBr3_Interior[]={0, 0, 0, 0, 0, 0};
	G4double rOuter_LaBr3_Interior[]={2.2*25.4/2*mm - 0.5*mm, 2.2*25.4/2*mm - 0.5*mm, 3.15*25.4/2*mm - 0.5*mm, 3.15*25.4/2*mm - 0.5*mm, 2.31*25.4/2*mm - 0.5*mm, 2.31*25.4/2*mm - 0.5*mm};

	G4VSolid* LaBr3_Interior_Solid = new G4Polycone("LaBr3_Interior_Solid",
													startPhi_LaBr3_Interior,
													endPhi_LaBr3_Interior,
													nrRZ_LaBr3_Interior,
													zPlane_LaBr3_Interior,
													rInner_LaBr3_Interior,
													rOuter_LaBr3_Interior);

	G4LogicalVolume* LaBr3_Interior_Logical = new G4LogicalVolume(LaBr3_Interior_Solid, fMatLaBr3Interior, "LaBr3_Interior_Logical");

	new G4PVPlacement(G4Transform3D(G4RotationMatrix(0,0,0), 
									G4ThreeVector(0,0,0)),
					  LaBr3_Interior_Logical,			// Logical volume
					  "LaBr3_Interior_Physical",		// Name
					  LaBr3_Housing_Logical,			// Mother volume
					  false,							// Unused boolean parameter
					  0,								// Copy number
					  fCheckOverlaps);					// Overlap Check

	////////////////////////////////////////////////////////////////////////
	// LaBr3 - Detector Reflector
	//
	// Material: Teflon
	// The LaBr3(Ce) crystal is wrapped in a material acting as a Lambertian reflector surface
	// which is designed to increase light collection efficiency 

	G4VSolid* LaBr3_Reflector_Solid = new G4Tubs("LaBr3_Reflector_Solid",
												  0., 
												  50.8*mm/2 + 1.5*mm, 
												  (50.8*mm + 1.5*mm + 5.0*mm)/2,
												  0.,
												  360.*deg);

	G4LogicalVolume* LaBr3_Reflector_Logical = new G4LogicalVolume(LaBr3_Reflector_Solid, fMatLaBr3Reflector, "LaBr3_Reflector_Logical");

	new G4PVPlacement(G4Transform3D(G4RotationMatrix(0, 0, 0), 
									G4ThreeVector(0, 0, 0.5*mm + (50.8*mm + 1.5*mm + 5.0*mm)/2)),
					  LaBr3_Reflector_Logical,			// Logical volume
					  "LaBr3_Reflector_Physical",		// Name
					  LaBr3_Interior_Logical,			// Mother volume
					  false,							// Unused boolean parameter
					  0,								// Copy number
					  fCheckOverlaps);					// Overlap Check

	////////////////////////////////////////////////////////////////////////
	// LaBr3 - Detector Crystal
	//
	// Material: LaBr3(Ce)

	G4VSolid* LaBr3_Crystal_Solid = new G4Tubs("LaBr3_Crystal_Solid",
												0., 
												50.8*mm/2, 
												50.8*mm/2,
												0.,
												360.*deg);

	LaBr3_Crystal_Logical = new G4LogicalVolume(LaBr3_Crystal_Solid, fMatLaBr3Crystal, "LaBr3_Crystal_Logical");

	new G4PVPlacement(G4Transform3D(G4RotationMatrix(0, 0, 0), 
									G4ThreeVector(0, 0, -(5.0*mm-1.5*mm)/2)),
					  LaBr3_Crystal_Logical,			// Logical volume
					  "LaBr3_Crystal_Physical",			// Name
					  LaBr3_Reflector_Logical,			// Mother volume
					  false,							// Unused boolean parameter
					  0,								// Copy number
					  fCheckOverlaps);					// Overlap Check

	// Create a region for the LaBr3 Detector Crystal so we can apply the PAI model to it
	G4Region* reg_LaBr3_Crystal = new G4Region("Region_LaBr3_Crystal");
  	LaBr3_Crystal_Logical->SetRegion(reg_LaBr3_Crystal);
	reg_LaBr3_Crystal->AddRootLogicalVolume(LaBr3_Crystal_Logical);

	////////////////////////////////////////////////////////////////////////
	// LaBr3 - Detector Light Guide
	//
	// Material: Glass

	G4VSolid* LaBr3_LightGuide_Solid = new G4Tubs("LaBr3_LightGuide_Solid",
													0., 
													50.8*mm/2, 
													5.0*mm/2,
													0.,
													360.*deg);

	G4LogicalVolume* LaBr3_LightGuide_Logical = new G4LogicalVolume(LaBr3_LightGuide_Solid, fMatLaBr3LightGuide, "LaBr3_LightGuide_Logical");

	new G4PVPlacement(G4Transform3D(G4RotationMatrix(0, 0, 0), 
									G4ThreeVector(0, 0, +(50.8*mm+1.5*mm)/2)),
					  LaBr3_LightGuide_Logical,			// Logical volume
					  "LaBr3_LightGuide_Physical",		// Name
					  LaBr3_Reflector_Logical,			// Mother volume
					  false,							// Unused boolean parameter
					  0,								// Copy number
					  fCheckOverlaps);					// Overlap Check

	////////////////////////////////////////////////////////////////////////
	// LaBr3 - Detector PMT
	//
	// Material: Glass (0.8 mm)

	G4double startPhi_LaBr3_PMT = 0.0*deg;
   	G4double endPhi_LaBr3_PMT = 360.0*deg;
   	G4int nrRZ_LaBr3_PMT = 4;
   	G4double zPlane_LaBr3_PMT[]={2.29*25.4*mm, 4.11*25.4*mm, 5.21*25.4*mm, 6.06*25.4*mm - 0.5*mm};
   	G4double rInner_LaBr3_PMT[]={0, 0, 0, 0};
	G4double rOuter_LaBr3_PMT[]={3.0*25.4/2*mm - 0.5*mm, 3.0*25.4/2*mm - 0.5*mm, 2.31*25.4/2*mm - 0.5*mm, 2.31*25.4/2*mm - 0.5*mm};

	G4VSolid* LaBr3_PMT_Solid = new G4Polycone("LaBr3_PMT_Solid",
													startPhi_LaBr3_PMT,
													endPhi_LaBr3_PMT,
													nrRZ_LaBr3_PMT,
													zPlane_LaBr3_PMT,
													rInner_LaBr3_PMT,
													rOuter_LaBr3_PMT);

	G4LogicalVolume* LaBr3_PMT_Logical = new G4LogicalVolume(LaBr3_PMT_Solid, fMatLaBr3PMT, "LaBr3_PMT_Logical");

	new G4PVPlacement(G4Transform3D(G4RotationMatrix(0, 0, 0), 
									G4ThreeVector(0, 0, 0)),
					  LaBr3_PMT_Logical,			// Logical volume
					  "LaBr3_PMT_Physical",		// Name
					  LaBr3_Interior_Logical,			// Mother volume
					  false,							// Unused boolean parameter
					  0,								// Copy number
					  fCheckOverlaps);					// Overlap Check

	////////////////////////////////////////////////////////////////////////
	// LaBr3 - Detector PMT Interior
	//
	// Material: Vacuum

	G4double startPhi_LaBr3_PMTInterior = 0.0*deg;
   	G4double endPhi_LaBr3_PMTInterior = 360.0*deg;
   	G4int nrRZ_LaBr3_PMTInterior = 4;
   	G4double zPlane_LaBr3_PMTInterior[]={2.29*25.4*mm + 1*mm, 4.11*25.4*mm, 5.21*25.4*mm, 6.06*25.4*mm - 0.5*mm - 1*mm};
   	G4double rInner_LaBr3_PMTInterior[]={0, 0, 0, 0};
	G4double rOuter_LaBr3_PMTInterior[]={3.0*25.4/2*mm - 0.5*mm - 1*mm, 3.0*25.4/2*mm - 0.5*mm - 1*mm, 2.31*25.4/2*mm - 0.5*mm - 1*mm, 2.31*25.4/2*mm - 0.5*mm - 1*mm};

	G4VSolid* LaBr3_PMTInterior_Solid = new G4Polycone("LaBr3_PMTInterior_Solid",
													startPhi_LaBr3_PMTInterior,
													endPhi_LaBr3_PMTInterior,
													nrRZ_LaBr3_PMTInterior,
													zPlane_LaBr3_PMTInterior,
													rInner_LaBr3_PMTInterior,
													rOuter_LaBr3_PMTInterior);

	G4LogicalVolume* LaBr3_PMTInterior_Logical = new G4LogicalVolume(LaBr3_PMTInterior_Solid, fMatLaBr3PMTInterior, "LaBr3_PMTInterior_Logical");

	new G4PVPlacement(G4Transform3D(G4RotationMatrix(0, 0, 0), 
									G4ThreeVector(0, 0, 0)),
					  LaBr3_PMTInterior_Logical,			// Logical volume
					  "LaBr3_PMTInterior_Physical",		// Name
					  LaBr3_PMT_Logical,			// Mother volume
					  false,							// Unused boolean parameter
					  0,								// Copy number
					  fCheckOverlaps);					// Overlap Check
							
	////////////////////////////////////////////////////////////////////////
	// Plastic - Detector Housing	
	//
	// Material: Aluminium

	G4double startPhi_Plastic_Housing = 0.0*deg;
   	G4double endPhi_Plastic_Housing = 360.0*deg;
   	G4int nrRZ_Plastic_Housing = 8;
   	G4double zPlane_Plastic_Housing[]={0, 11.05*mm, 11.05*mm, 37.77*mm, 37.77*mm, 58.34*mm, 58.34*mm, 240.28*mm};
   	G4double rInner_Plastic_Housing[]={0, 0, 0, 0, 0, 0, 0, 0};
	G4double rOuter_Plastic_Housing[]={76.2/2*mm, 76.2/2*mm, 54.36/2*mm, 54.36/2*mm, 76.2/2*mm, 76.2/2*mm, 60.40/2*mm, 60.40/2*mm};

   	G4VSolid* Plastic_Housing_Solid = new G4Polycone("Plastic_Housing_Solid",
													startPhi_Plastic_Housing,
													endPhi_Plastic_Housing,
													nrRZ_Plastic_Housing,
													zPlane_Plastic_Housing,
													rInner_Plastic_Housing,
													rOuter_Plastic_Housing);
	
	G4LogicalVolume* Plastic_Housing_Logical = new G4LogicalVolume(Plastic_Housing_Solid, fMatPlasticHousing, "Plastic_Housing_Logical");
	
	new G4PVPlacement(G4Transform3D(G4RotationMatrix(0, 0, 0), 
									G4ThreeVector(0, 2*25.4*mm,-2.5*25.4*mm)),
							Plastic_Housing_Logical,					// Logical volume
							"Plastic_Housing_Physical",	// Name
							WorldLogical,					// Mother volume
							false,							// Unused boolean parameter
							0,								// Copy number
							fCheckOverlaps);				// Overlap Check

	////////////////////////////////////////////////////////////////////////
	// Plastic - Detector Interior	
	//
	// Material: Air

	G4double startPhi_Plastic_Interior = 0.0*deg;
   	G4double endPhi_Plastic_Interior = 360.0*deg; 
	G4int nrRZ_Plastic_Interior = 6;
   	G4double zPlane_Plastic_Interior[]={0, 4.7*mm, 4.7*mm, 44.12*mm, 44.12*mm, 236.28*mm};
   	G4double rInner_Plastic_Interior[]={0, 0, 0, 0, 0, 0};
	G4double rOuter_Plastic_Interior[]={55.245/2*mm, 55.245/2*mm, 50.8/2*mm, 50.8/2*mm, 58.4/2*mm, 58.4/2*mm};
   	

   	G4VSolid* Plastic_Interior_Solid = new G4Polycone("Plastic_Interior_Solid",
														startPhi_Plastic_Interior,
														endPhi_Plastic_Interior,
														nrRZ_Plastic_Interior,
														zPlane_Plastic_Interior,
														rInner_Plastic_Interior,
														rOuter_Plastic_Interior);
	
	G4LogicalVolume* Plastic_Interior_Logical = new G4LogicalVolume(Plastic_Interior_Solid, fMatPlasticInterior, "Plastic_Interior_Logical");
			
	new G4PVPlacement(G4Transform3D(G4RotationMatrix(0, 0, 0), 
									G4ThreeVector(0, 0, 0)),
						Plastic_Interior_Logical,					// Logical volume
						"Plastic_Interior_Physical",	// Name
						Plastic_Housing_Logical,					// Mother volume
						false,							// Unused boolean parameter
						0,								// Copy number
						fCheckOverlaps);				// Overlap Check

	////////////////////////////////////////////////////////////////////////
	// Plastic - Detector Entrance Window
	//
	// Material: Mylar

   	G4VSolid* Plastic_EntranceWindow_Solid = new G4Tubs("Plastic_EntranceWindow_Solid",
														0.,
														50./2*mm,
														8.69/2*um,
														0.,
														360.*deg);

	G4LogicalVolume* Plastic_EntranceWindow_Logical = new G4LogicalVolume(Plastic_EntranceWindow_Solid, fMatPlasticEntranceWindow, "Plastic_EntranceWindow_Logical");

	new G4PVPlacement(G4Transform3D(G4RotationMatrix(0, 0, 0), 
									G4ThreeVector(0, 0, 4.7*mm - 8.69/2*um)),
						Plastic_EntranceWindow_Logical,	// Logical volume
						"Plastic_EntranceWindow_Physical",	// Name
						Plastic_Interior_Logical,		// Mother volume
						false,							// Unused boolean parameter
						0,								// Copy number
						fCheckOverlaps);				// Overlap Check

	////////////////////////////////////////////////////////////////////////
	// Plastic - Detector Crystal
	//
	// Material: EJ-204 Scintillator

   	G4VSolid* Plastic_Crystal_Solid = new G4Tubs("Plastic_Crystal_Solid",
												0.,
												50./2*mm,
												20./2*mm,
												0.,
												360.*deg);

	Plastic_Crystal_Logical = new G4LogicalVolume(Plastic_Crystal_Solid, fMatPlasticCrystal, "Plastic_Crystal_Logical");

	new G4PVPlacement(G4Transform3D(G4RotationMatrix(0, 0, 0), 
									G4ThreeVector(0,0,4.7*mm + 20./2*mm)),
						Plastic_Crystal_Logical,		// Logical volume
						"Plastic_Crystal_Physical",		// Name
						Plastic_Interior_Logical,		// Mother volume
						false,							// Unused boolean parameter
						0,								// Copy number
						fCheckOverlaps);				// Overlap Check

	// Create a region for the Plastic Detector Crystal so we can apply the PAI model to it
	G4Region* reg_Plastic_Crystal = new G4Region("Region_Plastic_Crystal");
  	Plastic_Crystal_Logical->SetRegion(reg_Plastic_Crystal);
	reg_Plastic_Crystal->AddRootLogicalVolume(Plastic_Crystal_Logical);
	
	////////////////////////////////////////////////////////////////////////
	// Plastic - Detector Light Guide
	//
	// Material: Plexiglass

   	G4VSolid* Plastic_LightGuide_Solid = new G4Tubs("Plastic_LightGuide_Solid",
													0.,
													50./2*mm,
													20./2*mm,
													0.,
													360.*deg);

	G4LogicalVolume* Plastic_LightGuide_Logical = new G4LogicalVolume(Plastic_LightGuide_Solid, fMatPlasticLightGuide, "Plastic_LightGuide_Logical");

	new G4PVPlacement(G4Transform3D(G4RotationMatrix(0, 0, 0), 
									G4ThreeVector(0,0,4.7*mm + 20.*mm + 20./2*mm)),
						Plastic_LightGuide_Logical,	// Logical volume
						"Plastic_LightGuide_Physical",	// Name
						Plastic_Interior_Logical,		// Mother volume
						false,							// Unused boolean parameter
						0,								// Copy number
						fCheckOverlaps);				// Overlap Check

	////////////////////////////////////////////////////////////////////////
	// Plastic - Detector PMT
	//
	// Material: Glass (0.8 mm)

   	G4VSolid* Plastic_PMT_Solid = new G4Tubs("Plastic_PMT_Solid",
											0.,
											52./2*mm,
											112./2*mm,
											0.,
											360.*deg);

	G4LogicalVolume* Plastic_PMT_Logical = new G4LogicalVolume(Plastic_PMT_Solid, fMatPlasticPMT, "Plastic_PMT_Logical");

	new G4PVPlacement(G4Transform3D(G4RotationMatrix(0, 0, 0), 
									G4ThreeVector(0,0,4.7*mm + 20.*mm + 20.*mm + 112./2*mm)),
						Plastic_PMT_Logical,			// Logical volume
						"Plastic_PMT_Physical",			// Name
						Plastic_Interior_Logical,		// Mother volume
						false,							// Unused boolean parameter
						0,								// Copy number
						fCheckOverlaps);				// Overlap Check
						
	////////////////////////////////////////////////////////////////////////
	// Plastic - Detector PMT Interior
	//
	// Material: Vacuum

   	G4VSolid* Plastic_PMTInterior_Solid = new G4Tubs("Plastic_PMTInterior_Solid",
													0.,
													52./2*mm - 1.*mm,
													112./2*mm - 1.*mm,
													0.,
													360.*deg);

	G4LogicalVolume* Plastic_PMTInterior_Logical = new G4LogicalVolume(Plastic_PMTInterior_Solid, fMatPlasticPMTInterior, "Plastic_PMTInterior_Logical");

	new G4PVPlacement(G4Transform3D(G4RotationMatrix(0, 0, 0), 
									G4ThreeVector(0, 0, 0)),
						Plastic_PMTInterior_Logical,	// Logical volume
						"Plastic_PMTInterior_Physical",	// Name
						Plastic_PMT_Logical,			// Mother volume
						false,							// Unused boolean parameter
						0,								// Copy number
						fCheckOverlaps);				// Overlap Check

	////////////////////////////////////////////////////////////////////////
	// PIPS - Detector Housing	
	//
	// Material: Stainless Steel

	G4double startPhi_PIPS_Housing = 0.0*deg;
   	G4double endPhi_PIPS_Housing = 360.0*deg;
   	G4int nrRZ_PIPS_Housing = 4;
   	G4double zPlane_PIPS_Housing[]={0, 12.3*mm, 12.3*mm, (12.3+7.)*mm};
   	G4double rInner_PIPS_Housing[]={0, 0, 0, 0};
	G4double rOuter_PIPS_Housing[]={32.0/2*mm, 32.0/2*mm, 7./2*mm, 7./2*mm};

   	G4VSolid* PIPS_Housing_Solid = new G4Polycone("PIPS_Housing_Solid",
	   											   startPhi_PIPS_Housing,
												   endPhi_PIPS_Housing,
												   nrRZ_PIPS_Housing,
												   zPlane_PIPS_Housing,
												   rInner_PIPS_Housing,
												   rOuter_PIPS_Housing);
	
	G4LogicalVolume* PIPS_Housing_Logical = new G4LogicalVolume(PIPS_Housing_Solid, fMatPIPSHousing, "PIPS_Housing_Logical");
	
	new G4PVPlacement(G4Transform3D(G4RotationMatrix(0, 0, 0), 
									G4ThreeVector(-2*25.4*mm, 0.25*25.4*mm, 0)),
							PIPS_Housing_Logical,					// Logical volume
							"PIPS_Housing_Physical",	// Name
							WorldLogical,					// Mother volume
							false,							// Unused boolean parameter
							0,								// Copy number
							fCheckOverlaps);				// Overlap Check

	////////////////////////////////////////////////////////////////////////
	// PIPS - Detector Interior	
	//
	// Material: Air

	G4double startPhi_PIPS_Interior = 0.0*deg;
   	G4double endPhi_PIPS_Interior = 360.0*deg; 
	G4int nrRZ_PIPS_Interior = 6;
   	G4double zPlane_PIPS_Interior[]={0, 0.5*mm, 0.5*mm, (0.5+12.3-4.0)*mm, (0.5+12.3-4.0)*mm, (12.3+7)*mm};
   	G4double rInner_PIPS_Interior[]={0, 0, 0, 0, 0, 0};
	G4double rOuter_PIPS_Interior[]={23.9/2*mm, 23.9/2*mm, (32.0/2 - 1.0)*mm, (32.0/2 - 1.0)*mm, 6.25/2*mm, 6.25/2*mm};
   	

   	G4VSolid* PIPS_Interior_Solid = new G4Polycone("PIPS_Interior_Solid",
														startPhi_PIPS_Interior,
														endPhi_PIPS_Interior,
														nrRZ_PIPS_Interior,
														zPlane_PIPS_Interior,
														rInner_PIPS_Interior,
														rOuter_PIPS_Interior);
	
	G4LogicalVolume* PIPS_Interior_Logical = new G4LogicalVolume(PIPS_Interior_Solid, fMatPIPSInterior, "PIPS_Interior_Logical");
			
	new G4PVPlacement(G4Transform3D(G4RotationMatrix(0, 0, 0), 
									G4ThreeVector(0, 0, 0)),
						PIPS_Interior_Logical,					// Logical volume
						"PIPS_Interior_Physical",	// Name
						PIPS_Housing_Logical,					// Mother volume
						false,							// Unused boolean parameter
						0,								// Copy number
						fCheckOverlaps);				// Overlap Check

	////////////////////////////////////////////////////////////////////////
	// PIPS - Elastomeric Ring (Front)
	//
	// Material: Polyethylene

	G4VSolid* PIPS_Elastomeric_Ring_Solid = new G4Tubs("PIPS_Elastomeric_Ring_Solid", 
														23.9/2*mm, 
														23.9/2*mm + 2*0.5*mm, 
														0.5/2*mm, 
														0., 360.*deg);

	G4LogicalVolume* PIPS_Elastomeric_Ring_Logical = new G4LogicalVolume(PIPS_Elastomeric_Ring_Solid, fMatPIPSElastomer, "PIPS_Elastomeric_Ring_Logical");

	new G4PVPlacement(G4Transform3D(G4RotationMatrix(0, 0, 0), 
									G4ThreeVector(0, 0, 0.5*mm + 0.5/2*mm)),
						PIPS_Elastomeric_Ring_Logical,	// Logical volume
						"PIPS_Elastomeric_Ring_Physical",		// Name
						PIPS_Interior_Logical,			// Mother volume
						false,							// Unused boolean parameter
						0,								// Copy number
						fCheckOverlaps);				// Overlap Check

	////////////////////////////////////////////////////////////////////////
	// PIPS - LBI polymer
	//
	// Material: Polyethylene
	
	G4VSolid* PIPS_LBI_Solid = new G4SubtractionSolid("PIPS_LBI_Solid", 
														new G4Tubs("LBI_Cyl", 23.9/2*mm + 2*0.5*mm, (32.0/2 - 1.0)*mm, 2.0/2*mm, 0., 360.*deg), 
														new G4Tubs("LBI_Cyl_Inner", 0., (32.0/2 - 1.0)*mm - 0.3, (2.0 - 0.5)/2*mm, 0., 360.*deg), 
														0, G4ThreeVector(0,0,0.5/2*mm));

	G4LogicalVolume* PIPS_LBI_Logical = new G4LogicalVolume(PIPS_LBI_Solid, fMatPIPSLBI, "PIPS_LBI_Logical");

	new G4PVPlacement(G4Transform3D(G4RotationMatrix(0, 0, 0), 
									G4ThreeVector(0, 0, 0.5*mm + 2.0/2*mm)),
						PIPS_LBI_Logical,	// Logical volume
						"PIPS_LBI_Physical",		// Name
						PIPS_Interior_Logical,			// Mother volume
						false,							// Unused boolean parameter
						0,								// Copy number
						fCheckOverlaps);				// Overlap Check

	////////////////////////////////////////////////////////////////////////
	// PIPS - Detector
	//
	// Material: Si Chip

   	G4VSolid* PIPS_Detector_Solid = new G4Tubs("PIPS_Detector_Solid",
												0.,
												28./2*mm,
												500./2*um,
												0.,
												360.*deg);

	PIPS_Detector_Logical = new G4LogicalVolume(PIPS_Detector_Solid, fMatPIPSDetector, "PIPS_Detector_Logical");

	new G4PVPlacement(G4Transform3D(G4RotationMatrix(0, 0, 0), 
									G4ThreeVector(0, 0, 0.5*mm + 0.5*mm + 500./2*um)),
						PIPS_Detector_Logical,		// Logical volume
						"PIPS_Detector_Physical",		// Name
						PIPS_Interior_Logical,		// Mother volume
						false,							// Unused boolean parameter
						0,								// Copy number
						fCheckOverlaps);				// Overlap Check

	// Create a region for the PIPS Detector so we can apply the PAI model to it
	G4Region* reg_PIPS_Detector = new G4Region("Region_PIPS_Detector");
  	PIPS_Detector_Logical->SetRegion(reg_PIPS_Detector);
	reg_PIPS_Detector->AddRootLogicalVolume(PIPS_Detector_Logical);

	////////////////////////////////////////////////////////////////////////
	// PIPS - Elastomer pad on the back of the Si chip
	//
	// Material: Polyethylene
	
	G4VSolid* PIPS_Elastomeric_Pad_Solid = new G4Tubs("PIPS_Elastomeric_Pad_Solid", 0., 23.9/2*mm, 0.5/2*mm, 0., 360.*deg);

	G4LogicalVolume* PIPS_Elastomeric_Pad_Logical = new G4LogicalVolume(PIPS_Elastomeric_Pad_Solid, fMatPIPSElastomer, "PIPS_Elastomeric_Pad_Logical");

	new G4PVPlacement(G4Transform3D(G4RotationMatrix(0, 0, 0), 
									G4ThreeVector(0, 0, 0.5*mm + 0.5*mm + 500.*um + 0.5/2*mm)),
						PIPS_Elastomeric_Pad_Logical,	// Logical volume
						"PIPS_Elastomeric_Pad_Physical",		// Name
						PIPS_Interior_Logical,			// Mother volume
						false,							// Unused boolean parameter
						0,								// Copy number
						fCheckOverlaps);				// Overlap Check

	////////////////////////////////////////////////////////////////////////
	// PIPS - Rear Contact (RC) on the back of Elastomer pad
	//
	// Material: Brass
	
	G4VSolid* PIPS_RC_Solid = new G4UnionSolid("PIPS_RC_Solid",
													  new G4Tubs("RC_Cyl1", 0., 26./2*mm, 1./2*mm, 0., 360.*deg), 
													  new G4Tubs("RC_Cyl2", 0., 1./2*mm, ((12.3+7.)*mm - 2.*mm)/2, 0., 360.*deg),
													  0, G4ThreeVector(0, 0, ((12.3+7.)*mm - 2.*mm - 1.*mm)/2));

	G4LogicalVolume* PIPS_RC_Logical = new G4LogicalVolume(PIPS_RC_Solid, fMatPIPSRC, "PIPS_RC_Logical");

	new G4PVPlacement(G4Transform3D(G4RotationMatrix(0, 0, 0), 
									G4ThreeVector(0, 0, 0.5*mm + 0.5*mm + 500.*um + 0.5*mm + 1./2*mm)),
						PIPS_RC_Logical,	// Logical volume
						"PIPS_RC_Physical",		// Name
						PIPS_Interior_Logical,			// Mother volume
						false,							// Unused boolean parameter
						0,								// Copy number
						fCheckOverlaps);				// Overlap Check

	////////////////////////////////////////////////////////////////////////
	// PIPS - Rear Insulator (RI)
	//
	// Material: Polyethylene
	
	G4VSolid* PIPS_RI_Solid = new G4Tubs("PIPS_RI_Solid",
										  4./2*mm, 
										  26/2*mm, 
										  4./2*mm, 
										  0., 360.*deg);

	G4LogicalVolume* PIPS_RI_Logical = new G4LogicalVolume(PIPS_RI_Solid, fMatPIPSRI, "PIPS_RI_Logical");

	new G4PVPlacement(G4Transform3D(G4RotationMatrix(0, 0, 0), 
									G4ThreeVector(0, 0, 0.5*mm + 0.5*mm + 500.*um + 0.5*mm + 1.*mm + 4./2*mm)),
						PIPS_RI_Logical,	// Logical volume
						"PIPS_RI_Physical",		// Name
						PIPS_Interior_Logical,			// Mother volume
						false,							// Unused boolean parameter
						0,								// Copy number
						fCheckOverlaps);				// Overlap Check

	////////////////////////////////////////////////////////////////////////
	// PIPS - BNC Insulator
	//
	// Material: Polyethylene
	
	G4VSolid* PIPS_BNC_Ins_Solid = new G4Tubs("PIPS_BNC_Ins_Solid",
										  1./2*mm, 
										  6.25/2*mm, 
										  7./2*mm, 
										  0., 360.*deg);

	G4LogicalVolume* PIPS_BNC_Ins_Logical = new G4LogicalVolume(PIPS_BNC_Ins_Solid, fMatPIPSBNCIns, "PIPS_BNC_Ins_Logical");

	new G4PVPlacement(G4Transform3D(G4RotationMatrix(0, 0, 0), 
									G4ThreeVector(0, 0, 12.3*mm + 7./2*mm)),
						PIPS_BNC_Ins_Logical,	// Logical volume
						"PIPS_BNC_Ins_Physical",		// Name
						PIPS_Interior_Logical,			// Mother volume
						false,							// Unused boolean parameter
						0,								// Copy number
						fCheckOverlaps);				// Overlap Check


	////////////////////////////////////////////////////////////////////////
  	// Visualisation attributes
  	
  	// World Volume (White)
  	G4VisAttributes* Vis_World = new G4VisAttributes(G4Colour(0.1,0.1,0.1,0.4));
  	Vis_World->SetForceWireframe(false);
  	//WorldLogical->SetVisAttributes(Vis_World);
	WorldLogical->SetVisAttributes(G4VisAttributes::GetInvisible());
	  
	// Source Volume (Light Yellow)
    G4VisAttributes* Vis_Source = new G4VisAttributes(G4Colour(1.,1.,1.,0.3));
    Vis_Source->SetForceWireframe(false);
    SourceLogical->SetVisAttributes(Vis_Source);

	// Detector Mounting Bracket (Gray)
    G4VisAttributes* Vis_DMB = new G4VisAttributes(G4Colour(0.1,0.1,0.1,.5));
    Vis_DMB->SetForceWireframe(false);
    DMB_Base_Logical->SetVisAttributes(Vis_DMB);

	// LaBr3 - Housing (Gray)
	G4VisAttributes* Vis_LaBr3_Housing = new G4VisAttributes(G4Colour(0.5,0.5,0.5,0.4));
    Vis_LaBr3_Housing->SetForceWireframe(true);
    LaBr3_Housing_Logical->SetVisAttributes(Vis_LaBr3_Housing);

	// LaBr3 - Interior (Green)
	G4VisAttributes* Vis_LaBr3_Interior = new G4VisAttributes(G4Colour(0.,1.0,0.,0.1));
    Vis_LaBr3_Interior->SetForceWireframe(false);
    LaBr3_Interior_Logical->SetVisAttributes(Vis_LaBr3_Interior);

	// LaBr3 - Reflector (White)
	G4VisAttributes* Vis_LaBr3_Reflector = new G4VisAttributes(G4Colour(1.,1.,1.,0.5));
    Vis_LaBr3_Reflector->SetForceWireframe(false);
    LaBr3_Reflector_Logical->SetVisAttributes(Vis_LaBr3_Reflector);

	// LaBr3 - Crystal (Blue)
	G4VisAttributes* Vis_LaBr3_Crystal = new G4VisAttributes(G4Colour(0.,0.,1.,0.5));
    Vis_LaBr3_Crystal->SetForceWireframe(false);
    LaBr3_Crystal_Logical->SetVisAttributes(Vis_LaBr3_Crystal);

	// LaBr3 - Light Guide (Cyan)
	G4VisAttributes* Vis_LaBr3_LightGuide = new G4VisAttributes(G4Colour(0.,1.0,1.0,0.4));
    Vis_LaBr3_LightGuide->SetForceWireframe(false);
    LaBr3_LightGuide_Logical->SetVisAttributes(Vis_LaBr3_LightGuide);

	// LaBr3 - PMT (Magenta)
	G4VisAttributes* Vis_LaBr3_PMT = new G4VisAttributes(G4Colour(1.,0.,1.,0.5));
    Vis_LaBr3_PMT->SetForceWireframe(false);
    LaBr3_PMT_Logical->SetVisAttributes(Vis_LaBr3_PMT);

	// LaBr3 - PMT Interior (White)
	G4VisAttributes* Vis_LaBr3_PMTInterior = new G4VisAttributes(G4Colour(1.,1.,1.,0.5));
    Vis_LaBr3_PMTInterior->SetForceWireframe(false);
    LaBr3_PMTInterior_Logical->SetVisAttributes(Vis_LaBr3_PMTInterior);

    // Plastic - Housing (Gray)
    G4VisAttributes* Vis_Plastic_Housing = new G4VisAttributes(G4Colour(0.5,0.5,0.5,0.4));
    Vis_Plastic_Housing->SetForceWireframe(true);
    Plastic_Housing_Logical->SetVisAttributes(Vis_Plastic_Housing);

	// Plastic - Interior (Green)
    G4VisAttributes* Vis_Plastic_Interior = new G4VisAttributes(G4Colour(0.,1.0,0.,0.1));
    Vis_Plastic_Interior->SetForceWireframe(false);
    Plastic_Interior_Logical->SetVisAttributes(Vis_Plastic_Interior);

	// Plastic - Entrance Window (Green)
    G4VisAttributes* Vis_Plastic_EntranceWindow = new G4VisAttributes(G4Colour(1.,0.,1.0,0.4));
    Vis_Plastic_EntranceWindow->SetForceWireframe(false);
    Plastic_EntranceWindow_Logical->SetVisAttributes(Vis_Plastic_EntranceWindow);

	// Plastic - Crystal (Blue)
    G4VisAttributes* Vis_Plastic_Crystal = new G4VisAttributes(G4Colour(0.,0.,1.,0.5));
    Vis_Plastic_Crystal->SetForceWireframe(false);
    Plastic_Crystal_Logical->SetVisAttributes(Vis_Plastic_Crystal);

	// Plastic - Light Guide (Cyan)
    G4VisAttributes* Vis_Plastic_LightGuide = new G4VisAttributes(G4Colour(0.,1.0,1.0,0.4));
    Vis_Plastic_LightGuide->SetForceWireframe(false);
    Plastic_LightGuide_Logical->SetVisAttributes(Vis_Plastic_LightGuide);
	
	// Plastic - PMT Glass (Magenta)
    G4VisAttributes* Vis_Plastic_PMT = new G4VisAttributes(G4Colour(1.,0.,1.,0.5));
    Vis_Plastic_PMT->SetForceWireframe(false);
    Plastic_PMT_Logical->SetVisAttributes(Vis_Plastic_PMT);

	// Plastic - PMT Interior (White)
    G4VisAttributes* Vis_Plastic_PMTInterior = new G4VisAttributes(G4Colour(1.,1.,1.,0.5));
    Vis_Plastic_PMTInterior->SetForceWireframe(false);
    Plastic_PMTInterior_Logical->SetVisAttributes(Vis_Plastic_PMTInterior);

	// PIPS - Housing (Gray)
    G4VisAttributes* Vis_PIPS_Housing = new G4VisAttributes(G4Colour(0.5,0.5,0.5,0.4));
    Vis_PIPS_Housing->SetForceWireframe(true);
    PIPS_Housing_Logical->SetVisAttributes(Vis_PIPS_Housing);

	// PIPS - Interior (Green)
    G4VisAttributes* Vis_PIPS_Interior = new G4VisAttributes(G4Colour(0.,1.0,0.,0.1));
    Vis_PIPS_Interior->SetForceWireframe(true);
    PIPS_Interior_Logical->SetVisAttributes(Vis_PIPS_Interior);

	// PIPS - Elastomeric Ring (Magenta)
    G4VisAttributes* Vis_PIPS_ElastomericRing = new G4VisAttributes(G4Colour(1.,0.,1.,1.));
    Vis_PIPS_ElastomericRing->SetForceWireframe(true);
    PIPS_Elastomeric_Ring_Logical->SetVisAttributes(Vis_PIPS_ElastomericRing);

	// PIPS - Elastomeric Ring (Blue)
    G4VisAttributes* Vis_PIPS_LBI = new G4VisAttributes(G4Colour(0.,0.,1.,0.3));
    Vis_PIPS_LBI->SetForceWireframe(true);
    PIPS_LBI_Logical->SetVisAttributes(Vis_PIPS_LBI);

	// PIPS - Detector (Orange)
    G4VisAttributes* Vis_PIPS_Detector = new G4VisAttributes(G4Colour(1.,1.,0.,1.));
    Vis_PIPS_Detector->SetForceWireframe(true);
    PIPS_Detector_Logical->SetVisAttributes(Vis_PIPS_Detector);

	// PIPS - Elastomeric Pad (Magenta)
    G4VisAttributes* Vis_PIPS_ElastomericPad = new G4VisAttributes(G4Colour(1.,0.,1.,1.));
    Vis_PIPS_ElastomericPad->SetForceWireframe(true);
    PIPS_Elastomeric_Pad_Logical->SetVisAttributes(Vis_PIPS_ElastomericPad);

	// PIPS - Rear Contact (Gold)
    G4VisAttributes* Vis_PIPS_RC = new G4VisAttributes(G4Colour(1.,.2,.0,1.));
    Vis_PIPS_RC->SetForceWireframe(true);
    PIPS_RC_Logical->SetVisAttributes(Vis_PIPS_RC);

	// PIPS - Rear Insulator (Dark Gray)
    G4VisAttributes* Vis_PIPS_RI = new G4VisAttributes(G4Colour(0.1,0.1,0.1,1.));
    Vis_PIPS_RI->SetForceWireframe(true);
    PIPS_RI_Logical->SetVisAttributes(Vis_PIPS_RI);

	// PIPS - BNC Insulator (Light Yellow)
    G4VisAttributes* Vis_PIPS_BNSIns = new G4VisAttributes(G4Colour(0.3,0.3,0.,.5));
    Vis_PIPS_BNSIns->SetForceWireframe(true);
    PIPS_BNC_Ins_Logical->SetVisAttributes(Vis_PIPS_BNSIns);

	////////////////////////////////////////////////////////////////////////
	// Return world volume
	return WorldPhysical; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
	G4String filterName, particleName;
  
  	G4SDParticleFilter* gammaFilter = new G4SDParticleFilter(filterName="gammaFilter",particleName="gamma");
  	G4SDParticleFilter* electronFilter = new G4SDParticleFilter(filterName="electronFilter",particleName="e-");
	
	////////////////////////////////////////////////////////////////////////
	// Construct the Multi Functional Detector (MDF) for the particle source scorer

	G4MultiFunctionalDetector* SourceScorer = new G4MultiFunctionalDetector("Source");
	G4SDManager::GetSDMpointer()->AddNewDetector(SourceScorer);	
	G4SDManager::GetSDMpointer()->SetVerboseLevel(0);
	SourceLogical->SetSensitiveDetector(SourceScorer);

	G4VPrimitiveScorer* kinEGamma = new G4PSIncidentKineticEnergy("kinEGamma", fCurrent_Out);
	kinEGamma->SetFilter(gammaFilter);
    SourceScorer->RegisterPrimitive(kinEGamma);

	G4VPrimitiveScorer* kinEElectron = new G4PSIncidentKineticEnergy("kinEElectron", fCurrent_Out);
	kinEElectron->SetFilter(electronFilter);
    SourceScorer->RegisterPrimitive(kinEElectron);

	////////////////////////////////////////////////////////////////////////
	// Construct the MFD for the LaBr3 Detector

	G4MultiFunctionalDetector* LaBr3_Scorer = new G4MultiFunctionalDetector("LaBr3_Detector");
	G4SDManager::GetSDMpointer()->AddNewDetector(LaBr3_Scorer);	
	G4SDManager::GetSDMpointer()->SetVerboseLevel(0);
	LaBr3_Crystal_Logical->SetSensitiveDetector(LaBr3_Scorer);
 	
    LaBr3_Scorer->RegisterPrimitive(new G4PSEnergyDeposit("eDep"));

	////////////////////////////////////////////////////////////////////////
	// Construct the MFD for the Plastic Detector

	G4MultiFunctionalDetector* Plastic_Scorer = new G4MultiFunctionalDetector("Plastic_Detector");
	G4SDManager::GetSDMpointer()->AddNewDetector(Plastic_Scorer);	
	G4SDManager::GetSDMpointer()->SetVerboseLevel(0);
	Plastic_Crystal_Logical->SetSensitiveDetector(Plastic_Scorer);
 	
    Plastic_Scorer->RegisterPrimitive(new G4PSEnergyDeposit("eDep"));

	////////////////////////////////////////////////////////////////////////
	// Construct the MFD for the PIPS Detector

	G4MultiFunctionalDetector* PIPS_Scorer = new G4MultiFunctionalDetector("PIPS_Detector");
	G4SDManager::GetSDMpointer()->AddNewDetector(PIPS_Scorer);	
	G4SDManager::GetSDMpointer()->SetVerboseLevel(0);
	PIPS_Detector_Logical->SetSensitiveDetector(PIPS_Scorer);
 	
    PIPS_Scorer->RegisterPrimitive(new G4PSEnergyDeposit("eDep"));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSourceRadius(G4double val)
{
	if(WorldPhysical) {    
    	G4Exception ("DetectorConstruction::SetSourceRadius()", "G4McMasterTDS", 
                 	JustWarning, 
                 	"Attempt to change already constructed geometry is ignored");
  	} else {
   		sourceRadius = val;
  	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSourceInnerRadius(G4double val)
{
	sourceInnerRadius = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetSourceRadius()
{
	// Return the source radius
	return sourceRadius;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetSourceInnerRadius()
{
	// Return the source inner radius
	return sourceInnerRadius;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineCommands()
{
    // Define command directory using generic messenger class
    fMessenger = new G4GenericMessenger(this, "/G4McMasterTDS/", "Geometry control");

	// Source Radius Command
	G4GenericMessenger::Command& SourceRadiusCmd
      = fMessenger->DeclareMethodWithUnit("SourceRadius","cm",
                                  &DetectorConstruction::SetSourceRadius, 
                                  "Set the radius of the source volume within the world volume.");
    SourceRadiusCmd.SetParameterName("radius", true);
    SourceRadiusCmd.SetDefaultValue("20.0");
}
