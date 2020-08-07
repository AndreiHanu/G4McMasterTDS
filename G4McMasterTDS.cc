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
// G4McMasterTDS.cc
//
// Description: 
//
// Main Geant4 simulation program for the McMaster Triple Detector System (TDS).
//
// The McMaster TDS consists of a Saint Gobain BrilLanCe 380 LaBr3(Ce) detector, 
// Eljen Technologies M550-20x8-1 Plastic Scintillator (EJ-204) detector, and a 
// Canberra PD450-15-500AM Passivated Implanted Planar Silicon (PIPS) detector. 
// To enable simultaneous and co-located measurements, all three detectors are mounted 
// in a common detector braket.
//
// The TDS is a spectroscopic system used by McMaster University 
// to measure the beta- and gamma-ray source term in CANDU nuclear power plants 
// in order to perform dosimetry assessments for the lens of the eye.
// ********************************************************************

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "Randomize.hh"

// Simulation Files
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
	// Detect interactive mode (if no arguments) and define UI session
  	//
  	G4UIExecutive* ui = 0;
  	if ( argc == 1 ) {
    		ui = new G4UIExecutive(argc, argv);
  	}
  
	// Choose the Random engine
  	G4Random::setTheEngine(new CLHEP::RanecuEngine);
     
  	// Construct the default run manager
	#ifdef G4MULTITHREADED
		G4cout << "Simulation is running on " << G4Threading::G4GetNumberOfCores() << " threads";
        G4MTRunManager* runManager = new G4MTRunManager;
        G4int nCores = G4Threading::G4GetNumberOfCores();
        runManager->SetNumberOfThreads(nCores);
	#else
  		G4RunManager* runManager = new G4RunManager;
	#endif

	//G4RunManager* runManager = new G4RunManager;

	// Set mandatory initialization classes
  	DetectorConstruction* detector = new DetectorConstruction();
  	runManager->SetUserInitialization(detector);
	  
	// Select a physics list
	PhysicsList* physicsList = new PhysicsList();
  	runManager->SetUserInitialization(physicsList);

	// Set user action classes
	runManager->SetUserInitialization(new ActionInitialization(detector)); 
  
	// Initialize visualization
	G4VisManager* visManager = new G4VisExecutive;
	// G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
	// G4VisManager* visManager = new G4VisExecutive("Quiet");
	visManager->Initialize();

  	// Get the pointer to the User Interface manager
  	G4UImanager* UImanager = G4UImanager::GetUIpointer();

  	// Process macro or start UI session
	if ( ! ui ) 
	{ 
		// batch mode
		G4String command = "/control/execute ";
		G4String fileName = argv[1];
		UImanager->ApplyCommand(command+fileName);
	} 
	else
	{  
		// interactive mode : define UI session
		ui->SessionStart();
		delete ui;
	}

	// Job termination
	// Free the store: user actions, physics_list and detector_description are
	// owned and deleted by the run manager, so they should not be deleted 
	// in the main() program !

	delete visManager;
	delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
