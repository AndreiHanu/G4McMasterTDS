#include "Run.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "G4Event.hh"
#include "G4Run.hh"
#include "G4Track.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4THitsMap.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

// Select output format for Analysis Manager
#include "Analysis.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Run::Run(DetectorConstruction* det, PrimaryGeneratorAction* primary):G4Run(),
detector(det), particleGun(primary)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::RecordEvent(const G4Event* event)
{ 	
  	// Get hits collections
  	G4HCofThisEvent* HCE = event->GetHCofThisEvent();
  	if(!HCE) { 
    	G4ExceptionDescription msg; 
    	msg << "No hits collection of this event found.\n"; 
    	G4Exception("Run::RecordEvent()","Code001", JustWarning, msg); 
    	return; 
  	} 
  	
	// Zero out the variables
	G4double eDep_LaBr3 = 0.;
	G4double eDep_Plastic = 0.;
	G4double eDep_PIPS = 0.;
	G4double kinEGamma = 0.;
	G4double kinEElectron = 0.;
	
	// Get the HitMaps for this event
	G4THitsMap<G4double>* event_eDep_LaBr3 = (G4THitsMap<G4double>*)(HCE->GetHC(G4SDManager::GetSDMpointer()->GetCollectionID("LaBr3_Detector/eDep")));
	G4THitsMap<G4double>* event_eDep_Plastic = (G4THitsMap<G4double>*)(HCE->GetHC(G4SDManager::GetSDMpointer()->GetCollectionID("Plastic_Detector/eDep")));
	G4THitsMap<G4double>* event_eDep_PIPS = (G4THitsMap<G4double>*)(HCE->GetHC(G4SDManager::GetSDMpointer()->GetCollectionID("PIPS_Detector/eDep")));
	G4THitsMap<G4double>* event_kinEGamma = (G4THitsMap<G4double>*)(HCE->GetHC(G4SDManager::GetSDMpointer()->GetCollectionID("Source/kinEGamma")));
	G4THitsMap<G4double>* event_kinEElectron = (G4THitsMap<G4double>*)(HCE->GetHC(G4SDManager::GetSDMpointer()->GetCollectionID("Source/kinEElectron")));
	
	std::map<G4int,G4double*>::iterator itr;

	// Get analysis manager
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

	// Calculate the fluence for this event taking into account any angular biasing
	G4double fluence = 1/(3.14159*std::pow(detector->GetSourceInnerRadius()/cm, 2)*(std::pow(std::sin(particleGun->GetGPS()->GetCurrentSource()->GetAngDist()->GetMaxTheta()), 2)-std::pow(std::sin(particleGun->GetGPS()->GetCurrentSource()->GetAngDist()->GetMinTheta()), 2)));
	//G4double fluence = 1.;
	
	// Get the total energy deposited for this event
	for (itr = event_eDep_LaBr3->GetMap()->begin(); itr != event_eDep_LaBr3->GetMap()->end(); itr++) {
		eDep_LaBr3 += *(itr->second);
	}

	for (itr = event_eDep_Plastic->GetMap()->begin(); itr != event_eDep_Plastic->GetMap()->end(); itr++) {
		eDep_Plastic += *(itr->second);
	}

	for (itr = event_eDep_PIPS->GetMap()->begin(); itr != event_eDep_PIPS->GetMap()->end(); itr++) {
		eDep_PIPS += *(itr->second);
	}

	// Get the incident kinetic energy for gammas in this event
	for (itr = event_kinEGamma->GetMap()->begin(); itr != event_kinEGamma->GetMap()->end(); itr++) {
		kinEGamma = *(itr->second);

		if (kinEGamma > 0){ 
			analysisManager->FillH1(analysisManager->GetH1Id("Source Spectrum (Gamma)"), kinEGamma/keV, fluence);

			// Add it to the energy migration matrices for each detector
			
			if (eDep_LaBr3 > 0 && eDep_LaBr3 < kinEGamma){
				analysisManager->FillH2(analysisManager->GetH2Id("LaBr3 Detector - Energy Migration Matrix (Gamma)"), kinEGamma/keV, eDep_LaBr3/keV);
			}

			if (eDep_Plastic > 0 && eDep_Plastic < kinEGamma){
				analysisManager->FillH2(analysisManager->GetH2Id("Plastic Detector - Energy Migration Matrix (Gamma)"), kinEGamma/keV, eDep_Plastic/keV);
			}

			if (eDep_PIPS > 0 && eDep_PIPS < kinEGamma){
				analysisManager->FillH2(analysisManager->GetH2Id("PIPS Detector - Energy Migration Matrix (Gamma)"), kinEGamma/keV, eDep_PIPS/keV);
			}
		}
	}

	// Get the incident kinetic energy for electrons in this event
	for (itr = event_kinEElectron->GetMap()->begin(); itr != event_kinEElectron->GetMap()->end(); itr++) {
		kinEElectron = *(itr->second);

		if (kinEElectron > 0){
			analysisManager->FillH1(analysisManager->GetH1Id("Source Spectrum (Electron)"), kinEElectron/keV, fluence);

			// Add it to the energy migration matrices for each detector

			if (eDep_LaBr3 > 0 && eDep_LaBr3 < kinEElectron){
				analysisManager->FillH2(analysisManager->GetH2Id("LaBr3 Detector - Energy Migration Matrix (Electron)"), kinEElectron/keV, eDep_LaBr3/keV);
			}

			if (eDep_Plastic > 0 && eDep_Plastic < kinEElectron){
				analysisManager->FillH2(analysisManager->GetH2Id("Plastic Detector - Energy Migration Matrix (Electron)"), kinEElectron/keV, eDep_Plastic/keV);
			}

			if (eDep_PIPS > 0 && eDep_PIPS < kinEElectron){
				analysisManager->FillH2(analysisManager->GetH2Id("PIPS Detector - Energy Migration Matrix (Electron)"), kinEElectron/keV, eDep_PIPS/keV);
			}
		}
	}

	// Record events with non-zero deposited energy

	if (eDep_LaBr3 > 0) {
		analysisManager->FillH1(analysisManager->GetH1Id("LaBr3 Detector - Measured Spectrum"), eDep_LaBr3/keV);
	}

	if (eDep_Plastic > 0) {
		analysisManager->FillH1(analysisManager->GetH1Id("Plastic Detector - Measured Spectrum"), eDep_Plastic/keV);
	}

	if (eDep_PIPS > 0) {
		analysisManager->FillH1(analysisManager->GetH1Id("PIPS Detector - Measured Spectrum"), eDep_PIPS/keV);
	}
	
	// Invoke base class method
  	G4Run::RecordEvent(event); 
}