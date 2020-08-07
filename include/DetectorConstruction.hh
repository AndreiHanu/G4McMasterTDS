#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class G4GenericMessenger;

class DetectorConstruction : public G4VUserDetectorConstruction
{
    public:
  	    // Constructor
        DetectorConstruction();
        // Destructor
        virtual ~DetectorConstruction();

        // Defines the detector geometry and returns a pointer to the physical World Volume
        virtual G4VPhysicalVolume* Construct();
    
        // Sensitive Detector 
	    virtual void ConstructSDandField();
    
        // Set Methods
        void SetSourceRadius(G4double val);
    
        // Get Methods
        G4double GetSourceRadius();
        G4double GetSourceInnerRadius();
    
    private:
        // Set Methods
        void SetSourceInnerRadius(G4double val);
        
        // Defines all the detector materials
        void DefineMaterials();
    
        // Define commands to change the geometry
        void DefineCommands();
    
        G4GenericMessenger* fMessenger;
        G4bool fCheckOverlaps;
    
        // Standard Materials
        G4Material* fMatWorld;
        G4Material* fMatDMB;

        G4Material* fMatLaBr3Housing;
        G4Material* fMatLaBr3Interior;
        G4Material* fMatLaBr3Reflector;
        G4Material* fMatLaBr3Crystal;
        G4Material* fMatLaBr3LightGuide;
        G4Material* fMatLaBr3PMT;
        G4Material* fMatLaBr3PMTInterior;

        G4Material* fMatPlasticHousing;
        G4Material* fMatPlasticInterior;
        G4Material* fMatPlasticEntranceWindow;
        G4Material* fMatPlasticCrystal;
        G4Material* fMatPlasticLightGuide;
        G4Material* fMatPlasticPMT;
        G4Material* fMatPlasticPMTInterior;

        G4Material* fMatPIPSHousing;
        G4Material* fMatPIPSInterior;
        G4Material* fMatPIPSDetector;
        G4Material* fMatPIPSElastomer;
        G4Material* fMatPIPSLBI;
        G4Material* fMatPIPSRC;
        G4Material* fMatPIPSRI;
        G4Material* fMatPIPSBNCIns;
    
        // Logical Volumes
        G4LogicalVolume* WorldLogical;
        G4LogicalVolume* SourceLogical;
        G4LogicalVolume* LaBr3_Crystal_Logical;
        G4LogicalVolume* Plastic_Crystal_Logical;
        G4LogicalVolume* PIPS_Detector_Logical;
    
        // Physical Volumes
        G4VPhysicalVolume* WorldPhysical;
        G4VPhysicalVolume* SourcePhysical;

	    // Rotation Angles
        G4double sourceRadius;
        G4double sourceInnerRadius;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

