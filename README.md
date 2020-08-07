## Welcome to G4McMasterTDS
G4McMasterTDS is an application based on the Geant4 Monte Carlo toolkit for simulating the response (i.e. geometric factor, efficiency, energy migration matrix, etc) of the McMaster University triple detector system (TDS) that was used for measuring the beta- and gamma-ray source term within CANDU nuclear power plants.

The McMaster TDS consists of a Saint Gobain BrilLanCe 380 LaBr3(Ce) detector, Eljen Technologies M550-20x8-1 Plastic Scintillator (EJ-204) detector, and a Canberra PD450-15-500AM Passivated Implanted Planar Silicon (PIPS) detector. To enable simultaneous and co-located measurements, all three detectors are mounted in a common detector braket.

## Build Notes
To build this application, you should have a working version of CMake
and Geant 4.10.3.p01 or higher.

Geant4: https://geant4.web.cern.ch/geant4/ CMake: https://cmake.org/

I recommend setting up your directory structure as follows:

- $G4WORKDIR/G4McMasterTDS : This is the source directory and is a clone of the code from this GitHub 
- $G4WORKDIR/G4McMasterTDS-build : Contains directories and files for each build. This directory also contains the executable files and example macros.

The simulation works both in sequential and in multi-threaded (MT) modes of Geant4. However, a significant speedup can be achieved by running the simulation in MT mode.

## Steps to compile:
### Step 1 - Source the Geant4 environment setup script
### NOTE: You'll need to change the version of Geant4 in the path

    source /opt/Geant4/geant4.10.03.p01-install/bin/geant4.sh

### Step 2 - Create the build directory and navigate to it
    
    mkdir G4McMasterTDS-build && cd $_

### Step 3 - Setup CMake, make the build, and run the build

    sudo cmake -DGeant4_DIR=/opt/Geant4/geant4.10.03.p01-install/lib/Geant4-10.3.1/ ~/G4WORK/G4McMasterTDS; sudo make -j8; ./G4McMasterTDS


To recompile, I typically just re-run the command in Step 3

## Scoring Physical Quantities
Scoring of physical quantities is carried out through the G4MultiFunctionalDetector which allows multiple primitive scorers to be assigned to a single volume. In this simulation, a G4MultiFunctionalDetector is assigned to the sensitive volume of Silicon chip and the source volume. The value recorded by the primitive scorers is collected on event-by-event basis enabling calculation of the detector response.

The primitive scorer(s) registered for the sensitive volume of Si chip are:
* G4PSEnergyDeposit

The primitive scorer(s) registered for the source volume are custom kinetic energy scorers:
* G4PSIncidentKineticEnergy

During a run, the data collected on an event-by-event basis is histogrammed into five logarithmically binned ROOT histograms. 

* The first histogram, labeled "Source Fluence (Gamma)", records fluence and kinetic energy of incident gamma-ray source particles.
* The second histogram, labeled "Source Fluence (Electron)", records fluence and kinetic energy of incident electron source particles.
* The third histogram, labeled "LaBr3 Detector - Measured Spectrum", records the measured or deposited energy in the Saint Gobain BrilLanCe 380 LaBr3(Ce) detector.
* The fourth histogram, labeled "LaBr3 Detector - Energy Migration Matrix (Gamma)", records a 2D histogram of the true kinetic energy versus measured/deposited energy in the Saint Gobain BrilLanCe 380 LaBr3(Ce) detector from gamma-rays.
* The fifth histogram, labeled "LaBr3 Detector - Energy Migration Matrix (Electron)", records a 2D histogram of the true kinetic energy versus measured/deposited energy in the Saint Gobain BrilLanCe 380 LaBr3(Ce) detector from electrons.
* The sixth histogram, labeled "Plastic Detector -  Measured Spectrum", records the measured or deposited energy in the Eljen Technologies M550-20x8-1 Plastic Scintillator (EJ-204) detector.
* The seventh histogram, labeled "Plastic Detector -  Energy Migration Matrix (Gamma)", records a 2D histogram of the true kinetic energy versus measured/deposited energy in the Eljen Technologies M550-20x8-1 Plastic Scintillator (EJ-204) detector from gamma-rays.
* The eight histogram, labeled "Plastic Detector -  Energy Migration Matrix (Electron)", records a 2D histogram of the true kinetic energy versus measured/deposited energy in the Eljen Technologies M550-20x8-1 Plastic Scintillator (EJ-204) detector from electrons.
* The ninth histogram, labeled "PIPS Detector -  Measured Spectrum", records the measured or deposited energy in the Canberra PD450-15-500AM PIPS detector.
* The tenth histogram, labeled "PIPS Detector -  Energy Migration Matrix (Gamma)", records a 2D histogram of the true kinetic energy versus measured/deposited energy in the Canberra PD450-15-500AM PIPS detector from gamma-rays.
* The eleventh histogram, labeled "PIPS Detector -  Energy Migration Matrix (Electron)", records a 2D histogram of the true kinetic energy versus measured/deposited energy in the Canberra PD450-15-500AM PIPS detector from electrons.

## Simulation Pictures

