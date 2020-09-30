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
/// \file Run.hh
/// \brief Definition of the Run class
//
// $Id: Run.hh 71375 2013-06-14 07:39:33Z maire $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Run_h
#define Run_h 1

#include "DetectorConstruction.hh"

#include "G4Run.hh"
#include "G4VProcess.hh"
#include "globals.hh"
#include <map>

class DetectorConstruction;
class G4ParticleDefinition;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Run : public G4Run
{
  public:
    Run(DetectorConstruction* detector);
   ~Run();

  public:
    void SetPrimary(G4ParticleDefinition* particle, G4double energy);
    void CountProcesses(const G4VProcess* process);
    void ParticleCount(G4int, G4String, G4double);
    void AddEdep (G4int i, G4double e);
    void AddTotEdep     (G4double e);
    void AddScintEdep (G4int j, G4double e2);
    void Add2TotEdep     (G4double e2);
    void AddTrackStatus (G4int i);
    //void Add2TrackStatus (G4int j)

    virtual void Merge(const G4Run*);
    void EndOfRun();
    
private:
  struct ParticleData {
   ParticleData()
     : fCount(0), fEmean(0.), fEmin(0.), fEmax(0.) {}
   ParticleData(G4int count, G4double ekin, G4double emin, G4double emax)
     : fCount(count), fEmean(ekin), fEmin(emin), fEmax(emax) {}
   G4int     fCount;
   G4double  fEmean;
   G4double  fEmin;
   G4double  fEmax;
  };
 struct ParticlesData {
   ParticlesData()
     : f2Count(0), f2Emean(0.), f2Emin(0.), f2Emax(0.) {}
   ParticlesData(G4int counts, G4double e2kin, G4double e2min, G4double e2max)
     : f2Count(counts), f2Emean(e2kin), f2Emin(e2min), f2Emax(e2max) {}
   G4int     f2Count;
   G4double  f2Emean;
   G4double  f2Emin;
   G4double  f2Emax;
  };
  private:
    DetectorConstruction*  fDetector;
    G4ParticleDefinition*  fParticle;
    G4double  fEkin; 

    G4int      fStatus[3];

    G4double   fEdeposit[kMaxAbsor],  fEmin[kMaxAbsor], fEmax[kMaxAbsor];
    G4double   f2Edeposit[mMaxScint], f2Emin[mMaxScint], f2Emax[mMaxScint];

    G4double   fTotEdep[3];
    G4double   f2TotEdep[3];

    std::map<G4String,G4int>        fProcCounter;
    std::map<G4String,ParticleData> fParticleDataMap[kMaxAbsor];
    std::map<G4String,ParticlesData> f2ParticleDataMap[mMaxScint];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

