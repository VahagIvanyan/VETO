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
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
// $Id: DetectorConstruction.hh 78560 2014-01-07 10:06:52Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Cache.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class DetectorMessenger;

const G4int kMaxAbsor = 43;                        
const G4int mMaxScint = 193;                       
class G4GlobalMagFieldMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  
  DetectorConstruction();
 ~DetectorConstruction();

public:

  G4Material* 
  MaterialWithSingleIsotope(G4String, G4String, G4double, G4int, G4int);
  
  void SetNbOfAbsor     (G4int);      
  void SetAbsorMaterial (G4int,const G4String&);     
  void SetAbsorThickness(G4int,G4double);
                
  void SetAbsorSizeX   (G4double);          
  void SetAbsorSizeZ   (G4double); 

  void SetNbOfScint     (G4int);      
  void SetScintMaterial (G4int,const G4String&);     
  void SetScintThickness(G4int,G4double);
                
  void SetScintSizeX   (G4double);          
  void SetScintSizeZ   (G4double); 
         
  virtual G4VPhysicalVolume* Construct();
  virtual void ConstructSDandField();
     
public:

  G4int       GetNbOfAbsor()             {return nb_plastic;}  
  G4int	      GetNbOfScint()   		 {return nb_cryst;} 

  G4Material* GetAbsorMaterial (G4int i) {return fAbsorMaterial[i];};
  G4double    GetAbsorThickness(G4int i) {return fAbsorThickness[i];}; 

  G4Material* GetScintMaterial (G4int j) {return fScintMaterial[j];};
  G4double    GetScintThickness(G4int j) {return fScintThickness[j];}; 
   
  G4double    GetXfront        (G4int i) {return fXfront[i];};
            
  G4double GetAbsorSizeX()              {return fAbsorSizeX;}; 
  G4double GetAbsorSizeY()              {return fAbsorSizeY;};
  G4double GetAbsorSizeZ()              {return fAbsorSizeZ;};

  G4double GetScintSizeX()              {return fScintSizeX;}; 
  G4double GetScintSizeY()              {return fScintSizeY;};
  G4double GetScintSizeZ()              {return fScintSizeZ;};

  G4double GetWorldSizeX()               {return fWorldSizeXY;}; 
  G4double GetWorldSizeYZ()              {return fWorldSizeZ;}; 
  
  void PrintParameters();
   
private:

  G4int		     nb_plastic;
  G4int              nb_cryst;
  G4Material*        fAbsorMaterial [kMaxAbsor];
  G4double           fAbsorThickness[kMaxAbsor];
  G4Material*        fScintMaterial[mMaxScint];
  G4double           fScintThickness[mMaxScint];

  G4double           fXfront[kMaxAbsor];  

  G4double           fAbsorSizeX;
  G4double           fAbsorSizeY;
  G4double           fAbsorSizeZ;

  G4double           fScintSizeX;
  G4double           fScintSizeY;
  G4double           fScintSizeZ;

  G4double           fWorldSizeXY;
  G4double           fWorldSizeZ; 
  
  G4Material*        fDefaultMaterial;
  
  G4Material*		fWindowMaterial;
  G4Material*		fTargetMaterial;
  
  G4Material*		fDetectorMaterial;

 //G4LogicalVolume**  fLogicCrist;
  //G4LogicalVolume** flogicAbsor;
  G4VPhysicalVolume* fPhysiWorld;
 
  DetectorMessenger* fDetectorMessenger;
  G4Cache<G4GlobalMagFieldMessenger*> fFieldMessenger;

private:

  void DefineMaterials();
  void ComputeParameters();
  G4VPhysicalVolume* ConstructVolumes();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

