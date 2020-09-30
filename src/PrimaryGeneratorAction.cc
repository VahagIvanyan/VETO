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
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
// $Id: PrimaryGeneratorAction.cc 95740 2016-02-23 09:34:37Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorMessenger.hh"
#include "HistoManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* det)
:G4VUserPrimaryGeneratorAction(),
 fDetector(det),fRndmBeam(0.),fGunMessenger(0)
{
  fParticleGun  = new G4ParticleGun(1);
  SetDefaultKinematic();
  //create a messenger for this class
  fGunMessenger = new PrimaryGeneratorMessenger(this); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fGunMessenger;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void PrimaryGeneratorAction::SetDefaultKinematic()
{
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*---------------------------------------------------------**/
//G4AnalysisManager* analysisManager1 = G4AnalysisManager::Instance();
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

G4double FluxMu;
G4ParticleDefinition* particle;

 FluxMu = G4UniformRand()*100;
if(FluxMu > 44){
  particle = G4ParticleTable::GetParticleTable()->FindParticle("mu+");
}
else{
  particle = G4ParticleTable::GetParticleTable()->FindParticle("mu-");
}
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleEnergy(4*GeV);  

 G4double thetas,phi,R_Sphere = 5*m;

G4double ThetaStep = twopi/96,ThetaStep2 = twopi/48,ThetaStep3 = twopi/48;
G4int steps = (G4int)(G4UniformRand()*96),steps2 = (G4int)(G4UniformRand()*48), steps3 = (G4int)(G4UniformRand()*48);
//steps2 = steps3;
//G4AnalysisManager::Instance()->FillH1(32, steps);
G4double R = 56.55*cm,R_2 = 46.65*cm,R_3 = 42.4*cm,xCryst,yCryst,zCryst,xMidCryst,yMidCryst,xInnerCryst,yInnerCryst;
G4double FullLength = 1.7*cm;
G4double FullWidth = 5*mm;
G4double length = (G4UniformRand()*FullLength) - FullLength/2;
G4double width = (G4UniformRand()*FullWidth) - FullWidth/2;

G4double ScintCenter = steps*ThetaStep + 0.032724923;//Outer Layer
G4double ScintMidCenter = steps2*ThetaStep2 + 0.06544985;//Mid Layer
G4double ScintInnerCenter = steps3*ThetaStep3;//Inner Layer

//Selection of hit's position from scintillator. Each scintillator SHOULD BE hited. 

xCryst = (R+length)*std::cos(ScintCenter); 
yCryst = (R+width)*std::sin(ScintCenter);
//xCryst = R*std::cos(ScintCenter); 
//yCryst = R*std::sin(ScintCenter);
zCryst = (G4UniformRand()*49.8*cm) - 24.9*cm;
  //For Second Ring Crystals
xMidCryst = (R_2+length)*std::cos(ScintMidCenter); 
yMidCryst = (R_2+width)*std::sin(ScintMidCenter);
 //xMidCryst = R_2*std::cos(ScintMidCenter); 
 //yMidCryst = R_2*std::sin(ScintMidCenter);
  //For Third Ring Crystals
xInnerCryst = (R_3+length)*std::cos(ScintInnerCenter); 
yInnerCryst = (R_3+width)*std::sin(ScintInnerCenter);
 //xInnerCryst = R_3*std::cos(ScintInnerCenter); 
 //yInnerCryst = R_3*std::sin(ScintInnerCenter);
  G4AnalysisManager::Instance()->FillH2(2,xCryst,yCryst);
  G4AnalysisManager::Instance()->FillH2(3,yCryst,zCryst);
  G4AnalysisManager::Instance()->FillH2(4,xCryst,zCryst);
//G4AnalysisManager::Instance()->FillH2(5,xMidCryst,zCryst);
G4double X_pos,Y_pos,Z_pos;
G4double sinTheta,cosTheta,ThetaAngle; 
G4double u,v,w;
  phi = G4UniformRand()*twopi; 

//Outer Ring hits	
if((ScintCenter>=90)&&(ScintCenter<=270))
{
	do{
		ThetaAngle = -G4UniformRand()*twopi/4; 
		cosTheta = std::cos(ThetaAngle);
 		if(G4UniformRand() < cosTheta*cosTheta) 
 		break;
	}
	while(true);
	sinTheta = std::sqrt(1-pow(cosTheta,2));
	thetas = std::acos(cosTheta);
	X_pos =xCryst + R_Sphere*cosTheta;
	Y_pos =yCryst + R_Sphere*sinTheta*std::sin(phi);
	Z_pos =zCryst + R_Sphere*std::cos(phi)*sinTheta;
	u = -cosTheta;
	v = -sinTheta*std::sin(phi);
	w = -std::cos(phi)*sinTheta;
}
else
{
	do{
		ThetaAngle = G4UniformRand()*twopi/4; 
		cosTheta = std::cos(ThetaAngle);	
		if(G4UniformRand() < cosTheta*cosTheta) 
 		break;
	}
	while(true);
	sinTheta = std::sqrt(1-pow(cosTheta,2));
	thetas = std::acos(cosTheta);
	X_pos =xCryst + R_Sphere*cosTheta;
	Y_pos =yCryst + R_Sphere*sinTheta*std::sin(phi);
	Z_pos =zCryst + R_Sphere*std::cos(phi)*sinTheta;
	u = -cosTheta;
	v = -sinTheta*std::sin(phi);
	w = -std::cos(phi)*sinTheta;
}
/// Middle Ring hits
if((ScintMidCenter>=90)&&(ScintMidCenter<=270))
{
	do{
		ThetaAngle = -G4UniformRand()*twopi/4; 
		cosTheta = std::cos(ThetaAngle);
 		if(G4UniformRand() < cosTheta*cosTheta) 
 		break;
	}
	while(true);
	sinTheta = std::sqrt(1-pow(cosTheta,2));
	thetas = std::acos(cosTheta);
	X_pos =xMidCryst + R_Sphere*cosTheta;
	Y_pos =yMidCryst + R_Sphere*sinTheta*std::sin(phi);
	Z_pos =zCryst + R_Sphere*std::cos(phi)*sinTheta;
	u = -cosTheta;
	v = -sinTheta*std::sin(phi);
	w = -std::cos(phi)*sinTheta;
}
else
{
	do{
		///// This is an incorrect way of cosTheta selection. Should be changed!!!!
		ThetaAngle = G4UniformRand()*twopi/4; 
		cosTheta = std::cos(ThetaAngle);	
		if(G4UniformRand() < cosTheta*cosTheta) 
 		break;
	}
	while(true);
	sinTheta = std::sqrt(1-pow(cosTheta,2));
	thetas = std::acos(cosTheta);
	X_pos =xMidCryst + R_Sphere*cosTheta;
	Y_pos =yMidCryst + R_Sphere*sinTheta*std::sin(phi);
	Z_pos =zCryst + R_Sphere*std::cos(phi)*sinTheta;
	u = -cosTheta;
	v = -sinTheta*std::sin(phi);
	w = -std::cos(phi)*sinTheta;
}
/// Inner Ring hits	
if((ScintInnerCenter>=90)&&(ScintInnerCenter<=270))
{
	do{
		ThetaAngle = -G4UniformRand()*twopi/4; 
		cosTheta = std::cos(ThetaAngle);
 		if(G4UniformRand() < cosTheta*cosTheta) 
 		break;
	}
	while(true);
	sinTheta = std::sqrt(1-pow(cosTheta,2));
	thetas = std::acos(cosTheta);
	X_pos =xInnerCryst + R_Sphere*cosTheta;
	Y_pos =yInnerCryst + R_Sphere*sinTheta*std::sin(phi);
	Z_pos =zCryst + R_Sphere*std::cos(phi)*sinTheta;
	u = -cosTheta;
	v = -sinTheta*std::sin(phi);
	w = -std::cos(phi)*sinTheta;
}
else
{
	do{
		ThetaAngle = G4UniformRand()*twopi/4; 
		cosTheta = std::cos(ThetaAngle);	
		if(G4UniformRand() < cosTheta*cosTheta) 
 		break;
	}
	while(true);
	sinTheta = std::sqrt(1-pow(cosTheta,2));
	thetas = std::acos(cosTheta);
	X_pos =xInnerCryst + R_Sphere*cosTheta;
	Y_pos =yInnerCryst + R_Sphere*sinTheta*std::sin(phi);
	Z_pos =zCryst + R_Sphere*std::cos(phi)*sinTheta;
	u = -cosTheta;
	v = -sinTheta*std::sin(phi);
	w = -std::cos(phi)*sinTheta;
}
	G4AnalysisManager::Instance()->FillH1(58, thetas);
	fParticleGun->SetParticlePosition(G4ThreeVector(X_pos,Y_pos,Z_pos));
	G4AnalysisManager::Instance()->FillH1(59, phi);
 	fParticleGun->SetParticleMomentumDirection(G4ThreeVector(u,v,w)); 

 //create vertex
 fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

