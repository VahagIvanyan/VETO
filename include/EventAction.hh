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
/// \file EventAction.hh
/// \brief Definition of the EventAction class
//
// $Id: EventAction.hh 95740 2016-02-23 09:34:37Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "DetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
  public:
    EventAction(DetectorConstruction*);
   ~EventAction();

  public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);
    
    void AddEdep(G4int k, G4double edep) { fEdepAbsor[k] += edep; }
    void AddCutEdep(G4int k, G4double edep) { fEdepAbsor[k] += edep; }
    void AddScintEdep(G4int numS, G4double scint_edep) { fEdepScint[numS] += scint_edep; }
  private:
    DetectorConstruction* fDetector;
	G4int Total,nVeto,nJPET;
	G4bool jpetExist, existInVeto;
    G4double  fEdepAbsor[kMaxAbsor],fEdepScint[mMaxScint];
};
/*

switch(nJPET){
	case 0:
	G4AnalysisManager::Instance()->FillH1(66,nJPET);
	break;
	case 1:
	G4AnalysisManager::Instance()->FillH1(67,nJPET);
	break;
	case 2:
	G4AnalysisManager::Instance()->FillH1(68,nJPET);
	break;
	case 3:
	G4AnalysisManager::Instance()->FillH1(69,nJPET);
	break;
	case 4:
	G4AnalysisManager::Instance()->FillH1(70,nJPET);
	break;
	case 5:
	G4AnalysisManager::Instance()->FillH1(71,nJPET);
	break;
	case 6:
	G4AnalysisManager::Instance()->FillH1(72,nJPET);
	break;
	case 7:
	G4AnalysisManager::Instance()->FillH1(73,nJPET);
	break;
	case 8:
	G4AnalysisManager::Instance()->FillH1(74,nJPET);
	break;
	case 9:
	G4AnalysisManager::Instance()->FillH1(75,nJPET);
	break;
	case 10:
	G4AnalysisManager::Instance()->FillH1(76,nJPET);
	break;
	case 11:
	G4AnalysisManager::Instance()->FillH1(77,nJPET);
	break;
	case 12:
	G4AnalysisManager::Instance()->FillH1(78,nJPET);
	break;
	case 13:
	G4AnalysisManager::Instance()->FillH1(79,nJPET);
	break;
	case 14:
	G4AnalysisManager::Instance()->FillH1(80,nJPET);
	break;
	default:
		G4cout<<nJPET<<G4endl;
	break;

}
if(nVeto == 1){
	switch(nJPET){
	case 0:
	G4AnalysisManager::Instance()->FillH1(81,nJPET);
	break;
	case 1:
	G4AnalysisManager::Instance()->FillH1(82,nJPET);
	break;
	case 2:
	G4AnalysisManager::Instance()->FillH1(83,nJPET);
	break;
	case 3:
	G4AnalysisManager::Instance()->FillH1(84,nJPET);
	break;
	case 4:
	G4AnalysisManager::Instance()->FillH1(85,nJPET);
	break;
	case 5:
	G4AnalysisManager::Instance()->FillH1(86,nJPET);
	break;
	case 6:
	G4AnalysisManager::Instance()->FillH1(87,nJPET);
	break;
	case 7:
	G4AnalysisManager::Instance()->FillH1(88,nJPET);
	break;
	case 8:
	G4AnalysisManager::Instance()->FillH1(89,nJPET);
	break;
	case 9:
	G4AnalysisManager::Instance()->FillH1(90,nJPET);
	break;
	case 10:
	G4AnalysisManager::Instance()->FillH1(91,nJPET);
	break;
	case 11:
	G4AnalysisManager::Instance()->FillH1(92,nJPET);
	break;
	case 12:
	G4AnalysisManager::Instance()->FillH1(93,nJPET);
	break;
	case 13:
	G4AnalysisManager::Instance()->FillH1(94,nJPET);
	break;
	case 14:
	G4AnalysisManager::Instance()->FillH1(95,nJPET);
	break;
	default:
		G4cout<<nJPET<<G4endl;
	break;

}
}
*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
