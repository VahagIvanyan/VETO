#include "EventAction.hh"

#include "Run.hh"
#include "HistoManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(DetectorConstruction* det)
:G4UserEventAction(), fDetector(det)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
  //energy deposited per event
Total = 0,nVeto = 0,nJPET = 0, existInVeto = false, jpetExist = false;
  for (G4int k=0; k<kMaxAbsor; k++) { fEdepAbsor[k] = 0.0; }
  for (G4int numS=0; numS<mMaxScint; numS++) { fEdepScint[numS] = 0.0; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* )
{
  //get Run
  Run* run = static_cast<Run*>(
             G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  G4int maxNbOfAbsor = fDetector->GetNbOfAbsor();
  G4int maxNbOfScint = fDetector->GetNbOfScint();
  //G4int NbEvt = evt->GetEventID();
  //plot energy deposited per event
  //

  G4double TotalEdep(0.);
  for (G4int k=1; k<=maxNbOfAbsor; k++) {
    if ((fEdepAbsor[k] > 500.*keV)){ 
//&&(fEdepAbsor[k] <= 500.*keV)) {
	existInVeto = true;	
	run->AddEdep(k,fEdepAbsor[k]);
	TotalEdep += fEdepAbsor[k];
	nVeto++;
     	  G4AnalysisManager::Instance()->FillH1(62, k);           
	  G4AnalysisManager::Instance()->FillH1(k, fEdepAbsor[k]);
     } 

  }
  if (TotalEdep > 0.) {
    run->AddTotEdep(TotalEdep);
  }
G4double TotalEdep2(0.);
	  	for (G4int numS=1,numhisto = 43; numS<=maxNbOfScint; numS++,numhisto++) {
	    	if ((fEdepScint[numS] > 0.*keV)&&(fEdepScint[numS] <= 1000.*keV)){
			jpetExist = true;
			run->AddScintEdep(numS,fEdepScint[numS]);
		      	TotalEdep2 += fEdepScint[numS];
		      	//G4AnalysisManager::Instance()->FillH1(63, numS);
			nJPET++;
		switch(nJPET){
				case 0:
				G4AnalysisManager::Instance()->FillH1(66,nJPET);
				break;
				default:
				G4AnalysisManager::Instance()->FillH1(67,nJPET);
				break;
	}
			G4AnalysisManager::Instance()->FillH1(64,nJPET);
		      if((numS>=1)&&(numS<=15)){
			   G4AnalysisManager::Instance()->FillH1(numhisto, fEdepScint[numS]);
		      }
	    	}			
	  	}
	  if (TotalEdep2 > 0) {
	    run->Add2TotEdep(TotalEdep2);   
	  }

	
	//if(nVeto > 0){
		if((existInVeto == true) || (jpetExist == true)){
			switch(nJPET){
				case 0:
				G4AnalysisManager::Instance()->FillH1(81,nJPET);
				break;
				default:
				G4AnalysisManager::Instance()->FillH1(82,nJPET);
				break;
			}
		}
	//}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

