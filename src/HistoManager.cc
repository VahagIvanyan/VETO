#include "HistoManager.hh"
#include "DetectorConstruction.hh"
#include "G4UnitsTable.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
HistoManager::HistoManager()
  : fFileName("Veto")
{
  Book();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
HistoManager::~HistoManager()
{
  delete G4AnalysisManager::Instance();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HistoManager::Book()
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetFileName(fFileName);
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);
  // Define histograms start values
  const G4int kMaxHisto = 96;
  const G4String id[] = { "0", "Absorber 1", "Absorber 2", "Absorber 3", "Absorber 4", "Absorber 5", "Absorber 6", "Absorber 7", "Absorber 8", "Absorber 9","Absorber 10","Absorber 11", "Absorber 12", "Absorber 13", "Absorber 14", "Absorber 15", "Absorber 16", "Absorber 17", "Absorber 18", "Absorber 19","Absorber 20","Absorber 21", "Absorber 22", "Absorber 23", "Absorber 24", "Absorber 25", "Absorber 26", "Absorber 27", "Absorber 28", "Absorber 29","Absorber 30","Absorber 31", "Absorber 32", "Absorber 33", "Absorber 34", "Absorber 35", "Absorber 36", "Absorber 37", "Absorber 38", "Absorber 39","Absorber 40","Absorber 41", "Absorber 42","Scintillator 1","Scintillator 2","Scintillator 3","Scintillator 4","Scintillator 5","Scintillator 6","Scintillator 7","Scintillator 8","Scintillator 9", "Scintillator 10","Scintillator 11","Scintillator 12","Scintillator 13","Scintillator 14","Scintillator 15","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95"};

  const G4String title[] = 
                { "dummy",                                        	//0
                  "total Energy deposited in absorber 1",   		
                  "total Energy deposited in absorber 2",       	
                  "total Energy deposited in absorber 3",         	
                  "total Energy deposited in absorber 4",         	
                  "total Energy deposited in absorber 5",         	
                  "total Energy deposited in absorber 6",         	
		  "total Energy deposited in absorber 7",         	
                  "total Energy deposited in absorber 8",         	
                  "total Energy deposited in absorber 9",         	
                  "total Energy deposited in absorber 10", 
                  "total Energy deposited in absorber 11",   		
                  "total Energy deposited in absorber 12",       	
                  "total Energy deposited in absorber 13",         	
                  "total Energy deposited in absorber 14",         	
                  "total Energy deposited in absorber 15",         	
                  "total Energy deposited in absorber 16",         	
		  "total Energy deposited in absorber 17",         	
                  "total Energy deposited in absorber 18",         	
                  "total Energy deposited in absorber 19",         	
                  "total Energy deposited in absorber 20",
                  "total Energy deposited in absorber 21",   		
                  "total Energy deposited in absorber 22",       	
                  "total Energy deposited in absorber 23",         	
                  "total Energy deposited in absorber 24",         	
                  "total Energy deposited in absorber 25",         	
                  "total Energy deposited in absorber 26",         	
		  "total Energy deposited in absorber 27",         	
                  "total Energy deposited in absorber 28",         	
                  "total Energy deposited in absorber 29",         	
                  "total Energy deposited in absorber 30",
                  "total Energy deposited in absorber 31",   		
                  "total Energy deposited in absorber 32",       	
                  "total Energy deposited in absorber 33",         	
                  "total Energy deposited in absorber 34",         	
                  "total Energy deposited in absorber 35",         	
                  "total Energy deposited in absorber 36",         	
		  "total Energy deposited in absorber 37",         	
                  "total Energy deposited in absorber 38",         	
                  "total Energy deposited in absorber 39",         	
                  "total Energy deposited in absorber 40", 
		  "total Energy deposited in absorber 41",   		
                  "total Energy deposited in absorber 42",         	
		  "total Energy deposited in Scintillator 1",         	
                  "total Energy deposited in Scintillator 2",         	
                  "total Energy deposited in Scintillator 3",         	
                  "total Energy deposited in Scintillator 4",         	
		  "total Energy deposited in Scintillator 5",         	
                  "total Energy deposited in Scintillator 6",         	
                  "total Energy deposited in Scintillator 7",         	
                  "total Energy deposited in Scintillator 8",         	
		  "total Energy deposited in Scintillator 9",         	
                  "total Energy deposited in Scintillator 10",         	
                  "total Energy deposited in Scintillator 11",         	
                  "total Energy deposited in Scintillator 12",  
                  "total Energy deposited in Scintillator 13",         	
                  "total Energy deposited in Scintillator 14",         	
                  "total Energy deposited in Scintillator 15",         	
       	          "Theta",           					
                  "Phi",     						
                  "Edep (MeV/mm) along absorbers",     			
		  "Angular distribution of Muons(Theta)",
                  "Muons on Scintillators of Veto",
		  "Muons on Scintillators of JPET",
		  "Muons on VetoAndJPET",
		  "Muons on JPET",
		  "Veto for JPET",
			"67",
			"68",
			"69",
			"70",
			"71",
			"72",
			"73",
			"74",
			"75",
			"76",
			"77",
			"78",
			"79",
			"80",
			"81",
			"82",
			"83",
			"84",
			"85",
			"86",
			"87",
			"88",
			"89",
			"90",
			"91",
			"92",
			"93",
			"94",
			"95"
                };
 //G4String name = "scin_"+G4UIcommand::ConvertToString(icopy);
 //G4AnalysisManager::Instance()->FillH2(1,nJPET,NbEvt);
  // Default values (to be reset via /analysis/h1/set command)
  G4int nbins = 100;
  G4double vmin = 0.;
  G4double vmax = 100.;

  // Create all histograms as inactivated 
  // as we have not yet set nbins, vmin, vmax
  for (G4int k=0; k<kMaxHisto; k++) {
    G4int ih = analysisManager->CreateH1(id[k], title[k], nbins, vmin, vmax);
    analysisManager->SetH1Activation(ih, false);
  }
  const G4int DMaxHisto = 5;
  const G4String ij[] = { "h2.0", "JPET","X and Y","Y and Z","X and Z"};
  const G4String title2D[] = 
                { "2D dummy",                  	//0
                  "Enteries per Scintillator",  //1
		  "X and Y",			//2
		  "Y and Z",			//3
		  "X and Z"			//4
                };
 for (G4int D=0; D<DMaxHisto; D++) {
    G4int il = analysisManager->CreateH2(ij[D], title2D[D], nbins, vmin, vmax, nbins,vmin,vmax);
    analysisManager->SetH2Activation(il, false);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
