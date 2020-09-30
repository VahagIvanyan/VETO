#include "DetectorMessenger.hh"

#include <sstream>

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:G4UImessenger(),fDetector(Det),
 fTestemDir(0),
 fDetDir(0), 
 fNbAbsorCmd(0), 
 fNbScintCmd(0),          
 fAbsorCmd(0),
 fScintCmd(0),
 fSizeXCmd(0),
 fSizeZCmd(0),
 fSizeSXCmd(0),
 fSizeSZCmd(0),
 fIsotopeCmd(0)
{ 
  fTestemDir = new G4UIdirectory("/testhadr/");
  fTestemDir->SetGuidance(" detector control.");
  
  fDetDir = new G4UIdirectory("/testhadr/det/");
  fDetDir->SetGuidance("detector construction commands");
  
  fNbAbsorCmd = new G4UIcmdWithAnInteger("/testhadr/det/setNbOfAbsor",this);
  fNbAbsorCmd->SetGuidance("Set number of Absorbers.");
  fNbAbsorCmd->SetParameterName("NbAbsor",false);
  fNbAbsorCmd->SetRange("NbAbsor>0");
  fNbAbsorCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fNbAbsorCmd->SetToBeBroadcasted(false);

  fNbScintCmd = new G4UIcmdWithAnInteger("/testhadr/det/setNbOfScint",this);
  fNbScintCmd->SetGuidance("Set number of Scintilators.");
  fNbScintCmd->SetParameterName("NbScint",false);
  fNbScintCmd->SetRange("NbScint>0");
  fNbScintCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fNbScintCmd->SetToBeBroadcasted(false);

   
  fAbsorCmd = new G4UIcommand("/testhadr/det/setAbsor",this);
  fAbsorCmd->SetGuidance("Set the absor nb, the material, the thickness.");
  fAbsorCmd->SetGuidance("  absor number : from 1 to NbOfAbsor");
  fAbsorCmd->SetGuidance("  material name");
  fAbsorCmd->SetGuidance("  thickness (with unit) : t>0.");

  fScintCmd = new G4UIcommand("/testhadr/det/setScint",this);
  fScintCmd->SetGuidance("Set the Scint nb, the material, the thickness.");
  fScintCmd->SetGuidance("  Scint number : from 1 to NbOfScint");
  fScintCmd->SetGuidance("  material name");
  fScintCmd->SetGuidance("  thickness (with unit) : t>0.");
  //
  G4UIparameter* AbsNbPrm = new G4UIparameter("AbsorNb",'i',false);
  AbsNbPrm->SetGuidance("absor number : from 1 to NbOfAbsor");
  AbsNbPrm->SetParameterRange("AbsorNb>0");
  fAbsorCmd->SetParameter(AbsNbPrm);
  //
  G4UIparameter* MatPrm = new G4UIparameter("material",'s',false);
  MatPrm->SetGuidance("material name");
  fAbsorCmd->SetParameter(MatPrm);
  //    
  G4UIparameter* ThickPrm = new G4UIparameter("thickness",'d',false);
  ThickPrm->SetGuidance("thickness of absorber");
  ThickPrm->SetParameterRange("thickness>0.");
  fAbsorCmd->SetParameter(ThickPrm);
  //
  G4UIparameter* unitPrm = new G4UIparameter("unit",'s',false);
  unitPrm->SetGuidance("unit of thickness");
  G4String unitList = G4UIcommand::UnitsList(G4UIcommand::CategoryOf("mm"));
  unitPrm->SetParameterCandidates(unitList);
  fAbsorCmd->SetParameter(unitPrm);
  //
  fAbsorCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fAbsorCmd->SetToBeBroadcasted(false);
  ///
  G4UIparameter* ScintNbPrm = new G4UIparameter("ScintNb",'i',false);
  ScintNbPrm->SetGuidance("Scint number : from 1 to NbOfScint");
  ScintNbPrm->SetParameterRange("ScintNb>0");
  fScintCmd->SetParameter(ScintNbPrm);
  //
  G4UIparameter* MatScintPrm = new G4UIparameter("material",'s',false);
  MatScintPrm->SetGuidance("material name");
  fScintCmd->SetParameter(MatScintPrm);
  //    
  G4UIparameter* ThickPrm2 = new G4UIparameter("thickness",'d',false);
  ThickPrm2->SetGuidance("thickness of Scint");
  ThickPrm2->SetParameterRange("thickness>0.");
  fScintCmd->SetParameter(ThickPrm2);
  //
  G4UIparameter* unitPrm2 = new G4UIparameter("unit",'s',false);
  unitPrm2->SetGuidance("unit of thickness");
  G4String unitList2 = G4UIcommand::UnitsList(G4UIcommand::CategoryOf("mm"));
  unitPrm2->SetParameterCandidates(unitList2);
  fScintCmd->SetParameter(unitPrm2);
  //
  fScintCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fScintCmd->SetToBeBroadcasted(false);
  ///
  fSizeXCmd = new G4UIcmdWithADoubleAndUnit("/testhadr/det/setSizeX",this);
  fSizeXCmd->SetGuidance("Set sizeX of the absorber");
  fSizeXCmd->SetParameterName("SizeX",false);
  fSizeXCmd->SetRange("SizeX>0.");
  fSizeXCmd->SetUnitCategory("Length");
  fSizeXCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fSizeXCmd->SetToBeBroadcasted(false);

  fSizeZCmd = new G4UIcmdWithADoubleAndUnit("/testhadr/det/setSizeZ",this);
  fSizeZCmd->SetGuidance("Set sizeZ of the absorber");
  fSizeZCmd->SetParameterName("SizeZ",false);
  fSizeZCmd->SetRange("SizeZ>0.");
  fSizeZCmd->SetUnitCategory("Length");
  fSizeZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fSizeZCmd->SetToBeBroadcasted(false);

  fSizeSXCmd = new G4UIcmdWithADoubleAndUnit("/testhadr/det/setSizeSX",this);
  fSizeSXCmd->SetGuidance("Set sizeX of the scintilator");
  fSizeSXCmd->SetParameterName("SizeSX",false);
  fSizeSXCmd->SetRange("SizeSX>0.");
  fSizeSXCmd->SetUnitCategory("Length");
  fSizeSXCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fSizeSXCmd->SetToBeBroadcasted(false);

  fSizeSZCmd = new G4UIcmdWithADoubleAndUnit("/testhadr/det/setSizeSZ",this);
  fSizeSZCmd->SetGuidance("Set sizeZ of the scintilator");
  fSizeSZCmd->SetParameterName("SizeSZ",false);
  fSizeSZCmd->SetRange("SizeSZ>0.");
  fSizeSZCmd->SetUnitCategory("Length");
  fSizeSZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fSizeSZCmd->SetToBeBroadcasted(false);

  fIsotopeCmd = new G4UIcommand("/testhadr/det/setIsotopeMat",this);
  fIsotopeCmd->SetGuidance("Build and select a material with single isotope");
  fIsotopeCmd->SetGuidance("  symbol of isotope, Z, A, density of material");
  //
  G4UIparameter* symbPrm = new G4UIparameter("isotope",'s',false);
  symbPrm->SetGuidance("isotope symbol");
  fIsotopeCmd->SetParameter(symbPrm);
  //      
  G4UIparameter* ZPrm = new G4UIparameter("Z",'i',false);
  ZPrm->SetGuidance("Z");
  ZPrm->SetParameterRange("Z>0");
  fIsotopeCmd->SetParameter(ZPrm);
  //      
  G4UIparameter* APrm = new G4UIparameter("A",'i',false);
  APrm->SetGuidance("A");
  APrm->SetParameterRange("A>0");
  fIsotopeCmd->SetParameter(APrm);  
  //    
  G4UIparameter* densityPrm = new G4UIparameter("density",'d',false);
  densityPrm->SetGuidance("density of material");
  densityPrm->SetParameterRange("density>0.");
  fIsotopeCmd->SetParameter(densityPrm);
  //
  G4UIparameter* unitPrm1 = new G4UIparameter("unit",'s',false);
  unitPrm1->SetGuidance("unit of density");
  G4String unitList1 = G4UIcommand::UnitsList(G4UIcommand::CategoryOf("g/cm3"));
  unitPrm1->SetParameterCandidates(unitList1);
  fIsotopeCmd->SetParameter(unitPrm1);
  //
  fIsotopeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fNbAbsorCmd;
  delete fNbScintCmd;
  delete fAbsorCmd;
  delete fScintCmd;
  delete fSizeXCmd;
  delete fSizeZCmd;
  delete fSizeSXCmd;
  delete fSizeSZCmd;
  delete fIsotopeCmd;
  delete fDetDir;
  delete fTestemDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{    
  if( command == fNbAbsorCmd )
   { fDetector->SetNbOfAbsor(fNbAbsorCmd->GetNewIntValue(newValue));}
   
 if( command == fNbScintCmd )
   { fDetector->SetNbOfScint(fNbScintCmd->GetNewIntValue(newValue));}


  if (command == fAbsorCmd)
   {
     G4int num; G4double tick;
     G4String unt, mat;
     std::istringstream is(newValue);
     is >> num >> mat >> tick >> unt;
     G4String material=mat;
     tick *= G4UIcommand::ValueOf(unt);
     fDetector->SetAbsorMaterial (num,material);
     fDetector->SetAbsorThickness(num,tick);
   }
   
  if (command == fScintCmd)
   {
     G4int num2; G4double tick2;
     G4String unt2, mat2;
     std::istringstream is(newValue);
     is >> num2 >> mat2 >> tick2 >> unt2;
     G4String material2=mat2;
     tick2 *= G4UIcommand::ValueOf(unt2);
     fDetector->SetScintMaterial (num2,material2);
     fDetector->SetScintThickness(num2,tick2);
   }
  if( command == fSizeXCmd )
   { fDetector->SetAbsorSizeX(fSizeXCmd->GetNewDoubleValue(newValue));}   
  if( command == fSizeZCmd )
   { fDetector->SetAbsorSizeZ(fSizeZCmd->GetNewDoubleValue(newValue));}   

  if( command == fSizeSXCmd )
   { fDetector->SetScintSizeX(fSizeSXCmd->GetNewDoubleValue(newValue));}   
  if( command == fSizeSZCmd )
   { fDetector->SetScintSizeZ(fSizeSZCmd->GetNewDoubleValue(newValue));}   

  if (command == fIsotopeCmd)
   {
     G4int Z; G4int A; G4double dens;
     G4String name, unt;
     std::istringstream is(newValue);
     is >> name >> Z >> A >> dens >> unt;
     dens *= G4UIcommand::ValueOf(unt);
     fDetector->MaterialWithSingleIsotope (name,name,dens,Z,A);
   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
