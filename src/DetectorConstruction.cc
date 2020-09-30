#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Trap.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4RunManager.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),fDefaultMaterial(0),fWindowMaterial(0),fTargetMaterial(0),fPhysiWorld(0),
 fDetectorMessenger(0)
{
  // default parameter values of the absorbers
  nb_plastic = 42;
  nb_cryst = 192;
 fAbsorThickness[0] = 0*mm;        //dummy, for initialization
 fScintThickness[0] = 0*mm;

  DefineMaterials();

for(G4int amount = 1; amount <= nb_plastic; amount++)  { 
  fAbsorThickness[amount] = 2*cm; 
  fAbsorSizeX = 1*m;
  fAbsorSizeZ = 5*cm;
  //ComputeParameters();
  // materials
  SetAbsorMaterial(amount,"PlasticScintillator");
}
for(G4int ScintNumber = 1; ScintNumber <= nb_cryst; ScintNumber++)  { 
  fScintThickness[ScintNumber] = 19*mm; 
  fScintSizeX = 0.5*m;
  fScintSizeZ = 7*mm;
  //ComputeParameters();
  // materials
  SetScintMaterial(ScintNumber,"G4_PLASTIC_SC_VINYLTOLUENE");
}
  ComputeParameters();
  // create commands for interactive definition of the calorimeter
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  G4NistManager* man = G4NistManager::Instance();
  man->FindOrBuildMaterial("G4_AIR");

	G4Element* H = man->FindOrBuildElement("H");
	G4Element* C = man->FindOrBuildElement("C");
	//G4Element* N = man->FindOrBuildElement("N");
	G4Element* O = man->FindOrBuildElement("O");
	//G4Element* F = man->FindOrBuildElement("F"); 
	G4Element* Na = man->FindOrBuildElement("Na");
	G4Element* Cl = man->FindOrBuildElement("Cl");
	G4Element* Fe = man->FindOrBuildElement("Fe");
	G4Element* Al = man->FindOrBuildElement("Al");
	G4Element* Pb = man->FindOrBuildElement("Pb");
	G4Element* Mo = man->FindOrBuildElement("Mo");
	G4Element* Si = man->FindOrBuildElement("Si");
    	G4Element* Bi = man->FindOrBuildElement("Bi");
	G4Element* W = man->FindOrBuildElement("W");

	///////////////////////////////////////////////////////////////////////////////////
	///				  AIR					        ///
	///////////////////////////////////////////////////////////////////////////////////
	
	G4double density,temperature,pressure; 
	/*
		density = 0.001225*g/cm3;
       		pressure = 98658.96*pascal;
		temperature = 273*kelvin; 
		G4Material* Air = new G4Material("Air", density, 2, kStateGas, temperature, pressure);
		Air->AddElement(N, 79.*perCent);
		Air->AddElement(O, 21.*perCent);
		G4NistManager *nist_man = G4NistManager::Instance();
		G4Material *AIR_mat = nist_man->FindOrBuildMaterial("Air");
	*/
	/////Materials

	G4Material* Graphite = new G4Material("Graphite", 2.266 *g/cm3, 1);
	Graphite->AddElement(C,1);

	G4Material* LiquidCarbon = new G4Material("LiquidCarbon", 1.2 *g/cm3, 1);
	LiquidCarbon->AddElement(C,1);
	
	G4Material* Hydrogen = new G4Material("Hydrogen", 0.071 *g/cm3, 1);
	Hydrogen->AddElement(H,1);
	
	  G4Material* H2O = 
	  new G4Material("Water", 1.000*g/cm3, 2);
	  H2O->AddElement(H, 2);
	  H2O->AddElement(O, 1);
	  H2O->GetIonisation()->SetMeanExcitationEnergy(78.0*eV);
 
	  //GetIonisation()->SetMeanExcitationEnergy(78.0*eV);
  
	  G4Material* NaCl = 
	  new G4Material("Salt", 2.165*g/cm3, 2);
	  NaCl->AddElement(Na, 1);
	  NaCl->AddElement(Cl, 1);

	  // Example of Vacuum				
	  density     = universe_mean_density;    //from PhysicalConstants.h
	  pressure    = 3.e-18*pascal;
	  temperature = 2.73*kelvin;
	  G4Material* Galactic =   
	  new G4Material("Galactic", 1., 1.008*g/mole, density,
		                     kStateGas,temperature,pressure);

    //  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
	
	//SiO2 
	G4Material* SiO2 = new G4Material("Silicon_Dioxide", 2.196 *g/cm3, 2);
	SiO2->AddElement(Si,1);
	SiO2->AddElement(O,2);

	//Fe 
	G4Material* Iron = new G4Material("Iron", 7.850 *g/cm3, 1);
	Iron->AddElement(Fe,1);

	//Lead
	G4Material* Lead = new G4Material("Lead", 11.34 *g/cm3, 1);
	Lead->AddElement(Pb,1);

	//Tungsten
	G4Material* Tungsten = new G4Material("Tungsten", 19.25 *g/cm3, 1);
	Tungsten->AddElement(W,1);

	//Molybdenium
	G4Material* Molybdenium = new G4Material("Molybdenium", 10.28 *g/cm3, 1);
	Molybdenium->AddElement(Mo,1);

	//Bismuth
	G4Material* Bismuth = new G4Material("Bismuth", 9.78 *g/cm3, 1);
	Bismuth->AddElement(Bi,1);
	
	//Aluminum
	G4Material* Aluminum = new G4Material("Al", 2.7 *g/cm3, 1);
	Aluminum->AddElement(Al,1);

	//Carbon(graphite)
	G4Material* Carbon = new G4Material("Carbon", 2.3 *g/cm3, 1);
	Carbon->AddElement(C,1);
       
	//Low_Carbon_Stainless_Steel
	G4Material* Low_Carbon_Steel = new G4Material("Low_Carbon_Steel", 7.8499*g/cm3, 2);
  	Low_Carbon_Steel->AddMaterial(Iron,99.75*perCent);
  	Low_Carbon_Steel->AddMaterial(Carbon,0.25*perCent);
	 //JPET PLASTIC
	 G4Material* BC408 = man->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
	//PlasticScintillator
	G4Material* BC418 = new G4Material("PlasticScintillator", 1.032*g/cm3, 2);
   	BC418->AddElement(H, 30);
  	BC418->AddElement(C, 27);

  fWindowMaterial = BC408;//EJ309B5;
  fDefaultMaterial = Galactic;//AIR_mat;
  fTargetMaterial = BC418;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::MaterialWithSingleIsotope( G4String name,
                           G4String symbol, G4double density, G4int Z, G4int A)
{
 // define a material from an isotope
 //
 G4int ncomponents;
 G4double abundance, massfraction;

 G4Isotope* isotope = new G4Isotope(symbol, Z, A);
 
 G4Element* element  = new G4Element(name, symbol, ncomponents=1);
 element->AddIsotope(isotope, abundance= 100.*perCent);
 
 G4Material* material = new G4Material(name, density, ncomponents=1);
 material->AddElement(element, massfraction=100.*perCent);

 return material;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ComputeParameters()
{
  // Compute total thickness of absorbers
  fAbsorSizeY = 0.;
  for (G4int iAbs = 1; iAbs < nb_plastic; iAbs++) {
    fAbsorSizeY += fAbsorThickness[iAbs];
  }
  fScintSizeY = 0.;
  for (G4int jScint = 1; jScint < nb_cryst; jScint++) {
    fScintSizeY += fScintThickness[jScint];
  }
  fWorldSizeXY = 11*m;
  fWorldSizeZ = 11*m;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  // complete the Calor parameters definition
  ComputeParameters();

  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();



 // rotation
	G4RotationMatrix* rotation = new G4RotationMatrix();
	 	rotation->rotateX(180*deg);
	G4RotationMatrix* rotation1 = new G4RotationMatrix();
     		rotation1->rotateY(90*deg);
	G4RotationMatrix* rotation2 = new G4RotationMatrix();
     		rotation2->rotateX(-90*deg);
	 	G4RotationMatrix* rotation3 = new G4RotationMatrix();
   	 	rotation3->rotateX(42*deg);
	 	rotation3->rotateY(90*deg);
	 
	
		
  //
  // World
  //
  G4Box* solidWorld =
    new G4Box("World",                                             //name
               5*cm+fWorldSizeXY/2,fWorldSizeXY/2,fWorldSizeZ/2);       //size

  G4LogicalVolume* logicWorld =
    new G4LogicalVolume(solidWorld,              //solid
                        fDefaultMaterial,        //material
                        "World");                //name

  fPhysiWorld = 
    new G4PVPlacement(0,                        //no rotation
                      G4ThreeVector(),          //at (0,0,0)
                      logicWorld,               //logical volume
                      "World",                  //name
                       0,                       //mother volume
                       false,                   //no boolean operation
                       0);                      //copy number

/*===================================================================*/
//				Veto Detector			     //
/*===================================================================*/
  G4double veto_dX = 2.1*m, veto_dY = 2.1*m, veto_dZ = 1*m;
  G4double Plastics_dX = 100*cm, Plastics_dY = 2*cm, Plastics_dZ = 5*cm;

  G4Box* solidVeto = new G4Box("Veto", veto_dX/2, veto_dY/2, veto_dZ/2);
                     
  G4LogicalVolume* logicVeto = 
    new G4LogicalVolume(solidVeto,          //its solid
                        fDefaultMaterial,   //its material
                        "VetoLV");          //its name
  for (G4int iplastic = 1; iplastic <= nb_plastic; iplastic++) {

 G4Box* solidPlastic = new G4Box("solidPlastic", Plastics_dX/2, Plastics_dY/2, Plastics_dZ/2);
    G4double posX = iplastic*(Plastics_dZ);
    G4RotationMatrix rotmat  = G4RotationMatrix();
	rotmat.rotateY(90*deg); 
	rotmat.rotateZ(90*deg); 
  G4LogicalVolume* logicAbsor= 
    new G4LogicalVolume(solidPlastic,          //its solid
                        fTargetMaterial,          //its material
                        "solidPlastic");


  G4double sizex1 = -veto_dX/2-20*mm, sizex2 = veto_dX/2, sizey1 = veto_dY/2, sizey2 = -veto_dY/2; 
  G4double y01 = sizex1+posX,y02 = sizex2-posX;
  G4double z01 = 1*cm, z02 = -1*cm;
  G4double x01 = sizey1-Plastics_dY/2,x02 = sizey2+Plastics_dY/2; 
    	G4ThreeVector posUp = G4ThreeVector(x01,y01,z01);     
    	G4Transform3D transform1 = G4Transform3D(rotmat,posUp);
    
	//G4ThreeVector posDown = G4ThreeVector(x02,y02,z02);     
    	//G4Transform3D transform2 = G4Transform3D(rotmat,posDown);
                                    
    new G4PVPlacement(transform1,            //rotation,position
                      logicAbsor,	//logicUpPlastics,         //its logical volume
                      "solidPlastic",             //its name
                      logicVeto,             //its mother  volume
                      true,                 //boolean operation
                      iplastic              //copy number
                      );      

  }

new G4PVPlacement(rotation,			   //rotation
		  G4ThreeVector(0.,0.,0.), //position
                  logicVeto,               //its logical volume
                      "Veto",              //its name
                      logicWorld,          //its mother  volume
                      true,               //no boolean operation
                      0);                //copy number) 
  //
  // JPET rings
  //
G4double dPhi = twopi/96,dPhi_Mid = twopi/48; // dPhi is filling the ring (2*Pi divided into numbers of crystals)
							
// Inner and Outer Radiuses of the rings
  G4double ring_R1 = 56.55*cm,ring2_R1 = 46.65*cm,ring3_R1 = 42.4*cm;
  //G4double ring_R2 = 58.5*cm,ring2_R2 = 47.6*cm,ring3_R2 = 43.35*cm
/////////////////////////////////////////////////////////

		  for (G4int icrys = 1; icrys <= 192; ++icrys) {
			if(icrys<=96){
			   G4double dX = 19*mm, dY = 7*mm, dZ = 500*mm;
			  G4Box* solidCryst = new G4Box("Scintil", dX/2, dY/2, dZ/2);
					     
			  G4LogicalVolume* logicCryst = 
			    new G4LogicalVolume(solidCryst,          //its solid
						fWindowMaterial,     //its material
						"CrystalLV");        //its name     
			  // place crystals within a ring 
			  //
			  //G4String name = "scin_"+G4UIcommand::ConvertToString(icopy);
			    G4double phi = icrys*dPhi + 0.032724923;//0.032724923 rad = 1.875 degree moved from its previous position 
			    G4RotationMatrix rotm1  = G4RotationMatrix();
			    rotm1.rotateZ(phi);

			    G4ThreeVector uz = G4ThreeVector(std::cos(phi),std::sin(phi),0.);     
			    G4ThreeVector position = (ring_R1)*uz;
			    G4Transform3D transform = G4Transform3D(rotm1,position);
						            
			    new G4PVPlacement(transform,             //rotation,position
					      logicCryst,            //its logical volume
					      "Scintil",             //its name
					      logicVeto,//logicRing,             //its mother  volume
					      true,                 //no boolean operation
					      icrys);             //copy number
				

			  G4VisAttributes* VisAttr = 0;
			  VisAttr = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
			  VisAttr->SetVisibility(true);
			  logicCryst->SetVisAttributes(VisAttr); 
	
			} if ((icrys > 96)&&(icrys<=144)){
			   G4double dX = 19*mm, dY = 7*mm, dZ = 500*mm;
			  G4Box* solidCryst = new G4Box("Scintil", dX/2, dY/2, dZ/2);
					     
			  G4LogicalVolume* logicCryst = 
			    new G4LogicalVolume(solidCryst,          //its solid
						fWindowMaterial,     //its material
						"CrystalLV");        //its name
				       
			  // place crystals within a ring 
			  //
			 G4double phi = icrys*dPhi_Mid + 0.06544985;//0.06544985 rad = 3.75 degree moved from its previous position 
			    G4RotationMatrix rotm2  = G4RotationMatrix();
			    rotm2.rotateZ(phi);

			    G4ThreeVector uz = G4ThreeVector(std::cos(phi),std::sin(phi),0.);   
			    G4ThreeVector position = (ring2_R1)*uz;
			    G4Transform3D transform = G4Transform3D(rotm2,position);
						            
			    new G4PVPlacement(transform,             //rotation,position
					      logicCryst,            //its logical volume
					      "Scintil",             //its name
					      logicVeto,//logicRing,             //its mother  volume
					      true,                 //no boolean operation
					      icrys);             //copy number
				
			  G4VisAttributes* VisAttr = 0;
			  VisAttr = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
			  VisAttr->SetVisibility(true);
			  logicCryst->SetVisAttributes(VisAttr); 
		
			}
			if(icrys > 144)
			{
			  G4double dX = 19*mm, dY = 7*mm, dZ = 500*mm;
			  G4Box* solidCryst = new G4Box("Scintil", dX/2, dY/2, dZ/2);   
			  G4LogicalVolume* logicCryst = 
			    new G4LogicalVolume(solidCryst,          //its solid
						fWindowMaterial,     //its material
						"CrystalLV");        //its name 
			  // place crystals within a ring 
			  //
			  //G4String name = "scin_"+G4UIcommand::ConvertToString(icopy);
			    G4double phi = icrys*dPhi_Mid;
			    G4RotationMatrix rotm3  = G4RotationMatrix();
			    rotm3.rotateZ(phi);
			    G4ThreeVector uz = G4ThreeVector(std::cos(phi),std::sin(phi),0.);  
			    G4ThreeVector position3 = (ring3_R1)*uz;
			    G4Transform3D transform = G4Transform3D(rotm3,position3);
						            
		       	    new G4PVPlacement(transform,             //rotation,position
					      logicCryst,            //its logical volume
					      "Scintil",             //its name
					      logicVeto, //logicRing,             //its mother  volume
					      true,                 //no boolean operation
					      icrys);             //copy number

			  G4VisAttributes* VisAttr = 0;
			  VisAttr = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
			  VisAttr->SetVisibility(true);
			  logicCryst->SetVisAttributes(VisAttr); 

			}

  
		}
		   /* G4RotationMatrix rotRing  = G4RotationMatrix();
		    rotRing.rotateX(180*deg); 
		    G4ThreeVector posRing = G4ThreeVector(0.,0.,0.);
		    G4Transform3D transformRing = G4Transform3D(rotRing,posRing);
			new G4PVPlacement(transformRing,         //rotation,position
				      logicRing,         //its logical volume
				      "Ring3",            //its name
				       logicWorld, //its mother  volume
				      false,             //no boolean operation
				      0);                //copy number
*/
  G4VisAttributes* VisAtt = 0;
  VisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  VisAtt->SetVisibility(true);
  logicWorld->SetVisAttributes(VisAtt);
  PrintParameters();

  //always return the physical World
  //
  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4cout << "\n-------------------------------------------------------------"
         << "\n ---> The Absorber is " << nb_plastic << " layers of:";
  for (G4int i=1; i <= nb_plastic; i++)
     {
      G4cout << "\n \t" << std::setw(12) << fAbsorMaterial[i]->GetName() <<": "
              << std::setw(6) << G4BestUnit(fAbsorThickness[i],"Length");
     }
  G4cout << "\n-------------------------------------------------------------\n"
         << G4endl;
  G4cout << "\n-------------------------------------------------------------"
         << "\n ---> The Scintilator is " << nb_cryst << " layers of:";
  for (G4int j=1; j <= nb_cryst; j++)
     {
      G4cout << "\n \t" << std::setw(12) << fScintMaterial[j]->GetName() <<": "
              << std::setw(6) << G4BestUnit(fScintThickness[j],"Length");
     }
  G4cout << "\n-------------------------------------------------------------\n"
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void DetectorConstruction::SetNbOfAbsor(G4int ival)
{
  // set the number of Absorbers
  //
  if (ival < 1 || ival > (kMaxAbsor-1))
    { G4cout << "\n ---> warning from SetfNbOfAbsor: "
             << ival << " must be at least 1 and and most " << kMaxAbsor-1
             << ". Command refused" << G4endl;
      return;
    }
  nb_plastic = ival;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
void DetectorConstruction::SetNbOfScint(G4int jval)
{
  // set the number of Absorbers
  //
  if (jval < 1 || jval > (mMaxScint-1))
    { G4cout << "\n ---> warning from SetfNbOfScint: "
             << jval << " must be at least 1 and and most " << mMaxScint-1
             << ". Command refused" << G4endl;
      return;
    }
  nb_cryst = jval;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorMaterial(G4int iabs,const G4String& material1)
{
  // search the material by its name
  //
  if (iabs > nb_plastic || iabs <= 0)
    { G4cout << "\n --->warning from SetfAbsorMaterial: absor number "
             << iabs << " out of range. Command refused" << G4endl;
      return;
    }

  G4Material* pttoMaterial1 = 
    G4NistManager::Instance()->FindOrBuildMaterial(material1);
  if (pttoMaterial1) {
      fAbsorMaterial[iabs] = pttoMaterial1;
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetScintMaterial(G4int jabs,const G4String& material2)
{
  // search the material by its name
  //
  if (jabs > nb_cryst || jabs <= 0)
    { G4cout << "\n --->warning from SetfScintMaterial: Scint number "
             << jabs << " out of range. Command refused" << G4endl;
      return;
    }

  G4Material* pttoMaterial2 = 
    G4NistManager::Instance()->FindOrBuildMaterial(material2);
  if (pttoMaterial2) {
      fScintMaterial[jabs] = pttoMaterial2;
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}
void DetectorConstruction::SetAbsorThickness(G4int iabs,G4double val)
{
  // change Absorber thickness
  //
  if (iabs > nb_plastic || iabs <= 0)
    { G4cout << "\n --->warning from SetfAbsorThickness: absor number "
             << iabs << " out of range. Command refused" << G4endl;
      return;
    }
  if (val <= DBL_MIN)
    { G4cout << "\n --->warning from SetfAbsorThickness: thickness "
             << val  << " out of range. Command refused" << G4endl;
      return;
    }
  fAbsorThickness[iabs] = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
void DetectorConstruction::SetScintThickness(G4int jabs,G4double jval)
{
  // change Scintilator thickness
  //
  if (jabs > nb_cryst || jabs <= 0)
    { G4cout << "\n --->warning from SetfScintThickness: Scint number "
             << jabs << " out of range. Command refused" << G4endl;
      return;
    }
  if (jval <= DBL_MIN)
    { G4cout << "\n --->warning from SetfScintThickness: thickness "
             << jval  << " out of range. Command refused" << G4endl;
      return;
    }
  fScintThickness[jabs] = jval;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorSizeX(G4double val)
{
  // change the transverse size
  //
  if (val <= DBL_MIN)
    { G4cout << "\n --->warning from SetfAbsorSizeX: thickness "
             << val  << " out of range. Command refused" << G4endl;
      return;
    }
  fAbsorSizeX = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
void DetectorConstruction::SetAbsorSizeZ(G4double val)
{
  // change the transverse size
  //
  if (val <= DBL_MIN)
    { G4cout << "\n --->warning from SetfAbsorSizeZ: thickness "
             << val  << " out of range. Command refused" << G4endl;
      return;
    }
  fAbsorSizeZ = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
void DetectorConstruction::SetScintSizeX(G4double jval)
{
  // change the transverse size
  //
  if (jval <= DBL_MIN)
    { G4cout << "\n --->warning from SetfScintSizeX: thickness "
             << jval  << " out of range. Command refused" << G4endl;
      return;
    }
  fScintSizeX = jval;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
void DetectorConstruction::SetScintSizeZ(G4double jval)
{
  // change the transverse size
  //
  if (jval <= DBL_MIN)
    { G4cout << "\n --->warning from SetfScintSizeZ: thickness "
             << jval  << " out of range. Command refused" << G4endl;
      return;
    }
  fScintSizeZ = jval;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
    if ( fFieldMessenger.Get() == 0 ) {
        // Create global magnetic field messenger.
        // Uniform magnetic field is then created automatically if
        // the field value is not zero.
        G4ThreeVector fieldValue = G4ThreeVector();
        G4GlobalMagFieldMessenger* msg =
        new G4GlobalMagFieldMessenger(fieldValue);
        //msg->SetVerboseLevel(1);
        G4AutoDelete::Register(msg);
        fFieldMessenger.Put( msg );
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


