typedef struct _drawopt {
  const char* volume;
  int         color;
} drawopt;

void icarus_geo(TString volName="volWorld"){

gSystem->Load("libGeom");
gSystem->Load("libGdml");

//TGeoManager::Import("icarus_complete_test_nowires_15jan2019.gdml"); //gdml da visualizzare
TGeoManager::Import("icarus_complete_test_no_ovbd_9jul2020.gdml"); //gdml da visualizzare
//TGeoManager::Import("icarus_complete_test_crt_03aug_nowires.gdml"); //gdml da visualizzare
//TGeoManager::Import("icarus_complete_test_24sep2019_nowires.gdml"); //gdml da visualizzare

//cambio colori per i volumi
drawopt optIcarus[] = {
   {"volWorld",     kAzure},
   {"volDetEnclosure", kBlue},
   {"volCryostat", kOrange},
   {"volCathode", kRed},
   {"volTPC0",	kYellow},
   {"volPMTPlane",	kRed},
   {"volGaseousArgon",	kRed},
   //{"volMech_Structure",  kGray},
   {"volThermIns", kBlue},
   {"volWarmVessel", kRed},
//pmtbig elements   
//   {"volPMTBig",	kCyan-10},
//   {"vol_PMT_TPBCoating",	kMagenta+4},
//   {"vol_PMT_AcrylicPlate",	kBlue-9},
 //  {"vol_PMT_Stalk",	kBlue},
 //  {"vol_PMT_SteelBase",	kGreen},
  // {"vol_PMT_Underside",	kGreen+3},
  // {"volOpDetSensitive",	kOrange},
//fine pmt big elements
//   {"volTPCWireUCommon", kAzure},
//  {"volTPCWireVCommon", kPink},
 //  {"volTPCWireYCommon", kBlack},
 //  {"volTPCWireW", kAzure},
//   {"volSteelShell", kGray},
//   {"volTPCPlaneU", kMagenta},
//   {"volTPCPlaneV", kMagenta-4},
//   {"volTPCPlaneY", kMagenta-7},
//  {"volTPCWireU1124", kRed},
//  {"volTPCWireV1124", kBlue},
{"volTPCActive",kGreen-6},
   {0, 0}
};

for (int i=0;; ++i) {
  if (optIcarus[i].volume==0) break;
    gGeoManager->FindVolumeFast(optIcarus[i].volume)->SetLineColor(optIcarus[i].color);
}

TList* mat = gGeoManager->GetListOfMaterials();
TIter next(mat);
TObject *obj;
while (obj == next()) {
 obj->Print();
}

//Check Overlaps
gGeoManager->CheckOverlaps(0.01);
gGeoManager->PrintOverlaps();
gGeoManager->SetMaxVisNodes(70000);

gGeoManager->GetTopVolume()->Draw();
 //gGeoManager->FindVolumeFast(volName)->Draw();

//Save in a root file
 TFile *tf = new TFile("icarus_geometry.root", "RECREATE");
 gGeoManager->Write();
 tf->Close();
}
