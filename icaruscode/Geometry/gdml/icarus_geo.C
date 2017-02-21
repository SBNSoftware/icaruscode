typedef struct _drawopt {
  const char* volume;
  int         color;
} drawopt;

//void icarus_geo(TString volName="volCryostat"){
void icarus_geo(TString volName="volWorld"){

gSystem->Load("libGeom");
gSystem->Load("libGdml");

TGeoManager::Import("double.gdml");

drawopt optIcarus[] = {
   {"volWorld",     0},
   {"volDetEnclosure", kBlue},
   {"volCryostat", kOrange},
   {"volCathode", kRed},
   {"volTPC",	kGreen},
 //  {"volPMT",	kAzure},
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
//   {"volTPCWireVert", kBlack},
  // {"volTPCWireX", kAzure},
   {"volSteelShell", kGray},
   {"volTPCPlaneU", kMagenta},
   {"volTPCPlaneV", kMagenta-4},
   {"volTPCPlaneX", kMagenta-7},
//  {"volTPCWireU1124", kRed},
//  {"volTPCWireV1124", kBlue},
//  {"volTPCWireVert", kOrange},
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

// gGeoManager->CheckOverlaps(0.01);
// gGeoManager->PrintOverlaps();
// gGeoManager->SetMaxVisNodes(70000);

 gGeoManager->GetTopVolume()->Draw();
 //gGeoManager->FindVolumeFast(volName)->Draw();

 TFile *tf = new TFile("double.root", "RECREATE");
 gGeoManager->Write();
 tf->Close();
}
