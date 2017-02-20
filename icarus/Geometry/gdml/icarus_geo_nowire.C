typedef struct _drawopt {
  const char* volume;
  int         color;
} drawopt;

void icarus_geo_nowire(TString volName="volCryostat"){

gSystem->Load("libGeom");
gSystem->Load("libGdml");

//TGeoManager::Import("icarus_nowire.gdml");
TGeoManager::Import("icarus_nowireModif.gdml");

drawopt optIcarus[] = {
   {"volWorld",     kYellow},
   {"volDetEnclosure", kBlue},
   {"volCryostat", kOrange},
   {"volCathode", kRed},
   {"volTPC",	kGreen},
   {"volTPCPlaneU", kAzure},
   {"volTPCPlaneV", kAzure},
   {"volTPCPlaneX", kAzure},
   {"volSteelShell", kGray},
//   {"volTPCWirePlaneWidthSide", kRed},
//  {"volTPCWireU1124", kRed},
//  {"volTPCWireV1124", kBlue},
//  {"volTPCWireVert", kOrange},
{"volTPCActive",kGreen-6},
  {0, 0}
};

for (int i=0;; ++i) {
  if (optIcarus[i].volume==0) break;
    gGeoManager->FindVolumeFast(optIcarus[i].volume)->SetLineColor(optIcarus[i].color);
    //gGeoManager->FindVolumeFast(optIcarus[i].volume)->SetVisibility(kTRUE);
}

   /* non porduce risultati
    gGeoManager->FindVolumeFast("volTPCPlaneU")->SetFillColor(kAzure);
    gGeoManager->FindVolumeFast("volTPCPlaneV")->SetFillColor(kAzure);
    gGeoManager->FindVolumeFast("volTPCPlaneX")->SetFillColor(kAzure);
    gGeoManager->FindVolumeFast("volCathode")->SetFillColor(kRed); */

TList* mat = gGeoManager->GetListOfMaterials();
TIter next(mat);
TObject *obj;
while (obj == next()) {
 obj->Print();
}

 gGeoManager->CheckOverlaps(0.01);
// gGeoManager->PrintOverlaps();
 gGeoManager->SetMaxVisNodes(70000);

 gGeoManager->GetTopVolume()->Draw();
 //gGeoManager->FindVolumeFast(volName)->Draw();

 TFile *tf = new TFile("icarus_nowireModif.root", "RECREATE");
 gGeoManager->Write();
 tf->Close();
}
