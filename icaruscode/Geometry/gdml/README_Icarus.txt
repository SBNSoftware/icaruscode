How to generate ICARUS gdml



To create the icarus.gdml file:

# ./generate_icarus.pl | ./make_gdml.pl -o icarus.gdml



To visulize the geometry using ROOT:

# root
# gSystem->Load("libGeom");
# gSystem->Load("libGdml");
# TGeoManager *geom = TGeoManager::Import("icarus.gdml");
# TGeoVolume *top = gGeoManager->GetTopVolume();
# top->Draw();

Macro used to visualize the geometry: icarus_geo.C


