#ifndef SBN_TrackTreeStoreObj
#define SBN_TrackTreeStoreObj

namespace sbn {

  struct selTrackInfo {
    int trackID;
    float t0;
    float start_x;
    float start_y;
    float start_z;
    float end_x;
    float end_y;
    float end_z;
    float length;
    float dir_x;
    float dir_y;
    float dir_z;
  selTrackInfo():
    trackID(-1),
      t0(-1),
      start_x(-1),
      start_y(-1),
      start_z(-1),
      end_x(-1),
      end_y(-1),
      end_z(-1),
      length(-1),
      dir_x(-1),
      dir_y(-1),
      dir_z(-1)
    {}
  };
}

#endif
