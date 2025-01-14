/**
 * @file   icaruscode/CRT/CRTUtils/CRTMatchingUtils.cxx
 * @author Francesco Poppi (poppi@bo.infn.it)
 * @date   January 2025
 */

#include "icaruscode/CRT/CRTUtils/CRTMatchingUtils.h"

#include <string>
#include <limits>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include <cetlib/search_path.h>
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "TMatrixD.h"
#include "TMatrixDEigen.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"

namespace icarus::crt{


    TopCRTCentersMap LoadTopCRTCenters()
    {
        TopCRTCentersMap TopCRTCenters;

        //  The follwing numbers have been extracted from the Top CRT modules geometry.
        //  The numbers are respectively moduleID and X,Y,Z coordinates of the
        //  module center (X markes the spot in the sketch below).
        //  The CRT modules are perfect square with 8 bars per side.
        //  The fixed coordinate (e.g. Y for Top CRT Horizontal) is returned from geometry.
        //  The transverce coordinate is returned from the average position of P1 and P2
        //  the average position of P1 and P3.
        //  The transformed CRT Hits are in cm.
        //
        //          ---------------------------------------------------------
        //          |      |      |      |      |      |      |      |      |       
        //          |      |      |      |      |      |      |      |      |       
        //          |      |      |      |      |      |      |      |      |       
        //          ---------------------------------------------------------
        //          |      |      |      |      |      |      |      |      |       
        //          |      |      |      |      |      |      |      |      |       
        //          |      |      |      |      |      |      |      |      |       
        //          ---------------------------------------------------------
        //          |      |      |      |      |      |      |      |      |       
        //          |      |      |      |      |      |      |      |      |       
        //          |      |      |      |      |      |      |      |      |       
        //          ---------------------------------------------------------
        //          |      |      |      |      |      |      |      |      |       
        //          |      |      |      |  P1  |  P2  |      |      |      |       
        //          |      |      |      |      |      |      |      |      |       
        //          --------------------------- X ---------------------------
        //          |      |      |      |      |      |      |      |      |       
        //          |      |      |      |  P3  |      |      |      |      |       
        //          |      |      |      |      |      |      |      |      |       
        //          ---------------------------------------------------------
        //          |      |      |      |      |      |      |      |      |       
        //          |      |      |      |      |      |      |      |      |       
        //          |      |      |      |      |      |      |      |      |       
        //          ---------------------------------------------------------
        //          |      |      |      |      |      |      |      |      |       
        //          |      |      |      |      |      |      |      |      |       
        //          |      |      |      |      |      |      |      |      |       
        //          ---------------------------------------------------------
        //          |      |      |      |      |      |      |      |      |       
        //          |      |      |      |      |      |      |      |      |       
        //          |      |      |      |      |      |      |      |      |       
        //          ---------------------------------------------------------
        //
        TopCRTCenters = {{108, { -460.975, 617.388,-1050.61}},
                        {109, { -460.975,617.388,-866.215}},
                        {110, { -460.975,617.388,-681.825}},
                        {111, { -460.975,617.388,-497.435}},
                        {112, { -460.975,617.388,-313.045}},
                        {113, { -460.975,617.388,-128.655}},
                        {114, { -460.975,617.388,55.735}},
                        {115, { -460.975,617.388,240.125}},
                        {116, { -460.975,617.388,424.515}},
                        {117, { -460.975,617.388,608.905}},
                        {118, { -460.975,617.388,793.295}},
                        {119, { -460.975,617.388,977.685}},
                        {120, { -460.975,617.388,1162.07}},
                        {121, { -460.975,617.388,1346.46}},
                        {122, { -276.585,617.388,-1050.61}},
                        {123, { -276.585,617.388,-866.215}},
                        {124, { -276.585,617.388,-681.825}},
                        {125, { -276.585,617.388,-497.435}},
                        {126, { -276.585,617.388,-313.045}},
                        {127, { -276.585,617.388,-128.655}},
                        {128, { -276.585,617.388, 55.735}},
                        {129, { -276.585,617.388, 240.125}},
                        {130, { -276.585,617.388, 424.515}},
                        {131, { -276.585,617.388, 608.905}},
                        {132, { -276.585,617.388, 793.295}},
                        {133, { -276.585,617.388, 977.685}},
                        {134, { -276.585,617.388, 1162.07}},
                        {135, { -276.585,617.388, 1346.46}},
                        {136, { -92.195 ,617.388,-1050.61}},
                        {137, { -92.195 ,617.388,-866.215}},
                        {138, { -92.195 ,617.388,-681.825}},
                        {139, { -92.195 ,617.388,-497.435}},
                        {140, { -92.195 ,617.388,-313.045}},
                        {141, { -92.195 ,617.388,-128.655}},
                        {142, { -92.195 ,617.388, 55.735}},
                        {143, { -92.195 ,617.388, 240.125}},
                        {144, { -92.195 ,617.388, 424.515}},
                        {145, { -92.195 ,617.388, 608.905}},
                        {146, { -92.195 ,617.388, 793.295}},
                        {147, { -92.195 ,617.388, 977.685}},
                        {148, { -92.195 ,617.388, 1162.07}},
                        {149, { -92.195 ,617.388, 1346.46}},
                        {150, { 92.195 , 617.388, -1050.61}},
                        {151, { 92.195 , 617.388, -866.215}},
                        {152, { 92.195 , 617.388, -681.825}},
                        {153, { 92.195 , 617.388, -497.435}},
                        {154, { 92.195 , 617.388, -313.045}},
                        {155, { 92.195 , 617.388, -128.655}},
                        {156, { 92.195 ,617.388, 55.735}},
                        {157, { 92.195 ,617.388, 240.125}},
                        {158, { 92.195 ,617.388, 424.515}},
                        {159, { 92.195 ,617.388, 608.905}},
                        {160, { 92.195 ,617.388, 793.295}},
                        {161, { 92.195 ,617.388, 977.685}},
                        {162, { 92.195 ,617.388, 1162.07}},
                        {163, { 92.195 ,617.388, 1346.46}},
                        {164, { 276.585,617.388, -1050.61}},
                        {165, { 276.585,617.388, -866.215}},
                        {166, { 276.585,617.388, -681.825}},
                        {167, { 276.585,617.388, -497.435}},
                        {168, { 276.585,617.388, -313.045}},
                        {169, { 276.585,617.388, -128.655}},
                        {170, { 276.585,617.388, 55.735}},
                        {171, { 276.585,617.388, 240.125}},
                        {172, { 276.585,617.388, 424.515}},
                        {173, { 276.585,617.388, 608.905}},
                        {174, { 276.585,617.388, 793.295}},
                        {175, { 276.585,617.388, 977.685}},
                        {176, { 276.585,617.388, 1162.07}},
                        {177, { 276.585,617.388, 1346.46}},
                        {178, { 460.975,617.388, -1050.61}},
                        {179, { 460.975,617.388, -866.215}},
                        {180, { 460.975,617.388, -681.825}},
                        {181, { 460.975,617.388, -497.435}},
                        {182, { 460.975,617.388, -313.045}},
                        {183, { 460.975,617.388, -128.655}},
                        {184, { 460.975,617.388, 55.735}},
                        {185, { 460.975,617.388, 240.125}},
                        {186, { 460.975,617.388, 424.515}},
                        {187, { 460.975,617.388, 608.905}},
                        {188, { 460.975,617.388, 793.295}},
                        {189, { 460.975,617.388, 977.685}},
                        {190, { 460.975,617.388, 1162.07}},
                        {191, { 460.975,617.388, 1346.46}},
                        {192, { 555.265, 496.038, -1050.61}},
                        {193, { 555.265, 496.038, -866.215}},
                        {194, { 555.265, 496.038, -681.825}},
                        {195, { 555.265, 496.038, -497.435}},
                        {196, { 555.265, 496.038, -313.045}},
                        {197, { 555.265, 496.038, -128.655}},
                        {198, { 555.265, 496.038, 55.735}},
                        {199, { 555.265, 496.038, 240.125}},
                        {200, { 555.265, 496.038, 424.515}},
                        {201, { 555.265, 496.038, 608.905}},
                        {202, { 555.265, 496.038, 793.295}},
                        {203, { 555.265, 496.038, 977.685}},
                        {204, { 555.265, 496.038, 1162.07}},
                        {205, { 555.265, 496.038, 1346.46}},
                        {206, { -555.265, 496.038, -1050.61}},
                        {207, { -555.265, 496.038, -866.215}},
                        {208, { -555.265, 496.038, -681.825}},
                        {209, { -555.265, 496.038, -497.435}},
                        {210, { -555.265, 496.038, -313.045}},
                        {211, { -555.265, 496.038, -128.655}},
                        {212, { -555.265, 496.038, 55.735}},
                        {213, { -555.265, 496.038, 240.125}},
                        {214, { -555.265, 496.038, 424.515}},
                        {215, { -555.265, 496.038, 608.905}},
                        {216, { -555.265, 496.038, 793.295}},
                        {217, { -555.265, 496.038, 977.685}},
                        {218, { -555.265, 496.038, 1162.07}},
                        {219, { -555.265, 496.038, 1346.46}},
                        {220, { -460.975, 496.038, -1143.4}},
                        {221, { -276.585, 496.038, -1143.4}},
                        {222, { -92.195, 496.038, -1143.4}},
                        {223, { 92.195, 496.038, -1143.4}},
                        {224, { 276.585, 496.038, -1143.4}},
                        // This module does not exist in reality, but exists in simulation
                        {225, { 0, 0, 0}},
                        {226, { -460.975, 525.038, 1533.608}},
                        {227, { -276.585, 525.038, 1533.608}},
                        {228, { -92.195, 525.038, 1533.608}},
                        {229, { 92.195, 525.038, 1533.608}},
                        {230, { 276.585, 525.038, 1533.608}},
                        {231, { 460.975, 525.038, 1533.608}}};
        return TopCRTCenters;
    }

    TransformedCRTHit AffineTransformation(double DX, double DZ,AffineTrans affine)
    {
        double CRTX=affine.B1+DZ*affine.A12+DX*affine.A11;
        double CRTZ=affine.B2+DZ*affine.A22+DX*affine.A21;
        return std::make_pair(CRTX, CRTZ);
    }

    TopCRTTransformations LoadTopCRTTransformations()
    {
        std::string fullFileName;
        cet::search_path searchPath("FW_SEARCH_PATH");
        searchPath.find_file("TopCrtCorrectionByTPC.txt", fullFileName);
        /*if (!searchPath.find_file("TopCrtCorrectionByTPC.txt", fullFileName)) {
            mf::LogError("CRTMatchingUtils_LoadTopCrtTransformation")
            << "Top CRT Correction transformation file not found in FW_SEARCH_PATH.";
        }*/
        std::ifstream corrFile(fullFileName, std::ios::in);

        std::map<int, AffineTrans> Type0Corr;
        std::map<int, AffineTrans> Type1Corr;
        std::map<int, AffineTrans> Type2Corr;
        std::map<int, AffineTrans> Type3Corr;
        std::map<int, AffineTrans> Type4Corr;
        std::map<int, AffineTrans> Type5Corr;

        if (!corrFile.is_open()) {
            mf::LogError("CRTMatchingUtils_LoadTopCRTTransformation")
            << "Failed to open Top CRT Correction transformation file: " << fullFileName;
        }
        std::string line;
        while (std::getline(corrFile, line)) {
            std::istringstream iss(line);
            std::vector<double> variables;
            double var;
            while (iss >> var) {
                variables.push_back(var);
            }

            int pos=0;

            AffineTrans Type0 = {variables.at(1+pos),variables.at(2+pos),variables.at(3+pos),variables.at(4+pos),variables.at(5+pos),variables.at(6+pos), variables.at(7+pos),variables.at(8+pos)};
            pos=pos+8;
            AffineTrans Type1 = {variables.at(1+pos),variables.at(2+pos),variables.at(3+pos),variables.at(4+pos),variables.at(5+pos),variables.at(6+pos), variables.at(7+pos),variables.at(8+pos)};
            pos=pos+8;
            AffineTrans Type2 = {variables.at(1+pos),variables.at(2+pos),variables.at(3+pos),variables.at(4+pos),variables.at(5+pos),variables.at(6+pos), variables.at(7+pos),variables.at(8+pos)};
            pos=pos+8;
            AffineTrans Type3 = {variables.at(1+pos),variables.at(2+pos),variables.at(3+pos),variables.at(4+pos),variables.at(5+pos),variables.at(6+pos), variables.at(7+pos),variables.at(8+pos)};
            pos=pos+8;
            AffineTrans Type4 = {variables.at(1+pos),variables.at(2+pos),variables.at(3+pos),variables.at(4+pos),variables.at(5+pos),variables.at(6+pos), variables.at(7+pos),variables.at(8+pos)};
            pos=pos+8;
            AffineTrans Type5 = {variables.at(1+pos),variables.at(2+pos),variables.at(3+pos),variables.at(4+pos),variables.at(5+pos),variables.at(6+pos), variables.at(7+pos),variables.at(8+pos)};
            Type0Corr.insert({variables.at(0), Type0});
            Type1Corr.insert({variables.at(0), Type1});
            Type2Corr.insert({variables.at(0), Type2});
            Type3Corr.insert({variables.at(0), Type3});
            Type4Corr.insert({variables.at(0), Type4});
            Type5Corr.insert({variables.at(0), Type5});
        }
        TopCRTTransformations LoadedTransformations={Type2Corr, Type1Corr, Type0Corr, Type4Corr, Type5Corr, Type3Corr};
        return LoadedTransformations;
    }

    CRTMatchingAlg::CRTMatchingAlg(const fhicl::ParameterSet& pset)
    {
        this->reconfigure(pset);
        return;
    }

    CRTMatchingAlg::CRTMatchingAlg() = default;

    void CRTMatchingAlg::reconfigure(const fhicl::ParameterSet& pset){
        fTickPeriod = pset.get<double>("TickPeriod", 0.4);
        fTickAtAnode = pset.get<double>("TickAtAnode", 850.);
        fAllowedOffsetCM = pset.get<double>("AllowedOffsetCM", 1.57);
        return;
    }

    Direction CRTMatchingAlg::PCAfit (std::vector<double> const& x, std::vector<double> const& y, std::vector<double> const& z)
    {
        int min=0;
        int max=x.size();
        int size=max-min;
        double xavg=0, yavg=0, zavg=0;
        for(int k=min; k<max; k++) {
            xavg+=x[k]; yavg+=y[k]; zavg+=z[k];
        }
        xavg=xavg/size;
        yavg=yavg/size;
        zavg=zavg/size;
        double x2=0, y2=0, z2=0, xiy=0, xiz=0, yiz=0;
        for(int k=min; k<max; k++){
            double xadj = x[k] - xavg;
            double yadj = y[k] - yavg;
            double zadj = z[k] - zavg;
            x2+=xadj*xadj;
            y2+=yadj*yadj;
            z2+=zadj*zadj;
            xiy+=xadj*yadj;
            xiz+=xadj*zadj;
            yiz+= yadj*zadj;
        }//end loop calculating covariance matrix elements
        x2=x2/size;
        y2=y2/size;
        z2=z2/size; 
        xiy=xiy/size;
        xiz=xiz/size;
        yiz=yiz/size;

        //Now covariance elements are ready, make covariance matrix and do PCA analysis
        TMatrixD covmat(3,3);
        covmat[0][0] = x2; covmat[1][1]=y2; covmat[2][2]=z2;
        covmat[0][1] = xiy; covmat[1][0] = xiy;
        covmat[0][2] = xiz; covmat[2][0] = xiz;
        covmat[2][1] = yiz; covmat[1][2] = yiz;
        TMatrixDEigen thiscov(covmat);
        TMatrixD this_eval = thiscov.GetEigenValues();
        TMatrixD this_evec = thiscov.GetEigenVectors();
        if(this_eval.GetNrows()!=3 || this_eval.GetNcols()!=3 || this_evec.GetNrows()!=3 || this_evec.GetNcols()!=3){ 
            throw std::logic_error("CRTMatchingAlg::PCAfit: evals/evects wrong size!");
        }//end if evals/evects aren't 3x3 matrices as expected
        double max_eval = std::numeric_limits<double>::lowest(); int maxevalpos = -1;
        for(int k = 0; k < 3; k++){
            if(this_eval[k][k]>max_eval){
                max_eval = this_eval[k][k];
                maxevalpos = k;
            }//end if(this_eval[k][k]>max_eval)
        }//end loop looking for best eval

        Direction thDirection= {this_evec[0][maxevalpos] , this_evec[1][maxevalpos], this_evec[2][maxevalpos], xavg, yavg, zavg};
        return thDirection;
    }

    CRTPlane CRTMatchingAlg::DeterminePlane(sbn::crt::CRTHit const& CRThit)
    {
        int Plane;
        double Pos;
        switch (CRThit.plane) {
          case 30:
            Plane=0;
            Pos=CRThit.y_pos;
            break;
          case 31: case 32:
          case 40: case 41: case 42: case 43: case 44: case 45:
            Plane=1;
            Pos=CRThit.x_pos;
            break;
          default:
            Plane=2;
            Pos=CRThit.z_pos;
            break;
        } // switch
        return std::make_pair(Plane, Pos);
    }
    
    ProjectionPoint CRTMatchingAlg::TranslatePointTo(double dirx, double diry, double dirz, double x0, double y0, double z0, double position)
    {
        if(dirx==0) throw std::invalid_argument("CRTMatchingAlg::TranslatePointTo: dirx is 0");
        double Lambda = (position-x0)/dirx;
        double PosAtY = (Lambda*diry+y0);
        double PosAtZ = (Lambda*dirz+z0);

        CrossPoint CrossingPoint = {position, PosAtY, PosAtZ};
        return CrossingPoint;
    }

    CrossPoint CRTMatchingAlg::CalculateForPlane(const Direction& dir, int plane, double position)
    {
        double dirX, dirY, dirZ;
        double meanX, meanY, meanZ;

        switch(plane) {
            case 0: // Fixed Coordinate: Y
                dirX = dir.diry; dirY = dir.dirx; dirZ = dir.dirz;
                meanX = dir.meany; meanY = dir.meanx; meanZ = dir.meanz;
                break;
            case 1: // Fixed Coordinate: X
                dirX = dir.dirx; dirY = dir.diry; dirZ = dir.dirz;
                meanX = dir.meanx; meanY = dir.meany; meanZ = dir.meanz;
                break;
            case 2: // Fixed Coordinate: Z
                dirX = dir.dirz; dirY = dir.diry; dirZ = dir.dirx;
                meanX = dir.meanz; meanY = dir.meany; meanZ = dir.meanx;
                break;
            default:
                throw std::invalid_argument("CRTMatchingAlg::CalculateForPlane(): invalid plane");
        }
        return CRTMatchingAlg::TranslatePointTo(dirX, dirY, dirZ, meanX, meanY, meanZ, position);
    }

    CrossPoint CRTMatchingAlg::DetermineProjection(const Direction& dir, CRTPlane plane)
    {
        CrossPoint thisCase = CRTMatchingAlg::CalculateForPlane(dir, plane.first, plane.second);

        switch(plane.first) {
            case 0: return {thisCase.Y, thisCase.X, thisCase.Z}; // Plane at Y e.g. Top CRT Horizontal Plane
            case 1: return {thisCase.X, thisCase.Y, thisCase.Z}; // Plane at X e.g. Side CRT West, East, Top CRT
            case 2: return {thisCase.Z, thisCase.Y, thisCase.X}; // Plane at Z e.g. Side CRT South, North, Top CRT 
            default:
                throw std::invalid_argument("CRTMatchingAlg::DetermineProjection(): invalid plane");
        }
    }

    TrackBarycenter CRTMatchingAlg::GetTrackBarycenter (std::vector<double> hx, std::vector<double> hy, std::vector<double> hz, std::vector<double> hw)
    {
        double average_x_charge=0, average_y_charge=0, average_z_charge=0, total_charge=0;
        bool isGood;
        for (unsigned i = 0; i < hz.size(); i++) {
            if(isnan(hz[i])) continue;

            if(hz[i]>-1000 && hz[i]<1000){
                average_z_charge+=hz[i]*hw[i];
                average_y_charge+=hy[i]*hw[i];
                average_x_charge+=hx[i]*hw[i];
                total_charge+=hw[i];
            }
        }
        average_z_charge=average_z_charge/total_charge;
        average_y_charge=average_y_charge/total_charge;
        average_x_charge=average_x_charge/total_charge;
        if(total_charge==0) isGood=false;
        else isGood=true;
        TrackBarycenter ThisTrackBary = {average_x_charge, average_y_charge, average_z_charge, isGood};
        return ThisTrackBary;
    }

    DriftedTrack CRTMatchingAlg::DriftTrack(const std::vector<art::Ptr<recob::Hit>>& trkHits, const std::vector<const recob::TrackHitMeta*>& trkHitMetas, const geo::GeometryCore *GeometryService, detinfo::DetectorPropertiesData const& detProp, double time, const recob::Track& tpcTrack) const
    {
        int outBound=0;
        int outBound1=0;
        int count=0;
        std::vector<double> recX, recY, recZ, recI;
        for(size_t i=0; i<trkHits.size(); i++){
            bool badhit = !trkHitMetas[i] || (trkHitMetas[i]->Index() == std::numeric_limits<unsigned int>::max()) ||
                        (!tpcTrack.HasValidPoint(trkHitMetas[i]->Index()));
            if(badhit) continue;
            geo::Point_t loc = tpcTrack.LocationAtPoint(trkHitMetas[i]->Index());
            //if(loc.X()==-999) continue;
            const geo::TPCGeo& tpcGeo = GeometryService->TPC(trkHits[i]->WireID());
            int const cryo = trkHits[i]->WireID().Cryostat;
            int tpc=trkHits[i]->WireID().TPC;
            double vDrift=detProp.DriftVelocity();
            // recoX: distance of the hit charge from the plane
            //  * `trkHits[i]->PeakTime()-fTickAtAnode`: TPC ticks from the trigger to when charge gets to plane
            //  * `time`: t0 [CRT or Flash] (with respect to trigger) in TPC ticks
            //  * difference is how many ticks passed from the track to when charge gets to plane: drift ticks
            double recoX=(trkHits[i]->PeakTime()-fTickAtAnode-time/fTickPeriod)*fTickPeriod*vDrift;

            double plane=tpcGeo.PlanePtr(trkHits[i]->WireID().Plane)->GetCenter().X();
            double cathode=tpcGeo.GetCathodeCenter().X();
            double X=plane-tpcGeo.DriftDir().X()*recoX;
            if(cryo==0 && (tpc==0 || tpc==1) && (X>(cathode+fAllowedOffsetCM)||X<(plane-fAllowedOffsetCM))) outBound++;
            else if(cryo==0 && (tpc==2 || tpc==3)&& (X<(cathode-fAllowedOffsetCM)||X>(plane+fAllowedOffsetCM))) outBound++;
            else if(cryo==1 && (tpc==0 || tpc==1)&& (X>(cathode+fAllowedOffsetCM)||X<(plane-fAllowedOffsetCM))) outBound++;
            else if(cryo==1 && (tpc==2 || tpc==3)&& (X<(cathode-fAllowedOffsetCM)||X>(plane+fAllowedOffsetCM))) outBound++;

            double fAllowedRel = 1+fAllowedOffsetCM / tpcGeo.DriftDistance();
            if(!tpcGeo.ActiveBoundingBox().ContainsX(X, fAllowedRel)) outBound1++;

            if(outBound!=outBound1 && count<=10){
                std::cout<<"Active Bounding Box MinX "<<tpcGeo.ActiveBoundingBox().MinX()<<" MaxX "<<tpcGeo.ActiveBoundingBox().MaxX()<<std::endl;
                std::cout<<"OutBound "<<outBound<<" outBound1 "<<outBound1<<" cryo "<<cryo<<" TPC "<<tpc<<" recoX "<<recoX<<" plane "<<plane<<" view "<<trkHits[i]->WireID().Plane<<" X "<<X<<" allowed cm "<<fAllowedOffsetCM<<" allowed % "<<fAllowedRel<<std::endl;
                count++;
            }
            recX.push_back(X);
            recY.push_back(loc.Y());
            recZ.push_back(loc.Z());
            recI.push_back(trkHits[i]->Integral());
        }
        DriftedTrack thisDriftedTrack = {recX, recY, recZ, recI, outBound};
        std::cout<<"--> Outbound cm "<<outBound<<" ; Outbound % "<<outBound1<<std::endl;
        return thisDriftedTrack;
    }

}
