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
#include "TMatrixDSym.h"
#include "TMatrixDEigen.h"
#include "TMatrixDSymEigen.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
namespace icarus::crt{
    
    CRTMatchingAlg::CRTMatchingAlg(const fhicl::ParameterSet& pset)
    {
        this->reconfigure(pset);
        return;
    }

    CRTMatchingAlg::CRTMatchingAlg() = default;

    void CRTMatchingAlg::reconfigure(const fhicl::ParameterSet& pset)
    {
        fAllowedOffsetCM = pset.get<double>("AllowedOffsetCM", 1.57);
        return;
    }

    PCAResults CRTMatchingAlg::PCAfit (std::vector<geo::Point_t> const& sp)
    {
        int min=0;
        int max=sp.size();
        int size=max-min;
        double xavg=0, yavg=0, zavg=0;
        for(int k=min; k<max; k++) {
            xavg+=sp[k].X(); yavg+=sp[k].Y(); zavg+=sp[k].Z();
        }
        xavg=xavg/size; yavg=yavg/size; zavg=zavg/size;
        double x2=0, y2=0, z2=0, xiy=0, xiz=0, yiz=0;
        for(int k=min; k<max; k++){
            double xadj = sp[k].X() - xavg;
            double yadj = sp[k].Y() - yavg;
            double zadj = sp[k].Z() - zavg;
            x2+=xadj*xadj; y2+=yadj*yadj; z2+=zadj*zadj;
            xiy+=xadj*yadj; xiz+=xadj*zadj; yiz+= yadj*zadj;
        }//end loop calculating covariance matrix elements
        x2=x2/size; y2=y2/size; z2=z2/size; 
        xiy=xiy/size; xiz=xiz/size; yiz=yiz/size;
        //Now covariance elements are ready, make covariance matrix and do PCA analysis
        TMatrixDSym covmat(3);
        covmat(0, 0) = x2; covmat(1, 1) = y2; covmat(2, 2) = z2;
        covmat(0, 1) = xiy; covmat(1, 0) = xiy;
        covmat(0, 2) = xiz; covmat(2, 0) = xiz;
        covmat(1, 2) = yiz; covmat(2, 1) = yiz;
        TMatrixDSymEigen thiscov(covmat);
        TVectorD this_eval = thiscov.GetEigenValues();
        TMatrixD this_evec = thiscov.GetEigenVectors();

        if(this_eval.NonZeros()!=3) throw std::logic_error("CRTMatchingAlg::PCAfit: evals/evects wrong size!");
        geo::Point_t mean = {xavg, yavg, zavg};        
        geo::Vector_t first = {this_evec[0][0], this_evec[1][0], this_evec[2][0]};
        geo::Vector_t second = {this_evec[0][1], this_evec[1][1], this_evec[2][1]};
        geo::Vector_t third = {this_evec[0][2], this_evec[1][2], this_evec[2][2]};

        PCAResults thPCAResult = {first, second, third, this_eval[0], this_eval[1], this_eval[2], mean};
        return thPCAResult;
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

    CrossingPoint CRTMatchingAlg::TranslatePointTo(geo::Vector_t dirVector, geo::Point_t meanPos, CRTPlane CRTwall)
    {
        
        if(dirVector.X()==0) throw std::invalid_argument("CRTMatchingAlg::TranslatePointTo: Cosine Director is null");
        return meanPos + dirVector * (CRTwall.second - meanPos.X()) / dirVector.X();
    }

    TranslationVector CRTMatchingAlg::RotateToLocalCRTPlane(const TranslationVector& transl, CRTPlane CRTwall)
    {
        geo::Vector_t rotatedDir;
        geo::Point_t rotatedMean;

        switch(CRTwall.first) {
            case 0: // Fixed Coordinate: Y
                rotatedDir={transl.dir.Y(), transl.dir.X(), transl.dir.Z()};
                rotatedMean={transl.mean.Y(), transl.mean.X(), transl.mean.Z()};
                break;
            case 1: // Fixed Coordinate: X
                rotatedDir={transl.dir.X(), transl.dir.Y(), transl.dir.Z()};
                rotatedMean={transl.mean.X(), transl.mean.Y(), transl.mean.Z()};
                break;
            case 2: // Fixed Coordinate: Z
                rotatedDir={transl.dir.Z(), transl.dir.Y(), transl.dir.X()};
                rotatedMean={transl.mean.Z(), transl.mean.Y(), transl.mean.X()};
                break;
            default:
                throw std::invalid_argument("CRTMatchingAlg::RotateToLocalCRTPlane(): invalid plane");
        }
        return {rotatedDir, rotatedMean};
    }

    CrossingPoint CRTMatchingAlg::RotateFromLocalCRTPlane(CrossingPoint crossPointCRT, CRTPlane CRTWall)
    {
        switch(CRTWall.first) {
            case 0: return {crossPointCRT.Y(), crossPointCRT.X(), crossPointCRT.Z()}; // Plane at Y e.g. Top CRT Horizontal Plane
            case 1: return {crossPointCRT.X(), crossPointCRT.Y(), crossPointCRT.Z()}; // Plane at X e.g. Side CRT West, East, Top CRT
            case 2: return {crossPointCRT.Z(), crossPointCRT.Y(), crossPointCRT.X()}; // Plane at Z e.g. Side CRT South, North, Top CRT 
            default:
                throw std::invalid_argument("CRTMatchingAlg::RotateFromLocalCRTPlane(): invalid plane");
        }
    }

    CrossingPoint CRTMatchingAlg::DetermineProjection(const TranslationVector& dir, CRTPlane CRTWall)
    {
        TranslationVector directionPlaneCoordinate = CRTMatchingAlg::RotateToLocalCRTPlane(dir, CRTWall);
        CrossingPoint crossingPointPlaneCoordinate = CRTMatchingAlg::TranslatePointTo(directionPlaneCoordinate.dir, directionPlaneCoordinate.mean, CRTWall);
        return CRTMatchingAlg::RotateFromLocalCRTPlane(crossingPointPlaneCoordinate, CRTWall);
    }

    TrackBarycenter CRTMatchingAlg::GetTrackBarycenter (std::vector<geo::Point_t> sp, std::vector<double> hw)
    {
        double average_x_charge=0, average_y_charge=0, average_z_charge=0, total_charge=0;
        bool isGood;
        for (unsigned i = 0; i < sp.size(); i++) {
            if(isnan(sp[i].Z())) continue;

            if(sp[i].Z()>-1000 && sp[i].Z()<1000){
                average_z_charge+=sp[i].Z()*hw[i];
                average_y_charge+=sp[i].Y()*hw[i];
                average_x_charge+=sp[i].X()*hw[i];
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

    DriftedTrack CRTMatchingAlg::DriftTrack(const std::vector<art::Ptr<recob::Hit>>& trkHits, const std::vector<const recob::TrackHitMeta*>& trkHitMetas, const geo::GeometryCore *GeometryService, detinfo::DetectorPropertiesData const& detProp, detinfo::DetectorClocksData  const& detClock, double time, const recob::Track& tpcTrack, int maxOutBoundPoints) const
    {
        int outBound=0;
        std::vector<double> recI;
        std::vector<geo::Point_t> driftedPositionVector;
        for(size_t i=0; i<trkHits.size(); i++){
            bool badhit = !trkHitMetas[i] || (trkHitMetas[i]->Index() == std::numeric_limits<unsigned int>::max()) ||
                        (!tpcTrack.HasValidPoint(trkHitMetas[i]->Index()));
            if(badhit) continue;
            geo::Point_t loc = tpcTrack.LocationAtPoint(trkHitMetas[i]->Index());
            const geo::TPCGeo& tpcGeo = GeometryService->TPC(trkHits[i]->WireID());
            double vDrift=detProp.DriftVelocity();

            // recoX: distance of the hit charge from the plane
            //  * `trkHits[i]->PeakTime()-fTickAtAnode`: TPC ticks from the trigger to when charge gets to plane
            //  * `time`: t0 [CRT or Flash] (with respect to trigger) in TPC ticks
            //  * difference is how many ticks passed from the track to when charge gets to plane: drift ticks
            double recoX=(trkHits[i]->PeakTime()-detClock.Time2Tick(detClock.TriggerTime())-time/detClock.TPCClock().TickPeriod())*detClock.TPCClock().TickPeriod()*vDrift;
            double plane=tpcGeo.PlanePtr(trkHits[i]->WireID().Plane)->GetCenter().X();
            double X=plane-tpcGeo.DriftDir().X()*recoX;
            double AllowedRel = 1+fAllowedOffsetCM / abs(tpcGeo.GetCathodeCenter().X()-tpcGeo.PlanePtr(trkHits[i]->WireID().Plane)->GetCenter().X());
            if(!tpcGeo.ActiveBoundingBox().ContainsX((X), AllowedRel)) outBound++;
            if(outBound>maxOutBoundPoints) {
                return {{},{},maxOutBoundPoints+1};
            }
            driftedPositionVector.push_back({X, loc.Y(), loc.Z()});
            recI.push_back(trkHits[i]->Integral());
        }
        DriftedTrack thisDriftedTrack = {driftedPositionVector, recI, outBound};
        return thisDriftedTrack;
    }
}