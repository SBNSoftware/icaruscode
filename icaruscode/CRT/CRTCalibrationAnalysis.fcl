#include "services_common_icarus.fcl"
#include "services_icarus_simulation.fcl"

process_name: CRTAnalysis

services:
{

  TFileService:           { fileName: "CRTCalibrationAnalysis.root" }

  TimeTracker:            {}

  message:                @local::icarus_message_services_prod_debug
                          @table::icarus_common_services

} 

source:
{
  module_type: RootInput

  maxEvents:  -1 

}

outputs:{}

physics:
{
  analyzers:
  {
    CRTCalibrationAnalysis: 
    {

      module_type:     "CRTCalibrationAnalysis"

      CRTDAQLabel:  "daqCRT"


    }
  }
  stream1:  [ out1 ]


  analysis: [ CRTCalibrationAnalysis ]


  end_paths: [ analysis ]
}
outputs: {
   out1: {
      compressionLevel: 1
      dataTier: "analysis"
      fastCloning: false
      fileName: "%ifb_%tc_ana.root"
      module_type: "RootOutput"
      saveMemoryObjectThreshold: 0
   }   
}

