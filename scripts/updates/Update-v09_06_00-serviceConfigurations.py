#!/usr/bin/env python

import sys, re

import SerialSubstitution
from SerialSubstitution import AddProcessor, RunSubstitutor


################################################################################
if __name__ == "__main__":
  
  #############################################################################
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # C++ source code
  #
  Subst = AddProcessor(SerialSubstitution.ProcessorClass("configuration"))
  
  Subst.AddFileType("fcl")
  
  # include files
  Subst.AddWord         ("icarus_basic_services",  "icarus_common_services")
  Subst.AddWord         ("services_icarus.fcl",  "services_common_icarus.fcl")
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #############################################################################
  
  sys.exit(RunSubstitutor())
# 
