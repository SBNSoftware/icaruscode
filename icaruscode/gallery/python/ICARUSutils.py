#!/usr/bin/env python


__doc__ = """
Collection of utilities to interface ICARUS with python, gallery and LArSoft.

This module requires ROOT.
"""

__all__ = [
  'loadICARUSgeometry',
  'justLoadICARUSgeometry',
]

import LArSoftUtils
import ROOT


################################################################################
### Geometry
###
def loadICARUSgeometry(config = None, registry = None):
  """Loads and returns ICARUS geometry with the standard ICARUS channel mapping.
  
  See `loadGeometry()` for the meaning of the arguments.
  """
  SourceCode = LArSoftUtils.SourceCode # alias
  
  SourceCode.loadHeaderFromUPS('icaruscode/Geometry/ChannelMapIcarusAlg.h')
  SourceCode.loadLibrary('icaruscode_Geometry')
  return LArSoftUtils.loadGeometry \
    (config=config, registry=registry, mapping=ROOT.geo.ChannelMapIcarusAlg)
# loadICARUSgeometry()


def justLoadICARUSgeometry(configFile, mapping = None):
  """Loads and returns ICARUS geometry from the specified configuration file.
  
  This is a one-stop procedure recommended only when running interactively.
  """
  if mapping is not None:
    raise NotImplementedError("Support for non-standard mapping not implemented yet.")
  return loadICARUSgeometry(config=LArSoftUtils.ConfigurationClass(configFile))
# justLoadICARUSgeometry()


################################################################################
