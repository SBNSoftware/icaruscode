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
import ROOTutils
from ROOTutils import ROOT


################################################################################
ICARUSchannelMappings = {
  'ICARUSsplitInductionChannelMapSetupTool': {
    'tool_type':       'ICARUSsplitInductionChannelMapSetupTool',
    'mapperClassName': 'icarus::ICARUSChannelMapAlg',
    'load':          [
                       'larcorealg_Geometry',
                       'icaruscode/Geometry/ICARUSChannelMapAlg.h',
                       'icaruscode/Geometry/ICARUSstandaloneGeometrySetup.h',
                       'icaruscode_Geometry',
                     ],
  }, # 'ICARUSsplitInductionChannelMapSetupTool'
  'ICARUSsingleInductionChannelMapSetupTool': {
    'tool_type':       'ICARUSsingleInductionChannelMapSetupTool',
    'mapperClassName': 'geo::ChannelMapIcarusAlg',
    'load':          [
                       'larcorealg_Geometry',
                       'icaruscode/Geometry/ChannelMapIcarusAlg.h',
                       'icaruscode/Geometry/ICARUSstandaloneGeometrySetup.h',
                       'icaruscode_Geometry',
                     ],
  }, # 'ICARUSsingleInductionChannelMapSetupTool'
}
DefaultChannelMapping = 'ICARUSsplitInductionChannelMapSetupTool'

################################################################################
### Geometry
###
def loadICARUSchannelMappingClass(
  config: "ConfigurationClass object with complete job configuration",
  ) -> "Class object for the proper channel mapping":
  
  #
  # we need to:
  # 1. find out which mapping is required
  # 2. load the proper libraries
  # 3. return the Python class object for the mapping class we want
  #
  
  #
  # 1. find out which mapping is required: known configurations
  #
  helperConfig = config.service('ExptGeoHelperInterface')
  if helperConfig.get(str)('service_provider') != 'IcarusGeometryHelper':
    raise RuntimeError(
     "ExptGeoHelperInterface in configuration is '{}', not IcarusGeometryHelper"
     .format(helperConfig['service_provider'])
     )
  # if
  if helperConfig.has_key('Mapper'):
    try:
      mapperConfig = helperConfig.get(ROOT.fhicl.ParameterSet)('Mapper')
      plugin_type = mapperConfig.get(str)('tool_type')
    except:
      raise RuntimeError(
        "Geometry helper configuration is missing the tool_type:\n"
        + helperConfig.to_indented_string()
        )
    # try ... except
  else: plugin_type = DefaultChannelMapping
  
  
  #
  # 2. load the proper libraries
  #
  
  # get the specification record
  try: mappingInfo = ICARUSchannelMappings[plugin_type]
  except KeyError:
    # when you get to this error, check that the tool name in the configuration
    # is actually spelled correctly first...
    raise RuntimeError(
     "Mapping plug in not supported: '{}': Python library needs to be updated."
     .format(plugin_type)
     )
  # try ... except
  
  # 
  # load the library which defines ROOT.lar.standalone.SetupGeometry as template;
  # otherwise, specializations will actually be interpreted as methods
  # (which blows up everything); on top of that, ROOT must be forced to
  # instantiate ROOT.lar.standalone.SetupGeometry, because once the
  # specialization is seen in a header before its template, it's over
  # (that's the difference with C++).
  #
  LArSoftUtils.SourceCode.loadHeaderFromUPS('larcorealg/Geometry/StandaloneGeometrySetup.h')
  # remember: the side effects of asking for ROOT.lar.standalone.SetupGeometry are essential
  if not isinstance(ROOT.lar.standalone.SetupGeometry, ROOT.TemplateProxy):
    try:
      raise RuntimeError(
        "Internal error: ROOT.lar.standalone.SetupGeometry should have been a template, it's a '{}'"
        .format(ROOT.lar.standalone.SetupGeometry.__class__)
        )
    except AttributeError:
      raise RuntimeError(
        "Internal error: ROOT.lar.standalone.SetupGeometry should have been a template, it's nothing."
        )
  # if
  
  # load the libraries
  for codeObj in mappingInfo.get('load', []):
    LArSoftUtils.SourceCode.load(codeObj)
  
  # get the class object
  try: mapperClass = ROOTutils.getROOTclass(mappingInfo['mapperClassName'])
  except AttributeError:
    # this needs investigation, as the code above should be sufficient to it
    raise RuntimeError(
      "The library with '{}' has not been correctly loaded!"
      .format(mappingInfo['mapperClassName'])
      )
  # try ... except
  
  #
  # 3. return the Python class object for the mapping class we want
  #
  return mapperClass
  
# loadICARUSchannelMappingClass()


def loadICARUSgeometry(
  config = None, registry = None, mappingClass = None,
  ):
  """Loads and returns ICARUS geometry with the standard ICARUS channel mapping.
  
  See `loadGeometry()` for the meaning of the arguments.
  """
  
  if mappingClass is None: mappingClass = loadICARUSchannelMappingClass(config)
  return LArSoftUtils.loadGeometry \
    (config=config, registry=registry, mapping=mappingClass)
# loadICARUSgeometry()


def justLoadICARUSgeometry(configFile, mappingClass = None):
  """Loads and returns ICARUS geometry from the specified configuration file.
  
  This is a one-stop procedure recommended only when running interactively.
  """
  return loadICARUSgeometry(config=LArSoftUtils.ConfigurationClass(configFile))
# justLoadICARUSgeometry()


################################################################################
