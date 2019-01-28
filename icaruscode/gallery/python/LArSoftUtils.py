#!/usr/bin/env python


__doc__ = """
Collection of utilities to interface LArSoft with python and gallery.

This module requires ROOT.
"""

__all__ = [
  'readHeader',
  'SourceCode',
  'make_getValidHandle',
  'makeFileList',
  'forEach',
  'eventLoop',
  'findFHiCL',
  'loadConfiguration',
  'ConfigurationClass',
  'startMessageFacility',
  'ServiceRegistry',
  'ConfigurationHelper',
  'loadGeometry',
  'justLoadGeometry',
  'loadSimpleService',
  ]

import sys, os
import ROOT
import galleryUtils
from galleryUtils import (
  make_getValidHandle, makeFileList, forEach, eventLoop,
  findFHiCL, loadConfiguration, ConfigurationHelper, ConfigurationClass,
  startMessageFacility, ServiceRegistry, 
  )


################################################################################
### Pass-through
### 
SourceCode = galleryUtils.SourceCode
readHeader = galleryUtils.readHeader


################################################################################
### LArSoft
################################################################################
def loadGeometry(config=None, registry=None, mapping=None):
  """The argument `config` is an instance of `ConfigurationClass`.
  
  If a config object is provided, configurations will be read from there.
  Otherwise, they will be read from the registry.
  If a registry is provided, the services will be registered in there.
  """
  assert(config or registry)
  serviceName = 'Geometry'
  
  # deal with messagefacility first
  if not (registry and registry.has("message")):
    messageConfig = config.service("message") if config else registry.config("message")
    startMessageFacility(messageConfig) # use default name
    if registry: registry.register("message", None) # there is no direct access, sorry
  # if need to load message facility
    
  geometryConfig = config.service(serviceName) if config else registry.config(serviceName)
  if geometryConfig is None:
    raise RuntimeError("Failed to retrieve the configuration for %s service" % serviceName)
  
  if not mapping:
    SourceCode.loadHeaderFromUPS('larcorealg/Geometry/ChannelMapStandardAlg.h')
    mapping = ROOT.geo.ChannelMapStandardAlg
  SourceCode.loadHeaderFromUPS("larcorealg/Geometry/StandaloneGeometrySetup.h")
  SourceCode.loadLibrary("larcorealg_Geometry")
  service = ROOT.lar.standalone.SetupGeometry(mapping)(geometryConfig)
  if registry: registry.register(serviceName, service)
  return service
# loadGeometry()


################################################################################
def justLoadGeometry(configFile, mapping=None):
  """Loads and returns the geometry from the specified configuration file.
  
  This is a one-stop procedure recommended only when running interactively.
  """
  return loadGeometry(config=ConfigurationClass(configFile), mapping=mapping)
# justLoadGeometry()


################################################################################
def loadSimpleService \
  (serviceClass, config=None, registry=None, interfaceClass=None):
  """Loads a service assuming some simple requirements:
   * no dependency from other services
   * constructor accepts a FHiCL parameter set
  
  If a config object is provided, configurations will be read from there.
  Otherwise, they will be read from the registry.
  If a registry is provided, the services will be registered in there.
  
  The service configuration is read from an item called as `interfaceClass`,
  or `serviceClass` itself if `interfaceClass` is None, with "Service" appended.
  """
  
  assert(config or registry)
  serviceName = (interfaceClass if interfaceClass else serviceClass).__name__
  configKey = serviceName + "Service"
  
  config = config.service(configKey) if config else registry.config(configKey)
  if config is None:
    raise RuntimeError("Failed to retrieve the configuration for %s service" % serviceName)
  
  service = serviceClass(config)
  if registry: registry.register(serviceName, service)
  return service
# loadSimpleService()


################################################################################
