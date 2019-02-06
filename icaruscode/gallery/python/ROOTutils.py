#!/usr/bin/env python


__doc__ = """
Collection of utilities to ease interaction with ROOT.

Unsurprisingly, this module requires ROOT.
"""

__all__ = [
  "splitROOTpath", "createROOTpath",
  ]

import ROOT


################################################################################
### Print vectors easily
def TVector2ToString(v):
  return "( %g, %g )" % (v.X(), v.Y())
def TVector3ToString(v):
  return "( %g, %g, %g )" % (v.X(), v.Y(), v.Z())
def TLorentzVectorToString(v):
  return "( %g, %g, %g; %g )" % (v.X(), v.Y(), v.Z(), v.T())

ROOT.TVector2.__str__ = TVector2ToString
ROOT.TVector3.__str__ = TVector3ToString
ROOT.TLorentzVector.__str__ = TLorentzVectorToString


################################################################################
###  File management
###  
def splitROOTpath(path):
  """
  Returns the specified path split into file path and ROOT directory path.
  
  The `path` is in the form:
  "/UNIX/path/to/file.root:internal/ROOT/directory/and/object".
  The returned value is a pair `(filePath, dirPath)`: in the example, that
  would be `("/UNIX/path/to/file.root", "internal/ROOT/directory/and/object")`.
  
  Note: for compatibility with some ROOT tradition, the separator ':' can be
  replaced by '/'
  """
  
  # this implementation is not robust, but I am in a hurry :-P
  try:
    filePath, ROOTpath = path.rsplit('.root')
  except ValueError:
    raise RuntimeError("Path '{}' does not include a ROOT file.".format(path))
  filePath += '.root'
  ROOTpath.lstrip(':')
  ROOTpath.lstrip('/')
  return filePath, ROOTpath
  
# splitROOTpath()


def createROOTpath(path, fileMode = "UPDATE"):
  """
  Creates a complete ROOT directory path.
  
  The `path` is in the form:
  "/UNIX/path/to/file.root:internal/ROOT/directory/structure".
  The ROOT file `/UNIX/path/to/file.root` will be created with the specified
  `fileMode` (path may be relative), then the `TDirectoryFile` hierarchy
  `internal/ROOT/directory/structure` will be created under it.
  The return value is a pair `(file, dir)`, where `file` is a open `TFile`
  for `/UNIX/path/to/file.root` and `dir` is the `TDirectory` object of
  `structure`.
  
  Remember to keep track of `file`, or else python may close it compromising
  `dir` as well.
  """
  
  filePath, ROOTpath = splitROOTpath(path)
  
  ROOTfile = ROOT.TFile(filePath, fileMode)
  if not ROOTfile.IsOpen():
    raise RuntimeError \
      ("Can't open ROOT file '{}' in '{}' mode".format(filePath, fileMode))
  
  # instead of using `TDirectory.mkdir()`, we do that manually
  ROOTpathElements = ROOTpath.split('/')
  ROOTdir = ROOTfile
  for ROOTdirName in ROOTpathElements:
    if not ROOTdirName: continue # empty name does nothing
    daughterDir = ROOTdir.GetDirectory(ROOTdirName)
    if not daughterDir:
      daughterDir = ROOTdir.CreateDirectory(ROOTdirName)
    if not daughterDir:
      raise RuntimeError("Can't access directory '{}' under '{}'".format
       (ROOTdirName, ROOTdir.GetPath()))
    ROOTdir = daughterDir
  # for
  
  return ROOTfile, ROOTdir
  
# createROOTpath()


class DirectoryChanger:
  """
  Object changing ROOT directory while on scope.
  
  The purpose is to make a ROOT directory current only as long as it is needed.
  The most typical uses of this objects include the automatic restoration of
  the previous directory as the object falls out of scope.
  Two methods are supported:
  1. function scope:
        
        def writeEverythingInto(dir, everything):
          dirChanger = ROOTutils.DirectoryChanger(dir)
          for item in everything: item.Write()
        # writeEverythingInto()
        
  2. local scope (equivalent to using `activateDirectory()`):
        
        with DirectoryChanger(dir):
          for item in everything: item.Write()
        # with
        
  
  """
  def __init__(self, newDir = None, saveDir = None):
    if saveDir: self.saveDir(saveDir)
    else:       self.saveCurrentDir()
    if newDir: newDir.cd()
  # __init__()
  
  def saveCurrentDir(self): self.saveDir(ROOT.gDirectory)
  
  def saveDir(self, ROOTdir): self.oldDir = ROOTdir
  
  def restoreDir(self):
    if self.oldDir: self.oldDir.cd()
  
  def forget(self): self.oldDir = None
  
  def __del__(self):
    if self.oldDir: self.oldDir.cd()
  
  def __enter__(self):
    
  def __exit__(self, exc_type, exc_value, traceback):
    self.restoreDir()
    self.forget()
  
# DirectoryChanger()


def activateDirectory(ROOTdir):
  """
  Sets a directory with `DirectoryChanger`.
  
  Example:
        
        dir = outputFile.GetDirectory("plots")
        with activateDirectory(dir):
          for plot in plots: item.Write()
        # with
        
  """
  return DirectoryChanger(ROOTdir)


################################################################################
