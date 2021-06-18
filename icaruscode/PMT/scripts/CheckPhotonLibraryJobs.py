#!//usr/bin/env python
#
# Run with `--help` for terse help.
#
# Changes:
# 20201216 (petrillo@slac.stanford.edu) [v2.1]
#   way faster if not collecting output files and using existing good job list;
#   keyboard interruption informs that the file lists are unchanged.
#

__doc__ = """Checks the output of the jobs specified by their XML configuration file.

The jobs must have been submitted by `project.py` and must have completed already.

"""
__version__ = "2.1"

import sys, os
import re
import time
import logging


def removeSuffixes(
 s: "string to be processed",
 *suffixes: "strings to be removed from the end of `s`; better longer first",
 ) -> "a copy of `s` with all suffixes removed from its end":
  while True:
    for suffix in suffixes:
      if not s.endswith(suffix): continue
      s = s[:-len(suffix)]
      break # bash would say `continue 2`
    else: return s
  # while
# removeSuffixes()


class CachedValue:
  """On the first call, it creates and returns an object.
  On next calls, it returns the same object.
  """
  def __init__(self, fetchProc): self.fetchProc = fetchProc
  def __call__(self, *args, **kwargs) -> "Returns the cached value.":
    try: return self.cachedValue
    except AttributeError: return self._fetchValue(*args, **kwargs)
  # __call__()
  def __nonzero__(self) -> "Returns whether the object is cached":
    return hasattr(self, "cachedValue")
  def _fetchValue(self, *args, **kwargs):
    self.cachedValue = self.fetchProc(*args, **kwargs)
    del self.fetchProc # not needed any more
    return self.cachedValue
  # _fetchValue()
# CachedValue


class JobIDclass:
  """Class parsing a job ID."""
  
  class InvalidJobID(RuntimeError): pass
  
  JobIDpattern = re.compile(r'^([0-9]+)\.([0-9]+)@(.*)$')
  
  def __init__(self, jobIDstring): self.parse(jobIDstring)
  
  def parse(self, jobIDstring):
    res = JobIDclass.JobIDpattern.match(jobIDstring)
    if not res: raise InvalidJobID(jobIDstring)
    self._jobNo = int(res.group(1))
    self._subjobIndex = int(res.group(2))
    self._server = res.group(3)
    return self
  # parse()
  
  def jobNo(self): return self._jobNo
  def subjobIndex(self): return self._subjobIndex
  def server(self): return self._server
  def subjobID(self): return str(self.jobNo()) + "." + str(self.subjobIndex())
  def jobID(self): return self.subjobID() + "@" + self.server()
  def subjobTag(self): return str(self.jobNo()) + "_" + str(self.subjobIndex())
  
  def __str__(self): return self.jobID()
  def __repr__(self): return "%s(%r)" % (self.__class__.__name__, self.jobID())
  
# class JobIDclass


class FileListIterator:
  """Iterator: returns one file name at each iteration.
  
  The text file is expected to have one file name per line. Lines are stripped
  (`str.strip()`) of spaces.
  Empty lines and lines which are a comment are skipped.
  A comment line is a line whose first non-blank character is a hash character
  ('#').
  
  The value returned at each iteration is (lineNo, fileNo, fileName); if
  `withLineNo` or `withFileNo` are `False`, the respective element is omitted.
  If only the file name is requested (all other option `False`) the return value
  is not a tuple but just the file name.
  """
  def __init__(self,
  #  /, # Python 3.8 only
    listFile: "text file containing the list" = None,
    listName: "path of the text file containing the list" = None,
    withLineNo: "whether to return also the number of line in file" = False,
    withFileNo: "whether to return also the number of file" = False,
   ):
    assert listFile or listName
    
    if not listFile: listFile = open(listName, 'r')
    self.fileIter = iter(listFile)
    self.lineNo = 0 if withLineNo else None
    self.fileNo = 0 if withFileNo else None
  # __init__()
  
  def __iter__(self): return self
  
  def __next__(self):
    while True:
      fileName = next(self.fileIter).strip()
      if self.lineNo is not None: self.lineNo += 1
      if not fileName or fileName[0] == '#': continue
    
      if self.fileNo is not None: self.fileNo += 1
      if not self.fileNo and not self.lineNo: return fileName
      else: return tuple(filter(None, ( self.lineNo, self.fileNo, fileName, )))
    # while
  # next()
  
# class FileListIterator


class JobChecker:
  """Performs checks on job output and collects output files.
  
  """
  
  def __init__(self,
   baseName: "name used for the check (defaults are shaped after it)",
   goodList: "name of good job list (None: automatic; False: no list)" = None,
   badList: "name of bad job list (None: automatic; False: no list)" = None,
   fileList: "name of output file list (None: automatic; False: no list)" = None,
   skipKnownGoodJobs: "do not check the jobs already in the good list" = False,
   skipKnownBadJobs: "do not check the jobs already in the bad list" = False,
   ):
    
    self.checkBaseDir = os.path.dirname(baseName)
    # remove suffix:
    self.checkBaseName = removeSuffixes(
     os.path.splitext(os.path.basename(baseName))[0], # start with basename
     '_xml', '-xml'
     )
    
    self.goodListName, self.goodList, self.knownGoodJobs \
      = self.setupList(goodList, 'goodxml', mustExist=skipKnownGoodJobs)
    self.badListName, self.badList, self.knownBadJobs \
      = self.setupList(badList, 'badxml', mustExist=skipKnownBadJobs)
    self.outputFileListName, self.outputFileList, _ \
      = self.setupList(fileList, 'outputfile')
    
    if not skipKnownGoodJobs: self.knownGoodJobs = set()
    if not skipKnownBadJobs: self.knownBadJobs = set()
    
  # __init__()
  
  
  def setupList(self, listName, listTag, mustExist = False):
    """Returns the name of the list, an empty list to be filled and the existing
    content.
    """
    
    if listName is False: # set up for no list at all
      return None, [], set()
    
    if listName is None:
      listPath = os.path.join \
        (self.checkBaseDir, self.checkBaseName + "-" + listTag + ".list")
    elif not os.path.dirname(listName):
      listPath = os.path.join(self.checkBaseDir, listName)
    else: listPath = listName
    
    if os.path.isfile(listPath):
      listContent = set(FileListIterator(listName=listPath))
      (logging.info if mustExist else logging.debug) \
        ("File list %s contains already %d entries.", listPath, len(listContent))
    elif mustExist:
      raise RuntimeError("File list '{}' ({} list) is required to exist."
       .format(listPath, listTag))
    else: listContent = set()
    
    return listPath, [], listContent
  # setupList()
  
  
  def reset(self):
    """Resets all the counters and records."""
    
    self.goodList = []
    self.badList = []
    self.outputFileList = []
    
  # reset()
  
  def isCollectingOutput(self) -> "if the output file list is being filled":
    return self.outputFileListName is not None
  
  def checkFromFile(self,
   XMLfilePath: "path to the file containing XML job configuration",
   projectName: "target the specified project name in the file" = "",
   stageName: "target the specified stage name in the file" = "",
   maxJobs: "if not None, process at most this many jobs" = None
   ):
    
    nJobs = 0
    for lineNo, fileName in FileListIterator(listName=XMLfilePath, withLineNo=True):
      
      if maxJobs is not None and nJobs >= maxJobs:
        logging.info("Maximum number of jobs checked (%d).", maxJobs)
        break
      # if
      nJobs += 1
      
      try:
        self.checkJob(fileName, projectName=projectName, stageName=stageName)
      except KeyboardInterrupt: raise
      except Exception as e:
        logging.error("Error processing job '%s' (file list line %d): %s",
          fileName, lineNo, e)
      # try ... except
      
    # for file line
  # checkFromFile()
  
  
  def checkJob(self,
   XMLfilePath: "path to the file containing XML job configuration",
   projectName: "target the specified project name in the file" = "",
   stageName: "target the specified stage name in the file" = "",
   ):
    
    if not os.path.isfile(XMLfilePath):
      raise RuntimeError("Can't open file '{}'.".format(XMLfilePath))
    
    # jobInfo is a lazy callable returning the job information:
    # information will be extracted only the first time `jobInfo()` is executed
    jobInfo = CachedValue(lambda: self.getJobInfo \
      (XMLfilePath, projectName=projectName, stageName=stageName))
    
    XMLfileDir, XMLfileName = os.path.split(XMLfilePath)
    
    if self.knownGoodJobs and (XMLfilePath in self.knownGoodJobs):
      logging.info("%s: known as good, check skipped.", XMLfileName)
      good = True
    elif self.knownBadJobs and (XMLfilePath in self.knownBadJobs):
      logging.info("%s: known as bad, check skipped.", XMLfileName)
      good = False
    else:
      good = self.checkJobGoodness(jobInfo(), XMLfilePath)
    #
    
    if good:
      self.goodList.append(XMLfilePath)
      if self.isCollectingOutput(): self.collectJobOutputFiles(jobInfo())
      #
    else:
      self.badList.append(XMLfilePath)
    
  # checkJob()
  
  
  def getJobInfo(self,
   jobConfigFile: "path to the file containing XML job configuration",
   projectName: "target the specified project name in the file" = "",
   stageName: "target the specified stage name in the file" = "",
   ):
    
    #
    # get the project parsed by project.py
    #
    projInfo = project.get_project \
     (jobConfigFile, projectname=projectName, stagename=stageName)
    
    if not projInfo: # this message should be improved...
      raise RuntimeError("Job '{}' does not have project {} stage {}".format(
        jobConfigFile, repr(projectName), repr(stageName)
        ))
    #
    
    stageInfo = \
     next(filter(lambda stage: stage.name == stageName, projInfo.stages), None) \
     if stageName else projInfo.stages[0]
    
    if not stageInfo:
      raise RuntimeError("Job '{}' project {} does not have a stage {}".format(
        jobConfigFile, repr(project.name), repr(stageName)
        ))
    #
    
    return stageInfo
  # getJobInfo()
  
  
  def checkJobGoodness(self, jobInfo, jobName):
    
    class JobCheckError(RuntimeError): pass
    
    try:
      outputDir = jobInfo.outdir
      logging.debug("Job '%s' output directory: '%s'", jobName, outputDir)
      
      if not os.path.isdir(outputDir):
        raise JobCheckError("no output directory present ('{}')".format(outputDir))
      
      if not os.path.exists(os.path.join(outputDir, 'checked')):
        raise JobCheckError("not checked (run `project.py --checkana` first)")
      
      for jobID in map(JobIDclass, open(os.path.join(outputDir, 'jobids.list'), 'r')):
        logging.debug("Checking subjob '%s'", jobID)
        
        subjobDir = os.path.join(outputDir, jobID.subjobTag())
        logging.debug("Subjob '%s' output directory: '%s'", jobID, subjobDir)
        if not os.path.isdir(subjobDir):
          raise JobCheckError("job %s missing output directory" % jobID)
        
        statusFile = os.path.join(subjobDir, 'larStage0.stat')
        if not os.path.isfile(statusFile):
          raise JobCheckError("job %s missing status file" % jobID)
        
        try:
          status = int(open(statusFile, 'r').readline().strip())
        except KeyboardInterrupt: raise
        except Exception as e:
          raise JobCheckError("job %s failed reading status file '%s': %s"
           % (jobID, statusFile, e))
        #
        
        if status != 0:
          raise JobCheckError("job %s exited with error code %d" % (jobID, status))
        
      # for subjob
      
      expectedOutputFileList = os.path.join(outputDir, 'filesana.list')
      if not os.path.exists(expectedOutputFileList):
        raise JobCheckError("no output file list ('%s')" % expectedOutputFileList)
      
      expectedOutputFiles = list(FileListIterator(listName=expectedOutputFileList))
      if len(expectedOutputFiles) == 0:
        raise JobCheckError("job has no output file")
      
      foundOutputFiles = list(filter(os.path.isfile, expectedOutputFiles))
      if len(foundOutputFiles) != len(expectedOutputFiles):
        raise JobCheckError("only %d/%d output files still present"
          % (len(foundOutputFiles), len(expectedOutputFiles)))
      # if
      
    except JobCheckError as e:
      logging.error("%s: %s", jobName, e)
      return False
    else:
      logging.info("%s succeeded.", jobName)
      return True
    
  # checkJobGoodness()
  
  
  def collectJobOutputFiles(self, jobInfo):
    
    outputFileList = os.path.join(jobInfo.outdir, 'filesana.list')
    self.outputFileList.extend(FileListIterator(listName=outputFileList))
    
  # collectJobOutputFiles()
  
  
  def writeList(self, content, fileName, tag) -> "Whether the file list was written":
    
    if not fileName: return False # we are not asked to write the list
    
    # some file systems do not support overwriting
    if os.path.exists(fileName):
      try: os.remove(fileName)
      except IOError as e:
        logging.warning("Could not delete the old %s file list '%s': %s.",
         tag, fileName, e)
      # try ... except
    # if
    
    # we do not write the list if it would be empty
    if not len(content): return False
    
    listDir = os.path.dirname(fileName)
    os.makedirs(os.path.normpath(listDir), exist_ok=True)
    
    with open(fileName, 'w') as listFile:
      
      print("# {} file list created on {}: {:d} entries".format(
        tag, time.ctime(), len(content)
       ), file=listFile, 
       )
      listFile.write("\n".join(content))
      listFile.write("\n")
    # with
    
    logging.info("File list for %s created as '%s' with %d entries.",
      tag, fileName, len(content))
    
    return True
  # writeList()
  
  
  def writeSummary(self):
    
    nJobs = len(self.badList) + len(self.goodList)
    
    # save file lists
    if len(self.badList) > 0:
      
      if len(self.goodList) == 0:
        logging.info("None of the %d jobs was successful!!", nJobs)
      else:
        logging.info("%d/%d jobs were not successful.", len(self.badList), nJobs)
      
    elif nJobs > 0:
      
      logging.info("All %d jobs were successful.", nJobs)
      
    else:
      logging.error("No jobs checked.")
    
    try:
      self.writeList(self.goodList, self.goodListName, "successful jobs")
    except IOError as e:
      logging.critical("Could not write good job file list '%s': %s",
        self.goodListName, e)
    # try
    
    try:
      self.writeList(self.badList, self.badListName, "non-successful jobs")
    except IOError as e:
      logging.critical("Could not write bad job file list '%s': %s",
        self.badListName, e)
    # try
    
    try:
      self.writeList(self.outputFileList, self.outputFileListName, "output files")
    except IOError as e:
      logging.critical("Could not write output file list '%s': %s",
        self.outputFileListName, e)
    # try
    
    
  # writeSummary()
  
# class JobChecker



if __name__ == "__main__":
  
  logging.basicConfig()
  
  import argparse
  
  Parser = argparse.ArgumentParser(description=__doc__)
  
  Parser.add_argument \
    ("XMLfileList", help="list of XML configuration of the jobs to check")
  
  Parser.add_argument("--maxjobs", "-n", dest="MaxJobs", default=None, type=int,
    help="if specified, process at most this number of jobs")
  
  Parser.add_argument("--debug", "-d", action="store_true",
    help="enable debugging messages")
  
  Parser.add_argument \
   ('--version', '-V', action='version', version="%(prog)s v" + __version__)
  
  jobListGroup = Parser.add_argument_group("Job lists")
  
  jobListGroup.add_argument("--goodlist", "-g", dest="GoodJobList",
    default=None, help="name of the list to be created with all good jobs")
  jobListGroup.add_argument("--badlist", "-b", dest="BadJobList",
    default=None, help="name of the list to be created with all bad jobs")
  jobListGroup.add_argument("--outputlist", "-o", dest="OutputFileList",
    default=None, help="name of the list to be created with ROOT output files")
  
  jobListGroup.add_argument("--nogoodlist", "-G", dest="NoGoodJobList",
    action="store_true", help="do not create a list with all good jobs")
  jobListGroup.add_argument("--nobadlist", "-B", dest="NoBadJobList",
    action="store_true", help="do not create a list with all bad jobs")
  jobListGroup.add_argument("--nooutputlist", "-O", dest="NoOutputFileList",
    action="store_true", help="do not create a list with all output files")
  
  jobListGroup.add_argument("--skipgood", dest="SkipGoodJobs",
    action="store_true",
    help="do not check jobs that are already in the good list"
    )
  jobListGroup.add_argument("--skipbad", dest="SkipBadJobs",
    action="store_true",
    help="do not check jobs that are already in the bad list"
    )
  
  args = Parser.parse_args()
  
  logging.getLogger().setLevel(logging.DEBUG if args.debug else logging.INFO)
  
  try:
    import project
  except ImportError as e:
    logging.error("Could not load `project.py` as Python module: %s", e)
  
  jobChecker = JobChecker(
    args.XMLfileList,
    goodList=(False if args.NoGoodJobList else args.GoodJobList),
    badList=(False if args.NoBadJobList else args.BadJobList),
    fileList=(False if args.NoOutputFileList else args.OutputFileList),
    skipKnownGoodJobs=args.SkipGoodJobs,
    skipKnownBadJobs=args.SkipBadJobs,
    )
  
  try:
    jobChecker.checkFromFile(args.XMLfileList, maxJobs=args.MaxJobs)
  except KeyboardInterrupt:
    logging.warning("\nCheck interrupted; file lists will not be changed.")
    sys.exit(1)
  #
  
  jobChecker.writeSummary()
  
  sys.exit(0)
# __main__

