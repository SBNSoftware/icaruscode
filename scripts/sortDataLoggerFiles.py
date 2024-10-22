#!/usr/bin/env python
#
# Documentation can be rendered in python interactive session with
# `import sortDataLoggerFiles; help(sortDataLoggerFiles)` and on command line
# with `sortDataLoggerFiles --help`.
# 
# It may require ROOT.
# 
# Changes:
# 20210411 (petrillo@slac.fnal.gov) [1.0]
#   first public version
# 20210518 (petrillo@slac.fnal.gov) [1.1]
#   added options for duplicate events
# 20210602 (petrillo@slac.fnal.gov) [1.2]
#   added optional stream name to the file name pattern;
#   fixed a bug where first logger option value would be ignored
# 20220222 (petrillo@slac.fnal.gov) [1.3]
#   added support for a new file name format, and for multiple formats
# 20230822 (petrillo@slac.fnal.gov) [1.4]
#   including the first timestamp in the file stats, if available
# 20240326 (petrillo@slac.fnal.gov) [1.5]
#   added support for another a new file name format
#

import sys, os
import re
import logging

__doc__ = """Sorts a list of data logger output files.

File paths are read from all the specified file lists in sequence, or from
standard input if no file list is specified.

If a line is encountered that does not match the typical file name pattern,
that line is ignored and a warning is printed.

Comments and empty lines at the beginning of the first file list are printed
at the top of the output as they are. All other comments and empty lines are
printed at the end of the output.

Note that it is possible to sort "in place" by specifying the same file list as
input and output.

Duplicate files are files on the same run, data logger cycle and data logger
number. By default, only the first of the duplicate files is written into the
output list, and only the number of duplicates is printed. Options allow to
write a detailed list of duplicate files on screen and on disk, or not to check
for duplication altogether.

"""

__author__ = 'Gianluca Petrillo (petrillo@slac.stanford.edu)'
__date__ = 'February 22, 2022'
__version__ = '1.5'


class CycleCompareClass:
  """Provides less() to compare numbers with a single offset cycle.
  
  For example, with offset 3 the order of [0:20] would be, from the lowest:
  [ 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 0, 1, 2, ]
  """
  def __init__(self, first): self.first = first
  def less(self, a, b): return (a < b) == ((a < self.first) == (b < self.first))
# class CycleCompareClass


class DAQIndex:
  def __init__(self, builder, thread = 1):
    if isinstance(builder, str): builder = eval(builder)
    try:
      self.builder, self.thread = map(int, builder) # support first parameter being a pair
      assert thread == 1, "When builder is a pair, thread should not be specified."
    except TypeError:
      self.builder = int(builder)
      self.thread = int(thread)
  def key(self): return ( self.builder, self.thread )
  def __lt__(self, other): return self.key() < other.key()
  def __str__(self): return f"EventBuilder{self.builder}_art{self.thread}"
# class DAQIndex


class FileNameParser:
  """Static object (namespace?) containing file name parsing utilities.
  
  All supported file name patterns are in the `Patterns` class variable.
  The static method `match()` tries all of them, in order.
  """
  Patterns = [
    {
      'Name': 'multibuilder',
      # pattern parameters:           run            builder  thread  pass  stream name    timestamp      filler (timestamp)
      #                               <1>               <2>      <3>  <4>         <5>    . <7>             <8>
      'Pattern': re.compile(r"data_run(\d+)_EventBuilder(\d+)_art(\d+)_(\d+)_fstrm([^_]*)(_(\d{8}T\d{6}))?(_.*)\.root"),
      'Parameters': {
        'DataLogger': ( (2, 3), DAQIndex ),
        'StreamName': ( 5, str ),
        'RunNumber' : ( 1, int ),
        'PassCount' : ( 4, int ),
        'Timestamp' : ( 7, str ),
      },
    }, # multibuilder
    {
      'Name': 'multistream',
      # pattern parameters:   data logger stream stream name run  pass     timestamp      filler (timestamp)
      #                              <1>  <2>    <3>         <4>   <5>   . <7>            <8>
      'Pattern': re.compile(r"data_dl(\d+)(_fstrm([^_]*))?_run(\d+)_(\d+)(_(\d{8}T\d{6}))?(_.*)\.root"),
      'Parameters': {
        'DataLogger': ( 1, DAQIndex ),
        'StreamName': ( 3, str ),
        'RunNumber' : ( 4, int ),
        'PassCount' : ( 5, int ),
        'Timestamp' : ( 7, str ),
      },
    }, # multistream
    {
      'Name': 'general',
      # pattern parameters:   stream name data logger run  pass    timestamp      filler (timestamp)
      #                       <1>            <2>      <3>   <4>  . <6>             <7>
      'Pattern': re.compile(r"([^_]+)_data_dl(\d+)_run(\d+)_(\d+)(_(\d{8}T\d{6}))?(_.*)?\.root"),
      'Parameters': {
        'DataLogger': ( 2, DAQIndex ),
        'StreamName': ( 1, str ),
        'RunNumber' : ( 3, int ),
        'PassCount' : ( 4, int ),
        'Timestamp' : ( 6, str ),
      },
    }, # multistream
  ] # Patterns
  
  class ParsedNameClass:
    def __init__(self, name, fields):
      self.name = name
      self.fields = fields
    # __init__()
    def __bool__(self): return len(self.fields) > 0
    def get(self, *fieldNames):
      return tuple(self.fields.get(fieldName, None) for fieldName in fieldNames)
  # class ParsedNameClass
  
  def __init__(self): pass

  @staticmethod
  def match(name):
    # the first successful pattern is used
    mkl = lambda v: (( v, ) if isinstance(v, int) else v) # make list
    for patternInfo in FileNameParser.Patterns:
      match = patternInfo['Pattern'].match(name)
      if match is None: continue
      d = {
        name: ( type_(value) if (value := match.group(*mkl(index))) else None )
        for name, ( index, type_ ) in patternInfo['Parameters'].items()
        }
      d['Name'] = patternInfo['Name']
      return FileNameParser.ParsedNameClass(name, d)
    else: return FileNameParser.ParsedNameClass(name, {})
  # match()
# class FileNameParser


class FileInfoClass:
  """This class collects information about a input file, including a sorting
  criterium.
  """
  
  DefaultXRootDprotocolHead = 'root://fndca1.fnal.gov:1094/'
  GenericXRootDprotocolHead = r'root://fn\w+\.fnal\.gov:[0-9]{1,5}/'
  XRootDprotocolDir = 'pnfs/fnal.gov/usr'
  POSIXprotocolHead = '/'
  POSIXprotocolDir = 'pnfs'
  
  POSIXPattern = re.compile(
    POSIXprotocolHead.replace('.', r'\.')
    + POSIXprotocolDir.replace('.', r'\.')
    + r"/([^/]+)/(.*)")
  XRootDPattern = re.compile(
    GenericXRootDprotocolHead
    + XRootDprotocolDir.replace('.', r'\.')
    + r"/(.*)"
    )
  _DataLoggerSorter = CycleCompareClass(first=DAQIndex(4))
  
  @staticmethod
  def getFirstDataLogger(index): return FileInfoClass._DataLoggerSorter.first
  @staticmethod
  def setFirstDataLogger(index):
    FileInfoClass._DataLoggerSorter = CycleCompareClass(first=index)
  
  def __init__(self,
               line: "input file line (should include endline)",
               source: "an arbitrary identifier to track the origin of the line" = None,
               ):
    """Constructor: use and parse the specified input file line."""
    self.line = line
    self.source = source
    self.path = line.strip()
    self.protocolAndDir, self.name = os.path.split(self.path)
    parsedName = FileNameParser.match(self.name)
    self.is_file = bool(parsedName)
    if self.is_file:
      self.dataLogger, self.run, self.stream, self.pass_, self.timestamp \
        = parsedName.get('DataLogger', 'RunNumber', 'StreamName', 'PassCount', 'Timestamp')
  # __init__()
  
  def __lt__(self, other):
    """Comparison: run, then pass, then (offset cycled) data logger number."""
    if not self.is_file:
      raise RuntimeError \
        ("Sorting not supported for non-file objects ('%s')" % self.path)
    # if
    if self.run < other.run: return True
    if self.run > other.run: return False
    
    if self.pass_ < other.pass_: return True
    if self.pass_ > other.pass_: return False
    
    if self.dataLogger != other.dataLogger:
      return \
        FileInfoClass._DataLoggerSorter.less(self.dataLogger, other.dataLogger)
    
    if self.timestamp < other.timestamp: return True
    if self.timestamp > other.timestamp: return False
    
    assert (self.stream is None) == (other.stream is None)
    return False if self.stream is None else self.stream < other.stream
    
  # __lt__()
  
  def __str__(self):
    s = f"Run {self.run} cycle {self.pass_} data logger {self.dataLogger}"
    if self.timestamp: s += f" time {self.timestamp}"
    if self.stream: s += f" stream {self.stream}"
    return s
  
  def pathToXRootD(self) -> "stored file path in XRootD format":
    if not self.is_file:
      raise RuntimeError(
       "XRootD conversion not supported for non-file objects ('%s')" % self.path
       )
    # if not file
    match = FileInfoClass.POSIXPattern.match(self.path)
    return os.path.join(
      FileInfoClass.DefaultXRootDprotocolHead, FileInfoClass.XRootDprotocolDir,
      *match.group(1, 2)
      ) if match else self.path
    
  # pathToXRootD()
  
  def pathToPOSIX(self) -> "stored file path in POSIX (PNFS local) format":
    if not self.is_file:
      raise RuntimeError(
       "XRootD conversion not supported for non-file objects ('%s')" % self.path
       )
    # if not file
    match = FileInfoClass.XRootDPattern.match(self.path)
    return os.path.join(
        FileInfoClass.POSIXprotocolHead, FileInfoClass.POSIXprotocolDir,
        match.group(1)
      ) if match else self.path
    
  # pathToXRootD()
  
# class FileInfoClass


class MinimumAccumulator:
  def add(self, data, key = None):
    if key is None: key = data
    try:
      if key >= self.minKey: return False
    except AttributeError: pass # no self.minKey yet?
    self.minKey = key
    self.minData = data
    return True
  # add()
  def min(self): return self.minData
# class MinimumAccumulator


def findFirstCycle(files, stream):
  firstLogger = None
  firstPassFiles = []
  wrapped = False
  for info in files:
    if info.stream != stream: continue
    if firstLogger == info.dataLogger: break # cycle completed
    if wrapped and info.dataLogger > firstLogger: break # cycle completed
    
    if firstLogger is None: firstLogger = info.dataLogger
    elif not wrapped and info.dataLogger < firstLogger: wrapped = True
    
    firstPassFiles.append(info)
    logging.debug("Added cycle %d logger %s stream %s time %s to first cycle list",
      info.pass_, info.dataLogger, info.stream,
      (info.timestamp if info.timestamp else "unknown"),
      )
  # for
  return firstPassFiles
# findFirstCycle()


def extractFirstEvent(filePath):
  try: import ROOT
  except ImportError:
    raise RuntimeError("""ROOT python module could not be loaded.
      In this condition, you'll have to skip the autodetection of the first logger
      by explicitly specifying its number as option to the script."""
      )
  # try ... except
  logging.debug("Opening '%s' for event number check...", filePath)
  srcFile = ROOT.TFile.Open(filePath, "READ")
  if not srcFile:
    raise RuntimeError \
      ("Failed to open '%s' for event number extraction." % filePath)
  #
  try: firstEvent = next(iter(srcFile.Events)) # go PyROOT
  except StopIteration:
    logging.debug("File '%s' appears to contain no events.", filePath)
    return None
  firstEventNumber = firstEvent.EventAuxiliary.event() # keep going PyROOT
  
  logging.debug("First event from '%s': %d", filePath, firstEventNumber)
  return firstEventNumber
# extractFirstEvent()


def detectFirstLogger(fileInfo):
  # in the end, we don't need a stream-aware algorithm to determine which
  # data logger received the first event, as long as we have all relevant
  # streams represented
  lowestEvent = MinimumAccumulator()
  for stream, files in fileInfo.items():
    if not len(files): continue
    for info in files:
      
      
      firstEvent = extractFirstEvent(info.pathToXRootD())
      if firstEvent is not None:
        lowestEvent.add(info, key=firstEvent)
        if firstEvent == 1: # can't get lower than this!
          firstLogger = info.dataLogger
          logging.debug("Definitively detected first logger: %s", firstLogger)
          return firstLogger
    # for files
  # for
  try: firstLogger = lowestEvent.min().dataLogger
  except AttributeError:
    # this is in general a problem because it implies that we are failing to
    # correctly parse the list of input files
    raise RuntimeError("No data found for the first data logger pass.")
  logging.debug("Detected first logger: %s", firstLogger)
  return firstLogger
# detectFirstLogger()


def buildFileIndex(
  fileInfo: "list with information from all files",
  ) -> "a dictionary: { key -> list of files }":
  
  fileKey = lambda info: ( info.run, info.pass_, info.dataLogger.key(), info.stream, info.timestamp )
  index = {}
  for info in fileInfo:
    index.setdefault(fileKey(info), []).append(info)
  return index
# buildFileIndex()


if __name__ == "__main__":
  
  logging.basicConfig(level=logging.INFO)
  
  import argparse
  
  parser = argparse.ArgumentParser(description=__doc__)
  parser.set_defaults(skipDuplicates=True)
  
  parser.add_argument('inputFiles', nargs="*", metavar='inputFileNames',
    help='input file lists [one from stdin by default]')
  parser.add_argument('--firstlogger',
    help='index of the first data logger in the cycle')
  parser.add_argument('--output', '-o', default=None,
    help=
     'name of the file to write the resulting list into (overwritten!) [stdout]'
    )
  parser.add_argument('--nooutput', action="store_true",
    help='do not print on screen nor write to file the files in input')
  
  duplGroup = parser.add_argument_group(title="duplicate file options")
  duplGroup.add_argument('--printduplicates', '-d', action="store_true",
    help='print duplicate files on screen')
  duplGroup.add_argument('--skipduplicates', '-S', dest='skipDuplicates',
    action="store_true",
    help='do not include duplicate files in the list (default)'
    )
  duplGroup.add_argument('--keepduplicates', '-K', dest='skipDuplicates',
    action="store_false",
    help='include also duplicate files in the list (default)'
    )
  duplGroup.add_argument('--duplicatelist', '-D', type=str, default=None,
    help='name of a file list to be created with duplicate entries')
  
  parser.add_argument('--xrootd', '--root', '-X', action="store_true",
    help='convert the paths to XRootD URL')
  parser.add_argument('--posix', '-P', action="store_true",
    help='convert the paths to local POSIX path')
  parser.add_argument('--debug', action="store_true",
    help='prints out debugging messages')
  parser.add_argument \
    ('--version', '-V', action='version', version='%(prog)s ' + __version__)

  args = parser.parse_args()
  
  if args.debug: logging.getLogger().setLevel(logging.DEBUG)
  
  if args.xrootd and args.posix:
    raise RuntimeError("XRootD and POSIX output format options are exclusive.")
  
  printDuplicates = args.printduplicates
  skipDuplicates = args.skipDuplicates
  makeDuplicateList = args.duplicatelist
  
  # "sources" are given directly as input (None = sys.stdin)
  sources = args.inputFiles if args.inputFiles else [ "<stdin>" ]
  
  # "inputFiles" are all the files found in the sources
  inputFiles = (
      [ file_ ] if file_.endswith('.root') else open(file_, 'r')
        for file_ in args.inputFiles 
    ) if args.inputFiles else [ sys.stdin, ]
  
  preComments = []
  postComments = []
  fileInfo = []
  sourceNames = []
  for iSource, file_ in enumerate(inputFiles):
    isSingleFile = isinstance(file_, list) and len(file_) <= 1
    for iLine, line in enumerate(file_):
      info = FileInfoClass(line, source=( iSource, None if isSingleFile else iLine + 1 ))
      if not info.is_file:
        if not info.path or info.path.startswith('#'):
          (postComments if fileInfo else preComments).append(info.line)
          continue
        else:
          logging.warning \
           ("Line %d ('%s') does not match file pattern." % (iLine, info.path))
          continue
      # if not file
      fileInfo.append(info)
    # for line in file
  # for input files
  
  Streams = list(set( info.stream for info in fileInfo ))
  logging.debug("%d data files in %d streams: %s",
    len(fileInfo), len(Streams),
    ", ".join(stream if stream else "<none>" for stream in Streams)
    )
  
  if fileInfo and (args.firstlogger is None):
    # uses internal FileInfoClass ordering (firstLogger not set: any will do)
    fileInfo.sort()
    firstPassFiles = dict( ( stream, findFirstCycle(fileInfo, stream) )
      for stream in Streams )
    assert firstPassFiles
    firstLogger = detectFirstLogger(firstPassFiles)
  else:
    firstLogger = DAQIndex(args.firstlogger if args.firstlogger is not None else 4)
  
  FileInfoClass.setFirstDataLogger(firstLogger)
  
  fileInfo.sort() # uses internal FileInfoClass ordering
  
  #
  # deal with duplicates
  #
  if printDuplicates or makeDuplicateList or skipDuplicates:
    nDuplicates = 0
    fileIndex = buildFileIndex(fileInfo)
    uniqueFiles = [] if skipDuplicates else None
    duplicateFiles = [] if makeDuplicateList else None
    # we rely on insertion-ordered dictionary guarantee of Python 3.7
    for fileList in fileIndex.values():
      mainInfo = fileList[0]
      if uniqueFiles is not None: uniqueFiles.append(mainInfo)
      if len(fileList) > 1:
        nDuplicates += len(fileList) - 1
        if duplicateFiles is not None: duplicateFiles.extend(fileList[1:])
        if printDuplicates:
          firstSource = mainInfo.source[0]
          msg = f"{mainInfo} with {len(fileList) - 1} duplicates of"
          
          if len(sources) > 1: msg += f" {sources[mainInfo.source[0]]}"
          if mainInfo.source[1] is not None: msg += f" line {mainInfo.source[1]}"
          msg += ":"
          for info in fileList[1:]:
            if info.source[0] != firstSource: msg += f"{sources[info.source[0]]}"
            if info.source[1] is not None: msg += f" line {info.source[1]}"
            msg += ";"
          # for
          logging.info(msg)
        # if print duplicates
      # if duplicates
    # for
    if nDuplicates: logging.info(f"Found {nDuplicates} duplicate files.")
    if duplicateFiles:
      with open(makeDuplicateList, 'w') as DuplicateListFile:
        for info in duplicateFiles: # lines still have their <CR>
          print(info.line, file=DuplicateListFile, end='')
      logging.info(f"{nDuplicates} duplicate file names written in '{makeDuplicateList}'.")
    # if we have duplicates and we write them
  # if print or store duplicates
  
  fileListContent = uniqueFiles if skipDuplicates else fileInfo
  
  
  #
  # print everything
  #
  
  # NOTE: keep this after all the input has been read,
  #       so that input files can be safely overwritten
  if not args.nooutput:
    outputFile = open(args.output, 'w') if args.output else sys.stdout
    
    # <CR> were not removed from `line`
    for line in preComments: outputFile.write(line)
    for info in fileListContent:
      if args.posix: line = info.pathToPOSIX() + '\n'
      elif args.xrootd: line = info.pathToXRootD() + '\n'
      else: line = info.line
      outputFile.write(line)
    for line in postComments: outputFile.write(line)
    
    if outputFile is not sys.stdout:
      logging.info \
        (f"{len(fileListContent)} file entries written into '{outputFile.name}'.")
      del outputFile
    # if
  else:
    logging.info(f"Found {len(fileListContent)} file entries.")
  # if ... else
  
  sys.exit(0)
  
# main
