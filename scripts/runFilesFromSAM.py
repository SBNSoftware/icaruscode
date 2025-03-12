#!/usr/bin/env python
#
# Changes:
# 20230922 (petrillo@slac.fnal.gov) [1.0]
#   first public version
#

import sys, os
import time
import logging

__doc__ = f"""Create file lists from ICARUS data runs.

For each query or definition, a SAM query is executed and the file locations
are written all on screen or each into their own input file, depending on the
arguments.
Only one output pattern can be specified. That pattern can hold placeholders:
  %n for the name as specified (may contain unacceptable characters!);
  %s for a sanitized version of the name;
  %t for a tag suitable for a file name;
  %m for the requested maximum limit of entries;
  %# for the number of input on the command line.

If the automatic naming option is specified instead (equivalent to specifying
an empty output pattern), a name is chosen with suffix '.filelist', that
includes the verbatim definition or a sanitized version of the query, and the
maximum number of entries requested, if any.

If the output name ends with a directory separator ({os.sep}), it is used as the
output directory, while the name is generated as if the output option were empty
(that is, as the automatic naming were requested).

By default, existing file lists are not overwritten.


"""

__author__ = 'Gianluca Petrillo (petrillo@slac.stanford.edu)'
__date__ = 'September 22, 2023'
__version__ = '1.0'

from samweb_client.client import SAMWebClient
import samweb_client.exceptions as samexcpt

logging.basicConfig()
logger = logging.getLogger(__name__)


# ------------------------------------------------------------------------------
class FatalError(RuntimeError): pass

class QueryElement: pass
class QueryDirect(QueryElement):
  def __init__(self, query):
    self.query = query
  def __str__(self): return self.query
# class QueryDirect

class QueryAtom(QueryElement):
  def __init__(self, key, *values):
    self.key = str(key)
    self.values = values[:]
  def add(self, *values): self.values.extend(values)
  def __str__(self):
    if not self.values: return ""
    elif len(self.values) == 1: return f"{self.key}={self.values[0]}"
    else: return f"{self.key} in ( {', '.join(self.values)} )"
  # __str__()
# class QueryAtom

class QueryWrappedElement(QueryElement):
  def __init__(self, content): self.value = content
  def __str__(self): return f"({self.value})"
# class QueryWrappedElement

class Query:
  def __init__(self, *elements):
    self.elements = list(map(Query.makeElement, elements))
  def __str__(self): return " AND ".join(map(str, self.elements))
  def __bool__(self): return bool(self.elements)
  def add(self, *elements):
    self.elements.extend(map(Query.makeElement, elements))
  def addAtom(self, key, *values):
    self.add(atom := QueryAtom(key, *values))
    return atom
  @staticmethod
  def makeElement(value):
    return value if isinstance(value, QueryElement) else QueryWrappedElement(value)
# class Query


class InputSpecification:
  Type = "run"
  
  def __init__(self, inputSpec):
    self.name = inputSpec
    self.inputType = self.Type
  
  def isType(self, type_): return self.inputType == type_

  def sanitized(self): return self.Sanitize(self.name)
  def __str__(self): return self.name
  def tag(self): return str(self)
  def describe(self): return str(self)
  def queryElement(self):
      raise NotImplementedError(f"Input specification of type {self.inputType} not known.")
  #
  @staticmethod
  def Sanitize(query):
    return "".join \
      (map(lambda c: c if c.isalnum() or c in "-_" else "_", str(query)))
  # Sanitize()
  
  @staticmethod
  def make(inputSpec):
    if isinstance(inputSpec, InputSpecification): return inputSpec
    try:
      run = int(inputSpec)
    except ValueError: pass
    else: return InputRun(run)
    if inputSpec.endswith('.root'): return InputROOTfile(inputSpec)
    return (InputSAMquery if ' ' in inputSpec else InputSAMdef)(inputSpec)
  # make()
# class InputSpecification

class InputRun(InputSpecification):
  Type = "run"
  
  def __init__(self, inputSpec):
    self.name = inputSpec
    self.inputType = InputRun.Type
  def __str__(self): return format(self.name, "05d")
  def tag(self): return f"run{self.name:05d}"
  def describe(self): return f"run {self.name}"
  def queryElement(self): return QueryAtom("run_number", self.name)
# class InputRun

class InputSAMquery(InputSpecification):
  Type = "query"
  
  def __init__(self, inputSpec):
    super().__init__(inputSpec)
    self.inputType = self.Type
  
  def tag(self): return self.sanitized()
  def describe(self): return f"SAM query '{self.sanitized()}'"
  def queryElement(self): return QueryDirect(self.name)
# class InputSAMquery

class InputSAMdef(InputSpecification):
  Type = "def"
  
  def __init__(self, inputSpec):
    super().__init__(inputSpec)
    self.inputType = self.Type
  
  def describe(self): return f"SAM definition '{self.name}'"
  def queryElement(self): return QueryAtom("dataset.tag", f'"{self.name}"') # ??
# class InputSAMdef

class InputROOTfile(InputSpecification):
  Type = "ROOT"
  
  def __init__(self, inputSpec):
    super().__init__(inputSpec)
    self.inputType = self.Type
    self.basename = os.path.basename(self.name)
    self.stem = os.path.splitext(self.basename)[0]
  
  def tag(self): return self.stem
  def describe(self): return f"ROOT file '{self.name}'"
  def queryElement(self): return QueryAtom("file_name", f'"{self.basename}"')
# class InputROOTfile


class FileListMakerClass:
  
  Stages = {
    'raw': [
      'raw',
      'decoded', # ancient
      ],
    'stage1': [
      'stage1',
      'stage1_caf_larcv',
      ],
    } # stages
  
  def __init__(self, samweb = None, options = {}):
    
    self.options = dict(
      schema=options.schema,
      output=getattr(options, 'output', None),
      experiment=getattr(options, 'SAMexperiment', None),
      maxFiles=getattr(options, 'maxFiles', None),
      header=getattr(options, 'header', False),
      overwrite=getattr(options, 'overwrite', False),
      stage=getattr(options, 'stage', None),
      streams=getattr(options, 'streams', []),
      project=getattr(options, 'project', None),
      location=getattr(options, 'location', None),
      )
    if not samweb:
      samweb = SAMWebClient(experiment=self.options['SAMexperiment'])
    self.samweb = samweb
  # __init__

  def process(self, inputSpec, index=None):
    
    inputSpec = InputSpecification.make(inputSpec)
    logger.debug("Input spec %s: %s", inputSpec.inputType, inputSpec)
    fileURLs = self.queryFileURLs(inputSpec)
    
    outputLines, nEntries = self.prepareOutput(inputSpec, fileURLs)
    
    if self.options['schema'].upper() == 'POSIX':
      outputLines = list(map(FileListMakerClass.XRootDtoPOSIX, outputLines))
    
    outputPath = self.writeOutput(inputSpec, index, outputLines)
    
    return outputPath, nEntries
  # process()
  
  def prepareQuery(self, inputSpec):
    if isinstance(inputSpec, InputSAMdef):
      return dict(dimensions=Query(), defname=inputSpec.name)
    else:
      return dict(dimensions=Query(inputSpec.queryElement()))
  # prepareQuery()
  
  def queryFileURLs(self, inputSpec):
    
    # get the list of URL
    locationArg = self.prepareQuery(inputSpec)
  
    if self.options['stage']:
      stages = FileListMakerClass.Stages.get(self.options['stage'], [ self.options['stage'] ])
      if self.options['stage'].lower() == "raw":
        locationArg['dimensions'].addAtom('data_tier', *stages)
      else:
        locationArg['dimensions'].addAtom('icarus_project.stage', *stages)
    elif isinstance(inputSpec, InputRun):
      raise RuntimeError("For run queries a stage must be specified")
    
    if self.options['streams']:
      locationArg['dimensions'].addAtom('data_stream', *self.options['streams'])
    
    if self.options['project']:
      locationArg['dimensions'].addAtom('icarus_project.version', self.options['project'])
    
    
    # replace the Query object by the query dimensions
    if locationArg['dimensions']:
      locationArg['dimensions'] = str(locationArg['dimensions'])
    else: del locationArg['dimensions']
    
    schema = self.options['schema']
    if schema.upper() == 'POSIX': schema = 'root'
  
    logger.debug("Query: %s", locationArg)
    try:
      filesAndLocations = self.samweb.listFilesAndLocations(
        schema=schema,
        **locationArg,
        )
    except samexcpt.DefinitionNotFound as e:
      raise FatalError(f"SAM definition not found: '{inputSpec}'.")
    except samexcpt.DimensionError as e:
      raise FatalError(f"Error in SAM query '{inputSpec}': {e}.")
    
    fileURLs = list(fileInfo[0] for fileInfo in filesAndLocations)
    logger.debug(" => %d entries", len(fileURLs))
    return fileURLs
  # queryFileURLs()
  
  def prepareOutput(self, inputSpec, fileURLs) -> "output lines and # entries":
    nFiles = len(fileURLs)
    fileURLs = fileURLs[:self.options['maxFiles']]
    nEntries = len(fileURLs)
    
    outputLines = []
    if self.options['header']:
      msg = str(nEntries)
      if nEntries != nFiles: msg += "/" + str(nFiles)
      msg += " files from " + inputSpec.describe() + " expanded on " \
        + time.ctime()
      outputLines.append("# " + msg)
    # if header
    
    outputLines.extend(fileURLs)
    return outputLines, nEntries
  # prepareOutput()
  
  def writeOutput(self, inputSpec, index, outputLines):
    
    outputPath, outputFile = self.createOutputFile(inputSpec, index)
    with outputFile:
      outputFile.write('\n'.join(outputLines))
      outputFile.write('\n')
    # with
    return outputPath
  # writeOutput()
  
  def createOutputFile(self, inputSpec, index) -> "name and file":
    
    if self.options['output'] is None: return '', sys.stdout
    outputFileName \
      = self.createOutputFileName(self.options['output'], inputSpec, index)
    logger.debug("Output file: %s", outputFileName)
    if os.path.exists(outputFileName) and not self.options['overwrite']:
      raise FatalError(f"File list '{outputFileName}' already exists,"
        f" not writing result of {inputSpec.describe()}")
    dirName = os.path.dirname(outputFileName)
    if dirName:
      try:
        os.makedirs(dirName, exist_ok=True)
      except:
        raise FatalError(f"Error creating output directory '{dirName}': {e}")
    # if there is a directory involved
    return outputFileName, open(outputFileName, 'w')
  
  # createOutputFile()
  
  def createOutputFileName(self, outputFilePattern, inputSpec, index):
    
    if outputFilePattern.endswith(os.sep):
      logger.debug("Output pattern '%s' represents a directory.")
      dirPath = outputFilePattern
      outputFilePattern = ""
    else:
      dirPath = ""
    
    if not outputFilePattern:
      outputFilePattern = dirPath
      
      patternComponents = [ "%t" ]
      
      if self.options['stage']: patternComponents.append(str(self.options['stage']))
      
      if self.options['streams']:
        sanitizeStream = lambda s: s.replace("%", "")
        patternComponents.extend(map(self.removeSAMwildcards, self.options['streams']))
      # if
      
      if self.options['project']:
        patternComponents.append(self.removeSAMwildcards(self.options['project']))
      
      if 'root' in self.options['schema'].lower(): patternComponents.append('XRootD')
      elif 'posix' in self.options['schema'].lower(): patternComponents.append('POSIX')
      
      outputFilePattern += "-".join(patternComponents)
      
      if self.options['maxFiles'] is not None: outputFilePattern += "_max%m"
      
      outputFilePattern += ".filelist"
      
      logger.debug("Default pattern built: %s", outputFilePattern)
    # if
    
    
    return outputFilePattern \
      .replace('%n', str(inputSpec)) \
      .replace('%s', inputSpec.sanitized()) \
      .replace('%t', inputSpec.tag()) \
      .replace('%m', str(self.options['maxFiles'])) \
      .replace('%#', str(0 if index is None else index))
    
  # createOutputFileName()
  
  PrefixKey = '/pnfs/fnal.gov/usr/'
  
  @staticmethod
  def XRootDtoPOSIX(path):
    if not path.startswith('root:'): return path
    try:
      prefixEnd = path.find(FileListMakerClass.PrefixKey)
    except ValueError: return path
    return '/pnfs/' + path[prefixEnd+len(FileListMakerClass.PrefixKey):]
  # XRootDtoPOSIX()
  
  @staticmethod
  def removeSAMwildcards(s): return s.replace("%", "")
  
# class FileListMakerClass
FileListMakerClass.__call__ = FileListMakerClass.process


# ------------------------------------------------------------------------------
if __name__ == "__main__":
  
  import argparse
  
  parser = argparse.ArgumentParser(description=__doc__)
  parser.set_defaults(
    header=False, schema='root', stage=None, location=None,
    )
  
  parser.add_argument("inputSpecs", nargs="+", metavar="inputSpecs",
    help="run number, SAM definition name, query or file name")
  
  inputGroup = parser.add_argument_group("Input options")
  
  inputGroup.add_argument("--stage", dest="stage",
    help="query for files at the specified stage")
  inputGroup.add_argument("--raw", dest="stage", action="store_const",
    const="raw", help="query for files at raw stage")
  inputGroup.add_argument("--decoded", "--stage0", dest="stage",
    action="store_const", const="stage0", help="query for files at stage0")
  inputGroup.add_argument("--reco", "--stage1", dest="stage",
    action="store_const", const="stage1", help="query for files at stage1")
  
  inputGroup.add_argument("--streams", "--stream", action="append",
    default=[], help="query for files in any of the specified streams")
  
  inputGroup.add_argument("--project", "--prjversion", default=None,
    help="query for files processed with the specified project version")
  
  # inputGroup.add_argument("--location", dest="location",
  #   help="query for files stored in the specified location type")
  # inputGroup.add_argument("--dcache", dest="location", action="store_const",
  #   const="dcache", help="query for files in dCache")
  # inputGroup.add_argument("--tape", "--enstore", dest="location",
  #   action="store_const", const="tape", help="query for files on tape")
  
  outputGroup = parser.add_argument_group("Output options")
  outputGroup.add_argument("--output", "-o", default=None,
    help="pattern of the file(s) to write the resulting list into [stdout]",
    )
  outputGroup.add_argument("--autooutput", "-O", action="store_const", const="",
    dest='output', help="create an output file name for each query",
    )
  outputGroup.add_argument("--force", "-f", action="store_true", dest="overwrite",
    help="if output file list already exist, overwrite it [%(default)s]"
    )
  
  outputGroup.add_argument("--schema", dest="schema", default="root",
    help="which type of URL to produce [%(default)s]")
  outputGroup.add_argument("--root", dest="schema", action="store_const",
    const="root", help="output file paths in XRootD format")
  outputGroup.add_argument("--posix", dest="schema", action="store_const",
    const="posix", help="output file paths in POSIX format")
  
  outputGroup.add_argument("--header", action="store_true", dest="header",
    help="prepend a header [%(default)s]")
  outputGroup.add_argument("--no-header", action="store_false", dest="header",
    help="do not prepend a header")
  outputGroup.add_argument("--max", dest="maxFiles", default=None, type=int,
    help="limit of the maximum number of files printed [%(default)s]")
  
  
  parser.add_argument("--experiment", "-e", dest="SAMexperiment",
    default=None, help="override the SAM station name to use")
  parser.add_argument("--debug", action="store_true",
    help="prints out debugging messages")
  parser.add_argument \
    ("--version", "-V", action="version", version="%(prog)s " + __version__)

  args = parser.parse_args()
  
  logger.setLevel(logging.DEBUG if args.debug else logging.INFO)
  
  samweb = SAMWebClient(experiment=args.SAMexperiment)
  processor = FileListMakerClass(samweb=samweb, options=args)
  
  inputSpecs = list(args.inputSpecs)
  nSpecs = len(inputSpecs)
  
  nErrors = 0
  for iInputSpec, inputSpec in enumerate(map(InputSpecification.make, inputSpecs)):
    try:
      outputPath, nEntries \
        = processor.process(inputSpec, index=iInputSpec if nSpecs > 1 else None)
      if outputPath:
        logger.info("%s expanded into '%s' (%d entries)",
          inputSpec.describe(), outputPath, nEntries)
      # if
    except FatalError as e:
      logger.error("Error encountered processing '%s':\n%s", inputSpec, e)
      nErrors += 1
    # try ... except
  # for
  
  if nErrors > 0:
    logger.critical("%d/%d requests processed with errors.", nErrors, nSpecs)
    sys.exit(1)
  # if errors
  
  if nSpecs > 1:
    logger.info("%d requests successfully processed.", nSpecs)
  
  sys.exit(0)
  
# main
