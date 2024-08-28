#!/usr/bin/env python

import logging
import time

__author__ = "Gianluca Petrillo (petrillo@slac.stanford.edu)"
__date__ = time.strptime("July 7, 2021", "%B %d, %Y")
__version__ = "1.1"
__doc__ = """
Manages SAM definitions for ICARUS data run.

This script can query SAM database for the files in a run, or create and manage
dataset definitions. The convenience lies in the uniform format for the
definition names, so that they are predictable and easy to encode in programs.

The program selects specific characteristics (stream, stage, project version) of
the file metadata. If no selection is provided, it will have one definition per
stage and per project version, while the rest of the metadata will be ignored.

Examples:

1. See how many BNB majority raw files (both off-beam and on-beam) are available
   for runs 11803 and 11806:
    
    %(prog)s  --query --stream=%%bnbmajority --stage=raw 11803 11806
    
2. Create two new datasets with all stage0 files of runs 11803 and 11806:
    
    %(prog)s  --create --stage=stage0 11803 11806
    
3. Have I already created a definition for stage1 for those runs?
    
    %(prog)s  --query --stage=stage1 11803 11806
    

"""

from samweb_client.client import SAMWebClient
import samweb_client.exceptions as samexcpt

logging.basicConfig(format="[%(levelname)s] %(message)s")

ExperimentName = "ICARUS"
SAMExperimentName = ExperimentName.lower()
DefaultStages = [ 'raw', 'stage0', 'stage1', ]

StageDimensions = {
  'raw':    "data_tier raw",
  'stage0': "icarus_project.stage stage0",
  'stage1': "icarus_project.stage in ( stage1, stage1_caf_larcv )",
} # StageDimensions

# this is a modal flag that is dangerous enough to be not user-controlled:
# if `True`, queries without specifying a run number will be allowed.
AllowAllRuns = False


# ------------------------------------------------------------------------------
class SampleInfo:
  def __init__(self,
    run=None, stage=None, stream=None, projectVersion=None,
    minSize: "MiB" = None,
    ):
    self.run = SampleInfo._copyList(run)
    self.stage = SampleInfo._copyList(stage)
    self.stream = SampleInfo._copyList(stream)
    self.projectVersion = SampleInfo._copyList(projectVersion)
    self.minSize = minSize if minSize else 0
  # __init__()
  
  def isRunDefined(self) -> "returns if run is collapsed to a value":
    return not SampleInfo.hasOptions(self.run) \
      and (AllowAllRuns or self.run is not None)
  # isRunDefined()
  
  def isStageComplete(self) -> "returns if stage is collapsed to a value":
    return not SampleInfo.hasOptions(self.stage)
  
  def isStreamComplete(self) -> "returns if stream is collapsed to a value":
    return not SampleInfo.hasOptions(self.stream)
  
  def isProjectVersionComplete(self) -> "returns if project version is collapsed to a value":
    return not SampleInfo.hasOptions(self.projectVersion)
  
  def isComplete(self) -> "returns if there are no ambiguities left in the info":
    return self.isRunDefined() and self.isStageComplete() \
      and self.isStreamComplete() and self.isProjectVersionComplete()
  # isComplete()
  
  def defName(self):
    components = [ ExperimentName, 'data' ]
    if self.isRunDefined() and self.run is not None:
      components.append(f"run{self.run}")
    if self.isStageComplete() and self.stage is not None:
      components.append(self.stage)
    if self.isStreamComplete() and self.stream is not None:
      components.append(self.stream.replace("%", "")) # remove SAM wildcards
    if self.isProjectVersionComplete() and self.projectVersion is not None:
      components.append(self.projectVersion)
    return "_".join(filter(None, components))
  # defName()
  
  def copy(self, **kwargs):
    kwargs.setdefault('run', self.run)
    kwargs.setdefault('stage', self.stage),
    kwargs.setdefault('stream', self.stream),
    kwargs.setdefault('projectVersion', self.projectVersion),
    kwargs.setdefault('minSize', self.minSize),
    return SampleInfo(**kwargs)
  # copy()
  
  def __str__(self):
    s = f"run {self.run if self.run else 'unspecified'}"
    if self.stage: s += f", stage {self.stage}"
    if self.stream: s += f", stream {self.stream}"
    if self.projectVersion: s += f", project version {self.projectVersion}"
    return s
  # __str__()
  
  @staticmethod
  def makeOptionList(value):
    return value if SampleInfo.hasOptions(value) else [ value ]
  
  @staticmethod
  def hasOptions(value):
    return SampleInfo._isIterable(value) and not isinstance(value, str)
  
  @staticmethod
  def _isIterable(value): return hasattr(value, "__iter__")
  
  @staticmethod
  def _copyList(value):
    return value[:] if SampleInfo._isIterable(value) else value
  
# class SampleInfo


# ------------------------------------------------------------------------------
class DimensionQueryMaker:
  """Callable object making a SAM query out of a `SampleInfo` object."""
  
  def __call__(self,
   info: "SampleInfo object with all the information for the query",
   minimum: "throws an exception if fewer than these elements are included" = 0,
   ) -> "a SAM dimensions query string":
    
    dims = [
      DimensionQueryMaker.simpleItem('run_number', info.run),
      DimensionQueryMaker.multiItem(StageDimensions, info.stage, 'stage'),
      DimensionQueryMaker.simpleItem('data_stream', info.stream),
      DimensionQueryMaker.simpleItem('icarus_project.version', info.projectVersion),
      ]
    if info.minSize > 0:
      dims.append(
        DimensionQueryMaker.comparedItem('file_size', info.minSize << 20, ">=")
        )
    query = " and ".join(filter(None, dims))
    if len(dims) < minimum:
      raise RuntimeError(f"Query resulted in only {len(dims)} constraints: '{query}'")
    return query
  # __call__()
  
  @staticmethod
  def addParentheses(s): return "( " + s + " )" if s else ""

  @staticmethod
  def simpleItem(key, value, typeName=None):
    # if typeName is None: typeName = key
    if SampleInfo.hasOptions(value):
      return f"{key} in ( { ' '.join(map(str, value)) } )"
    elif value is not None: return f"{key} {value}"
    else: return ""
  # simpleItem()


  @staticmethod
  def comparedItem(key, value, compType=None):
    if value is not None: return f"({key} {compType} {value})"
    else: return ""
  # comparedItem()


  @staticmethod
  def multiItem(queries, values, typeName):
    values = SampleInfo.makeOptionList(values)
    if not values or None in values: return ""
    dims = []
    for value in values:
      try:
        dims.append(queries[value])
      except KeyError:
        raise RuntimeError(
          f"{typeName} '{value}' not supported (expected one of: {', '.join(map(str, queries))})."
          )
      # try ... except
    # for
    if len(dims) > 1: dims = list(map(DimensionQueryMaker.addParentheses, dims))
    query = " or ".join(dims) if dims else ""
    if len(dims) > 1: query = DimensionQueryMaker.addParentheses(query)
    return query
  # multiItem()

# class DimensionQueryMaker


# ------------------------------------------------------------------------------
class SampleBrowser:
  
  def __init__(self, samweb):
    self.samweb = samweb if samweb else SAMWebClient()
  
  
  def iterateProjectVersions(self, info: "SampleInfo object defining iteration ranges"):
    assert self.samweb, "SAM web client not initialized. We do not go anywhere."
    
    forceSeparateVersions =  (info.projectVersion == '%')
    
    # expand the project versions: each version will be treated separately
    if info.projectVersion is None or forceSeparateVersions:
      if info.stage != 'raw' or forceSeparateVersions:
        versionQuery = DimensionQueryMaker()(info, minimum=2)
        projectVersions = self._discoverProjectVersions(versionQuery)
      else: projectVersions = [ None ] # no version constraint is ok for raw files
    else:
      projectVersions = SampleInfo.makeOptionList(info.projectVersion)
    # if ... else
    
    for prjVer in projectVersions:
      prjInfo = info.copy(projectVersion=prjVer)
      yield prjInfo
    # for
  # iterateProjectVersions()


  def iterateStreams(self, info: "SampleInfo object defining iteration ranges"):
    assert info.isStageComplete(), "Stage must be set."
    for stream in SampleInfo.makeOptionList(info.stream):
      streamInfo = info.copy(stream=stream)
      yield from self.iterateProjectVersions(streamInfo)
    # for streams
  # iterateStreams()


  def iterateStages(self, info: "SampleInfo object defining iteration ranges"):
    assert info.stage is not None, "Stages needs to be explicitly specified."
    for stage in SampleInfo.makeOptionList(info.stage):
      assert stage, "Stage must be specified."
      stageInfo = info.copy(stage=stage)
      logging.debug(f"Processing {stageInfo}")
      yield from self.iterateStreams(stageInfo)
    # for
  # iterateStages()

  
  def iterate(self, info: "SampleInfo object defining iteration ranges"):
    """Iterates through all elements in `info`."""
    for run in SampleInfo.makeOptionList(info.run):
      runInfo = info.copy(run=run)
      yield from self.iterateStages(runInfo)
    # for
  # iterate()
  def __call__(self, info): return self.iterate(info)
  
  
  def _discoverProjectVersions(self, dims):
    """Returns the project versions from the all files matching `dims`."""
    logging.debug("Discovering project versions for %s", dims)
    
    try:
      files = self.samweb.listFiles(dims, fileinfo=True)
    except samexcpt.Error:
      logging.error(f"SAM exception while translating query: '{dims}'")
      raise
    
    logging.debug(f"  => querying metadata for {len(files)} files")
    
    versions = set(
      meta['icarus_project.version']
      for meta in self.samweb.getMetadataIterator(fInfo.file_id for fInfo in files)
      )
    logging.debug(f"  => extracted {len(versions)} versions: %s", versions)
    return list(versions)
  # _discoverProjectVersions()

# class SampleBrowser


# ------------------------------------------------------------------------------
class SampleProcessClass:
  
  def __init__(self, samweb, create=False, query=False, printDefs=False,
   describe=False, check=False, delete=False,
   fake=False, force=False, prependUser=True, minSize=None,
   ):
    # action collection
    self.actions = []
    if check: self.actions.append("check")
    if create: self.actions.append("create")
    if query: self.actions.append("query")
    if printDefs: self.actions.append("printDefs")
    if describe: self.actions.append("describe")
    if delete: self.actions.append("delete")
    if not self.actions:
      raise RuntimeError("At least one action needs to be enabled.")
    self.fake = fake
    self.force = force
    self.prependUser = prependUser
    self.minSize = minSize if minSize else 0
    
    self.samweb = samweb if samweb else SAMWebClient()
    self.buildQuery = DimensionQueryMaker()
    
    try: self.SAMuser = samweb.get_user()
    except samexcpt.Error as e:
      if self.prependUser:
        logging.error("Could not find out your name! %s", e)
        raise
      self.SAMuser = None
    #
    
  # __init__()
  
  
  def __call__(self, info):
    assert self.samweb, "SAM web client not initialized. We do not go anywhere."
    
    if not info.isComplete():
      raise RuntimeError(f"Can't process incomplete specification: {info}")
    
    dim = self.buildQuery(info)
    
    logging.debug("Info: %s => dim='%s'", info, dim)
    defName = self.buildDefName(info)
    
    for action in self.actions:
    
      if action == "printDefs":
        print(defName)
      
      elif action == "check":
        self.doCheck(defName)
      
      elif action == "describe":
        self.doDescribe(defName)
      
      elif action == "query":
        self.doQuery(info=info, defName=defName, dim=dim)
    
      elif action == "create":
        self.doCreateDef(info, defName=defName, dim=dim)
      
      elif action == "delete":
        self.doDeleteDef(info, defName=defName, dim=dim)
      
      else:
        raise RuntimeError(f"LOGIC ERROR: action {action} not implemented.")
      
    # for
    
    return defName
  # __call__()
  
  
  def doDescribe(self, defName):
    if self.fake:
      print(f"DRYRUN> descDefinition({defName!r})")
      return
    try:
      print(self.samweb.descDefinition(defName))
    except samexcpt.DefinitionNotFound:
      print(f"Definition Name: {defName}  => NOT FOUND")
  # doDescribe()
  
  def doQuery(self, info, defName, dim):
    if self.fake:
      print(f"DRYRUN> {dim}")
      return None
    try:
      summary = self.getSummary(info=info, defName=defName, dims=dim)
      queryError = None
    except samexcpt.Error as e:
      logging.error(f"Query of definition {defName} (query: '{dim}') failed: %s", e)
      summary = None
      queryError = e
    else:
      msg = str(info) + ":"
      msg +=  " " + (
        str(summary['total_event_count']) if summary['total_event_count']
          else "unknown number of"
        )
      msg += f" events in {summary['file_count']} files"
      msg += (
        " (0 GiB)" if summary['total_file_size'] is None
        else f" ({summary['total_file_size']/(1 << 30):g} GiB)"
        )
      print(msg)
      e = None
    # try ... else
    return summary if summary else queryError
  # doQuery()
  
  def doCheck(self, defName):
    assert defName
    if self.fake:
      print(f"DRYRUN> countFiles(defname={defName!r})")
      return True
    count = self.defNameCount(defName)
    if count is None:
      print(f"{defName} not available")
      return False
    else:
      print(f"{defName} available ({count} files)")
      return True
  # doCheck()
  
  def doCreateDef(self, info, defName=None, dim=None):
    assert info
    if not defName: defName = self.buildDefName(info)
    if not dim: dim = self.buildQuery(info)
    descr = self.describeSample(info)
    
    count = self.defNameCount(defName)
    if count is not None:
      logging.error(f"Definition {defName!r} already exists (and matches {count} files).")
      return None
    
    try:
      count = self.samweb.countFiles(dimensions=dim)
    except samexcpt.Error as e:
      logging.error(f"Attempt to count matches with {dim!r} failed: %s", e)
      return None
    if count == 0:
      print(f"Definition {defName} NOT created as it would match no file (query: {dim!r})")
      return None
    logging.debug(f"Creating {defName!r}, now matching {count} files")
    if self.fake:
      print(f"DRYRUN> createDefinition(defname={defName!r}, dims={dim!r}, description={descr!r})")
    else:
      try:
        self.samweb.createDefinition(defname=defName, dims=dim, description=descr)
        print(f"{defName} created ({count} files)")
      except samexcpt.Error as e:
        logging.error \
          (f"Failed to create definition {defName} from query='{dim}': %s", e)
        return None
    # if
    return defName
  # doCreateDef()
  
  def doDeleteDef(self, info, defName=None, dim=None):
    assert self.samweb
    assert info
    if not defName: defName = self.buildDefName(info)
    if not dim: dim = self.buildQuery(info)
    descr = self.describeSample(info)
    
    try:
      defInfo = self.samweb.descDefinitionDict(defName)
    except samexcpt.DefinitionNotFound:
      print(f"Definition Name: {defName}  => NOT FOUND") # we are ok... I guess
      return None
    
    # we perform many checks to make sure this deletion is proper
    ForcedMsg = { True: "forced to delete it anyway", False: "won't delete unless forced to", }
    checksOk = True
    
    if not self.SAMuser:
      try: self.SAMuser = self.samweb.get_user()
      except samexcpt.Error as e:
        logging.error("Could not find out your name! %s", e)
    # if not cached already
    try: SAMgroup = self.samweb.get_group()
    except samexcpt.Error as e:
      logging.error("Could not find out the name of your group! %s", e)
    logging.debug(f"You appear to be {SAMuser!r} of group {SAMgroup!r}")
    
    if defInfo['username'] != SAMuser:
      logging.warning(
        f"Definition {defName!r} was created on {defInfo['create_time']}"
        f" by {defInfo['username']}/{defInfo['group']}, not by you ({SAMuser})"
        f": won't delete."
        )
      checksOk = False
    # if
    
    if defInfo['group'] != SAMgroup:
      logging.warning(
        f"Definition {defName!r} was created on {defInfo['create_time']}"
        f" by {defInfo['username']}/{defInfo['group']}, not by your group ({SAMgroup})"
        f": won't delete."
        )
      checksOk = False
    # if
    
    if defInfo['dimensions'] != dim:
      logging.warning(f"Definition {defName!r} has unexpected query:"
        f" ({defInfo['dimensions']!r}, expected: {dim!r}); {ForcedMsg[self.force]}."
        )
      if not self.force: checksOk = False
    #
    
    if defInfo['description'] != descr:
      logging.warning(
        f"Definition {defName!r} appears not to be created with this program:"
        f" description mismatch"
        f" ({defInfo['description']!r}, expected: {descr!r}); {ForcedMsg[self.force]}."
        )
      if not self.force: checksOk = False
    # if
    
    if not checksOk:
      logging.error(f"Definition {defName!r} will NOT be deleted.")
      return None
    
    if self.fake:
      print(f"DRYRUN> deleteDefinition({defName!r})")
      return defName
    try:
      self.samweb.deleteDefinition(defName)
    except samexcpt.DefinitionNotFound:
      logging.error(f"Definition {defName} NOT FOUND (can't be deleted)")
    # except samexcpt.DefinitionNotFound: # which exception to match?
    #   logging.error(f"Definition {defName!r} has already been used and can't be deleted.")
    except samexcpt.NoAuthorizationCredentials as e:
      logging.error(f"Failed to delete definition {defName!r} for lack of credentials: %s", e)
    except samexcpt.Error as e:
      logging.error(f"Failed to delete definition {defName!r}: %s", e)
    
    count = self.defNameCount(defName)
    if count is not None:
      logging.error(f"Deletion of definition {defName!r} silently FAILED"
                    f" (still there with its own {count} files).")
      return None
    # if still there
    print(f"Definition {defName} successfully deleted.")
    return defName
  # doDeleteDef()
  
  
  def buildDefName(self, info):
    if self.prependUser: return self.SAMuser + "_" + info.defName()
    else: return info.defName()
  # buildDefName()
  
  def getSummary(self, info=None, defName=None, dims=None):
    assert info or defName or dims
    e = RuntimeError("Insufficient parameters to get summary") # should not happen
    
    if not defName and info: defName = self.buildDefName(info)
    if defName:
      try: return self.samweb.listFilesSummary(defname=defName)
      except samexcpt.DefinitionNotFound as e: queryError = e
    
    if not dims and info: dims = self.buildQuery(info)
    if dims:
      try: return self.samweb.listFilesSummary(dims)
      except samexcpt.Error as e: queryError = e
    
    raise queryError
  # getSummary()
  
  
  def defNameCount(self, defName: "SAM definition name"):
    """Returns the count of files of `defName`, `None` if not found.
    
    Throws exception in all other error situations.
    """
    try: return self.samweb.countFiles(defname=defName)
    except samexcpt.DefinitionNotFound: return None
  # defNameCount()
  
  def describeSample(self, info): return "ICARUS data " + str(info)
  
# class SampleProcessClass


# ------------------------------------------------------------------------------
def collapseList(l):
  try:
    if not isinstance(l, str) and len(l) == 1: return next(iter(l))
  except TypeError: pass
  return l
# collapseList()


if __name__ == "__main__":
  import sys
  import argparse
  
  parser = argparse.ArgumentParser(description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter)
  
  SampleGroup = parser.add_argument_group(title="Sample selection")
  SampleGroup.add_argument("runs", nargs="*" if AllowAllRuns else "+", type=int,
                      help="runs to process")
  SampleGroup.add_argument("--stage", "-s", action="append",
    help=f"stages to include {DefaultStages}")
  SampleGroup.add_argument("--prjversion", "-p", action="append",
    help="project versions to include [autodetect (resource-intensive!)]")
  SampleGroup.add_argument("--stream", "-f", action="append",
    help="data streams to include (use 'any' for... any) [any]")
  SampleGroup.add_argument("--global", "-g", dest='globalDef',
    action="store_true", help="do not prepend SAM user name to definitions")
  SampleGroup.add_argument("--minsize", "-S", dest='minimumSize',
    action='store', nargs='?', default=None, const=8, type=int,
    help="include only files larger than this size (MiB) [8 MiB if no value]")
  
  ActionGroup = parser.add_argument_group(title="Actions")
  ActionGroup.add_argument("--check", action="store_true",
    help="prints whether the definition for the sample is available")
  ActionGroup.add_argument("--describe", action="store_true",
    help="describes an existing definition for the sample")
  ActionGroup.add_argument("--query", action="store_true",
    help="queries the definitions related to the samples")
  ActionGroup.add_argument("--defname", action="store_true",
    help="prints the name of the definitions related to the samples")
  ActionGroup.add_argument("--create", action="store_true",
    help="creates one definition per sample (use --defname to see their names)")
  ActionGroup.add_argument("--delete", action="store_true",
    help="attempts to remove one definition per sample")
  
  GeneralOptGroup = parser.add_argument_group(title="General options")
  GeneralOptGroup.add_argument("--experiment", "-e", default=ExperimentName,
    help="sets the experiment name (for definitions) [%(default)s]")
  GeneralOptGroup.add_argument("--samexperiment", "-E",
    help="sets the experiment name (chooses SAM database) [same as --experiment]")
  GeneralOptGroup.add_argument("--auth", dest="AuthMode", choices=[ 'token', 'cert' ],
    default='token', help="choose authentication method for SAM [%(default)s]")
  GeneralOptGroup.add_argument("--force", "-F", action="store_true",
    help="skips safety checks of some operations")
  GeneralOptGroup.add_argument("--fake", "--dryrun", "-n", action="store_true",
    help="does not perform actual creation and query actions")
  GeneralOptGroup.add_argument("--debug", action="store_true",
    help="enable verbose debugging output")
  GeneralOptGroup.add_argument("--version", "-V", action="version",
    version="%(prog)s v" + __version__, help="prints the version number")
  
  args = parser.parse_args()
  
  logging.getLogger().setLevel(logging.DEBUG if args.debug else logging.INFO)
  if args.stage is None: args.stage = DefaultStages
  if args.stream:
    args.stream = [ None if s == "any" else s for s in args.stream ]
  ExperimentName = args.experiment
  SAMexperiment = args.samexperiment if args.samexperiment else ExperimentName.lower()

  logging.debug("Using SAM database '%s'.", SAMexperiment)
  samweb = SAMWebClient(
    experiment=SAMexperiment,
    disable_cert_auth=(args.AuthMode != 'cert'),
    disable_token_auth=(args.AuthMode != 'token'),
    )
  iterateSamples = SampleBrowser(samweb)
  processInfo = SampleProcessClass(samweb,
    create=args.create, query=args.query, printDefs=args.defname,
    describe=args.describe, check=args.check, delete=args.delete,
    fake=args.fake, force=args.force,
    prependUser=not args.globalDef,
    )
  
  baseSampleInfo = SampleInfo(
    run=collapseList(args.runs if args.runs else None),
    stage=collapseList(args.stage),
    stream=collapseList(args.stream), # None means any
    projectVersion=collapseList(args.prjversion), # None for autodetect
    minSize=args.minimumSize,
    )
  
  for sampleInfo in iterateSamples(baseSampleInfo):
    processInfo(sampleInfo)
  
  sys.exit(0)
# main
