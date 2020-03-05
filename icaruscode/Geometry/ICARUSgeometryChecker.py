#!/usr/bin/env python

__doc__ = """
Performs simple checks on ICARUS geometry.
"""
__version__ = "%(prog)s 1.0"

import itertools
import logging


# ##############################################################################
class NestedIteration:
  class BadIterator:
    def __iter__(self): return self
    def next(self): raise StopIteration
  # class BadIterator
    
  def __init__(self, iterable):
    self.outerIter = iter(iterable)
    self.innerIter = NestedIteration.BadIterator()
  # __init__()
  
  def __iter__(self): return self
  
  def next(self):
    while True:
      try: return next(self.innerIter)
      except StopIteration: pass
      self.innerIter = iter(next(self.outerIter)) # let StopIteration through
    # while
  # next()
  
# class NestedIteration


# ##############################################################################
def GeoPairToString(pair):
  return " - ".join(map(lambda obj: str(obj.ID()), pair))

def BoxToString(box): return "%s -- %s" % (box.Min(), box.Max())


# ------------------------------------------------------------------------------
def boxID(box, default):
  try: return box.ID()
  except AttributeError: return default
# boxID()


# ------------------------------------------------------------------------------
def CheckBoxOverlaps(boxes, objName = None):
  """Returns a list of pairs of overlapping cryostats (empty if no overlap)."""
  if objName is None: objName = boxes.__class__.__name__
  overlaps = []
  for (iBox1, box1), (iBox2, box2) in itertools.combinations(enumerate(boxes), 2):
    if not box2.Overlaps(box1): continue
    logging.error("%s %s (%s to %s) and %s (%s to %s) overlap!", objName,
     boxID(box1, iBox1), box1.Min(), box1.Max(),
     boxID(box2, iBox2), box2.Min(), box2.Max(),
     )
    overlaps.append( ( box1, box2, ) )
  # for boxes
  return overlaps
# CheckBoxOverlaps()


# ------------------------------------------------------------------------------
def checkPlaneAlignment(planes, tolerance = 0.1):
  """Check alignment in z position of planes on the same x.
  
  Returns triples (distance, first plane, second plane) for each misaligned
  pair of planes.
  """
  assert tolerance >= 0.0
  PlanesByX = groupPlanesByX(planes, sortBy = 'z')
  inconsistentPlanes = {}
  for planesOnX in PlanesByX:
    if len(planesOnX) < 2: continue
    planeIter = iter(planesOnX)
    refPlane = next(planeIter)
    for nextPlane in planeIter:
      refBox = refPlane.BoundingBox()
      nextBox = nextPlane.BoundingBox()
      #
      # 1. check plane alignment on z direction
      #
      distance = nextBox.MinZ() - refBox.MaxZ()
      if abs(distance) > tolerance:
        logging.error(
         "Planes on x=%s are misaligned along z: ( %g -- %g ) (%s) vs. ( %g -- %g ) (%s)",
         refPlane.GetCenter().X(), refBox.MinZ(), refBox.MaxZ(), refPlane.ID(),
         nextBox.MinZ(), nextBox.MaxZ(), nextPlane.ID(),
         )
        inconsistentPlanes.setdefault('zalign', []) \
         .append( ( distance, refPlane, nextPlane, ) )
      # if misaligned on z
      
      #
      # 2. check wire view
      #
      if refPlane.View() != nextPlane.View():
        logging.error("Plane %s is on view '%s', %s is on view '%s'.",
         refPlane.ID(), ROOT.geo.PlaneGeo.ViewName(refPlane.View()),
         nextPlane.ID(), ROOT.geo.PlaneGeo.ViewName(nextPlane.View()),
         )
        inconsistentPlanes.setdefault('view', []) \
         .append( ( refPlane, nextPlane, ) )
      # if views do not match
      
      #
      # 3. check wire orientation
      #
      if 1.0 - refPlane.GetIncreasingWireDirection().Dot(nextPlane.GetIncreasingWireDirection()) > tolerance:
        logging.error("Plane %s measures direction %s, %s measures %s.",
         refPlane.ID(), refPlane.GetIncreasingWireDirection(),
         nextPlane.ID(), nextPlane.GetIncreasingWireDirection(),
         )
        inconsistentPlanes.setdefault('direction', []) \
         .append( ( refPlane, nextPlane, ) )
      # if views do not match
      
      refPlane = nextPlane
    # for all following planes on the same x
  # for groups of planes along the same x
  return inconsistentPlanes
# checkPlaneAlignment()

  
# ------------------------------------------------------------------------------
def checkPlaneWireAlignment(planeA, planeB, tolerance = 0.01):
  """
  Wires in the plane A are verified to have a continuation in plane B as another
  wire. Only wires that reach the side (extreme z coordinate) are tested.
  
  Wires are assumed to have to go top to bottom on y direction.
  
  The list of misaligned wires is returned in the 5-uple form:
      
      ( shift, left wire ID, left wire, right wire ID, right wire )
      
  if the wire on the left plane was matched to one on the right plane, and
      
      ( None, left wire ID, left wire, None, None )
      
  if the wire on the left plane was not matched to any wire on the right plane.
  """
  # throughout this code, left refers to lower _z_, right to higher _z_
  
  misalignedWires = []
  
  # define a left and a right plane
  if planeA.GetCenter().Z() < planeB.GetCenter().Z():
    leftPlane, rightPlane = planeA, planeB
  else:
    leftPlane, rightPlane = planeB, planeA
  
   # wire pitch should be the same for the two planes
  rightWirePitch = rightPlane.WirePitch()
  assert abs(rightWirePitch - leftPlane.WirePitch()) < tolerance
  
  #
  # 1. determine the sides
  #
  leftMaxZ = max(
   max(wire.GetStart().Z(), wire.GetEnd().Z())
   for wire in leftPlane.IterateWires()
   )
  # RminZ = min( min(wire.GetStart().Z(), wire.GetEnd().Z()) for wire in rightPlane.IterateWires() )
  
  #
  # 2. for each wire in left plane ending close to right side, pick and test
  #    the matching wire of the right plane
  #
  stats = StatCollector()
  wireNoDiff = None # offset in wire number only between matching wires
  for leftWireNo, leftWire in enumerate(leftPlane.IterateWires()):
    
    leftEnd \
      = leftWire.GetStart() if leftWire.GetStart().Z() > leftWire.GetEnd().Z() \
        else leftWire.GetEnd()
    
    # 
    # 2.1. if the wire does not end on the right edge, move on
    # 
    if abs(leftEnd.Z() - leftMaxZ) > 0.01: continue
    
    #
    # 2.2. find the closest wire on the left 
    #
    leftWireID = ROOT.geo.WireID(leftPlane.ID(), leftWireNo)
    try:
      rightWireID = rightPlane.NearestWireID(leftEnd)
    except ROOT.geo.InvalidWireError, e:
      logging.error(
       "No wire on %s is close enough to %s (closest is %s, but would have been %s)", 
       rightPlane.ID(), leftWire.ID(),
       (e.suggestedWireID() if e.hasSuggestedWire() else "unknown"),
       (e.badWireID() if e.hasBadWire() else "unknown"),
       )
      misalignedWires.append( ( None, leftWireID, leftWire, None, None, ) )
      continue
    # try ... except no wire matched
    rightWire = rightPlane.Wire(rightWireID)
    
    #
    # 2.3. check the distance of that wire from this one
    #
    shift = leftWire.DistanceFrom(rightWire)
    stats.add(shift)
    
    #
    # 2.4. compare wire orientation
    #
    leftWire.Direction().Dot(rightWire.Direction())
    if 1.0 - leftWire.Direction().Dot(rightWire.Direction()) > tolerance:
      logging.error(
       "Wire %s has direction %s, the matched %s has %s",
       leftWireID, leftWire.Direction(),
       rightWireID, rightWire.Direction(),
       )
    # if too far
    
  # for
  
  logging.debug("Shift for %d wires between %s and %s: %g +/- %g cm",
   stats.entries(), leftPlane.ID(), rightPlane.ID(),
   stats.average(), stats.RMS(),
   )
  
  return misalignedWires
# checkPlaneWireAlignment()


# ------------------------------------------------------------------------------
def checkWireAlignment(planes, tolerance = 0.1):
  """Check alignment in z position of planes on the same x.
  
  Returns triples (distance, first plane, second plane) for each misaligned
  pair of planes.
  """
  assert tolerance >= 0.0
  PlanesByX = groupPlanesByX(planes, sortBy = 'z')
  misalignedWires = [] # grouped by extended plane
  for planesOnX in PlanesByX:
    if len(planesOnX) < 2: continue
    misalignedWiresOnPlane = []
    planeIter = iter(planesOnX)
    leftPlane = next(planeIter)
    for rightPlane in planeIter:
      misalignedWiresOnPair \
       = checkPlaneWireAlignment(leftPlane, rightPlane, tolerance=tolerance)
      if misalignedWiresOnPair:
        misalignedWiresOnPlane.extend(misalignedWiresOnPair)
      leftPlane = rightPlane
    # for all following planes on the same x
    if misalignedWiresOnPlane: misalignedWires.append(misalignedWiresOnPlane)
  # for groups of planes along the same x
  return misalignedWires
# checkWireAlignment()

  
# ------------------------------------------------------------------------------
class SimpleProximityClusterer:
  def __init__(self, keyFunc = None, tolerance = None):
    self.keyFunc = (lambda obj: obj) if keyFunc is None else keyFunc
    self.tolerance = tolerance
  # __init__()
  def __call__(self, objects):
    groups = []
    refKey = None
    for obj in objects:
      key = self.keyFunc(obj)
      withPrevious = refKey is not None and self._compareKeys(key, refKey)
      if not withPrevious:
        refKey = key
        groups.append(list())
      groups[-1].append(obj)
    # for objects
    return groups
  # __call__()
  def _compareKeys(self, a, b):
    return (abs(b - a) < self.tolerance) if self.tolerance else (a == b)
  @staticmethod
  def cluster(objs, keyFunc = None, tolerance = None):
    return SimpleProximityClusterer(keyFunc=keyFunc, tolerance=tolerance)(objs)
# class SimpleProximityClusterer

  
def groupPlanesByX(planes, tolerance = 0.1, sortBy = None):
  xPos = lambda plane: plane.GetCenter().X()
  cluster = SimpleProximityClusterer(xPos, tolerance) # 1 mm
  groupedByX = cluster(sorted(planes, None, xPos))
  if sortBy:
    if   sortBy.lower() == 'x':
      sortKey = xPos
    elif sortBy.lower() == 'y':
      sortKey = lambda plane: plane.GetCenter().Y()
    elif sortBy.lower() == 'z':
      sortKey = lambda plane: plane.GetCenter().Z()
    else:
      raise RuntimeError("Unsupported sorting requested: '%s'" % sortBy)
    for planes in groupedByX: planes.sort(None, sortKey)
  # if sorting
  return groupedByX
# groupPlanesByX()


# ------------------------------------------------------------------------------
class AccumulateExtrema:
  def __init__(self): self.reset()
  def reset(self):
    self.min_ = None
    self.max_ = None
    self.n = 0
  # reset()
  def add(self, value):
    if self.min_ > value or self.min_ is None: self.min_ = value
    if self.max_ < value: self.max_ = value # None < any number
  # add()
  def min(self): return self.min_
  def max(self): return self.max_
  def minmax(self): return self.min(), self.max()
# class AccumulateExtrema


# ------------------------------------------------------------------------------
class StatCollector:
  def __init__(self): self.reset()
  def reset(self):
    self.n = 0
    self.w = 0.0
    self.wx = 0.0
    self.wx2 = 0.0
  # reset()
  def add(self, value, weight = 1.0):
    self.n += 1
    self.w += weight
    self.wx += weight * value
    self.wx2 += weight * value**2
  # add()
  def entries(self): return self.n
  def weightSum(self): return self.w
  def sum(self): return self.wx
  def sumSq(self): return self.wx2
  def average(self): return self.wx / self.w if self.w else None
  def averageSq(self): return self.wx2 / self.w if self.w else None
  def variance(self): return self.averageSq() - self.average()**2
  def RMS(self): return self.variance() ** 0.5
# class StatCollector


# ##############################################################################
def performGeometryChecks(argv):
  
  import argparse
  
  Parser = argparse.ArgumentParser(description=__doc__)
  
  Parser.set_defaults()
  
  Parser.add_argument("config", help="configuration file [FHiCL]")
  
  Parser.add_argument("--servicetable", "-T", default='services',
   help="name of the FHiCL table where all services are configured ['%(default)s']"
   )
  
  Parser.add_argument("--debug", "-d", action="count", default=0,
   help="be more verbose")
  
  Parser.add_argument("--version", "-V", action="version", version=__version__)
  
  #argGroup = Parser.add_argument_group(title="Numerical output format")
  #argGroup.add_argument("--integer", "-d", "-i", action="store_false",
    #dest="bFloat", help="sets the sum of integral numbers")
  
  args = Parser.parse_args(args=argv[1:])
  
  # `logging.DEBUG` is the standard debugging level; 
  # `logging.DEBUG + 1` (args.debug = 0) is a bit less verbose, but still well
  #    within `logging.INFO` level
  logging.basicConfig(level=logging.DEBUG - (args.debug - 1))
  
  from ICARUSservices import ServiceManager
  import ROOT
  import ROOTutils
  
  global ROOT

  ServiceManager.setConfiguration(args.config, args.servicetable)
  
  geom = ServiceManager('Geometry')
  
  FailureSummary = []
  
  Cryostats = list(geom.IterateCryostats())
  TPCs = list(geom.IterateTPCs())
  Planes = list(geom.IteratePlanes())
  
  #
  # cryostat overlap
  #
  overlappingCryostats = CheckBoxOverlaps(Cryostats, "cryostat")
  if overlappingCryostats:
    msg = "%s cryostat overlaps detected: %s." % (
     len(overlappingCryostats),
     ", ".join(map(GeoPairToString, overlappingCryostats)),
     )
    logging.error(msg)
    FailureSummary.append(msg)
  else: logging.info("No cryostat overlap detected.")
  
  #
  # TPC overlaps
  #
  overlappingTPCs = CheckBoxOverlaps(TPCs, "TPC")
  if overlappingTPCs:
    msg = "%s TPC overlaps detected: %s." % (
     len(overlappingTPCs), ", ".join(map(GeoPairToString, overlappingTPCs)),
     )
    logging.error(msg)
    FailureSummary.append(msg)
  else: logging.info("No TPC overlap detected.")
  
  #
  # TPC active volume overlaps
  #
  overlappingActiveVolTPCs = CheckBoxOverlaps \
    (map(ROOT.geo.TPCGeo.ActiveBoundingBox, TPCs), "active TPC volume")
  if overlappingActiveVolTPCs:
    msg = "%s TPC active volume overlaps detected: %s." % (
     len(overlappingActiveVolTPCs),
     ", ".join(map(GeoPairToString, overlappingActiveVolTPCs)),
     )
    logging.error(msg)
    FailureSummary.append(msg)
  else: logging.info("No TPC active volume overlap detected.")
  
  #
  # plane box overlaps
  #
  overlappingPlanes = CheckBoxOverlaps \
   (map(ROOT.geo.PlaneGeo.BoundingBox, Planes), "wire planes")
  if overlappingPlanes:
    logging.error("%s wire plane overlaps detected: %s.",
     len(overlappingPlanes),
     ", ".join(map(GeoPairToString, overlappingPlanes))
     )
    FailureSummary.append \
      ("%s wire plane overlaps detected." % len(overlappingPlanes))
  else: logging.info("No wire plane overlap detected.")
  
  #
  # check plane alignment
  #
  inconsistentPlanes = checkPlaneAlignment(Planes, tolerance=0.0001) # 1 um (!!)
  
  try: planesWithIssues = inconsistentPlanes['zalign']
  except KeyError:
    logging.info("All %d planes are correctly aligned.", len(Planes))
  else:
    logging.error("%s wire planes present misalignment: ",
      len(planesWithIssues), ", ".join(
        "%s vs. %s (%g mm)" % (planeA.ID(), planeB.ID(), distance)
        for distance, planeA, planeB in planesWithIssues
      ),
      )
    FailureSummary.append \
      ("%s wire planes present misalignment: " % len(misalignedPlanes))
  # if error
  
  try: planesWithIssues = inconsistentPlanes['view']
  except KeyError:
    logging.info("All %d planes have matching views.", len(Planes))
  else:
    logging.error("%s wire planes present direction inconsistencies:",
      len(planesWithIssues), ", ".join(
        "%s vs. %s" % (planeA.ID(), planeB.ID())
        for planeA, planeB in planesWithIssues
      ),
      )
    FailureSummary.append(
     "%s wire planes present direction inconsistencies." % len(planesWithIssues)
     )
  # if error
  
  try: planesWithIssues = inconsistentPlanes['direction']
  except KeyError:
    logging.info("All %d planes have consistent directions.", len(Planes))
  else:
    logging.error("%s wire planes present view inconsistencies:",
      len(planesWithIssues), ", ".join(
        "%s vs. %s" % (planeA.ID(), planeB.ID())
        for planeA, planeB in planesWithIssues
      ),
      )
    FailureSummary.append(
     "%s wire planes present view inconsistencies." % len(planesWithIssues)
     )
  # if error
  
  #
  # check wire alignment
  #
  misalignedWires = checkWireAlignment(Planes, tolerance=0.0001) # 1 um required
  if misalignedWires:
    logging.error("Misaligned wires found on %d extended planes:",
      len(misalignedWires),
      )
    for misalignedWiresOnPlane in misalignedWires:
      logging.error("%d on wires on plane around x=%g cm:",
       len(misalignedWiresOnPlane),
       misalignedWiresOnPlane[0].GetCenter().X()
       )
      for shift, wireLid, wireL, wireRid, wireR in misalignedWiresOnPlane:
        if shift is None:
          logging.error("  %s did not match any wire", wireLid)
        else:
          logging.error("  %s and %s are misaligned by %g um",
           wireLid, wireRid, shift * 1.0e4,
           )
      # for all misaligned wire pairs
    # for
    FailureSummary.append \
      ("Misaligned wires found on %d extended planes:" % len(misalignedWires))
  else:
    logging.info("No misaligned wires detected.")
  
  if FailureSummary:
    logging.info("%d failure categories detected:\n - %s",
     len(FailureSummary), "\n - ".join(FailureSummary),
     )
  # if failures
  return 1 if FailureSummary else 0
# performGeometryChecks()


# ##############################################################################
if __name__ == "__main__":
  import sys
  sys.exit(performGeometryChecks(sys.argv))
# main

# ##############################################################################
