#!/usr/bin/env python
#
# Applies shifts on a specified pairing pattern and produces an output
# good for being pasted in a configuration file in JSON/FHiCL format.
# 
# The applied pattern is coded in MainPattern.
#

from __future__ import print_function

#
#   y
# /\
# ||
# ||   |                              |                             |
# ||   |        4    7             14 | 17             24   27      |
# ||   |   1              9   11      |      19   21             29 |
# ||   |        3    6             13 | 16             23   26      |
# ||   |   0              8   10      |      18   20             28 |
# ||   |        2    5             12 | 15             22   25      |
# ||   |                              |                             |
# ||
# ##====================================================================> z
#
# pattern in three levels:
#   * inner:  pairs (or single channels)
#   * middle: line-group
#   * outer:  pattern-group
# 
# line-groups are printed each on its own line; there is no other meaning
# in that middle layer beyond output eye candy.
# 

MainPattern = [
  [
    [  0,  2, ],
    [  1,  4, ],
    [  3,  6, ],
    [  5,  8, ],
    [  7,  9, ],
    [ 10, 12, ],
    [ 11, 14, ],
    [ 13      ],
  ],
  [
    [ 15, 18, ],
    [ 16      ],
    [ 17, 19, ],
    [ 20, 22, ],
    [ 21, 24, ],
    [ 23, 26, ],
    [ 25, 28, ],
    [ 27, 29, ],
  ],
] # MainPattern

NChannels = 360


# ------------------------------------------------------------------------------
class Pattern:
  
  def __init__(self, patter, offset):
    self.pairs = self.shiftPattern(patter, offset)
  
  @staticmethod
  def shiftPattern(pattern, offset):
    l = list()
    for orig in pattern:
      l.append(
        Pattern.shiftPattern(orig, offset) if isinstance(orig, list)
          else orig + offset
      )
    # for
    return l
  # shiftPattern()
  
  def __str__(self):
    return "\n".join(
      " ".join(
        "{:12}".format(
          "[ {} ],".format(
            ",".join(
              "{:3d}".format(c) for c in p # elements in pair
            )
          )
        )
        for p in g # pairs in group
      )
      for g in self.pairs # groups
    )
  # __str__()
  
# class Pattern


def countElements(pattern):
  c = 0
  for orig in pattern:
    c += countElements(orig) if isinstance(orig, list) else 1
  return c
# countElements()


# ------------------------------------------------------------------------------
if __name__ == "__main__":
  
  ChannelsPerPattern = countElements(MainPattern)
  print("Pattern is made of {} elements.".format(ChannelsPerPattern))
  
  for offset in range(0, NChannels, ChannelsPerPattern):
    pattern = Pattern(MainPattern, offset)
    print("# offset {}\n{}".format(offset, pattern))
  # for
  
# main

