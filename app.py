from numpy.linalg import inv
import sys
import math

class node:
  def __init__(self, x, y):
    self.x = x
    self.y = y

class member:
  def __init__(self, _node1, type1, _node2, type2, EI, EA):
    self.L = math.sqrt(math.pow(_node1.x-_node2.x,2),math.pow(_node1.y-_node2.y,2))
    self.S = (_node2.y - _node1.y)/self.L
    self.C = (_node2.x - _node1.x)/self.L
    self.node1 = _node1
    self.type1 = type1
    self.node2 = _node2
    self.type2 = type2
    self.EI = EI
    self.EA = EA


if __name__ == "__main__":
  data = [line.rstrip('\n').split(',') for line in open(sys.argv[1])]
  print data


""" Structure of data file
    node name x y
    member node1 type node2 type EI EA
"""