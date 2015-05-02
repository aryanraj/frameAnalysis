from numpy.linalg import inv
from numpy import transpose
import sys
import math

class node:
  def __init__(self, x, y):
    self.x = x
    self.y = y

class member:
  def __init__(self, _node1, type1, _node2, type2, EI, EA):
    self.node1 = _node1
    self.type1 = type1
    self.node2 = _node2
    self.type2 = type2
    self.EI = EI
    self.EA = EA
    self.build_T_matrix()

  def build_T_matrix(self):
    L = math.sqrt(math.pow(_node1.x-_node2.x,2),math.pow(_node1.y-_node2.y,2))
    S = (_node2.y - _node1.y)/L
    C = (_node2.x - _node1.x)/L
    self.Tm = [
      [ C,-S, 0],
      [ S, C, 0],
      [ 0, 0, 1]
    ]
    self.TmT = transpose(self.Tm)


if __name__ == "__main__":
  data = [line.rstrip('\n').split(',') for line in open(sys.argv[1])]
  print data


""" 
    Structure of data file
    node name x y
    member node1 type node2 type EI EA
"""