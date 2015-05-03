from numpy.linalg import inv
from numpy import transpose
from numpy import dot
from numpy import concatenate as concat
import sys
import math
from pprint import pprint as pp

nodes = []
members = []
K_matrix = None
Delta_matrix = []
R_matrix = []

def get_node_by_name(name):
  for _node in nodes:
    if _node.name == name:
      return _node

class node:
  def __init__(self, _data):
    self.name = _data[1]
    self.x = float(_data[2])
    self.y = float(_data[3])
    self.constraint = (bool(int(_data[4])),bool(int(_data[5])),bool(int(_data[6])))

class member:
  def __init__(self, _data):
    self.node1 = get_node_by_name(_data[1])
    self.type1 = _data[2]
    self.node2 = get_node_by_name(_data[3])
    self.type2 = _data[4]
    self.EI = float(_data[5])
    self.EA = float(_data[6])
    self.build_T_matrix()
    self.build_Km_matrix()
    self.build_K_matrix()

  def build_T_matrix(self):
    node1 = self.node1
    node2 = self.node2
    L = math.sqrt(math.pow(node1.x-node2.x,2)+math.pow(node1.y-node2.y,2))
    S = (node2.y - node1.y)/L
    C = (node2.x - node1.x)/L
    self.Tm = [
      [ C,-S, 0],
      [ S, C, 0],
      [ 0, 0, 1]
    ]
    self.TmT = transpose(self.Tm)
    self.L = L

  def build_Km_matrix(self):
    self.Km = []
    EA = self.EA
    EI = self.EI
    L = self.L
    self.Km.append([  EA/L    ,      0      ,     0       ,   -EA/L    ,      0      ,     0      ])
    self.Km.append([    0     , 12*EI/L**3  , 6*EI/L**2   ,     0      , -12*EI/L**3 ,  6*EI/L**2 ])
    self.Km.append([    0     , 6*EI/L**2   , 4*EI/L      ,     0      , -6*EI/L**2  ,  2*EI/L    ])
    self.Km.append([  -EA/L   ,      0      ,     0       ,    EA/L    ,      0      ,     0      ])
    self.Km.append([    0     , -12*EI/L**3 , -6*EI/L**2  ,     0      , 12*EI/L**3  , -6*EI/L**2 ])
    self.Km.append([    0     , 6*EI/L**2   , 2*EI/L      ,     0      , -6*EI/L**2  ,  4*EI/L    ])

  def build_K_matrix(self):
    zeros = [
      [0,0,0],
      [0,0,0],
      [0,0,0]
    ]
    mat1 = concat((concat((self.Tm,zeros),axis=1),concat((zeros,self.Tm),axis=1)),axis=0)
    mat2 = concat((concat((self.TmT,zeros),axis=1),concat((zeros,self.TmT),axis=1)),axis=0)
    self.K = dot(mat1,self.Km,mat2)

  def extract_K(self, K):
    s = (int(self.node1.name)-1)*3
    t = (int(self.node2.name)-1)*3
    for i in range(3):
      for j in range(3):
        K[s+i][s+j] += self.K[ i ][ j ]
        K[t+i][t+j] += self.K[i+3][j+3]
        K[s+i][t+j] += self.K[ i ][j+3]
        K[t+i][s+j] += self.K[i+3][ j ]

if __name__ == "__main__":
  data = [line.rstrip('\n').split() for line in open(sys.argv[1])]
  for _data in data:
    if _data[0] == "node":
      nodes.append(node(_data))
    elif _data[0] == "member":
      members.append(member(_data))
  len_K = len(nodes)
  K_matrix = [[0]*len_K*3]*len_K*3
  for _member in members:
    _member.extract_K(K_matrix)
  pp(K_matrix)
  # pp(members[0].K)


"""
    node name x y constraintX constraintY constraintR
    member node1 type node2 type EI EA
"""