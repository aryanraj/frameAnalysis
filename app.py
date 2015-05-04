from numpy.linalg import inv
from numpy import transpose
from numpy import concatenate as concat
import numpy as np
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
    self.fx = _data[4]
    self.fy = _data[5]
    self.fr = _data[6]

  def extract_R(self, R):
    start = (int(self.name)-1)*3
    if self.fx != 'u':
      R[start] = float(self.fx)
    else:
      R[start] = 'R'+str(start)
    if self.fy != 'u':
      R[start+1] = float(self.fy)
    else:
      R[start+1] = 'R'+str(start+1)
    if self.fr != 'u':
      R[start+2] = float(self.fr)
    else:
      R[start+2] = 'R'+str(start+2)

  def extract_r(self, r):
    start = (int(self.name)-1)*3
    if self.fx != 'u':
      r[start] = 'r'+str(start)
    else:
      r[start] = 0
    if self.fy != 'u':
      r[start+1] = 'r'+str(start+1)
    else:
      r[start+1] = 0
    if self.fr != 'u':
      r[start+2] = 'r'+str(start+2)
    else:
      r[start+2] = 0

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
    # self.build_Final_K()

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
    if self.type1 == "fixed" and self.type2 == "fixed":
      self.Km.append([  EA/L    ,      0      ,     0       ,   -EA/L    ,      0      ,     0      ])
      self.Km.append([    0     , 12*EI/L**3  , 6*EI/L**2   ,     0      , -12*EI/L**3 ,  6*EI/L**2 ])
      self.Km.append([    0     , 6*EI/L**2   , 4*EI/L      ,     0      , -6*EI/L**2  ,  2*EI/L    ])
      self.Km.append([  -EA/L   ,      0      ,     0       ,    EA/L    ,      0      ,     0      ])
      self.Km.append([    0     , -12*EI/L**3 , -6*EI/L**2  ,     0      , 12*EI/L**3  , -6*EI/L**2 ])
      self.Km.append([    0     , 6*EI/L**2   , 2*EI/L      ,     0      , -6*EI/L**2  ,  4*EI/L    ])
    elif self.type1 == "fixed" and self.type2 == "pinned":
      self.Km.append([  EA/L    ,      0      ,     0       ,   -EA/L    ,      0      ,     0      ])
      self.Km.append([    0     ,  3*EI/L**3  , 3*EI/L**2   ,     0      ,  -3*EI/L**3 ,     0      ])
      self.Km.append([    0     ,  3*EI/L**2  , 3*EI/L      ,     0      ,  -3*EI/L**2 ,     0      ])
      self.Km.append([  -EA/L   ,      0      ,     0       ,    EA/L    ,      0      ,     0      ])
      self.Km.append([    0     , -3*EI/L**3  , -3*EI/L**2  ,     0      ,  3*EI/L**3  ,     0      ])
      self.Km.append([    0     ,      0      ,     0       ,     0      ,      0      ,     0      ])
    elif self.type1 == "pinned" and self.type2 == "fixed":
      self.Km.append([  EA/L    ,      0      ,     0       ,   -EA/L    ,      0      ,     0      ])
      self.Km.append([    0     ,  3*EI/L**3  ,     0       ,     0      ,  -3*EI/L**3 ,  3*EI/L**2 ])
      self.Km.append([    0     ,      0      ,     0       ,     0      ,      0      ,     0      ])
      self.Km.append([  -EA/L   ,      0      ,     0       ,    EA/L    ,      0      ,     0      ])
      self.Km.append([    0     ,  -3*EI/L**3 ,     0       ,     0      ,  3*EI/L**3  , -3*EI/L**2 ])
      self.Km.append([    0     ,  3*EI/L**2  ,     0       ,     0      , -3*EI/L**2  ,  3*EI/L    ])
    elif self.type1 == "pinned" and self.type2 == "pinned":
      self.Km.append([  EA/L    ,      0      ,     0       ,   -EA/L    ,      0      ,     0      ])
      self.Km.append([    0     ,      0      ,     0       ,     0      ,      0      ,     0      ])
      self.Km.append([    0     ,      0      ,     0       ,     0      ,      0      ,     0      ])
      self.Km.append([  -EA/L   ,      0      ,     0       ,    EA/L    ,      0      ,     0      ])
      self.Km.append([    0     ,      0      ,     0       ,     0      ,      0      ,     0      ])
      self.Km.append([    0     ,      0      ,     0       ,     0      ,      0      ,     0      ])


  def build_K_matrix(self):
    zeros = [
      [0,0,0],
      [0,0,0],
      [0,0,0]
    ]
    mat1 = concat((concat((self.Tm,zeros),axis=1),concat((zeros,self.Tm),axis=1)),axis=0)
    mat2 = concat((concat((self.TmT,zeros),axis=1),concat((zeros,self.TmT),axis=1)),axis=0)
    self.K = np.dot(np.dot(mat1,self.Km),mat2)
    # print('\n'.join([''.join(['{:20}'.format(item) for item in row]) for row in self.K]))

  def build_Final_K(self):
    node1 = self.node1
    node2 = self.node2
    L = math.sqrt(math.pow(node1.x-node2.x,2)+math.pow(node1.y-node2.y,2))
    C = (node2.x - node1.x)/L
    S = (node2.y - node1.y)/L
    C2 = C*C
    S2 = S*S
    CS = C*S
    AE_L = self.EA/L
    EI_L3 = self.EI/L**3
    EI_L2 = self.EI/L**2
    EI_L = self.EI/L
    self.K = []
    self.K.append([ (AE_L*C2+12*EI_L3*S2)   , (AE_L-12*EI_L3)*CS      , -6*EI_L2*S    ,  -(AE_L*C2+12*EI_L3*S2) , -(AE_L-12*EI_L3)*CS     , -6*EI_L2*S  ])
    self.K.append([ (AE_L-12*EI_L3)*CS      , (AE_L*S2+12*EI_L3*C2)   , 6*EI_L2*C     , -(AE_L-12*EI_L3)*CS     , -(AE_L*S2+12*EI_L3*C2)  , 6*EI_L2*C   ])
    self.K.append([ -6*EI_L2*S              , 6*EI_L2*C               , 4*EI_L        , 6*EI_L2*S               , -6*EI_L2*C              , 2*EI_L      ])
    self.K.append([ -(AE_L*C2+12*EI_L3*S2)  , -(AE_L-12*EI_L3)*CS     , 6*EI_L2*S     , (AE_L*C2+12*EI_L3*S2)   , (AE_L-12*EI_L3)*CS      , 6*EI_L2*S   ])
    self.K.append([ -(AE_L-12*EI_L3)*CS     , -(AE_L*S2+12*EI_L3*C2)  , -6*EI_L2*C    , (AE_L-12*EI_L3)*CS      , (AE_L*S2+12*EI_L3*C2)   , -6*EI_L2*C  ])
    self.K.append([ -6*EI_L2*S              , 6*EI_L2*C               , 2*EI_L        , 6*EI_L2*S               , -6*EI_L2*C              , 4*EI_L      ])
    # pp(self.K)


  def extract_K(self, K):
    s = (int(self.node1.name)-1)*3
    t = (int(self.node2.name)-1)*3  
    for i in range(3):
      for j in range(3):
        K[s+i][s+j] += self.K[ i ][ j ]
        K[t+i][t+j] += self.K[i+3][j+3]
        K[s+i][t+j] += self.K[ i ][j+3]
        K[t+i][s+j] += self.K[i+3][ j ]
    # print('\n'.join([''.join(['{:20}'.format(item) for item in row]) for row in K]))
    # pp(K)

  def extract_R(self, R):
    pass

def shuffle_matrix(K,r,R):
  temp_K = []
  temp_R = []
  count = 0
  for i in range(len(R)):
    if isinstance(R[i], basestring):
      temp_K.insert(0,K[i])
      temp_R.insert(0,R[i])
      count += 1
    else:
      temp_K.append(K[i])
      temp_R.append(R[i])

  K = transpose(temp_K)
  temp_K = []
  temp_r = []
  for i in range(len(r)):
    if isinstance(r[i], basestring):
      temp_K.append(K[i])
      temp_r.append(r[i])
    else:
      temp_K.insert(0,K[i])
      temp_r.insert(0,r[i])

  temp_K = transpose(temp_K)

  A = np.array(temp_K)[0:count,0:count]
  B = np.array(temp_K)[0:count,count:len(R)]
  C = np.array(temp_K)[count:len(R),0:count]
  D = np.array(temp_K)[count:len(R),count:len(R)]
  E = np.array(temp_r)[0:count]
  X = np.array(temp_r)[count:len(R)]
  Y = np.array(temp_R)[0:count]
  F = np.array(temp_R)[count:len(R)]

  return (A,B,C,D,E,F,X,Y)


if __name__ == "__main__":
  data = [line.rstrip('\n').split() for line in open(sys.argv[1])]
  for _data in data:
    if len(_data) != 0:
      if _data[0] == "node":
        nodes.append(node(_data))
      elif _data[0] == "member":
        members.append(member(_data))
  
  len_K = len(nodes)
  K_matrix = [[ 0 for j in range(len_K*3) ] for i in range(len_K*3)]
  R_matrix = [0 for i in range(len_K*3)]
  r_matrix = [0 for i in range(len_K*3)]

  for _member in members:
    _member.extract_K(K_matrix)
    _member.extract_R(R_matrix)
  for _node in nodes:
    _node.extract_R(R_matrix)
    _node.extract_r(r_matrix)

  A,B,C,D,E,F,X,Y = shuffle_matrix(K_matrix,r_matrix,R_matrix)

  A = np.array(A,dtype=float)
  B = np.array(B,dtype=float)
  C = np.array(C,dtype=float)
  D = np.array(D,dtype=float)
  E = np.array(E,dtype=float)
  F = np.array(F,dtype=float)

  _X = np.dot(inv(D),np.subtract(F,np.dot(C,E)))
  _Y = np.add(np.dot(B,_X),np.dot(A,E))

  print "A"
  print('\n'.join([''.join(['{:15}'.format(item) for item in row]) for row in A]))
  print "B"
  print('\n'.join([''.join(['{:15}'.format(item) for item in row]) for row in B]))
  print "C"
  print('\n'.join([''.join(['{:15}'.format(item) for item in row]) for row in C]))
  print "D"
  print('\n'.join([''.join(['{:15}'.format(item) for item in row]) for row in D]))

  print "K"
  print('\n'.join([''.join(['{:15}'.format(item) for item in row]) for row in K_matrix]))

  print "r"
  print r_matrix
  
  print "R"
  print R_matrix

  print X,_X
  print Y,_Y

  # pp(X)
  # pp(_X)
  # pp(Y)
  # pp(_Y)


"""
    node name x y fx fy fr
    member node1 type node2 type EI EA
"""

"""
    [ A B ][E] = [Y]
    [ C D ][X] = [F]
    X = inv(D)*(F-C*E)
    Y = B*X+A*E
"""