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
    self.FEM1 = 0
    self.FEM2 = 0
    self.FX1 = 0
    self.FX2 = 0
    self.FY1 = 0
    self.FY2 = 0
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
    self.L2G = concat((concat((self.Tm,zeros),axis=1),concat((zeros,self.Tm),axis=1)),axis=0)
    self.G2L = concat((concat((self.TmT,zeros),axis=1),concat((zeros,self.TmT),axis=1)),axis=0)
    self.K = np.dot(np.dot(self.L2G,self.Km),self.G2L)
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

  def manage_forceX(self, data):
    a = float(data[1])
    b = self.L - a
    P = float(data[2])
    FX1 = -P*b/L
    FX2 = -P*a/L
    self.FX1 += FX1
    self.FX2 += FX2

  def force_node_type_M(self, M1, M2):
    if self.type1 == "pinned" and self.type2 == "fixed":
      M2 -= M1/2
      M1 = 0
    if self.type1 == "fixed" and self.type2 == "pinned":
      M1 -= M2/2
      M2 = 0
    if self.type1 == "pinned" and self.type2 == "pinned":
      M1 = 0
      M2 = 0
    return (M1,M2)

  def manage_forceY(self, data):
    a = float(data[1])
    b = self.L - a
    P = float(data[2])
    L = self.L
    M1 = -P*b**2*a/L**2
    M2 = P*a**2*b/L**2
    M1,M2 = self.force_node_type_M(M1,M2)
    FY1 = -P*b/L + (M1+M2)/L
    FY2 = -P*a/L - (M1+M2)/L
    self.FEM1 += M1
    self.FEM2 += M2
    self.FY1 += FY1
    self.FY2 += FY2

  def manage_LDL(self, data):
    pos = data[1]
    W = float(data[2])
    L = self.L
    if pos == "left":
      M1 = -W*L**2/20
      M2 = W*L**2/30
      M1,M2 = self.force_node_type_M(M1,M2)
      FY1 = -W*L/3 + (M1+M2)/L
      FY2 = -W*L/6 - (M1+M2)/L
    elif pos == "right":
      M1 = -W*L**2/30
      M2 = W*L**2/20
      M1,M2 = self.force_node_type_M(M1,M2)
      FY1 = -W*L/6 + (M1+M2)/L
      FY2 = -W*L/3 - (M1+M2)/L
    self.FEM1 += M1
    self.FEM2 += M2
    self.FY1 += FY1
    self.FY2 += FY2

  def manage_UDL(self, data):
    W = float(data[1])
    L = self.L
    M1 = -W*L**2/12
    M2 = W*L**2/12
    M1,M2 = self.force_node_type_M(M1,M2)
    FY1 = -W*L/2
    FY2 = -W*L/2
    self.FEM1 += M1
    self.FEM2 += M2
    self.FY1 += FY1
    self.FY2 += FY2    

  def manage_moment(self, data):
    a = float(data[1])
    b = self.L - a
    M = float(data[2])
    L = self.L
    M1 = -M*b/L**2*(2*a-b)
    M2 = -M*a/L**2*(2*b-a)
    M1,M2 = self.force_node_type_M(M1,M2)
    FY1 = (M1+M2)/L
    FY2 = -(M1+M2)/L
    self.FEM1 += M1
    self.FEM2 += M2
    self.FY1 += FY1
    self.FY2 += FY2

  def extract_C(self, C):
    FL = [ self.FX1 , self.FY1 , self.FEM1 , self.FX2 , self.FY2 , self.FEM2 ]
    FG = np.dot(self.L2G,FL)
    s = (int(self.node1.name)-1)*3
    t = (int(self.node2.name)-1)*3
    
    C[s+0] += FG[0]
    C[s+1] += FG[1]
    C[s+2] += FG[2]
    C[t+0] += FG[3]
    C[t+1] += FG[4]
    C[t+2] += FG[5]


def shuffle_matrix(K,r,R,C):
  temp_K = []
  temp_R = []
  temp_C = []
  count = 0
  for i in range(len(R)):
    if isinstance(R[i], basestring):
      temp_K.insert(0,K[i])
      temp_R.insert(0,R[i])
      temp_C.insert(0,C[i])
      count += 1
    else:
      temp_K.append(K[i])
      temp_R.append(R[i])
      temp_C.append(C[i])

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
  F = np.array(temp_R)[count:len(R)]
  Y = np.array(temp_R)[0:count]
  G = np.array(temp_C)[0:count]
  H = np.array(temp_C)[count:len(R)]

  return (A,B,C,D,E,F,G,H,X,Y)


if __name__ == "__main__":
  data = [line.rstrip('\n').split() for line in open(sys.argv[1])]
  for _data in data:
    if len(_data) != 0:
      if _data[0] == "node":
        nodes.append(node(_data))
      elif _data[0] == "member":
        members.append(member(_data))
      elif _data[0] == "forceX":
        members[len(members)-1].manage_forceX(_data)
      elif _data[0] == "forceY":
        members[len(members)-1].manage_forceY(_data)
      elif _data[0] == "UDL":
        members[len(members)-1].manage_UDL(_data)
      elif _data[0] == "LDL":
        members[len(members)-1].manage_LDL(_data)
      elif _data[0] == "moment":
        members[len(members)-1].manage_moment(_data)
  
  len_K = len(nodes)
  K_matrix = [[ 0 for j in range(len_K*3) ] for i in range(len_K*3)]
  R_matrix = [0 for i in range(len_K*3)]
  r_matrix = [0 for i in range(len_K*3)]
  C_matrix = [0 for i in range(len_K*3)]

  for _member in members:
    _member.extract_K(K_matrix)
    _member.extract_C(C_matrix)
  for _node in nodes:
    _node.extract_R(R_matrix)
    _node.extract_r(r_matrix)

  A,B,C,D,E,F,G,H,X,Y = shuffle_matrix(K_matrix,r_matrix,R_matrix,C_matrix)

  A = np.array(A,dtype=float)
  B = np.array(B,dtype=float)
  C = np.array(C,dtype=float)
  D = np.array(D,dtype=float)
  E = np.array(E,dtype=float)
  F = np.array(F,dtype=float)
  G = np.array(G,dtype=float)
  H = np.array(H,dtype=float)

  _X = np.dot(inv(D),np.subtract(np.subtract(F,np.dot(C,E)),H))
  _Y = np.add(np.add(np.dot(B,_X),np.dot(A,E)),G)

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


"""
    node name x y fx fy fr
    member node1 type node2 type EI EA
    forceX x` magnitude
    forceY x` magnitude
    UDL magnitude
    LDL pos magnitude
    moment x` magnitude 
"""

"""
    [ A B ][E] = [Y]-[G]
    [ C D ][X] = [F]-[H]
    X = inv(D)*(F-C*E-H)
    Y = B*X+A*E+G
"""