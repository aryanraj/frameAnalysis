from numpy.linalg import inv
import sys

class node:
  def __init__(self):
    pass

class member:
  def __init__(self):
    pass

if __name__ == "__main__":
  data = [line.rstrip('\n').split(',') for line in open(sys.argv[1])]
  print data