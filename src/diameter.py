#!/home/hlma/anaconda3/bin/python3

import sys, os

#This python script is to create a class to compute diameter.
#2020-07-21


# Using floyd algorithm

class Diameter():

   def __init__(self,edge,V,adj):
      self.edge=edge
      self.V=V
      self.adj=adj
   def floyd(self):
    #  flag=2 
      inf = 99999999
      dis = {}  # dictionary of the shortest distance
      path = {}  # record the shortest path
      for i in self.V:
          for j in self.V:
              if i == j:
                 dis[(i,j)]=0
              else:
                 dis[(i,j)]=inf
      for i in self.V:
              for j in self.V:
                 path[(i,j)]=-1
        


      # read weight information
      for i in range(self.edge):
           u, v, w = self.adj[i][0],self.adj[i][1],self.adj[i][2]
          # if flag == 1:
          #    dis[(u,v)] = w
          # elif flag == 2:
           dis[(u,v)] = w
           dis[(v,u)] = w


      # floyd algorithm
      for k in self.V:
          for i in self.V:
              for j in self.V:
                  if dis[(i,j)] > dis[(i,k)] + dis[(k,j)]:
                      dis[(i,j)] = dis[(i,k)] + dis[(k,j)]
                      path[(i,j)] = k


      for i in self.V:
          for j in self.V:
              if dis[(i,j)] == inf:
                 dis[(i,j)] = 0

      diameter=max(dis.values())
      return diameter


#A=[[1,2,0.5],[2,3,0.6],[3,4,1]]
#e=3
#V=[1,2,3,4]

#D=Diameter(e,V,A)
#print(D.floyd())


