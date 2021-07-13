#!/usr/bin/env python

from ase.io import read,write
#from ase.geometry import cell_to_cellpar
import sys
import numpy as np
import random 
from copy import deepcopy as cp
import os

mol = read(sys.argv[1])
imol = read(sys.argv[2])
nmol = int(sys.argv[3])
cutoff = float(sys.argv[4]) 

ntry = 5000
inat = imol.get_global_number_of_atoms()
disp = np.zeros((inat,3))
refpos = imol.positions[0,:]
for i in range(1,inat):
  disp[i,:] =imol.positions[i,:] - refpos

print(disp)

cell = mol.cell[:,:]

try:
  os.remove('out.extxyz')
except:
  pass

success = True
for i in range(nmol):
  isPlaced = False
  for j in range(ntry):
    #random.seed(i*j*1234) 
    xf = random.random()
    yf = random.random()
    zf = random.random()
    factor = np.array([xf,yf,zf])
    frac = np.dot(factor,cell)
    im2 = cp(imol)
    for k in range(inat):
      im2.positions[k] = frac+disp[k] 
      #print(frac,disp[k],frac+disp[k])
#    print(im2.positions)


    #random.seed(i*j*1234*random.random()) 
    im2.rotate(random.random()*90,(random.random(),random.random(),random.random()))
    im1 = cp(mol)
    im1.extend(im2)
    im1.wrap()

    dist = im1.get_all_distances(mic=True)
    gnat = im1.get_global_number_of_atoms()
    indices = [x for x in range(gnat-inat) ]
    mindist = 1e5 
    for k in range(inat):
      mindist = min(mindist,np.min(im1.get_distances(gnat-k-1,indices,mic=True)))

    if mindist > cutoff:
      mol = cp(im1)
      write('out.extxyz',mol,append=True)
      isPlaced = True
      print("mol %d placed after %d attempts"%(i,j))
      break   
  if not isPlaced:
    print("STOP: unable to place molecules. last attempt %10.5f at %d, insert: %d"%(mindist,j,i)) 
    success = False
    break

if success:
  write("final.extxyz",mol) 
  
