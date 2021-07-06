#!/usr/bin/env python3

##################################################
#### HELP ######
#
# run 
# anyconvert2.py --help
# to get started 
#
################
##################################################


flag = False
try:
  from scm.plams import *
except:
  flag = True
  print("\n Unable to find scm.plams. Installing it for you...")
  import pip
  pip.main(['install','--user','plams'])
try:  
  import ase.io as io            
except:
  print("\n Unable to find ase. Installing it for you...")
  if not flag:    import pip
  flag = True
  pip.main(['install','--user','ase'])
if flag:
  print("\n\n!!!! Had to install packages (locally). Please run the program again.\n")  
  exit()
import numpy as np
import argparse
import sys
import os



#
# print help and exit
#
def help_n_exit(message):
    print(message)
    sys.exit(1)

def print_support_exit(i=0):
  print("\n See https://wiki.fysik.dtu.dk/ase/ase/io/io.html for ASE supported formats")
  print(" As you scroll down the page, there will be a table with <format> <description> <capability>")
  print(" Your input/output files should have the extension name same as the <format> listed in this table.")
  print(" e.g. POSCAR : POSCAR.vasp")
  print("      cif    : graphite.cif")
  print("      xmolout: xmolout.xmol")
  print(" The name before the .format can be arbitrary")
  print("\n In addition, this program supports xmol and bgf\n")
  print(" version July 4th, 2021\n Written by Karthik Ganeshan (kug46.psu@gmail.com)\n")
  sys.exit(i)


def readfile(infile,index):
  # Returns ase mol list from input file
  filename, extension    = os.path.splitext(infile)
  extension = extension.split('.')[-1]
  if index == -2: index = ":"
  if extension == "xmol":
    return read_xmol(infile,index)
  if extension == "bgf":
    return read_bgf(infile,index)
  mollist = io.read(infile,format=extension,index=index)
  if index != ":": mollist = [mollist]
  for mol in mollist:
    if np.max(mol.cell) == 0.0: 
      mol.pbc = True
      mol.cell = [[100,0,0],[0,100,0],[0,0,100]]
  return mollist

def read_xmol(infile,index):
  extension = "xmol"
  mollist = io.read(infile,format="xyz",index=index) 
  if index != ":": mollist = [mollist]
  with open(infile,'r') as f:
    for  mol in mollist:
      nat = int(f.readline())
      rcell = [ float(x) for x in f.readline().split()[3:9]]
      cell = reax2ase(rcell)
      mol.cell = cell
      mol.pbc = True
      mol.wrap()
      for i in range(nat): next(f) 
  return mollist
      
def read_bgf(infile,index):
  extension = "bgf"
  if index != 0:
    print("\n I'm converting only the first frame of .bgf input. Use the .xmol format for frame control.\n")
  cryst, atoms = [], []
  with open(infile,'r') as f:
    while True:
      line = f.readline()
      if "END" in line:
        break
      else:
        if "BIOGRF" in line:
          cryst = [100,100,100,90,90,90]
        elif "CRYSTX" in line:
          cryst = [float(x) for x in line.split()[1:]]
        elif "HETATM" in line:
          l = line.split()
          atoms.append([l[2],float(l[3]),float(l[4]),float(l[5])])
  nat = len(atoms)
  from ase import Atoms as aseatoms
  mol = aseatoms("C%d"%nat)
  mol.cell = reax2ase(cryst)
  symbols = [ at[0] for at in atoms ]
  mol.set_chemical_symbols(symbols)
  for i in range(nat):
    mol.positions[i] = [atoms[i][1],atoms[i][2],atoms[i][3]]
  mol.pbc = True
  mol.wrap
  return [mol]


def writefile(mollist,outfile):
  filename,extension    = os.path.splitext(outfile)
  extension = extension.split('.')[-1]
  # Write ase mol to file
  try:
    os.remove(outfile)
  except:
    pass
  if extension == "bgf" or extension == "xmol": 
    ase2reax(mollist,outfile)
  else:
    for mol in mollist:
      mol.wrap()   # Wrap atoms 
      io.write(outfile,mol,format=extension)

def convert_lattice(lattice):
  # Takes in a 3x3 cell matrix as per standard convention
  # returns abcABC as per ReaxFF convention
  # COPIED FROM SCM's xyz2bgf.py 
  a, b, c = map(np.linalg.norm, lattice)
  al = np.dot(lattice[1], lattice[2])/(b*c)
  be = np.dot(lattice[0], lattice[2])/(a*c)
  ga = np.dot(lattice[0], lattice[1])/(a*b)
  al,be,ga = map(lambda x: Units.convert(np.arccos(x),'rad','deg'), [al,be,ga])
  return a, b, c, al, be, ga
  

def reax2ase(lattice):
  # Takes in abcABC 6-element list 
  # returns VEC{1,2,3}
  # convention followed by Reax will be retained, but is OK for ase 
  # COPIED FROM SCM's bgf2ams.awk
  axes = [lattice[0],lattice[1],lattice[2]]
  angles = [lattice[3],lattice[4],lattice[5]]

  deg2rad = 3.1415926535897932384626433832795 / 180;
  halfa           = angles[0]*deg2rad;
  hbeta           = angles[1]*deg2rad;
  hgamma          = angles[2]*deg2rad;
  sinalf          = np.sin(halfa);
  cosalf          = np.cos(halfa);
  sinbet          = np.sin(hbeta);
  cosbet          = np.cos(hbeta);
  cosphi          = (np.cos(hgamma)-cosalf*cosbet)/(sinalf*sinbet);
  if (cosphi >  1): cosphi =  1;
  elif (cosphi < -1): cosphi = -1;
  sinphi    = np.sqrt(1.0 - cosphi*cosphi);
  tm11 = axes[0]*sinbet*sinphi;
  tm21 = axes[0]*sinbet*cosphi;
  tm31 = axes[0]*cosbet;
  tm22 = axes[1]*sinalf;
  tm32 = axes[1]*cosalf;
  tm33 = axes[2];

  cell = [[tm11,tm21,tm31],[0,tm22,tm32],[0,0,tm33]]

  return cell

def write_geofile(molecule, filename, description):
    header = ['XTLGRF 200\n']
    header.append('DESCR ' + str(description) + '\n')

    if molecule.align_lattice(convention='reax'):
        log("The lattice supplied for job did not follow the convention required by ReaxFF. I rotated the whole system for you. You're welcome", 3)

    f = lambda x: tuple([100.0 * int(i==x) for i in range(3)])
    while len(molecule.lattice) < 3:
        molecule.lattice.append(f(len(molecule.lattice)))

    header.append('CRYSTX  {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f}\n'.format(*convert_lattice(molecule.lattice)))

    atoms = []
    for i,at in enumerate(molecule):
        newline = 'HETATM {:>5d} {:<2}               {: >10.5f}{: >10.5f}{: >10.5f} {:<2}     1 1  0.0\n'.format(i+1, at.symbol, *at.coords, at.symbol)
        atoms.append(newline)
    atoms.append('END\n')

    with open(filename,'a') as f:
        f.writelines(header)
        f.writelines(atoms)
        f.write("\n")

def write_xmol(molecule, filename, description,nat):
    outfile = open(filename,'a')

    if molecule.align_lattice(convention='reax'):
        log("The lattice supplied for job did not follow the convention required by ReaxFF. I rotated the whole system for you. You're welcome", 3)

    f = lambda x: tuple([100.0 * int(i==x) for i in range(3)])
    while len(molecule.lattice) < 3:
        molecule.lattice.append(f(len(molecule.lattice)))
    header = []
    header.append('%d\n'%nat)
    header.append('{}                      0 {:10.2f} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f}\n'.format(description,0.0,*convert_lattice(molecule.lattice)))

    atoms = []
    for i,at in enumerate(molecule):
        newline = '{:<2}{: >10.5f}{: >10.5f}{: >10.5f}\n'.format(at.symbol, *at.coords)
        atoms.append(newline)

    with open(filename,'a') as f:
        f.writelines(header)
        f.writelines(atoms)

def ase2reax(mollist,outfile,description="converted"):
  filename,extension    = os.path.splitext(outfile)
  extension = extension.split('.')[1]
  try:
    os.remove(outfile)
  except:
    pass
  for frame,mol in enumerate(mollist):
    amsmol = Molecule()
    nat = mol.get_number_of_atoms()
    for i in range(nat):
      amsmol.add_atom(Atom(symbol=mol.get_chemical_symbols()[i],coords=(mol.positions[i][0],mol.positions[i][1],mol.positions[i][2])))
    amsmol.lattice = mol.cell[:,:]
    description = description + "%s"%frame
    if extension == "bgf":
      write_geofile(amsmol,outfile,description)
    elif extension == "xmol":
      write_xmol(amsmol,outfile,description,nat)
    else:
      print("Error. Write to reax format needs to be .bgf or .xmol only")
      sys.exit(5)
  
  



if __name__ == "__main__":
  
  #------ Support ---------------
  if len(sys.argv) < 2:
      help_n_exit("\n USAGE: " + str(sys.argv[0]) + " -i filename.format -o filename.format\n        or, "+str(sys.argv[0])+" infile.format outfile.format [optional: frame #]\n For more help: %s --help\n"%sys.argv[0])     
  #------------------------------  

  if 3 <= len(sys.argv) <= 4:
    infile = sys.argv[1]
    outfile = sys.argv[2]
    try:
      index = int(sys.argv[3])
    except:
      index = -2
  else:
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input',help='input file with ase format(or xmol/bgf): filename.format, e.g. graphene.xmol')
    parser.add_argument('-o','--output',help='output file with ase format(or xmol/bgf): filename.format, e.g. graphene.bgf')
    parser.add_argument('-f','--frame',default=-2,help="""frame number to use as input starting from 0. Use ....     
                                               -1: last frame,
                                               -2: all frames (default)""")
    parser.add_argument('-s','--supported-formats',action="store_true",help='List help for supported formats')
    args = parser.parse_args() 

    if args.supported_formats:
      print(print_support_exit())

    infile = args.input    
    index = int(args.frame)
    outfile = args.output

  mollist = readfile(infile,index)
  writefile(mollist,outfile)


