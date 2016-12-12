#!/bin/bash
# reads in a polymer configuration, writes it to a cylinder with SRD particles
# Cmd line args: poly.lammpstrj Npolymer cylr cylh

import numpy as np
import math
import sys


print "lammpstrj l_box"
infile = sys.argv[1]
npoly = 27
l_box = float(sys.argv[2])

# calculate number of SRD particles:
v_box = l_box*l_box*l_box
nsrd = int(10*v_box)		# 10 is the SRD particle density, can change later
print "There are ", nsrd, " SRD particles."

# Read Protein file
proteinfile=infile
ctr = 1
proteinpos=[]
with open(proteinfile, "r") as f:
    for line in f:
	if (ctr == 6):
	    line = line.split()
	    xl = float(line[0])
	    xh = float(line[1])
	if (ctr == 7):
	    line = line.split()
	    yl = float(line[0])
	    yh = float(line[1])
	if (ctr == 8):
	    line = line.split()
	    zl = float(line[0])
	    zh = float(line[1])
	if (ctr > 9):
	    line = line.split()
	    x = (xh-xl)*float(line[2]) + xl
	    y = (yh-yl)*float(line[3]) + yl
	    z = (zh-zl)*float(line[4]) + zl
	    proteinpos.append([x,y,z])
	ctr = ctr + 1

# Unwrap protein:
for i in range(len(proteinpos)-1):
    first = proteinpos[i]
    second = proteinpos[i+1]
    dx = second[0] - first[0]
    dy = second[1] - first[1]
    dz = second[2] - first[2]
    second[0] = second[0] - (xh-xl)*round(dx/(xh-xl))
    second[1] = second[1] - (yh-yl)*round(dy/(yh-yl))
    second[2] = second[2] - (zh-zl)*round(dz/(zh-zl))
    proteinpos[i+1] = second


# Calculate center of mass:
xcm = 0
ycm = 0
zcm = 0
ctr = 0
for pos in proteinpos:
    xcm = xcm + pos[0]
    ycm = ycm + pos[1]
    zcm = zcm + pos[2]
    ctr = ctr + 1
xcm = xcm/float(ctr)
ycm = ycm/float(ctr)
zcm = zcm/float(ctr)
print "center of mass ctr: ", ctr

# Remap protein to center:
temp = []
for pos in proteinpos:
    x = pos[0] - xcm
    y = pos[1] - ycm
    z = pos[2] - zcm
    print x,y,z
    temp.append([x,y,z])

proteinpos = temp



halfbox = 0.5*l_box
xlo = -halfbox 
xhi = halfbox 
ylo = -halfbox 
yhi = halfbox 
zlo = -halfbox 
zhi = halfbox 

atoms = npoly + nsrd
bonds = npoly - 1
angles = npoly -2
dihedrals = npoly-3
impropers = 0
atomtypes = 4
bondtypes = 1
angletypes = 1
dihedraltypes = 3
impropertypes = 0

g=open("shea_srd.dat", "w")
g.write("LAMMPS Description\n\n")
g.write("\t"+str(atoms)+"\tatoms\n")
g.write("\t"+str(bonds)+"\tbonds\n")
g.write("\t"+str(angles)+"\tangles\n")
g.write("\t"+str(dihedrals)+"\tdihedrals\n")
g.write("\t"+str(impropers)+"\timpropers\n\n")
g.write("\t"+str(atomtypes)+"\tatom types\n")
g.write("\t"+str(bondtypes)+"\tbond types\n")
g.write("\t"+str(angletypes)+"\tangle types\n")
g.write("\t"+str(dihedraltypes)+"\tdihedral types\n")
g.write("\t"+str(impropertypes)+"\timproper types\n\n\n")

#Write Box Lengths:
g.write(str(xlo)+" "+str(xhi)+" xlo xhi\n")
g.write(str(ylo)+" "+str(yhi)+" ylo yhi\n")
g.write(str(zlo)+" "+str(zhi)+" zlo zhi\n\n\n")


#Write Atomic Masses:
g.write("Masses\n\n")
g.write("\t1 "+str(10.0)+"\n")
g.write("\t2 "+str(10.0)+"\n")
g.write("\t3 "+str(10.0)+"\n")
g.write("\t4 "+str(1.0)+"\n")

#Write Atoms:
ctr = 1
zpos = zlo
g.write("\nAtoms\n\n")
typelist = [4,3,4,3,4,3,5,5,3,4,3,4,3,4,5,5,5,3,3,4,4,3,3,4,4,3,4]
for pos in proteinpos:
#    g.write("\t"+str(ctr)+" "+str(ctr)+str(" 1 0 "))
    g.write("\t"+str(ctr)+" "+str(typelist[ctr-1]-2)+" "+str(pos[0])+" "+str(pos[1])+" "+str(pos[2]))
    g.write(" "+str(ctr)+" 0 1 19.098593171027442\n")
    ctr = ctr + 1

for i in xrange(nsrd):
#    g.write("\t"+str(ctr)+" "+str(ctr)+str(" 2 0 "))
    rand_x = np.random.random_sample()
    rand_y = np.random.random_sample()
    rand_z = np.random.random_sample()
    x = rand_x*l_box-0.5*l_box
    y = rand_y*l_box-0.5*l_box
    z = rand_z*l_box-0.5*l_box
#    g.write(str(x)+" "+str(y)+" "+str(z)+" 0 1\n")
    g.write("\t"+str(ctr)+" 4 "+str(x)+" "+str(y)+" "+str(z))
    g.write(" "+str(ctr)+" 0 0 1\n")

    ctr = ctr + 1


# Write Bonds:
g.write("\nBonds\n\n")
bondnum = 1
# Brushes
for i in range(npoly-1):
    g.write("\t"+str(bondnum)+" 1 "+str(1+i)+" "+str(2+i)+"\n")
    bondnum = bondnum + 1 

# Write Angles:
g.write("\nAngles\n\n")
anglenum = 1
firstpolymono=1
for i in range(npoly-2):
    g.write("\t"+str(anglenum)+" 1 "+str(firstpolymono+i)+" "+str(firstpolymono+i+1)+" "+str(firstpolymono+i+2)+"\n")
    anglenum = anglenum + 1
    

# Write Dihedrals:
g.write("\nDihedrals\n\n")
dihedralnum = 1
# alpha = 1, beta = 2, turn = 3
dihedraltypelist=[2,2,2,2,3,3,3,2,2,2,2,2,3,3,1,1,1,1,1,1,1,1,1,1]
for i in range(npoly-3):
    g.write("\t"+str(dihedralnum)+" "+str(dihedraltypelist[dihedralnum-1])+" "+str(firstpolymono+i)+" "+str(firstpolymono+i+1)+" "+str(firstpolymono+i+2)+" "+str(firstpolymono+i+3)+"\n")
    dihedralnum = dihedralnum + 1
