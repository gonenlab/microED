#! /usr/bin/env python

# calculate a*, b*, and c* for xtals based on spots picked from diffraction patterns
# w/o any major plane images.  Only set up to work if a = b != c and alpha=beta-=gamma = 90.
 
## !! ## !!!! ## Be very conscious of y & z coordinates!  program currently uses  ** Brent style ** coords as opposed to image style coords:
## top left corner = 0,4096 bottom-left corner - 0,0

# Matt Iadanza 2013-06-24

# user entered variable
##################
numtoreturn = 20 # how many vector pairs to return
lowthresh = 50   # a threshold for low magnitudes to eliminate conparing the same spot on different images

findmults = "N"    #  (Y/N) find multiples of unit cells? 2*a, 3*b, 2*c, 4*b 

eab = 54.8       # expected unitcell a and b 
ec = 114.0      # expected unitcell c
abthresh = 3.0  # +- when searching for a and b unit cell vectors
cthresh = 4.5   # +- when searching for c unit cell vectors

##################


### imports

import json
import os
import sys
import math
import numpy


# get some variables:
data = json.load(open('cellfind_params.json'))
numberofimages = len(data["allimages"])
output = open("output_cellfind.txt", "w")

#make the spolist dictionarys
spotlist= {}    
spotcount = 0
xydic = {}  # original x,y coords to image#,spot#
xyzdic = {} # calculated x,y,z coords to image#,spot#
xyzrevdic = {}  # input imag#,spot# returns calculated x,z,y as vector

output.write("spot name - x,y - x,y,z\n")

for eachimage in range(1, numberofimages+1):
    spotrange = range(1, len(data["allimages"]["image"+str(eachimage)]["spots"])+1)
    theta = data["allimages"]["image"+str(eachimage)]["theta"]
    for eachspot in spotrange:
        spotcount = spotcount+1
        ox = data["allimages"]["image"+str(eachimage)]["spots"][str(eachspot)]["x"]
        oy = data["allimages"]["image"+str(eachimage)]["spots"][str(eachspot)]["y"]
        xydic[ox,oy] = eachimage,eachspot
        x = data["allimages"]["image"+str(eachimage)]["spots"][str(eachspot)]["x"] - data["allimages"]["image"+str(eachimage)]["xcenter"]
        y = -(oy - data["allimages"]["image"+str(eachimage)]["ycenter"])*math.cos(theta*math.pi/180)
        z = -(oy - data["allimages"]["image"+str(eachimage)]["ycenter"])*math.sin(theta*math.pi/180)
        spotlist[spotcount] = x,y,z
        xyzdic[x,y,z] = eachimage,eachspot
        xyzrevdic[eachimage,eachspot] = x,y,z
        ###
        output.write(str(xyzdic[x,y,z])+" "+str(x)+" "+str(y)+" "+str(xyzrevdic[eachimage,eachspot])+"\n")
        ###

# subtract every vector from every other and determine the magnitude of result

magdic = {}    # input magnitiude return the two vectors whose difference gave it
vectors = []
magnitude = {}  # input difference vector returns magnitude

for i in spotlist:
    vector = spotlist[i]
    vectors.append(vector)
for n in vectors:
    for i in vectors:
        vdiff = numpy.subtract(n,i)
        magnitude[str(vdiff)] = float(math.sqrt((vdiff[0])**2+(vdiff[1])**2+(vdiff[2])**2))
        magdic[str(magnitude[str(vdiff)])] = [i,n]

# find the n vector pairs with the smallest magnitude and return which spots gave them

keys = []
mag2spot = {}   # input magnitude returns the two spots that gave it
for i in magdic:
    xyz1 = magdic[i][0]
    xyz2 = magdic[i][1]
    keys.append(float(i))
    mag2spot[float(i)] = xyzdic[xyz1], xyzdic[xyz2]

#print out the n shortest vectors
output.write("\n%s shortest nonzero vectors:\n" % numtoreturn)

returns = 0
shorties = []   #the difference vectors for the shortest vectors in keys
keys.sort()
for i in keys:
    if i > lowthresh and returns <= numtoreturn:
        returns = returns +1    
        shorties.append(mag2spot[i][0])
        if not mag2spot[i][1] in shorties: 
            shorties.append(mag2spot[i][1])
        output.write(str(mag2spot[i])+'\t'+str(i)+"\n")

# for all of  the n smallest vectors r(n) and r(m)  calculate the angle
# cos(theta)= r(n) dot r(m) / magnitude r(n) * magnitude r(m)
# between it and all others

output.write("\ndifference vectors\n")

allvsubs = {}       # Input two vectors output their difference vector
allvrevsubs = {}    # input a difference vector output the two vectors that gave it
for i in shorties:
    vector1 = xyzrevdic[i]  
    for n in shorties:
        vector2 = xyzrevdic[n]
        allvsubs[i,n]  = numpy.subtract(vector1,vector2)
        veckey = str(numpy.subtract(vector1,vector2))
        allvrevsubs[veckey] = i,n
####
for i in allvsubs:
    output.write(str(i)+' '+str(allvsubs[i])+"\n")
####


goodangs = {}
goodlines = []
magnitude = {}          # input vector two spot #'s return difference vector magnitide
for i in allvsubs:
    v1 = allvsubs[i]
    if allvsubs[i][0] != allvsubs[i][1]:
        for n in allvsubs:    
            v2 = allvsubs[n]
            if allvsubs[n][0] != allvsubs[n][1]:
                v12 = numpy.dot(v1,v2)
                v1mag = math.sqrt(v1[0]**2+v1[1]**2+v1[2]**2)
                v2mag = math.sqrt(v2[0]**2+v2[1]**2+v2[2]**2)
                cosphi = (v12)/(v1mag*v2mag)
                angle =(180/math.pi)*math.acos(round(cosphi,12))
                if findmults == 'Y' and  (eab-abthresh < v1mag < eab+abthresh or (2*eab)-abthresh < v1mag < (2*eab)+abthresh or (3*eab)-abthresh < v1mag < (3*eab)+abthresh or (4*eab)-abthresh < v1mag < (4*eab)+abthresh or ec-cthresh < v1mag < ec+cthresh or (2*ec)-cthresh < v1mag < (2*ec)+cthresh or (3*ec)-cthresh < v1mag < (3*ec)+cthresh or (4*ec)-cthresh < v1mag < (4*ec)+cthresh) and (eab-abthresh < v2mag < eab+abthresh or (2*eab)-abthresh < v2mag < (2*eab)+abthresh or (3*eab)-abthresh < v2mag < (3*eab)+abthresh or (4*eab)-abthresh < v2mag < (4*eab)+abthresh) and 87 < angle < 93:
                    goodlines.append([(allvrevsubs[str(v1)]), (allvrevsubs[str(v2)]), str(v12), str(v1mag), str(v2mag), str(cosphi), str(angle)])

                if findmults == 'N' and  (eab-abthresh < v1mag < eab+abthresh or ec-cthresh < v1mag < ec+cthresh) and (eab-abthresh < v2mag < eab+abthresh or ec-cthresh < v2mag < ec+cthresh) and (85 < angle < 95):
                    goodlines.append([allvrevsubs[str(v1)], allvrevsubs[str(v2)], str(v12), str(v1mag), str(v2mag), str(cosphi), str(angle)])
                    goodangs[((allvrevsubs[str(v1)], allvrevsubs[str(v2)]))] = angle   
                    magnitude[allvrevsubs[str(v1)]] = v1mag
                    magnitude[allvrevsubs[str(v2)]] = v2mag
output.write("\n")
output.write("possible unitcell vectors:\ncalculations:\ncos phi = (V12dotV34)/(V12mag*V34mag)\nV12, V34, V12dotV34, V12mag, V34mag, cos phi, angle\n")

for i in goodlines:
    output.write("%s\n" % i)

## find the orthogonal vectors
output.write("\nall orthogonal vector sets:\n")

orthog1 = []
orthog2 = []
orthog3 = []
for i in goodlines:
    if i[0][0] == i[1][0] and i[0][1] !=  i[1][1]:
        if [i[0],i[1]] or [i[1],i[0]] not in orthog1:
            orthog1.append((i[0],i[1]))
for i in orthog1:
    for n in orthog1:
        if i[0] == n[0] and i[1] != n[1]:
            if [i,n] or [n,i] not in orthog2:
                orthog2.append((i,n))

orthogonals = []
for i in orthog1:
    for n in orthog2:
        if (i != n[0]) and (i != n[1]) and (i[0] == n[0][0] or i[0] == n[0][1] or i[1] == n[1][0] or i[1] == n[1][1]) and (i[0] != n[0][1] and i[1] != n[0][0]) and (i[0] != n[1][1] and i[1] != n[1][0]):
            orthogonals.append((i, n[0],n[1]))


for i in orthogonals:
    output.write(str(i)+"\n")
    output.write("%s: angle: %s magnitudes: %s %s \n" % (i[0],goodangs[i[0]],magnitude[i[0][0]], magnitude[i[0][1]]))
    output.write("%s: angle: %s magnitudes: %s %s \n" % (i[1],goodangs[i[1]],magnitude[i[1][0]], magnitude[i[1][1]]))
    output.write("%s: angle: %s magnitudes: %s %s \n\n" % (i[2],goodangs[i[2]],magnitude[i[2][0]], magnitude[i[2][1]]))

output.write("\n**** Full a,b,c sets ****\n")
for i in orthogonals:
    if 5*eab +ec -5*abthresh+cthresh < magnitude[i[0][0]] + magnitude[i[0][1]] + magnitude[i[1][0]] + magnitude[i[1][1]] + magnitude[i[2][0]] + magnitude[i[2][1]] < 5*eab+ec+5*abthresh+cthresh:
        output.write(str(i)+"\n")
        output.write("%s: angle: %s magnitudes: %s %s \n" % (i[0],goodangs[i[0]],magnitude[i[0][0]], magnitude[i[0][1]]))
        output.write("%s: angle: %s magnitudes: %s %s \n" % (i[1],goodangs[i[1]],magnitude[i[1][0]], magnitude[i[1][1]]))
        output.write("%s: angle: %s magnitudes: %s %s \n\n" % (i[2],goodangs[i[2]],magnitude[i[2][0]], magnitude[i[2][1]]))





