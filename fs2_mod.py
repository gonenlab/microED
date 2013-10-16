#! /usr/bin/env python
## imports

## Finds the average unit cell vector using find cell_params.json which is created by user or from imageJ files and i2cfiles.py

# Matt Iadanza 2013-07-05

import json
import os
import sys
import math
import numpy


athresh = 0.1 # angle threshold
eab = 55.0 # expected a and b ucell dimensions
abthresh = 1 
ec = 112.0 # expected c dimension
cthresh = 2
eabg = 90 # expected alpha=beta=gamma angles

# get some variables:
data = json.load(open('cellfind_params.json'))
numberofimages = len(data["allimages"])
output = open("output_cellfind.txt", "w")


print ""
print "make the spolist dictionarys"
spotlist= {}   # input raw spot number returns numpy array with vector 
spotcount = 0
xydic = {}  # original x,y coords to image#,spot#
xyzdic = {} # calculated x,y,z coords to image#,spot#
xyzrevdic = {}  # input imag#,spot# returns calculated x,z,y as vector
spotdic  = {} # input raw spot number  (1 through  total spots) return image,spotno
revspotdic = {} # input image, spot no  return raw spot no
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
        spotlist[spotcount] = numpy.array([x,y,z])
        spotdic[spotcount] = eachimage,eachspot
        revspotdic[eachimage,eachspot] = spotcount
        xyzdic[x,y,z] = eachimage,eachspot
        xyzrevdic[eachimage,eachspot] = numpy.array([x,y,z])
        ###
        output.write(str(xyzdic[x,y,z])+" "+str(ox)+" "+str(oy)+" "+str(xyzrevdic[eachimage,eachspot])+"\n")
        ###

print "%s spots picked" % len(xyzrevdic)


vectors = []
magnitude = {}  # input difference vector returns magnitude
twovsdiffdic= {} # input the spot numbers for two vectors eturns their difference vector IE (1,2) (3,4),(5,6)(78)
done = []
vdiff = {}  ## input raw spot numbers return the difference vector for the two

count = []
print"subtract every vector from every other and determine the magnitude of result"
for n in spotlist:
    for i in spotlist:
        magnitude[spotdic[n],spotdic[i]] = numpy.linalg.norm(numpy.subtract(spotlist[n],spotlist[i]))
        if i != n and (eab-abthresh < numpy.linalg.norm(numpy.subtract(spotlist[n],spotlist[i])) < eab+abthresh or ec-cthresh < numpy.linalg.norm(numpy.subtract(spotlist[n],spotlist[i])) < ec+cthresh):
            vdiff[n,i] = numpy.subtract(spotlist[n],spotlist[i])
            

print "%s possible unit cell vectors found" % len(vdiff)
print "calculate the angles between possible unit cell vectors "
output.write("/nunit cell vectors orthoganl to each other\n")
popairs = []
for i in vdiff:
    for n in vdiff:
        if i != n:
            v12 = numpy.dot(vdiff[i],vdiff[n])
            v1mag = numpy.linalg.norm(i)
            v2mag = numpy.linalg.norm(n)
            cosphi = (v12)/(v1mag*v2mag)
            if 1 >= cosphi >= -1:
                angle =(180/math.pi)*math.acos(round(cosphi,12))
                if eabg-athresh < angle < eabg+athresh:
                    output.write("%s %s %s %s %s\n" % ((spotdic[i[0]],spotdic[i[1]]),(spotdic[n[0]],spotdic[n[1]]),angle,magnitude[(spotdic[i[0]],spotdic[i[1]])],magnitude[(spotdic[n[0]],spotdic[n[1]])]))
                    popairs.append(((spotdic[i[0]],spotdic[i[1]]),(spotdic[n[0]],spotdic[n[1]]),angle,magnitude[(spotdic[i[0]],spotdic[i[1]])],magnitude[(spotdic[n[0]],spotdic[n[1]])]))



print "%s %s-degree pairs found" % (len(popairs),eabg)
print "find orhogonal groups"
output.write("\northogdnal groups\n")
pogroups = {}  # input a vector pair and return all pairs that are orhogonal to it
for i in popairs:
    pohold = []
    for n in popairs:
        if i[0] == n[0] and i[1] != n[1]:
            pohold.append(n[1])
    pogroups[i[0]] = pohold

print "%s possible orthogonal triplets found " % len(pogroups)

allas = []
allbs = []
allcs = []
otriplets = {}      ## input number (arbitrary just based on count) return a unique orthogonal triplet (but does include combos)
tripcount = 1
for i in pogroups:
    idiff = numpy.subtract(xyzrevdic[i[0]],xyzrevdic[i[1]])
    for n in pogroups[i]:
        v1 = numpy.subtract(xyzrevdic[n[0]],xyzrevdic[n[1]]) 
        for k in pogroups[i]:
            v2 = numpy.subtract(xyzrevdic[k[0]],xyzrevdic[k[1]])
            v12 = numpy.dot (v1,v2)
            v1mag = numpy.linalg.norm(v1)
            v2mag = numpy.linalg.norm(v2)
            cosphi = (v12)/(v1mag*v2mag)
            if 1 >= cosphi >= -1:
                angle =(180/math.pi)*math.acos(round(cosphi,12))
                if eabg-athresh < angle < eabg + athresh and magnitude[i]+magnitude[k]+magnitude[n] < 2*eab+ec+2*abthresh+cthresh:
                    output.write('\n%s\t%s\t%s\t%s\t%s\t%s\n' %(i,n,k,magnitude[i],magnitude[n],magnitude[k]))
                    output.write('%s\t%s\t%s\n' % (idiff,v1,v2))
                    otriplets[tripcount] = ((i,n,k,magnitude[i],magnitude[n],magnitude[k],idiff,v1,v2))  
                    tripcount = tripcount+1

#### identify which is a b and c:

for i in otriplets:
    if otriplets[i][3] > otriplets[i][4] and otriplets[i][3] > otriplets[i][5]:
        c = otriplets[i][6]
        ab1 = otriplets[i][7]
        ab2 = otriplets[i][8]
    if otriplets[i][4] > otriplets[i][3] and otriplets[i][4] > otriplets[i][5]:
        c = otriplets[i][7]
        
        ab1 = otriplets[i][6]
        ab2 = otriplets[i][8]
    if otriplets[i][5] > otriplets[i][4] and otriplets[i][5] > otriplets[i][3]:
        c = otriplets[i][8]
        ab1 = otriplets[i][6]
        ab2 = otriplets[i][7]
    if c[2] < 0:
            c = -c
    if ab1[0] < 0:
        ab1 = -ab1
    if ab2[0] < 0:
        ab2 = -ab2
    cross = numpy.cross(ab1,ab2)
    if cross[2] > 0:
        b = ab1 
        a = ab2 
    if cross[2] < 0:
        a = ab1
        b = ab2
    allas.append(a)
    allbs.append(b)
    allcs.append(c)


isum = 0
jsum = 0
ksum = 0
print "Average unit cell vectors:"
output.write("\nprediced unitcell vectors\n")
if len(allas) >0:    
    for i in allas:
        isum = isum + i[0]
        jsum = jsum + i[1]
        ksum = ksum + i[2]
    print"a: [%s, %s, %s] Magnitude: %s" % (isum/len(allas),jsum/len(allas),ksum/len(allas),numpy.linalg.norm([isum/len(allas),jsum/len(allas),ksum/len(allas)]))
    output.write("a: [%s, %s, %s] Magnitude: %s\n" % (isum/len(allas),jsum/len(allas),ksum/len(allas),numpy.linalg.norm([isum/len(allas),jsum/len(allas),ksum/len(allas)])))

isum = 0
jsum = 0
ksum = 0
if len(allbs) >0:
    for i in allbs:
        isum = isum + i[0]
        jsum = jsum + i[1]
        ksum = ksum + i[2]
    print"b: [%s, %s, %s] Magnitude: %s" % (isum/len(allbs),jsum/len(allbs),ksum/len(allbs),numpy.linalg.norm([isum/len(allbs),jsum/len(allbs),ksum/len(allbs)]))
    output.write("b: [%s, %s, %s] Magnitude: %s\n" % (isum/len(allbs),jsum/len(allbs),ksum/len(allbs),numpy.linalg.norm([isum/len(allbs),jsum/len(allbs),ksum/len(allbs)])))
isum = 0
jsum = 0
ksum = 0
if len(allcs) > 0:
    for i in allcs:
        isum = isum + i[0]
        jsum = jsum + i[1]
        ksum = ksum + i[2]
    print"c: [%s, %s, %s] Magnitude: %s" % (isum/len(allcs),jsum/len(allcs),ksum/len(allcs),numpy.linalg.norm([isum/len(allcs),jsum/len(allcs),ksum/len(allcs)]))
    output.write("c: [%s, %s, %s] Magnitude: %s\n" % (isum/len(allcs),jsum/len(allcs),ksum/len(allcs),numpy.linalg.norm([isum/len(allcs),jsum/len(allcs),ksum/len(allcs)])))

if len(allas) < 0 or len(allbs) < 0 or len(allcs) < 0:
    print "no orthogonal triplets found. Add more spots and try again"
print "checksum"
print len(allas), len(allbs), len(allcs)
print len(otriplets)
