#! /usr/bin/env python



### imports

import json
import os
import sys
import math
import numpy


# user entered variable
##################
abthresh = 2.0  # +- when searching for a and b unit cell vectors
cthresh = 5.0   # +- when searching for c unit cell vectors
acompvec = numpy.array([3.7,-52.9,-16.2])
bcompvec = numpy.array([55.5,3.4,0.9])
ccompvec = numpy.array([0.3,-33.4,109.2])
angthresh = 1.0

eab = (numpy.linalg.norm(acompvec)+numpy.linalg.norm(bcompvec))/2
ec = numpy.linalg.norm(ccompvec)      # expected unitcell c

##################





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
        spotlist[spotcount] = numpy.array([x,y,z])
        xyzdic[x,y,z] = eachimage,eachspot
        xyzrevdic[eachimage,eachspot] = x,y,z
        ###
        output.write(str(xyzdic[x,y,z])+" "+str(x)+" "+str(y)+" "+str(xyzrevdic[eachimage,eachspot])+"\n")
        ###

# subtract every vector from every other and determine the magnitude of result
# compare these to the reference vector

####### Find the A's
vectors = []
twoas = []
threeas = []
fouras = []
aas = []
for i in spotlist:
    vector = spotlist[i]
    vectors.append(vector)
for n in vectors:
    for i in vectors:
        if n[0] != i[0] and n[1] != i[1] and n[2] != i[2]:
            vdiff = numpy.subtract(n,i)
            magnitude = numpy.linalg.norm(vdiff)
            if eab-abthresh< magnitude < eab +abthresh:
                aas.append((vdiff,magnitude))
            if 2*eab-abthresh< magnitude < 2*eab +abthresh:
                twoas.append((vdiff,magnitude))
            if 3*eab-abthresh< magnitude < 3*eab +abthresh:
                threeas.append((vdiff,magnitude))
            if 4*eab-abthresh< magnitude < 4*eab +abthresh:
                fouras.append((vdiff,magnitude))
        ## calculation vs reference
        ## cos(theta)= r(n) dot r(m) / magnitude r(n) * magnitude r(m)
onesxs = []
onesys = []
oneszs = []
twosxs = []
twosys = []
twoszs = []
threesxs = []
threesys = []
threeszs = []
foursxs = []
foursys = []
fourszs = []
for i in aas:
    cosphi = numpy.dot(i[0],acompvec)/(numpy.linalg.norm(i[0])*numpy.linalg.norm(acompvec))
    angle =(180/math.pi)*math.acos(round(cosphi,12))
    if angle < angthresh:
        onesxs.append(i[0][0])
        onesys.append(i[0][1])
        oneszs.append(i[0][2])
for i in twoas:
    cosphi = numpy.dot(i[0],acompvec)/(numpy.linalg.norm(i[0])*numpy.linalg.norm(acompvec))
    angle =(180/math.pi)*math.acos(round(cosphi,12))
    if angle < angthresh:
        twosxs.append(i[0][0]/2)
        twosys.append(i[0][1]/2)
        twoszs.append(i[0][2]/2)
for i in threeas:
    cosphi = numpy.dot(i[0],acompvec)/(numpy.linalg.norm(i[0])*numpy.linalg.norm(acompvec))
    angle =(180/math.pi)*math.acos(round(cosphi,12))
    if angle < angthresh:
        threesxs.append(i[0][0]/3)
        threesys.append(i[0][1]/3)
        threeszs.append(i[0][2]/3)
for i in fouras:
    cosphi = numpy.dot(i[0],acompvec)/(numpy.linalg.norm(i[0])*numpy.linalg.norm(acompvec))
    angle =(180/math.pi)*math.acos(round(cosphi,12))
    if angle < angthresh:
        foursxs.append(i[0][0]/4)
        foursys.append(i[0][1]/4)
        fourszs.append(i[0][2]/4)

print len(onesxs),len(twosxs),len(threesxs),len(foursxs)
if len(onesxs) > 0 and len(onesys) > 0 and len(oneszs) > 0:
    print "1a mean: [%s, %s, %s] %s" % (sum(onesxs)/len(onesxs),sum(onesys)/len(onesys),sum(oneszs)/len(oneszs),numpy.linalg.norm(numpy.array([sum(onesxs)/len(onesxs),sum(onesys)/len(onesys),sum(oneszs)/len(oneszs)])))
if len(twosxs) > 0 and len(twosys) > 0 and len(twoszs) > 0:
    print "2a mean: [%s, %s, %s] %s" % (sum(twosxs)/len(twosxs),sum(twosys)/len(twosys),sum(twoszs)/len(twoszs),numpy.linalg.norm(numpy.array([sum(twosxs)/len(twosxs),sum(twosys)/len(twosys),sum(twoszs)/len(twoszs)])))
if len(threesxs) > 0 and len(threesys) > 0 and len(threeszs) > 0:
    print "3a mean: [%s, %s, %s] %s" % (sum(threesxs)/len(threesxs),sum(threesys)/len(threesys),sum(threeszs)/len(threeszs), numpy.linalg.norm(numpy.array([sum(threesxs)/len(threesxs),sum(threesys)/len(threesys),sum(threeszs)/len(threeszs)])))
if len(foursxs) > 0 and len(foursys) > 0 and len(fourszs) > 0:
    print "4a mean: [%s, %s, %s] %s" % (sum(foursxs)/len(foursxs),sum(foursys)/len(foursys),sum(fourszs)/len(fourszs), numpy.linalg.norm(numpy.array([sum(foursxs)/len(foursxs),sum(foursys)/len(foursys),sum(fourszs)/len(fourszs)])))

if (len(onesxs)+len(twosxs)+len(threesxs)+len(foursxs)) > 0 and (len(onesys)+len(twosys)+len(threesys)+len(foursys)) and (len(oneszs)+len(twoszs)+len(threeszs)+len(fourszs)):
    averagevector = numpy.array([(sum(onesxs)+sum(twosxs)+sum(threesxs)+sum(foursxs))/(len(onesxs)+len(twosxs)+len(threesxs)+len(foursxs)),(sum(onesys)+sum(twosys)+sum(threesys)+sum(foursys))/(len(onesys)+len(twosys)+len(threesys)+len(foursys)),(sum(oneszs)+sum(twoszs)+sum(threeszs)+sum(fourszs))/(len(oneszs)+len(twoszs)+len(threeszs)+len(fourszs))])
    print "%s %s" % (averagevector, numpy.linalg.norm(averagevector))


###### FIND THE B's
twobs = []
threebs = []
fourbs = []
bbs = []
for n in vectors:
    for i in vectors:
        if n[0] != i[0] and n[1] != i[1] and n[2] != i[2]:
            vdiff = numpy.subtract(n,i)
            magnitude = numpy.linalg.norm(vdiff)
            if eab-abthresh< magnitude < eab +abthresh:
                bbs.append((vdiff,magnitude))
            if 2*eab-abthresh< magnitude < 2*eab +abthresh:
                twobs.append((vdiff,magnitude))
            if 3*eab-abthresh< magnitude < 3*eab +abthresh:
                threebs.append((vdiff,magnitude))
            if 4*eab-abthresh< magnitude < 4*eab +abthresh:
                fourbs.append((vdiff,magnitude))
## calculation vs reference
## cos(theta)= r(n) dot r(m) / magnitude r(n) * magnitude r(m)
onesxs = []
onesys = []
oneszs = []
twosxs = []
twosys = []
twoszs = []
threesxs = []
threesys = []
threeszs = []
foursxs = []
foursys = []
fourszs = []
for i in bbs:
    cosphi = numpy.dot(i[0],bcompvec)/(numpy.linalg.norm(i[0])*numpy.linalg.norm(bcompvec))
    angle =(180/math.pi)*math.acos(round(cosphi,12))
    if angle < angthresh:
        onesxs.append(i[0][0])
        onesys.append(i[0][1])
        oneszs.append(i[0][2])
for i in twobs:
    cosphi = numpy.dot(i[0],bcompvec)/(numpy.linalg.norm(i[0])*numpy.linalg.norm(bcompvec))
    angle =(180/math.pi)*math.acos(round(cosphi,12))
    if angle < angthresh:
        twosxs.append(i[0][0]/2)
        twosys.append(i[0][1]/2)
        twoszs.append(i[0][2]/2)
for i in threebs:
    cosphi = numpy.dot(i[0],bcompvec)/(numpy.linalg.norm(i[0])*numpy.linalg.norm(bcompvec))
    angle =(180/math.pi)*math.acos(round(cosphi,12))
    if angle < angthresh:
        threesxs.append(i[0][0]/3)
        threesys.append(i[0][1]/3)
        threeszs.append(i[0][2]/3)
for i in fourbs:
    cosphi = numpy.dot(i[0],bcompvec)/(numpy.linalg.norm(i[0])*numpy.linalg.norm(bcompvec))
    angle =(180/math.pi)*math.acos(round(cosphi,12))
    if angle < angthresh:
        foursxs.append(i[0][0]/4)
        foursys.append(i[0][1]/4)
        fourszs.append(i[0][2]/4)

print len(onesxs),len(twosxs),len(threesxs),len(foursxs)
if len(onesxs) > 0 and len(onesys) > 0 and len(oneszs) > 0:
    print "1b mean: [%s, %s, %s] %s" % (sum(onesxs)/len(onesxs),sum(onesys)/len(onesys),sum(oneszs)/len(oneszs),numpy.linalg.norm(numpy.array([sum(onesxs)/len(onesxs),sum(onesys)/len(onesys),sum(oneszs)/len(oneszs)])))
if len(twosxs) > 0 and len(twosys) > 0 and len(twoszs) > 0:
    print "2b mean: [%s, %s, %s] %s" % (sum(twosxs)/len(twosxs),sum(twosys)/len(twosys),sum(twoszs)/len(twoszs),numpy.linalg.norm(numpy.array([sum(twosxs)/len(twosxs),sum(twosys)/len(twosys),sum(twoszs)/len(twoszs)])))
if len(threesxs) > 0 and len(threesys) > 0 and len(threeszs) > 0:
    print "3b mean: [%s, %s, %s] %s" % (sum(threesxs)/len(threesxs),sum(threesys)/len(threesys),sum(threeszs)/len(threeszs), numpy.linalg.norm(numpy.array([sum(threesxs)/len(threesxs),sum(threesys)/len(threesys),sum(threeszs)/len(threeszs)])))
if len(foursxs) > 0 and len(foursys) > 0 and len(fourszs) > 0:
    print "4b mean: [%s, %s, %s] %s" % (sum(foursxs)/len(foursxs),sum(foursys)/len(foursys),sum(fourszs)/len(fourszs), numpy.linalg.norm(numpy.array([sum(foursxs)/len(foursxs),sum(foursys)/len(foursys),sum(fourszs)/len(fourszs)])))

if (len(onesxs)+len(twosxs)+len(threesxs)+len(foursxs)) > 0 and (len(onesys)+len(twosys)+len(threesys)+len(foursys)) and (len(oneszs)+len(twoszs)+len(threeszs)+len(fourszs)):
    averagevector = numpy.array([(sum(onesxs)+sum(twosxs)+sum(threesxs)+sum(foursxs))/(len(onesxs)+len(twosxs)+len(threesxs)+len(foursxs)),(sum(onesys)+sum(twosys)+sum(threesys)+sum(foursys))/(len(onesys)+len(twosys)+len(threesys)+len(foursys)),(sum(oneszs)+sum(twoszs)+sum(threeszs)+sum(fourszs))/(len(oneszs)+len(twoszs)+len(threeszs)+len(fourszs))])
    print "%s %s" % (averagevector, numpy.linalg.norm(averagevector))


###### FIND THE C's
twocs = []
threecs = []
fourcs = []
ccs = []
for n in vectors:
    for i in vectors:
        if n[0] != i[0] and n[1] != i[1] and n[2] != i[2]:
            vdiff = numpy.subtract(n,i)
            magnitude = numpy.linalg.norm(vdiff)
            if ec-cthresh< magnitude < ec +cthresh:
                ccs.append((vdiff,magnitude))
            if 2*ec-cthresh< magnitude < 2*ec +cthresh:
                twocs.append((vdiff,magnitude))
            if 3*ec-cthresh< magnitude < 3*ec +cthresh:
                threecs.append((vdiff,magnitude))
            if 4*ec-cthresh< magnitude < 4*ec +cthresh:
                fourcs.append((vdiff,magnitude))
        ## calculation vs reference
        ## cos(theta)= r(n) dot r(m) / magnitude r(n) * magnitude r(m)
onesxs = []
onesys = []
oneszs = []
twosxs = []
twosys = []
twoszs = []
threesxs = []
threesys = []
threeszs = []
foursxs = []
foursys = []
fourszs = []
for i in ccs:
    cosphi = numpy.dot(i[0],ccompvec)/(numpy.linalg.norm(i[0])*numpy.linalg.norm(ccompvec))
    angle =(180/math.pi)*math.acos(round(cosphi,12))
    if angle < angthresh:
        onesxs.append(i[0][0])
        onesys.append(i[0][1])
        oneszs.append(i[0][2])
for i in twocs:
    cosphi = numpy.dot(i[0],ccompvec)/(numpy.linalg.norm(i[0])*numpy.linalg.norm(ccompvec))
    angle =(180/math.pi)*math.acos(round(cosphi,12))
    if angle < angthresh:
        twosxs.append(i[0][0]/2)
        twosys.append(i[0][1]/2)
        twoszs.append(i[0][2]/2)
for i in threecs:
    cosphi = numpy.dot(i[0],ccompvec)/(numpy.linalg.norm(i[0])*numpy.linalg.norm(ccompvec))
    angle =(180/math.pi)*math.acos(round(cosphi,12))
    if angle < angthresh:
        threesxs.append(i[0][0]/3)
        threesys.append(i[0][1]/3)
        threeszs.append(i[0][2]/3)
for i in fourcs:
    cosphi = numpy.dot(i[0],ccompvec)/(numpy.linalg.norm(i[0])*numpy.linalg.norm(ccompvec))
    angle =(180/math.pi)*math.acos(round(cosphi,12))
    if angle < angthresh:
        foursxs.append(i[0][0]/4)
        foursys.append(i[0][1]/4)
        fourszs.append(i[0][2]/4)

print len(onesxs),len(twosxs),len(threesxs),len(foursxs)
if len(onesxs) > 0 and len(onesys) > 0 and len(oneszs) > 0:
    print "1c mean: [%s, %s, %s] %s" % (sum(onesxs)/len(onesxs),sum(onesys)/len(onesys),sum(oneszs)/len(oneszs),numpy.linalg.norm(numpy.array([sum(onesxs)/len(onesxs),sum(onesys)/len(onesys),sum(oneszs)/len(oneszs)])))
if len(twosxs) > 0 and len(twosys) > 0 and len(twoszs) > 0:
    print "2c mean: [%s, %s, %s] %s" % (sum(twosxs)/len(twosxs),sum(twosys)/len(twosys),sum(twoszs)/len(twoszs),numpy.linalg.norm(numpy.array([sum(twosxs)/len(twosxs),sum(twosys)/len(twosys),sum(twoszs)/len(twoszs)])))
if len(threesxs) > 0 and len(threesys) > 0 and len(threeszs) > 0:
    print "3c mean: [%s, %s, %s] %s" % (sum(threesxs)/len(threesxs),sum(threesys)/len(threesys),sum(threeszs)/len(threeszs), numpy.linalg.norm(numpy.array([sum(threesxs)/len(threesxs),sum(threesys)/len(threesys),sum(threeszs)/len(threeszs)])))
if len(foursxs) > 0 and len(foursys) > 0 and len(fourszs) > 0:
    print "4c mean: [%s, %s, %s] %s" % (sum(foursxs)/len(foursxs),sum(foursys)/len(foursys),sum(fourszs)/len(fourszs), numpy.linalg.norm(numpy.array([sum(foursxs)/len(foursxs),sum(foursys)/len(foursys),sum(fourszs)/len(fourszs)])))

if (len(onesxs)+len(twosxs)+len(threesxs)+len(foursxs)) > 0 and (len(onesys)+len(twosys)+len(threesys)+len(foursys)) and (len(oneszs)+len(twoszs)+len(threeszs)+len(fourszs)):
    averagevector = numpy.array([(sum(onesxs)+sum(twosxs)+sum(threesxs)+sum(foursxs))/(len(onesxs)+len(twosxs)+len(threesxs)+len(foursxs)),(sum(onesys)+sum(twosys)+sum(threesys)+sum(foursys))/(len(onesys)+len(twosys)+len(threesys)+len(foursys)),(sum(oneszs)+sum(twoszs)+sum(threeszs)+sum(fourszs))/(len(oneszs)+len(twoszs)+len(threeszs)+len(fourszs))])
    print "%s %s" % (averagevector, numpy.linalg.norm(averagevector))
