#! /usr/bin/env python

# predict spots in all images (including non major plain images) based on the unit cell vectors produced by findcell.py and refined using findparallel.py
# outputs batch#.csv that is used by revcoords_tf.py to produce an illustration of the spots.  Uses Brent's new (more accurate) matrix method
#
# Matt Iadanza 2013-07-08 

#### Imports:

import numpy
import math
import os
import json


########### Variables

data = json.load(open('tf_parameters.json')) 
imgsize = data["globals"]["imgsize"]      
boxsize =  data["globals"]["circrad"]
imgmaxres = data["globals"]["imgmaxres"]
reslimit = data["globals"]["reslimit"]
hvals = range(-data["globals"]["hrange"], data["globals"]["hrange"]+1)
kvals = range(-data["globals"]["krange"], data["globals"]["krange"]+1)
lvals = range(-data["globals"]["lrange"], data["globals"]["lrange"]+1)
a = numpy.array(data["globals"]["ucva"])
b = numpy.array(data["globals"]["ucvb"])
c = numpy.array(data["globals"]["ucvc"])
maxresrad = (imgmaxres/reslimit)*(0.5*imgsize)
circrad =  data["globals"]["circrad"]   


################ build the list of all the matrix of hk coords

hklindices = []
for eachhval in kvals:
    for eachkval in hvals:
        for eachlval in lvals:
            hklindices.append((eachhval, eachkval, eachlval))


###############  Get teh images to process from the list

with open('imagelist.txt') as pfile:
    imagestoprocess = pfile.read().splitlines()
print "processing "+str(len(imagestoprocess))+" images"



for eachimage in imagestoprocess:

## get the image specific parametrs

    batch = data[eachimage]["batchnumber"]
    output = file(str(batch)+".csv", "w")
###### calculate the vectors for the reference spots and their dot product:

    xcenter = data[eachimage]["center"][0]
    ycenter = data[eachimage]["center"][1]
    theta =data[eachimage]["theta"]
    refpoints = data[eachimage]["refpoints"]
    bscentx = data[eachimage]["bscenter"][0]
    bscenty = data[eachimage]["bscenter"][1]
    bsrtopy = data[eachimage]["bscenter"][1]-100
    bsrboty = data[eachimage]["bscenter"][1]+100
    bsrad = data[eachimage]["bsrad"]
    imagename = data[eachimage]["giffile"]

    refvectors = []
    for i in refpoints:
        x = i[0] - xcenter
        y = (ycenter - i[1])*math.cos(theta*math.pi/180)
        z = (ycenter - i[1])*math.sin(theta*math.pi/180)
        refvectors.append(numpy.array([x,y,z]).reshape(3,1))

## make the unit matrix

    unitmatrix = numpy.array([[a[0],b[0],c[0]],[a[1],b[1],c[1]],[a[2],b[2],c[2]]])
    uinv =  numpy.linalg.inv(unitmatrix)
    r1 = numpy.dot(uinv, refvectors[0])
    r2 = numpy.dot(uinv, refvectors[1])
    
    r1round = numpy.array([int(round(r1.item((0,0)),0)),int(round(r1.item((1,0)),0)),int(round(r1.item((2,0)),0))])
    r2round = numpy.array([int(round(r2.item((0,0)),0)),int(round(r2.item((1,0)),0)),int(round(r2.item((2,0)),0))])

    checkplane = numpy.cross(r1round,r2round)
    print r1round
    print r2round
    print checkplane
## test every hkl index
    spots = []
    for i in hklindices:
        hkltest = numpy.array([i[0]*checkplane[0],i[1]*checkplane[1],i[2]*checkplane[2]])
        #print i
        #print hkltest
        #print hkltest[0]+hkltest[1]+hkltest[2]
        if  -0.25*numpy.linalg.norm(checkplane) < int(hkltest[0])+int(hkltest[1])+int(hkltest[2]) < numpy.linalg.norm(checkplane)*0.25:
            spots.append((i[0],i[1],i[2])) 
######## calculate the xyz coordinates of all of miller indices
    xydic = {}          # miller index returns point xy
    xyzvectors = {}        # miller index returns x,y,z vector

    for i in spots:
        x = (i[0]*a[0])+(i[1]*b[0])+(i[2]*c[0])
        y = (i[0]*a[1])+(i[1]*b[1])+(i[2]*c[1])
        z = (i[0]*a[2])+(i[1]*b[2])+(i[2]*c[2])
        xyzvectors[i] = numpy.array([x,y,z])
    
        xproj = x + xcenter
        yproj = ycenter - (y/(math.cos(theta*(math.pi/180))))
        xydic[i] = (round(xproj,1),round(yproj,1))
        if xydic[i][1] > 0 and xydic[i][0] > 0:
            if math.sqrt((math.ceil(xydic[i][0]) - xcenter)**2 + (math.ceil(xydic[i][1]) - ycenter)**2) < maxresrad:
                if math.ceil(xydic[i][0]) < bscentx+bsrad or (math.ceil(xydic[i][0]) > bscentx+bsrad and (bsrtopy-circrad > math.ceil(xydic[i][1]) or math.ceil(xydic[i][1]) > bsrboty+circrad)):
                    if (math.ceil(xydic[i][0]) - bscentx)**2 + (math.ceil(xydic[i][1])-bscenty)**2 > (bsrad**2)+circrad:
                        output.write("%s,%s,%s\t%s,%s\t%s\n" % (i[0],i[1],i[2],xydic[i][0],xydic[i][1],math.ceil(xydic[i][0])+(4096*math.ceil(xydic[i][1]))))
## calculate the test vectors
