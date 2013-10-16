#!/usr/bin/env python

###################
# finds pixel values in circles around points specified in coords.csv generated
# by lauecalc.py or batchtiltspot.py

# Matt Iadanza 2012-05-30


###################imports
from PIL import Image
import os
import csv
import math
import sys
import numpy
import json

## get global parameters from parfile.

data = json.load(open('tf_parameters.json')) 
imgsize = data["globals"]["imgsize"]      
boxsize =  data["globals"]["circrad"]   
radius = range(1,boxsize+1)            

## get list of images to process

with open('imagelist.txt') as pfile:
    imagestoprocess = pfile.read().splitlines()
print "processing "+str(len(imagestoprocess))+" images"
################### Start processing images


for eachimage in imagestoprocess:

## get the image specific parametrs

    batch = data[eachimage]["batchnumber"]
    spotthreshold = data[eachimage]["integrationthresh"]
    
## open the image and prepare dicts for processing

    im = Image.open(eachimage)
    imgdata = list(im.getdata())
    imgarray = numpy.asarray(imgdata, dtype=numpy.uint16)



    indexdic = {}
    coordsdic = {}
    points = []
    inapix = {}
    values = csv.reader(open(str(batch)+'.csv', 'rb'), delimiter='\t')
    for row in values:
        points.append(row[2])
        indexdic[row[2]] = row[0]
        coordsdic[row[2]] = row[1]
        inapix[row[0]] = row[2]

## calculate the mask used for background subtraction  mc = mask calc.

    mcx = 0
    mcy = -boxsize
    subtracts = {}      # how much to subtract from either side of each line to make a circle

    for n in range(-boxsize,boxsize+1):
        mcx = math.sqrt(math.pow(mcx,2)-(2*mcy)-1)
        mcy = n
        subtracts[n] = (boxsize - int(mcx))
 

############ get yer DATA!
### calculate background mean intensity for each spot individually
# above the center line
    rawoutput = file('rawout.txt', 'w')
    ispotmbg = {}       # Individal SPOT Mean BackGround 
    bglinecountrunning = 0     
    bgintensityrunning = 0
    for eachpoint in points:
        intpoint = int(float(eachpoint))    # the point (spot center)
    
        for eachline in list(reversed(radius)):
            linecount= len(imgarray[intpoint-((imgsize*eachline)-boxsize):intpoint-((imgsize*eachline)-boxsize-subtracts[eachline])])+len(imgarray[intpoint+(imgsize*eachline)+boxsize-subtracts[eachline]:intpoint+(imgsize*eachline+boxsize)])
            rawoutput.write(str(imgarray[intpoint-((imgsize*eachline)-boxsize):intpoint-((imgsize*eachline)-boxsize-subtracts[eachline])])+str(imgarray[intpoint+(imgsize*eachline)+boxsize-subtracts[eachline]:intpoint+(imgsize*eachline+boxsize)])+"------ BGlinesum: "+str(sum(imgarray[intpoint-((imgsize*eachline)-boxsize):intpoint-((imgsize*eachline)-boxsize-subtracts[eachline])])+sum(imgarray[intpoint+(imgsize*eachline)+boxsize-subtracts[eachline]:intpoint+(imgsize*eachline+boxsize)]))+','+str(linecount)+'\n')
            bglinecountrunning = bglinecountrunning + linecount
            bgintensityrunning = bgintensityrunning + sum(imgarray[intpoint-((imgsize*eachline)-boxsize):intpoint-((imgsize*eachline)-boxsize-subtracts[eachline])])+sum(imgarray[intpoint+(imgsize*eachline)+boxsize-subtracts[eachline]:intpoint+(imgsize*eachline+boxsize)])

# don't have to worry about the center line because it always = 0 subtracts
    

# below the center line
        for eachline in radius:
            linecount= len(imgarray[intpoint-((imgsize*eachline)-boxsize):intpoint-((imgsize*eachline)-boxsize-subtracts[eachline])])+len(imgarray[intpoint+(imgsize*eachline)+boxsize-subtracts[eachline]:intpoint+(imgsize*eachline+boxsize)])
            rawoutput.write(str(imgarray[intpoint-((imgsize*eachline)-boxsize):intpoint-((imgsize*eachline)-boxsize-subtracts[eachline])])+str(imgarray[intpoint+(imgsize*eachline)+boxsize-subtracts[eachline]:intpoint+(imgsize*eachline+boxsize)])+"------ BGlinesum: "+str(sum(imgarray[intpoint-((imgsize*eachline)-boxsize):intpoint-((imgsize*eachline)-boxsize-subtracts[eachline])])+sum(imgarray[intpoint+(imgsize*eachline)+boxsize-subtracts[eachline]:intpoint+(imgsize*eachline+boxsize)]))+','+str(linecount)+'\n')
            bglinecountrunning = bglinecountrunning + linecount
            bgintensityrunning = bgintensityrunning + sum(imgarray[intpoint-((imgsize*eachline)-boxsize):intpoint-((imgsize*eachline)-boxsize-subtracts[eachline])])+sum(imgarray[intpoint+(imgsize*eachline)+boxsize-subtracts[eachline]:intpoint+(imgsize*eachline+boxsize)])
        rawoutput.write(indexdic[eachpoint]+' - '+coordsdic[eachpoint]+'\n')
    
# calc BG intensity for this spot and push to dic 
        ispotmbg[indexdic[eachpoint]] = math.ceil((100*float(bgintensityrunning)/float(bglinecountrunning)))/100
        bglinecountrunning = 0     
        bgintensityrunning = 0

### calculate the intensities of the spots themselves

    bgsubstracted = []
    meani ={}           # MEAN Intensitys for all spots
    inti = {}           # integrated intensitys for all spots
    for eachpoint in points:
        intpoint = int(float(eachpoint))    # the point (spot center)


# varable for tracking number of pixels chosen from each line
        linecountrunning = 0
        intensityrunning = 0

####### For each spot
## get pixels, BG subtract, and sum the n lines before:

        for eachline in list(reversed(radius)):
            linecount= len(imgarray[intpoint-(imgsize*eachline)-boxsize+subtracts[eachline]:intpoint-(imgsize*eachline)+boxsize-subtracts[eachline]])
            rawoutput.write(str(imgarray[intpoint-(imgsize*eachline)-boxsize+subtracts[eachline]:intpoint-(imgsize*eachline)+boxsize-subtracts[eachline]])+"------ linesum: "+str(sum(imgarray[intpoint-(imgsize*eachline)-boxsize+subtracts[eachline]:intpoint-(imgsize*eachline)+boxsize-subtracts[eachline]]))+','+str(linecount) +'  '+str(ispotmbg[indexdic[eachpoint]])+'\n')
        
            bgsubtracted = imgarray[intpoint-(imgsize*eachline)-boxsize+subtracts[eachline]:intpoint-(imgsize*eachline)+boxsize-subtracts[eachline]]
            bgsubtracted = [x - ispotmbg[indexdic[eachpoint]] for x in bgsubtracted]
            linecountrunning = linecountrunning + linecount
            intensityrunning = intensityrunning + sum(bgsubtracted[:])


    ## get pixels, BG subtract, and sum the line containing the center pixel
        linecount = len(imgarray[intpoint-(imgsize*eachline)-boxsize+subtracts[eachline]:intpoint-(imgsize*eachline)+boxsize-subtracts[eachline]])
        rawoutput.write(str(imgarray[intpoint-boxsize:intpoint+boxsize])+"------ linesum: "+str(sum(imgarray[intpoint-boxsize:intpoint+boxsize]))+','+str(linecount)+'     '+indexdic[eachpoint]+' - '+coordsdic[eachpoint]+'  '+str(eachpoint)+'\n')
        bgsubtracted = imgarray[intpoint-(imgsize*eachline)-boxsize+subtracts[eachline]:intpoint-(imgsize*eachline)+boxsize-subtracts[eachline]]
        bgsubtracted = [x - ispotmbg[indexdic[eachpoint]] for x in bgsubtracted]
        linecountrunning = linecountrunning + linecount
        intensityrunning = intensityrunning + sum(bgsubtracted[:])

## get pixels, BG subtract, and sum the n lines after 
        for eachline in radius:
            linecount = len(imgarray[intpoint-(imgsize*eachline)-boxsize+subtracts[eachline]:intpoint-(imgsize*eachline)+boxsize-subtracts[eachline]])
            rawoutput.write(str(imgarray[intpoint+(imgsize*eachline)-boxsize+subtracts[eachline]:intpoint+(imgsize*eachline)+boxsize-subtracts[eachline]])+"------ raw: "+str(sum(imgarray[intpoint+(imgsize*eachline)-boxsize+subtracts[eachline]:intpoint+(imgsize*eachline)+boxsize-subtracts[eachline]]))+','+str(linecount)+'\n')
            bgsubtracted = imgarray[intpoint+(imgsize*eachline)-boxsize+subtracts[eachline]:intpoint+(imgsize*eachline)+boxsize-subtracts[eachline]]
            bgsubtracted = [x - ispotmbg[indexdic[eachpoint]] for x in bgsubtracted]
            linecountrunning = linecountrunning + linecount
            intensityrunning = intensityrunning + sum(bgsubtracted[:])   
        if float(intensityrunning) > 0:
	    meani[indexdic[eachpoint]] = float(intensityrunning)/float(linecountrunning)
	    inti[indexdic[eachpoint]] = float(intensityrunning)

# split up the hkl indices in meani into h, k, and l, dictionaries
# then calculate isym

    h = {}
    k = {}
    l = {}
    isym = {}
    for eachcoords in meani:
        line = []
        line = eachcoords.split(',')
        h[eachcoords] = line[0]
        k[eachcoords] = line[1]
        l[eachcoords] = line[2]
        if int(l[eachcoords]) > 0:
            isym[eachcoords] = 1
        if int(l[eachcoords]) < 0:
            isym[eachcoords] = 2
        if int(l[eachcoords]) == 0 and int(h[eachcoords]) > 0:
            isym[eachcoords] = 1
        if int(l[eachcoords]) == 0 and int(h[eachcoords]) < 0:
            isym[eachcoords] = 2
        if int(l[eachcoords]) == 0 and int(h[eachcoords]) == 0 and int(k[eachcoords]) > 0:
            isym[eachcoords] = 1
        if int(l[eachcoords]) == 0 and int(h[eachcoords]) == 0 and int(k[eachcoords]) < 0:
            isym[eachcoords] = 2

###### Print the results for all spots

    passcount = 0
    output = file(str(batch)+".int", 'w')
    spotlist = []
    for i in inti:
        if inti[i] >= spotthreshold:
        	passcount = passcount+1
        	output.write(str(h[i])+'\t'+str(k[i])+'\t'+str(l[i])+'\t'+str(isym[i])+'\t'+str(batch)+'\t'+str(inti[i])+'\t'+str(math.sqrt(inti[i]))+'\n')
        	spotlist.append(inapix[i])
    print ""
    print eachimage
    print "total spots: "+str(len(points))
    print "spots w/ + intensity: "+str(len(inti))  
    print "> "+str(spotthreshold)+" threshold: "+str(passcount)
    
    
### make the liiustration
    drawlist = []
    for i in spotlist:
        drawlist.append("circle "+str(int(float(i)%4096)) +","+str(int(float(i)/4096))+" "+str(int(float(i)%4096)+boxsize) +","+str(int(float(i)/4096)+boxsize) )
    command  = ('convert -size 4096x4096 xc:Black -stroke Green -strokewidth 2 -fill none -transparent black -draw "%s" points_px.gif' % ' '.join(drawlist))
    os.system('%s' % command)
    os.system("convert %s %s -gravity center -compose over -composite %s" % (str(batch)+".gif", "points_px.gif",str(batch)+"_spots.gif"))

## cleanup

os.system('rm points_px.gif')