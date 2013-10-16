#!/usr/bin/env python

# uses the csv output from batchtiltspot.py to make an illustration showing the picked spots and their miller indices.
# Useful for troubleshooting.

# Matt Iadanza 2013-07-09

import os


#################
                #
imagename = "2.gif"      # image batch number
                #
#################

# Parameters from param file
with open('2.csv') as pfile:
    lines = pfile.read().splitlines()


spotlist = []
for i in lines:
    values = i.split("\t")
    xy = values[1].split(",")
    spotlist.append("text %s,%s '%s' " % (xy[0],xy[1],values[0]))

print len(spotlist)
command  = ('convert -size 4096x4096 xc:Black -font Arial -pointsize 10 -stroke Blue -strokewidth 2 -fill none -transparent black -draw "%s" +compress points_px.tif' % ' '.join(spotlist))
os.system('%s' % command)
os.system("convert %s %s -gravity center -compose over -composite %s" % (imagename, "points_px.tif","spots.gif"))
os.system('open spots.gif' )
    