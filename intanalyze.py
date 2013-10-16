#! /usr/bin/env python

## make an intensity file for brent to play wih in excel

with open('combint.txt') as pfile:
    lines = pfile.read().splitlines()
output = open("intensity_analysis.csv", "w")



for i in lines:
    data = i.split("\t")
    output.write("%s,%s,%s\t%s\t%s\n" % (data[0],data[1],data[2],data[4],data[5]))    