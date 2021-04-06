#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import io
import subprocess
import numpy as np
import pandas as pd
import re
import egglib

from optparse import OptionParser


print ''
print 'Takes a fasta file and calculates Theta W, Pi and TajimasD, using the EggLib (may need to install via conda)'
print ''
print 'Mat Beale, Wellcome Sanger Institute, September 2019'
print ''

# Get options
usage = "Usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-f", "--file",action="store", type="string", dest="infile",help="Specify filename of fasta file to be analysed")
parser.add_option("-p", "--path",action="store", type="string", dest="inpath",help="Specify path of input and output files [./]", default='./')
parser.add_option("-m", "--missing",action="store", type="float", dest="mymissing",help="Specify proportion of genomes which need to have the site for it to be included [0.05]",default=0.05)

(options, args) = parser.parse_args()

indir = options.inpath
infile = options.infile
allow_missing = float(options.mymissing)

# pull in alignment and define tests
test_alignment = egglib.io.from_fasta(indir + infile, cls=egglib.Align)
cs = egglib.stats.ComputeStats()
cs.add_stats('S', 'thetaW', 'Pi', 'D', 'lseff', 'nseff')

# function calling egglib to calculate diversity stats
def calc_dev_stats(myaln):
    stats = cs.process_align(myaln,max_missing=allow_missing)
    return [stats['nseff'],stats['lseff'],stats['S'], stats['thetaW'], stats['Pi'], stats['D']]


colnames = ['NSamples','NSites','VSites','ThetaW.count','Pi.count','TajimasD','Pi','ThetaW']
mystatsout = calc_dev_stats(test_alignment)

# egglib seems to output count rather than actual pi, so divide that by Nsites
# need to check that this is actually a number (since can get errors, None, NaN etc)
if isinstance(mystatsout[4],float):
   mystatsout.append(mystatsout[4]/mystatsout[1])
else:
   mystatsout.append(mystatsout[4])

# and same for waterson's estimator (thetaW)
if isinstance(mystatsout[3], float):
   mystatsout.append(mystatsout[3]/mystatsout[1])
else:
   mystatsout.append(mystatsout[3])



# print out to command line
print pd.DataFrame(mystatsout, index =colnames, columns =[infile])

outfile=indir + infile + '.nucleotide-diversity.tsv'
outfile = open(outfile, 'w')
outfile.write('File\tNSamples\tNSites\tVSites\tThetaW\tPi.count\tPi\tTajimasD')
outfile.write('\n')
outfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (infile, mystatsout[0], mystatsout[1], mystatsout[2], round(mystatsout[6],5), round(mystatsout[4],4), round(mystatsout[6],7), round(mystatsout[5],7)))
outfile.write('\n')
outfile.close()

