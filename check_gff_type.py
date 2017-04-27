# -*- coding: utf-8 -*-
# __author__ = fiona
# time: 2017/3/10 09:34

import re, os, Bio, argparse, sys, fileinput, urllib2

gff = 'F:\\temp\\Homo_sapiens.GRCh38.87.gff3'
out= 'F:\\temp\\Homo_sapiens.GRCh38.87.gff3.type.list'
fw = open(out,'wb')
with open(gff) as fr:
	for line in fr:
		if not re.match(r'^#.+',line.strip()):
			arr = line.strip().split('\t')
			if len(arr) ==9:
				
				fw.write(arr[8]+'\n')
fw.close()
				