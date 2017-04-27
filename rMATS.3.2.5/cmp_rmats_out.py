# -*- coding: utf-8 -*-
# __author__ = fiona
# time: 2017/4/17 08:51

import re, os, Bio, argparse, sys, fileinput, urllib2
'''
结论是all——id 记录所有的as events
'''

import pandas
from pandas import *
all_file = 'I:\\rmats-out-eg\\ASEvents\\fromGTF.A3SS.txt'
novel_file = 'I:\\rmats-out-eg\\ASEvents\\fromGTF.novelEvents.A3SS.txt'
jc_only_file = 'I:\\rmats-out-eg\\MATS_output\\A3SS.MATS.JunctionCountOnly.txt'
jc_on_target_file = 'I:\\rmats-out-eg\\MATS_output\\A3SS.MATS.JunctionCountOnly.txt'

all_table = read_table(all_file,sep='\t')
all_id = set(all_table['ID'].values)
novel_table = read_table(novel_file,sep='\t')
novel_id = set(novel_table['ID'].values)
jc_only_table = read_table(jc_only_file,sep='\t')
jc_only_id = set(jc_only_table['ID'].values)
jc_on_target_id = set(read_table(jc_on_target_file,sep='\t')['ID'].values)

print len(all_id)
print len(novel_id)
print len(jc_only_id)
print len(jc_on_target_id)
print len(jc_on_target_id &jc_only_id)
print jc_on_target_id == jc_only_id
print len(all_id - jc_only_id)
# print len(all_id-novel_id)
# print len(novel_id &jc_only_id)
# print len((novel_id | jc_only_id) - all_id)
print all_id - jc_only_id

