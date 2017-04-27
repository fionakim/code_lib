# -*- coding: utf-8 -*-
# __author__ = fiona
# time: 2017/3/10 18:04

import re, os, Bio, argparse, sys, fileinput, urllib2
import Bio
import subprocess
from bs4 import BeautifulSoup
from file import File
from gff3_file import Gff3File



class FastaFile(File):
	def __init__(self):
		self._seq_ids = []
		self._seq_obj = []
		self._seq_id_len_dic = []
		self.seq_type = ''  # 有 prot nucl两种类型
		self._anno_gff3 = Gff3File()
