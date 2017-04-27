# -*- coding: utf-8 -*-
# __author__ = fiona
# time: 2017/3/10 18:05

import re, os, Bio, argparse, sys, fileinput, urllib2
import unittest

class File(object):
	def __init__(self):
		self._path = ''
		self._name = ''
	
	def set_path(self, path):
		if not isinstance(path, str):
			raise Exception('文件路径不为字符串')
			self._path = path
		
			
	
	def path(self):
		return self._path
	
	def name(self):
		return self._name
	
	def check(self):
		if os.path.exists(self._path):
			return True
		else:
			return False