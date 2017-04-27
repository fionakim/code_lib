# -*- coding: utf-8 -*-
# __author__ = fiona
# time: 2017/3/31 00:12

import re, os, Bio, argparse, sys, fileinput, urllib2

opts = {'a':1,'b':2}

def get_value_from_dic(d,opts):
    man = opts
    print man
    return d.keys()


if __name__ == '__main__':
    opts = 56789
    opts_test = {1:'a','c':3}
    a = 678
    get_value_from_dic(opts_test,a)