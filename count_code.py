# -*- coding: utf-8 -*-
# __author__ = fiona
# time: 2017/3/29 12:21

import re, os, Bio, argparse, sys, fileinput, urllib2, subprocess
from collections import defaultdict

if __name__ == '__main__':
    gtf = '/mnt/ilustre/users/sanger-dev/workspace/20170214/Single_assembly_module_tophat_cufflinks/Assembly/output/Cuffmerge/cuffmerge_gffcompare/merged.gtf'
    out = '/mnt/ilustre/users/sanger-dev/workspace/20170214/Single_assembly_module_tophat_cufflinks/Assembly/output/Cuffmerge/cuffmerge_gffcompare/merged_count.gtf'
    txpt_cls_cmd = ' awk -F \';\' \'{printf $2\"\;\"$(NF-2)\"\\n\";}\'  %s | uniq ' % (gtf)
    cls_cmd = 'awk -F \';\' \'{printf $(NF-2)"\\n";}\' %s | awk -F \'\"\' \'{printf $2"\\n";}\'|uniq |sort |uniq' % (gtf)
    
    txpt_cls_content = str(subprocess.check_output(txpt_cls_cmd, shell=True)).split('\n')
    cls_content = subprocess.check_output(cls_cmd, shell=True).split('\n')
    print(cls_content)
    
    cls_txpt_set_dic = {}
    for cls in cls_content:
        cls = cls.strip()
        if cls:
            cls_txpt_set_dic[cls] = {'txpt_set': set(), 'count': 0}
    
    for record in txpt_cls_content:
        m = re.search(r'\s*transcript_id\s+\"(\S+)\";\s*class_code\s+\"(\S+)\"', record.strip())
        if m:
            cls = m.group(2)
            txpt = m.group(1)
            cls_txpt_set_dic[cls]['txpt_set'].add(txpt)
    
    fw = open(out, 'wb')
    for cls in cls_txpt_set_dic.keys():
        cls_txpt_set_dic[cls]['count'] = len(cls_txpt_set_dic[cls]['txpt_set'])
        newline = '{}\t{}\t{}\n'.format(cls, ','.join(cls_txpt_set_dic[cls]['txpt_set']),
                                        str(cls_txpt_set_dic[cls]['count']))
        fw.write(newline)
    fw.close()
