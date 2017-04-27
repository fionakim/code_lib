# -*- coding: utf-8 -*-
# __author__ = fiona
# time: 2017/3/27 16:00

from __future__ import print_function
import re, os, Bio, argparse, sys, fileinput, urllib2

if __name__ == '__main__':
    tmap = '/mnt/ilustre/users/sanger-dev/workspace/20170214/Single_assembly_module_tophat_cufflinks/Assembly/output/Cuffmerge/cuffmerge_gffcompare/gffcmp_cuffmerge_.merged.gtf.tmap'
    refmap = '/mnt/ilustre/users/sanger-dev/workspace/20170214/Single_assembly_module_tophat_cufflinks/Assembly/output/Cuffmerge/cuffmerge_gffcompare/gffcmp_cuffmerge_.merged.gtf.refmap'
    merged = '/mnt/ilustre/users/sanger-dev/workspace/20170214/Single_assembly_module_tophat_cufflinks/Assembly/output/Cuffmerge/cuffmerge_gffcompare/merged.gtf'
    combine = '/mnt/ilustre/users/sanger-dev/workspace/20170214/Single_assembly_module_tophat_cufflinks/Assembly/output/Cuffmerge/cuffmerge_gffcompare/gffcmp_cuffmerge_.annotated.gtf'
    # tmap = 'F:\\temp\\cuffmerge_gffcompare\\gffcmp_cuffmerge_.merged.gtf.tmap'
    # refmap = 'F:\\temp\\cuffmerge_gffcompare\\gffcmp_cuffmerge_.merged.gtf.refmap'
    # merged = 'F:\\temp\\cuffmerge_gffcompare\\merged.gtf'
    # combine = 'F:\\temp\\cuffmerge_gffcompare\\gffcmp_cuffmerge_.annotated.gtf'
    tmap_set = set()
    tmap_txpt_cls = set()
    refmap_set = set()
    refmap_txpt_cls = set()
    combine_set = set()
    combine_txpt_cls = set()
    merge_set = set()
    merge_txpt_cls = set()
    
    for tmap_line in open(tmap):
        tmap_arr = re.split(r'\t', tmap_line.strip())
        tmap_internal_txpt_id = ''
        tmap_ref_txpt_id = ''
        tmap_cls = ''
        tmap_internal_txpt_id = tmap_arr[4]
        tmap_ref_txpt_id = tmap_arr[1]
        tmap_cls = tmap_arr[2]
        
        tmap_set.add((tmap_internal_txpt_id, tmap_ref_txpt_id, tmap_cls))
        tmap_txpt_cls.add((tmap_internal_txpt_id, tmap_ref_txpt_id, tmap_cls))
    #
    # for refmap_line in open(refmap):
    #     refmap_arr = re.split(r'\t', refmap_line.strip())
    #     refmap_ref_txpt_id = refmap_arr[1]
    #     refmap_internal_txpt_id = refmap_arr[3].split('|')[1]
    #     refmap_cls = refmap_arr[2]
    #     refmap_set.add((refmap_internal_txpt_id, refmap_ref_txpt_id, refmap_cls))
    #     refmap_txpt_cls.add((refmap_internal_txpt_id, refmap_cls))
    
    for combine_line in open(combine):
        combine_internal_tid = ''
        combine_ref_txpt_id = ''
        combine_cls = ''
        if re.search('\t+transcript\t+', combine_line):
            combine_internal_tid_m = re.search(r'transcript_id\s+\"(\S+)\";', combine_line.strip())
            # print(internal_txpt_id)
            combine_cls_m = re.search(r'(\s*)(class_code)(\s+)\"(\S+)\";', combine_line.strip())
            combine_ref_txpt_id_m = re.search(r'cmp_ref\s+\"(\S+)\";', combine_line.strip())
            if combine_internal_tid_m and combine_cls_m:
                combine_internal_tid = combine_internal_tid_m.group(1)
                combine_cls = combine_cls_m.group(4)
            if combine_ref_txpt_id_m:
                combine_ref_txpt_id = combine_ref_txpt_id_m.group(1)
            else:
                combine_ref_txpt_id = ''
            # print((internal_txpt_id, ref_txpt_id, cls))
            combine_txpt_cls.add((combine_internal_tid, combine_cls))
            combine_set.add((combine_internal_tid, combine_ref_txpt_id, combine_cls))
    
    for merge_line in open(merged):
        merge_internal_txpt_id_m = re.search(r'transcript_id\s+\"(\S+)\";', merge_line)
        merge_ref_txpt_id_m = re.search(r'nearest_ref\s+\"(\S+)\";', merge_line)
        merge_cls_m = re.search(r'(\s*)(class_code)(\s+)\"(\S+)\";', merge_line.strip())
        merge_internal_txpt_id = ''
        merge_ref_txpt_id = ''
        merge_cls = ''
        if merge_internal_txpt_id_m and merge_cls_m:
            merge_internal_txpt_id = merge_internal_txpt_id_m.group(1)
            merge_cls = merge_cls_m.group(4)
        if merge_ref_txpt_id_m:
            merge_ref_txpt_id = merge_ref_txpt_id_m.group(1)
        else:
            merge_ref_txpt_id = ''
        # print((internal_txpt_id, ref_txpt_id, cls))
        merge_txpt_cls.add((merge_internal_txpt_id, merge_cls))
        merge_set.add((merge_internal_txpt_id, merge_ref_txpt_id, merge_cls))
    
    # print('the count of tmap txpts,ref_txpt_clprindiff_combine_merged)))
    
    print('==============tmap================')
    print(len(tmap_txpt_cls))
    print(len(tmap_set))
    print('=================refmap=============')
    print(len(refmap_set))
    print(len(refmap_txpt_cls))
    print('================merge==============')
    print(len(merge_txpt_cls))
    print(len(merge_set))
    print('============combine==================')
    print(len(combine_txpt_cls))
    print(len(combine_set))
    print('===========combine and merge===================')
    print(len((combine_set & merge_set)))
    print(len((combine_set & merge_set) - combine_set))
