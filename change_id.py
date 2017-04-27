# -*- coding: utf-8 -*-
# __author__ = fiona
# time: 2017/3/25 18:34

from __future__ import print_function
import re, os, Bio, argparse, sys, fileinput, urllib2, pickle

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-input_gtf", "--gtf", help="file_path", required=True)
    parser.add_argument("-out", "--out", help="out_file_path", required=True)
    # parser.add_argument("-n", "--number", help="number th file", required=True)
    parser.add_argument("-ref_dic", "--ref_dic_pk", help="ref_gtf_pickle", required=True)
    # parser.add_argument("-method", "--method", help="gtf source tool: cufflinks|stringtie", required=True)
    parser.add_argument("-combined_gtf_dic_pk", "--combined_gtf_dic_pk", help="combined_gtf_dic_pk_file_path",
                        default=None)
    args = vars(parser.parse_args())
    gtf = args['gtf']
    out = args['out']
    ref_dic_pk = args['ref_dic_pk']
    
    combined_gtf_pk = args['combined_gtf_dic_pk']
    
    with open(ref_dic_pk, 'rb') as handle:
        txpt_gname_gid_dic = pickle.load(handle)
    # print 'ref gene name 和gene id 对应dic 已经装载完毕，开始修饰第{}个chunk文件：{},结果文件为{}'.format(number,gtf,out)
    
    # if re.match(r'^\s*stringtie\s*$', method):
    with open(combined_gtf_pk, 'rb') as handle:
        combined_gtf_content_dic = pickle.load(handle)  # 包含 new_gene_id: ref_gene_name 的字典
    
    print('将两个dic 的pickle文件内容装进内存结束')
    fw = open(out, 'w')
    # if re.match(r'^\s*cufflinks\s*$', method):
    for line in open(gtf):
        internal_gid = ""
        internal_tid = ""
        cls = ""
        ref_gid = ""
        tid = ""
        ref_gname = ''
        ref_txpt_id = ''
        newline = ''
        # internal_tid_m = re.search(r'transcript_id\s+\"(\S+)\";', line)
        internal_tid = re.search(r'transcript_id\s+\"(\S+)\";', line).group(1)
        # if internal_tid_m:
        #     internal_tid = internal_tid_m.group(1)
        # else:
        #     raise Exception('no valid internal_tid in merged.gtf line: {}'.format(line.strip()))
        # if internal_tid in combined_gtf_content_dic.keys():
        ref_txpt_id = combined_gtf_content_dic[internal_tid]['ref_txpt_id']
        combined_cls = combined_gtf_content_dic[internal_tid]['cls']
        
        newline = re.sub(r'\s*class_code(.*?);', '', line.strip())  # 去除已有的class_code
        newline = re.sub(r'^(.*)$', '\g<1>; class_code \"%s\"' % (combined_cls),
                         newline.strip(';'))  # 在末尾加上annotate.gtf的class_code
        
        if ref_txpt_id.strip():
            # if ref_txpt_id.strip() in txpt_gname_gid_dic.keys():
            ref_gid = txpt_gname_gid_dic[ref_txpt_id]['gene_id']
            ref_gname = txpt_gname_gid_dic[ref_txpt_id]['gname']
            
            if re.match(r'^[^u]$', combined_cls.strip()):
                newline = re.sub(r'(\s*)(gene_id)(\s+)\"(\S+)\";',
                                 '\g<1>\g<2>\g<3>\"{}\";'.format(ref_gid), newline)
                newline = re.sub(r'(\s*)(transcript_id)(\s+)\"(\S+?)\";(\s+?)(\S+.*)',
                                 '\g<1>\g<2>\g<3>\"\g<4>\";\g<5>cmp_ref_gname \"{}\";\g<6>'.format(ref_gname), newline)
            
            # if re.match(r'^[^=u]$', combined_cls.strip()):
            #     newline = re.sub(r'(\s*)(transcript_id)(\s+)\"(\S+)\";',
            #                      '\g<1>\g<2>\g<3>\"{}\";'.format(internal_tid),
            #                      newline)
            
            
            if re.match(r'^=$', combined_cls.strip()):
                newline = re.sub(r'(\s*)(transcript_id)(\s+)\"(\S+)\";',
                                 '\g<1>\g<2>\g<3>\"{}\";'.format(ref_txpt_id),
                                 newline)
            fw.write(newline.strip(';\s') + ';\n')
            # print('{}文件写入行：{}'.format(gtf,newline))
            continue
        else:
            fw.write(newline.strip(';\s') + ';\n')
            # print('the line {} has no valid ref txpt id'.format(line.strip()))
            continue
            # else:
            #     raise Exception(
            #         'combined gtf file {} has no record for internal txpt id:{}'.format(combined_gtf_pk, internal_tid))
            #
    fw.close()
    
    
    
    
    #         if re.search(r'class_code\s+\"=\";', line) and oId:
    #             if oId:
    #                 tid = oId
    #             else:
    #                 tid = internal_tid
    #             if gname and (gname in txpt_gname_gid_dic.keys()):
    #                 gid = txpt_gname_gid_dic[gname]
    #             else:
    #                 gid = internal_gid
    #                 print('{}在ref.gtf中没有记录'.format(gname))
    #             newline = re.sub(r'(\s*)(gene_id)(\s+)\"(\S+)\";', '\g<1>\g<2>\g<3>\"{},{}\"'.format(gid, gname), line)
    #             newline = re.sub(r'(\s*)(transcript_id)(\s+)\"(\S+)\";', '\g<1>\g<2>\g<3>\"{},{}\"'.format(tid, gname),
    #                              newline)
    #             fw.write(newline)
    #             continue
    #         if not re.search(r'class_code\s+\"[=u]\";', line):
    #             tid = internal_tid
    #             if gname and (gname in txpt_gname_gid_dic.keys()):
    #                 gid = txpt_gname_gid_dic[gname]
    #             else:
    #                 gid = internal_gid
    #             newline = re.sub(r'(\s*)(gene_id)(\s+)\"(\S+)\";', '\g<1>\g<2>\g<3>\"{},{}\"'.format(gid, gname), line)
    #             newline = re.sub(r'(\s*)(transcript_id)(\s+)\"(\S+)\";', '\g<1>\g<2>\g<3>\"{},{}\"'.format(tid, gname),
    #                              newline)
    #             fw.write(newline)
    #             continue
    #         if re.search(r'class_code\s+\"u\";', line):
    #             fw.write(line)
    #             continue
    # elif re.match(r'^\s*stringtie\s*$', method):
    #     for line in open(gtf):
    #         stringtie_gid = ""
    #         stringtie_tid = ""
    #         gid = ""
    #         tid = ""
    #         stringtie_gid_m = re.search(r'gene_id\s+\"(\S+)\";', line)
    #         stringtie_tid_m = re.search(r'transcript_id\s+\"(\S+)\";', line)
    #         if stringtie_tid_m:
    #             stringtie_tid = stringtie_tid_m.group(1)
    #         else:
    #             raise Exception('invalid stringtie merged.gtf line: {}, no stringtie txpt id value'.format(line))
    #         if stringtie_gid_m:
    #             stringtie_gid = stringtie_gid_m.group(1)
    #         else:
    #             raise Exception('invalid stringtie merged.gtf line: {},no stringtie gene id value'.format(line))
    #
    #
    #         if re.search(r'class_code\s+\"[^u]\";', line):
    #             if stringtie_tid in tmap_content_dic.keys():
    #                 ref_txpt_id = tmap_content_dic[stringtie_tid]
    #             else:
    #                 raise Exception(
    #                     'the stringtie internal gene id {} of txpts(class_code is not = or u) in merged.gtf file {} has no referenced gene name in its tmap file'.format(
    #                         stringtie_gid, gtf))
    #
    #             tid = stringtie_tid
    #             if ref_txpt_id in txpt_gname_gid_dic.keys():
    #                 ref_gid = txpt_gname_gid_dic[ref_txpt_id]['gene_id']
    #                 ref_gname = txpt_gname_gid_dic[ref_txpt_id]['gname']
    #             else:
    #                 ref_gid = stringtie_gid
    #                 ref_gname = ''
    #                 # # print('no ref gene id for gene name {} of stringtie gene id {} in merged.gtf'.format(gname,
    #                 #                                                                                      stringtie_gid))
    #             newline = re.sub(r'(\s*)(gene_id)(\s+)\"(\S+)\";', '\g<1>\g<2>\g<3>\"{},{}\"'.format(ref_gid, ref_gname), line)
    #             if  re.search(r'class_code\s+\"=\";', line):
    #                 newline = re.sub(r'(\s*)(transcript_id)(\s+)\"(\S+)\";',
    #                                  '\g<1>\g<2>\g<3>\"{},{}\"'.format(ref_txpt_id, ref_gname),
    #                                  newline)
    #
    #             else:
    #                 newline = re.sub(r'(\s*)(transcript_id)(\s+)\"(\S+)\";', '\g<1>\g<2>\g<3>\"{},{}\"'.format(tid, ref_gname),
    #                              newline)
    #             fw.write(newline)
    #             continue
    #         if re.search(r'class_code\s+\"u\";', line):
    #             fw.write(line)
    #             continue
    #
    # else:
    #     raise Exception('您输入的数据来源错误，只可以是cufflinks or stringtie')
