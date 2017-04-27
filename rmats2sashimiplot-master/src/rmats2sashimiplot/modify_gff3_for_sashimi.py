# -*- coding: utf-8 -*-
# __author__ = fiona
# time: 2017/4/17 16:23

import re, os, Bio, argparse, sys, fileinput, urllib2,regex

parser = argparse.ArgumentParser(description="")
parser.add_argument("-gff3", "--gff3", help="file_path", required=True)
parser.add_argument("-out", "--out", help="out_file_path", required=True)

args = parser.parse_args()
in_gff3 = args.gff3
out_gff3 = args.out
out = open(out_gff3, 'wb')

for line in open(in_gff3):
    m = regex.match(r'(^[^#]\S*?\t(\S+?\t){7,7})((\S+?=(\S+?);)+(\S+?=(\S+?)))$', line.strip())
    if m:
        short_str = m.captures(1)[0]
        seq_type = m.captures(2)[1].strip()
        if seq_type == "mRNA" or seq_type == "exon":
            part_a_attr = m.captures(3)
            part_b_attr = m.captures(5)
            id_str = ''
            parent_str = ''
            enst_str = ''
            enst_id = ''
            all_lst = m.captures(3)[0].split(';')
            if all_lst:
                for attr in all_lst:
                    id_m = re.match(r'ID=\S+?;?', attr)
                    if id_m:
                        id_str = attr
                        continue
                    parent_m = re.match(r'[pP]arent=\S+?;?', attr)
                    if parent_m:
                        parent_str = attr
                        continue
                if seq_type == "mRNA":
                    enst_id = re.split(r'=', id_str.strip(';'))[1]
                else:
                    enst_id = re.split(r'=', parent_str.strip(';'))[1]
                enst_str = 'Name={}'.format(enst_id)
                newline = '{}{};{};{}\n'.format(m.captures(1)[0], id_str, parent_str, enst_str)
                out.write(newline)
            else:
                print m.captures(3), m.captures(5)
        else:
            out.write(line)

out.close()
