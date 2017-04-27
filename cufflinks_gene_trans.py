#!/usr/bin/python
# -*- coding: utf-8 -*-
import argparse
import subprocess
import re

parser = argparse.ArgumentParser(description="file")
parser.add_argument("-file_path", "--file_path", help="file_path", required = True)
parser.add_argument("-new_file_path", "--new_file_path", help="new_file_path", required = True)
parser.add_argument("-ref_gtf", "--ref_gtf", help="ref_gtf", required = True)
parser.add_argument("-method", "--method", help="cufflinks/stringtie", required = True)
args = vars(parser.parse_args())

file_path = args["file_path"]
new_file_path = args["new_file_path"]
ref_gtf = args["ref_gtf"]
method = args["method"]

"""建立转录本和基因一一对应的字典结构"""
a1=""" awk -F "\\t" "{print $9}" %s | awk -F ';' '{print $1$2}' | awk -F '"' '{print $2"\\t"$4}' | uniq -c | awk '{print $2"\\t"$3}' """%(ref_gtf)
print a1
a2=subprocess.check_output(a1, shell=True)
a3=[i.split("\t")for i in a2.split("\n")[:-1]]
trans2genes = {}
for ll in a3:
    if ll[0] not in trans2genes.keys():
        trans2genes[ll[0]]=ll[1]
print '参考基因组gtf转换提取成功！'
    
def cufflinks_gene_trans(file_path, new_file_path, trans2genes):
    with open(file_path,'r+') as f1, open(new_file_path, 'w+') as f2:
        for lines in f1:
            line1=lines.strip().split("\t")
            if 'oId' not in line1[8]:
                f2.write(lines)
            else:
                new_line2 = line1[8].split(";")
                for llm in new_line2:
                    ll=llm.replace(' ','')
                    if re.search('oId', ll):
                        transcript = ll.split("\"")[1]
                        if transcript in  trans2genes.keys():
                            _line = 'gene_id '+'"'+trans2genes[transcript]+'"'+"; "+"transcript_id "+'"'+transcript+'"'+";"
                            _line1 = "\t".join(line1[:8])+"\t"+_line+";".join(new_line2[2:])+"\n"
                            f2.write(_line1)
                            break
                        else:
                            f2.write(lines)
                            break
                    else:
                        pass

def stringtie_gene_trans(file_path, new_file_path):
    with open(file_path,"r+") as f1, open(new_file_path,'w+') as f2:
        for lines in f1:
            line1=lines.strip().split("\t")
            if """class_code "=";""" in line1[8] and re.search('ref_gene_id',line1[8]):
                line2=line1[8].split(";")
                for ll in line2:
                    if re.search('ref_gene_id',ll):
                        gene_id = ll.split("\"")[1]
                        _line = 'gene_id '+'"'+gene_id+'"'+"; "
                        _line1 = "\t".join(line1[:8])+"\t"+_line+";".join(line2[1:])+"\n"
                        f2.write(_line1)
                    else:
                        pass
            else:
                f2.write(lines)

if method == 'cufflinks':
    cufflinks_gene_trans(file_path, new_file_path,trans2genes)
elif method == 'stringtie':
    stringtie_gene_trans(file_path, new_file_path)