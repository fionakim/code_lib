#!/usr/bin/python
# -*- coding: utf-8 -*-
import argparse
import subprocess
import re

parser = argparse.ArgumentParser(description="file")
parser.add_argument("-file_path", "--file_path", help="file_path", required=True)
parser.add_argument("-new_file_path", "--new_file_path", help="new_file_path", required=True)
parser.add_argument("-ref_gtf", "--ref_gtf", help="ref_gtf", required=True)
parser.add_argument("-method", "--method", help="cufflinks/stringtie", required=True)
args = vars(parser.parse_args())

file_path = args["file_path"]
new_file_path = args["new_file_path"]
ref_gtf = args["ref_gtf"]
method = args["method"]


def get_dic_from_list(**kwargs):
	d = {}
	lst = kwargs['lst']
	sub_sep = kwargs['sub_sep']
	sep = kwargs['primary_sep']
	primary_key_col = kwargs['primary_key_index']
	primary_value_index = kwargs['primary_value_index']
	for record in lst:
		record_items = record.split(sep)
		record_eles = [re.split(r''.format(sub_sep), ele)[1].strip("\"") for ele in record_items]
		d[record_eles[primary_key_col - 1]] = record_eles[primary_value_index - 1]
	return d


def get_multi_dic_from_list(**kwargs):
	d = {}
	lst = kwargs['lst']
	sub_sep = kwargs['sub_sep']
	sep = kwargs['primary_sep']
	primary_key_index = kwargs['primary_key_index']
	second_key_index = kwargs['second_key_index']
	for record in lst:
		record_eles = [re.split(r''.format(sub_sep), ele)[1].strip("\"") for ele in record.split(sep)]
		gene_name = record_eles[2]
		sub_d = {}
		oId = ""
		if 'oId' in record:
			oId = record_eles[3]
		sub_d[record_eles[second_key_index - 1]] = {'gene_name': gene_name, 'oId': oId}
		d[record_eles[primary_key_index - 1]] = sub_d
	return d


"""建立转录本和基因一一对应的字典结构"""
merged_tmp_cmd = """ awk -F '\t|;'  '{printf $9";"$10 ; for (i=12; i <= NF-2; i++) printf ";"$i; printf "\n" }' %s  | uniq """ % (
	file_path)
merged_tmp_content = subprocess.check_output(merged_tmp_cmd, shell=True).split('\n')
ref_tmp_cmd = """ awk -F '\t|;' '{for (i=10; i <= NF-1; i++) printf ";"$i; printf "\n" }' %s  | uniq """ % (ref_gtf)
ref_tmp_content = subprocess.check_output(ref_tmp_cmd, shell=True).split('\n')

merged_dic = get_multi_dic_from_list(lst=merged_tmp_cmd, primary_sep=';', sub_sep='\s+', primary_key_col=1,
                                     second_key_index=-1)
gname_gid_dic = get_dic_from_list(lst=ref_tmp_content, primary_sep=';', sub_sep='\s+', primary_key_index=1,primary_value_index=0)

fw = open(new_file_path,'w')
for line in open(file_path):
	if  'class_code "="' in line:
		arr = re.split(r'[;\t]',line.strip())
		old_gid = re.split(r'\s+',arr[8])[1].strip("\"")
		old_tid = re.split(r'\s+',arr[9])[1].strip("\"")
		cls = re.split(r'\s+',arr[-2])[1].strip("\"")
		if merged_dic["="][old_tid]['gene_name'] in gname_gid_dic.keys():
			gid = gname_gid_dic[merged_dic["="][old_tid]['gene_name']]
		else:
			gid = old_gid
		if merged_dic["="][old_tid]['oId']:
			tid = merged_dic["="][old_tid]['oId']
		else:
			tid = old_tid
		arr[8] = 'gene_id "{}"'.format(gid)
		arr[9] = ' transcript_id "{}"'.format(tid)
		newline = '\t'.join(arr[0:8])+'\t'+';'.join(arr[8:])
	else:
		fw.write(line)
fw.close()


def cufflinks_gene_trans(file_path, new_file_path, trans2genes):
	with open(file_path, 'r+') as f1, open(new_file_path, 'w+') as f2:
		for lines in f1:
			line1 = lines.strip().split("\t")
			if 'oId' not in line1[8]:
				f2.write(lines)
			else:
				new_line2 = line1[8].split(";")
				for llm in new_line2:
					ll = llm.replace(' ', '')
					if re.search('oId', ll):
						transcript = ll.split("\"")[1]
						if transcript in trans2genes.keys():
							_line = 'gene_id ' + '"' + trans2genes[
								transcript] + '"' + "; " + "transcript_id " + '"' + transcript + '"' + ";"
							_line1 = "\t".join(line1[:8]) + "\t" + _line + ";".join(new_line2[2:]) + "\n"
							f2.write(_line1)
							break
						else:
							f2.write(lines)
							break
					else:
						pass


def stringtie_gene_trans(file_path, new_file_path):
	with open(file_path, "r+") as f1, open(new_file_path, 'w+') as f2:
		for lines in f1:
			line1 = lines.strip().split("\t")
			if """class_code "=";""" in line1[8] and re.search('ref_gene_id', line1[8]):
				line2 = line1[8].split(";")
				for ll in line2:
					if re.search('ref_gene_id', ll):
						gene_id = ll.split("\"")[1]
						_line = 'gene_id ' + '"' + gene_id + '"' + "; "
						_line1 = "\t".join(line1[:8]) + "\t" + _line + ";".join(line2[1:]) + "\n"
						f2.write(_line1)
					else:
						pass
			else:
				f2.write(lines)


if method == 'cufflinks':
	cufflinks_gene_trans(file_path, new_file_path, trans2genes)
elif method == 'stringtie':
	stringtie_gene_trans(file_path, new_file_path)
