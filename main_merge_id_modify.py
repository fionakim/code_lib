# -*- coding: utf-8 -*-
# __author__ = fiona
# time: 2017/3/24 14:42


from __future__ import print_function
from __future__ import print_function
from __future__ import print_function
from __future__ import print_function
import argparse
import subprocess
import re
import time
import os, pickle

'''
此脚本主要用作cufflinks的merged.gtf 转换ID使用。
原则如下：
1. class code 为=时，认为此转录本是ref.gtf中的已知转录本，则gene/transcript name/ID 都要替换为ref.gtf中的ID和name。其中gene id 来源于 记录中gene name 在ref.gtf 对应的gene ID，
txpt 的ID为 记录中的oId
2. class code 为u时,不改变记录
3. 其他情况只替换其gene name
使用说明：
因merged.gtf较大，因此将其劈成指定行数的小文件集合，同步进行修饰，增加速度
参数：-s 处理脚本 必须在一起合用 是change_id.py 的路径
  -tmp  name 暂存文件的文件夹名字  绝对路径为merged。gtf 的文件夹路径下的name文件夹
   -rm 是否结束程序后删除 暂存文件夹
   -lines merged gtf 被split后，每个小文件里有多少行内容
   
使用命令举例：
python main_merge_id_modify.py  -tmp tem_folder_0328_no_raise  -merge merged.gtf  -out_merge new_jin_cufflinks_0328_no_raise.gtf  -ref_gtf ref.gtf -s change_id.py   -combined  cuffcmp.annotated.gtf -batch_no 5 > new_jin_cufflinks_0328_no_raise.log
   
'''


# file_path = "F:\\code_lib\\merged.gtf"
# new_file_path = "F:\\code_lib\\new_merged.gtf"
# ref_gtf = "F:\\code_lib\\ref.gtf"

# =================函数区==================
def get_gname_gid_dic_from_ref_gtf_cuff(ref):
    d = {}
    for line in open(ref):
        gene_id_m = re.search(r'gene_id\s+\"(\S+)\"', line.strip())
        gname_m = re.search(r'gene_name\s+\"(\S+)\"', line.strip())
        # if gene_id_m and gname_m:
        gname = gname_m.group(1)
        gene_id = gene_id_m.group(1)
        # if gname in d.keys() and d[gname] != gene_id:
        #     raise Exception('ref.gtf 中的gene_name {} 有一个以上的gene_id: {} and {}'.format(gname, d[gname], gene_id))
        d[gname] = gene_id
        # else:
        #     raise Exception('ref gtf文件不合法的第九列: {},没有完整的gene id 和gene_name 记录'.format(line.strip()))
    return d


def get_gname_gid_dic_from_ref_gtf(ref_gtf):
    d = {}
    for line in open(ref_gtf):
        gene_id_m = re.search(r'gene_id\s+\"(\S+)\"', line.strip())
        gname_m = re.search(r'gene_name\s+\"(\S+)\"', line.strip())
        txpt_id_m = re.search(r'transcript_id\s+\"(\S+)\"', line.strip())
        if gene_id_m and gname_m and txpt_id_m:
            gname = gname_m.group(1)
            gene_id = gene_id_m.group(1)
            txpt_id = txpt_id_m.group(1)
            # if txpt_id in d.keys() and d[txpt_id] != {'gene_id': gene_id, 'gname': gname}:
            #     raise Exception('ref.gtf 中的transcript id {} 有一个以上的gene_info: {} and {}'.format(txpt_id, d[txpt_id],
            #                                                                                    {'gene_id': gene_id,
            #                                                                                    {'gene_id': gene_id,
            #                                                                                     'gname': gname}))
            d[txpt_id] = {'gene_id': gene_id, 'gname': gname}
        else:
            raise Exception('ref gtf文件不合法的第九列: {},没有完整的gene id 和gene_name 记录'.format(line.strip()))
    return d


def split_merged_gtf(f, line_num, tmp_dir_name, name_prefix):
    tmp_dir = os.path.join(os.path.dirname(f), tmp_dir_name)
    abs_prefix = os.path.join(tmp_dir, name_prefix)
    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir)
    else:
        subprocess.call('rm -r {}'.format(tmp_dir),shell=True)
        os.mkdir(tmp_dir)
    subprocess.call('split -l {} {} {}'.format(line_num, f, abs_prefix), shell=True)
    chunk_files = [os.path.join(tmp_dir, chunk) for chunk in os.listdir(tmp_dir)]
    return chunk_files, tmp_dir


def write_dic_to_pk(pk_file, dic):
    with open(pk_file, 'wb') as handle:
        pickle.dump(dic, handle, protocol=pickle.HIGHEST_PROTOCOL)


def get_dic_from_combined_gtf(combined_gtf):
    d = {}
    for line in open(combined_gtf):
        if re.search(r'\t+transcript\t+', line):
            cls = re.search(r'class_code\s+\"(\S+)\";', line).group(1)
            internal_txpt_id = re.search(r'transcript_id\s+\"(\S+)\";', line).group(1)
            ref_txpt_id_m = re.search(r'cmp_ref\s+\"(\S+)\";', line)
            ref_txpt_id = ''
            if ref_txpt_id_m:
                ref_txpt_id = ref_txpt_id_m.group(1)
            d[internal_txpt_id] = {'ref_txpt_id': ref_txpt_id, 'cls': cls}
            
            #
            # # if (re.match(r'^\-$', ref_gene_name) or re.match(r'^\-$', ref_txpt_id)) and re.match(r'^[^u]$', cls_code):
            # #     raise Exception('{} 文件中出现非u转录本无参考基因name和参考转录本id的情况')
            # if not re.match(r'^-$', ref_txpt_id) and re.match(r'^-$', internal_txpt_id):
            #     d[internal_txpt_id] = ref_txpt_id
    return d


# =============================主函数=====================
if __name__ == '__main__':
    time1 = time.time()
    parser = argparse.ArgumentParser(description="file")
    parser.add_argument("-merge", "--merged_gtf", help="input merged gtf ", required=True)
    parser.add_argument("-out_merge", "--new_merged_gtf", help="out new merged gtf", required=True)
    parser.add_argument("-ref_gtf", "--ref_gtf", help="ref_gtf", required=True)
    parser.add_argument("-s", "--trans_script", help="trans content script path", required=True)
    parser.add_argument("-tmp", "--tmp_dir", help="temp dir name", default='temp_folder')
    parser.add_argument("-rm", "--rm_tmp", help="rm temp dir or not,default is true(1),false is 0  ", default=0)
    parser.add_argument("-lines", "--split_line_number",
                        help="split big file to  how many lines per file default is 40000 ", default=40000)
    # parser.add_argument("-method", "--method", help="gtf source tool: cufflinks|stringtie", required=True)
    # parser.add_argument("-combined_gtf ", "--combined_gtf _file", help="combined_gtf _file, src 为combined_gtf 时需要设定有意义的值", default=None)
    parser.add_argument("-combined", "--combined_gtf", help="combined.gtf file", required=True)
    parser.add_argument("-batch_no", "--batch_no", help="batch_no", required=True)
    
    args = vars(parser.parse_args())
    file_path = args["merged_gtf"]
    new_file_path = args["new_merged_gtf"]
    ref_gtf = args["ref_gtf"]
    trans_script = args['trans_script']
    tmp_dir = args['tmp_dir']
    rm = args['rm_tmp']
    line_number = int(args['split_line_number'])
    # method = args['method']
    combined_gtf = args['combined_gtf']
    batch_no = int(args['batch_no'])
    
    print('开始劈开merged.gtf')
    chunk_files_lst, tmp_folder = split_merged_gtf(file_path, line_number, tmp_dir, 'chunk_file_')
    tmp_ref_gtf = os.path.join(tmp_folder, 'temp_ref.gtf')
    ref_tmp_cmd = """ awk -F '\t' 'NF>=9{print $9}' %s  | uniq > %s """ % (ref_gtf, tmp_ref_gtf)
    print('结束劈开merged.gtf')
    print('开始awk ref.gtf')
    ref_tmp_content = subprocess.check_output(ref_tmp_cmd, shell=True).split('\n')
    print('ref gtf 暂时文件读取完毕')
    txpt_gname_gid_dic = get_gname_gid_dic_from_ref_gtf(tmp_ref_gtf)
    print('ref gtf 信息装载完毕')
    ref_dic_pickle = os.path.join(tmp_folder, 'ref_gname_gid_dic.pk')
    write_dic_to_pk(ref_dic_pickle, txpt_gname_gid_dic)
    print('ref 的info dic字典对象已保存，保存地址为{}'.format(ref_dic_pickle))
    
    # if re.match(r'^\s*stringtie\s*$', method) and (not combined_gtf _file):
    #     raise Exception('当merged.gtf为stringtie结果时，combined_gtf  file文件路径必须设定')
    # if re.match(r'^\s*stringtie\s*$', method):
    #     combined_gtf _content_dic = get_dic_from_cmp_gtf (combined_gtf _file)
    #     combined_gtf _content_dic_pickle = os.path.join(tmp_folder, 'combined_gtf _content_dic.pk')
    #     write_dic_to_pk(combined_gtf _content_dic_pickle, combined_gtf_content_dic)
    print('开始装载combined_gtf信息')
    combined_gtf_content_dic = get_dic_from_combined_gtf(combined_gtf)
    combined_gtf_content_dic_pickle = os.path.join(tmp_folder, 'combined_gtf_content_dic.pk')
    write_dic_to_pk(combined_gtf_content_dic_pickle, combined_gtf_content_dic)
    print('combined_gtf 信息写入文件完毕')
    print('combined_gtf 的info dic字典对象已保存，保存地址为{}'.format(ref_dic_pickle))
    
    son_pro_lst = []
    out_file_lst = []
    limit = (len(chunk_files_lst) / batch_no) * batch_no
    batch_lst = chunk_files_lst[0:limit]
    single_lst = chunk_files_lst[limit:len(chunk_files_lst)]
    for chunk in chunk_files_lst:
        number = chunk_files_lst.index(chunk) + 1
        out_file = os.path.join(tmp_folder, 'modified_' + os.path.basename(chunk))
        out_file_lst.append(out_file)
        mod_gtf_cmd = 'python {}  -input_gtf  {} -out {} -ref_dic {} -combined_gtf_dic_pk {} '.format(trans_script,
                                                                                                      chunk,
                                                                                                      out_file,
                                                                                                      ref_dic_pickle,
                                                                                                      combined_gtf_content_dic_pickle
                                                                                                      )
        print('开始修饰第{}个文件'.format(number))
        child_pro = subprocess.Popen(mod_gtf_cmd, shell=True)
        son_pro_lst.append(child_pro)
        
        if len(son_pro_lst) >= 5:
            for pro in son_pro_lst:
                pro.communicate()
                son_pro_lst = []
            continue
        else:
            continue
    
    single_pro_lst = []
    for chunk_file in single_lst:
        out_file = os.path.join(tmp_folder, 'modified_' + os.path.basename(chunk_file))
        out_file_lst.append(out_file)
        mod_gtf_cmd = 'python {}  -input_gtf  {} -out {} -ref_dic {} -combined_gtf_dic_pk {} '.format(trans_script,
                                                                                                      chunk_file,
                                                                                                      out_file,
                                                                                                      ref_dic_pickle,
                                                                                                      combined_gtf_content_dic_pickle
                                                                                                      )
        extra = subprocess.Popen(mod_gtf_cmd, shell=True)
        single_pro_lst.append(extra)
    for extra_pro in single_pro_lst:
        extra_pro.communicate()
    
    out_file_str = "  ".join(out_file_lst)
    subprocess.call('cat {} > {}'.format(out_file_str, new_file_path), shell=True)
    print('已将最后结果汇聚为{}'.format(new_file_path))
    if rm:
        subprocess.call('rm -rf {}'.format(tmp_folder), shell=True)
    time2 = time.time()
    duration = time2 - time1
    m, s = divmod(duration, 60)
    h, m = divmod(m, 60)
    print('整个程序运行的时间为{}h:{}m:{}s'.format(h, m, s))
