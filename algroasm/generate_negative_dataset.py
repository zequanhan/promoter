import random
from Bio import SeqIO
import pandas as pd
import pandas as pd
import numpy as np
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.metrics import accuracy_score, classification_report
import re
import os
import sys
sys.path.append('/home/hanzequan/test_bectiral/rf_model/model1/algroasm')
from feature import com_seq_feature

def inter_group_disruption(seq, k=3):
    """随机交换序列中的子序列位置来生成负样本"""
    subsequences = [seq[i:i+k] for i in range(0, len(seq), k)]
    random.shuffle(subsequences)
    return ''.join(subsequences)
def intra_group_disruption(seq, k=3):
    """随机交换子序列内的碱基位置来生成负样本"""
    def shuffle_subsequence(subseq):
        subseq_list = list(subseq)
        random.shuffle(subseq_list)
        return ''.join(subseq_list)
    
    subsequences = [shuffle_subsequence(seq[i:i+k]) for i in range(0, len(seq), k)]
    return ''.join(subsequences)

def find_noncoding_regions(gbk_record, seq_length):
    """识别非编码区域，并返回一个包含(start, end)元组的列表"""
    noncoding_regions = []
    coding_regions = [feature for feature in gbk_record.features if feature.type == "CDS"]
    sorted_coding_regions = sorted(coding_regions, key=lambda x: x.location.start)

    last_end = 0
    for feature in sorted_coding_regions:
        start = int(feature.location.start)
        if last_end < start:  # 找到非编码区域
            noncoding_regions.append((last_end, start))
        last_end = int(feature.location.end)

    # 处理最后一个编码区域后的序列
    if last_end < len(gbk_record.seq):
        noncoding_regions.append((last_end, len(gbk_record.seq)))
    
    # 筛选出长度至少为seq_length的非编码区域
    suitable_regions = [region for region in noncoding_regions if region[1] - region[0] >= seq_length]
    return suitable_regions

def random_interception_noncoding_extended(accession, seq_length, fasta_dir, gbk_dir):
    gbk_path = f"{gbk_dir}/{accession}.gbk"

    try:
        gbk_record = next(SeqIO.parse(gbk_path, "genbank"))
        genome_seq = str(gbk_record.seq)
        noncoding_regions = find_noncoding_regions(gbk_record, seq_length)
        
        if noncoding_regions:
            selected_region = random.choice(noncoding_regions)
            start_pos = random.randint(selected_region[0], selected_region[1] - seq_length)
        else:
            # 如果没有足够长的非编码区域，随机选择起始位置并截取，可能涉及编码区
            start_pos = random.randint(0, len(genome_seq) - seq_length)
            
        end_pos = start_pos + seq_length
        neg_seq = genome_seq[start_pos:end_pos]

        # 处理序列中的 'N'
        nucleotides = ['A', 'T', 'C', 'G']
        neg_seq = ''.join([random.choice(nucleotides) if nucleotide == 'N' else nucleotide for nucleotide in neg_seq])

        return neg_seq, start_pos, end_pos
    except Exception as e:
        print(f"Error processing {accession}: {e}")
    
    return None, None, None




def generate_wgri_negative_sample_with_start(seq_length, fasta_dir, accessions_list):
    """从整个基因组中随机截取相应长度的序列作为一个WGRI负样本，同时返回起始位置"""
    accession = random.choice(accessions_list)
    phage_fasta_path = f"{fasta_dir}/{accession}.fasta"
    try:
        seq_record = next(SeqIO.parse(phage_fasta_path, "fasta"))
        genome_seq = str(seq_record.seq)
        start_pos = random.randint(0, len(genome_seq) - seq_length)
        neg_seq = genome_seq[start_pos:start_pos + seq_length]
        return neg_seq, start_pos
    except:
        return None, None



def generate_negative_samples_for_missing_accessions(missing_accessions_list, fasta_dir, seq_length):
    negative_samples = []
#     num_samples_list = [13, 13, 12]  # 分配每个accession的样本数量
    #共生成正负数据集的差个序列
    for accession, num_samples_per_accession in zip(missing_accessions_list, num_samples_list):
        for _ in range(num_samples_per_accession):  # 对每个缺失的accession生成指定数量的负样本
            neg_seq, start_pos = generate_wgri_negative_sample_with_start(seq_length, fasta_dir, [accession])
            if neg_seq is not None:
                end_pos = start_pos + seq_length - 1
                negative_samples.append({
                    'accession': accession,
                    'start': start_pos,
                    'end': end_pos,
                    'sequence': neg_seq
                })
    
    return pd.DataFrame(negative_samples)
