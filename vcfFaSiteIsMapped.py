'''
Description: 
version: 
Author: wenyuhao
Date: 2023-05-15 10:27:46
LastEditors: wenyuhao
LastEditTime: 2023-05-15 10:27:47
'''
import gzip
import pandas as pd
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import Counter
from tqdm import tqdm
from loguru import logger
import numpy as np
parseFile = lambda x: SeqIO.parse(gzip.open(x, "rt"), "fasta") if 'gz' in x else SeqIO.parse(x,'fasta')
import argparse
parser = argparse.ArgumentParser(description='create database for ',prog='mutatedSequence')
parser.add_argument('-m','--mutations',type=str,help='vcf file of mutations')
args = parser.parse_args()

vcfFile = args.mutations #'/data/wenyuhao/55/ensembl/ensembl109/db/all.vcf.gz'
dnaFa = ''
gtf = ''
#加载基因组序列文件
logger.info('加载基因组序列文件...')
h = {i.id:str(i.seq) for i in tqdm(parseFile(dnaFa))}

#加载基因组注释文件
logger.info('加载基因组注释文件...')
import gtfparse
from tqdm import tqdm
gtf = gtfparse.read_gtf(gtf)
transcript = gtf[gtf['feature']=='transcript']
cds =  gtf[gtf['feature']=='CDS']

#检查vcf中文件和基因组序列是否匹配上
logger.info('检查vcf中文件和基因组序列是否匹配上')
import vcf
vcfRead = vcf.Reader(filename = vcfFile)
muts = [(r.CHROM,r.POS,r.REF) for r in tqdm(vcfRead)]
#[t  for t in muts if h[t[0]][t[1]-1] != t[2][0]]
assert all([h[t[0]][t[1]-1:t[1]+len(t[2])-1] == t[2] for t in muts]),'some position not match'