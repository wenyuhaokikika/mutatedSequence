#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   mutatedProteins.py
@Time    :   2023/05/07 16:59:09
@Author  :   wenyuhao 
@Version :   1.0
@Contact :   925201392@qq.com
@Desc    :   None
'''

# here put the import lib

import argparse
import pandas as pd
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import Counter
from tqdm import tqdm
from loguru import logger
# 对每个转录本的最后一个cds区域的end注释
def parsePos(x):#解析类似14:22438547-22438554这样的字符串，解析为chr,start,end
	c,p = x.split(':')
	return (c,*p.split('-'))

parser = argparse.ArgumentParser(description='join all cds and translate to protein')
parser.add_argument('-b','--bed',type=str,help='bed file file')
parser.add_argument('-o','--out',type=str,help='output dir')
parser.add_argument('--use_old_bed', action='store_true', help='use old bed to annot')

args = parser.parse_args()
outPath = args.out
#mutatedBED = './db/mutated.bed'
mutatedBED = args.bed

def parse(x):
	#print(x)
	x = x.replace('addHead:','addHead:|') if x.endswith('addHead:') else x.strip(';')# .replace('addHead:\t','addHead:|\t')
	p,i = x.split(' ')
	c,s,e = parsePos(p)
	trans,rank,strand,stopWithCodon,addHead = [i.split(':')[1] for i in i.split(';')]
	return c,s,e,trans,rank,strand,stopWithCodon,addHead
if args.use_old_bed:
    mutatedSeq = pd.read_table(mutatedBED)
else:
    mutatedSeq = pd.read_table(mutatedBED,header=None)
    mutatedSeq.columns = ['info','seq']+list(mutatedSeq.columns)[2:]
    mutatedSeq['chr'],mutatedSeq['start'],mutatedSeq['end'],mutatedSeq['transcript'],mutatedSeq['rank'],mutatedSeq['strand'],mutatedSeq['stopWithCodon'],mutatedSeq['addHead'] = zip(*mutatedSeq['info'].apply(parse))
    mutatedSeq['stopWithCodon'] = mutatedSeq['stopWithCodon'].apply(eval)
    for c in ['start','end','rank','strand']:
        mutatedSeq[c] = mutatedSeq[c].astype(int)
    mutatedSeq = mutatedSeq[['chr', 'start', 'end', 'transcript', 'rank', 'strand','stopWithCodon','seq','addHead']]
    logger.info('reformat mutated bed files')
    mutatedSeq.to_csv(mutatedBED,sep='\t',index=False)

def genateFa(h):
	h = h.sort_values('rank')
	# if  h['stopWithCodon'].iloc[0] == False and h.shape[0]==1:
	# 	print(h)
	if h['stopWithCodon'].iloc[0] == False and h.shape[0]>1:#如果转录本不是因为终止子停止的，这种情况不要在转录结束位点加尾巴
		h = h.iloc[:-1,:]
	else: #如果是因为终止子停止的，这种情况要在最后一个cds（伪cds）区域把第一个碱基cut掉，这一步为了和bcftools consensus保持一致。
		h.loc[h.index[-1], 'seq'] = h.loc[h.index[-1], 'seq'][1:]
		# h['seq'].iloc[-1] = h['seq'].iloc[-1][1:]
	ss = h.apply(lambda x:x['seq'] if x['strand']==1 else str(Seq(x['seq']).reverse_complement()),axis=1).to_list()  #这里在生成时已经算反向互补，这里还要算吗？
	#ss = h['seq'].to_list()
	ss = ''.join(ss)
	addHead = h['addHead'].iloc[0] if not pd.isna(h['addHead'].iloc[0])  else '|'
	addHead = addHead if h['strand'].iloc[0]==1 else str(Seq(addHead).reverse_complement()) if addHead!='|' else '|'
	ss = addHead+ss
	return ss
#print([trans for trans,df in tqdm(mutatedSeq.groupby('transcript')) if df.shape[0]==0])
sequences = [SeqRecord(Seq(genateFa(df)),id=trans,description=f'strand:{df.iloc[0,:]["strand"]},start:{df["start"].min()},end:{df["start"].max()}') for trans,df in tqdm(mutatedSeq.groupby('transcript'))]
logger.info('output mutated cds with expand fatsa seq')
with open(os.path.join(outPath,'mutatedCDS.fasta'), "w") as output_handle:
    SeqIO.write(sequences, output_handle, "fasta")
def recSeq(x):
	if '|' in x or ';' in x:
		tmp = str(x.seq).replace('|','').replace(';','')
		return '|'+Seq(tmp).translate().split('*',1)[0]
	else:
		return x.seq.translate().split('*',1)[0]
logger.info('output mutated proteins cut with first "*"')
sequences = [SeqRecord(recSeq(i),id=i.id,description = i.description) for i in tqdm(sequences)]
with open(os.path.join(outPath,'mutatedProteins.fasta'), "w") as output_handle:
    SeqIO.write(sequences, output_handle, "fasta")
