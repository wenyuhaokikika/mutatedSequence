#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   createDB.py
@Time    :   2023/05/05 12:46:31
@Author  :   wenyuhao 
@Version :   1.0
@Contact :   925201392@qq.com
@Desc    :   None
'''

# here put the import lib
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
parser.add_argument('-d','--DB',type=str,help='database path for cds and transcripts')
parser.add_argument('-r','--resource',type=str,help='ensembl database path for resource')
parser.add_argument('-r','--resource',type=str,help='input VCF file for a sample')
parser.add_argument('-e','--expand',default=5000,type=int,help='extend 5000bp for last cds region,make the 5kbp region as the last "cds" region for a transcript')


__release__ == ''

def findResourceFiles(path,key,standard):
  if path
  
#参数
refDir = './'
EXPAND = 5000 #最后一个终止子向后延伸的长度
dnaFa = 'Homo_sapiens.GRCh38.dna.toplevel.fa.gz'
gtf = 'Homo_sapiens.GRCh38.109.chr_patch_hapl_scaff.gtf.gz'
cdsFa = 'Homo_sapiens.GRCh38.cds.all.fa.gz'
enspFa = 'Homo_sapiens.GRCh38.pep.all.fa.gz'
vcfFile = '/data/wenyuhao/55/ensembl/ensembl109/db/all.vcf.gz'



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

# 生成cds的注释文件
logger.info('生成cds的注释文件')
df = pd.DataFrame([{'id':p.id,'description':p.description,'seq':str(p.seq)} for p in parseFile(cdsFa)])
df['description'] = df['description'].apply(lambda x:x.split(' ',2)[-1])
df['V'],df['chr'],df['start'],df['end'],df['strand'] = zip(*df['description'].apply(lambda x:x.replace('chromosome:','').split(' ')[0].replace('scaffold:','').split(':')))

## 对cds文件进行校验
logger.info('对cds文件进行校验')
length = df['end'].astype(int)-df['start'].astype(int)+1
# Counter(length==df['seq'].apply(len))
enspDf = pd.DataFrame([{'id':p.id,'name':p.name,'description':p.description,'seq':str(p.seq)} for p in parseFile(enspFa)])
enspDf['description'] = enspDf['description'].apply(lambda x:x.split(' ',2)[-1])
enspDf['V'],enspDf['chr'],enspDf['start'],enspDf['end'],enspDf['strand'] = zip(*enspDf['description'].apply(lambda x:x.replace('chromosome:','').split(' ')[0].replace('scaffold:','').split(':')))
enspDf['ENST'] = enspDf['description'].apply(lambda x:x.split(' ')[2].replace('transcript:',''))
enspDf = enspDf.drop(['name'],axis=1)

tmpD = enspDf.set_index('ENST')['seq'].to_dict()
df['ENSP'] = df['id'].apply(lambda x:tmpD[x])

dic1 = enspDf.assign(length=lambda x:x['end'].astype(int)-x['start'].astype(int)+1).set_index('ENST')['length'].to_dict()
dic2 = df.assign(length=lambda x:x['end'].astype(int)-x['start'].astype(int)+1).set_index('id')['length'].to_dict()
assert (len(set(dic1)-set(dic2)) == len(set(dic2)-set(dic1))) and (len(set(dic2)-set(dic1))==0),'cds序列文件和ensp的序列文件有差异，请检查原始下载文件'
panL = pd.DataFrame([{'enst':k,'ENSP':v,'cds':dic2[k]} for k,v in dic1.items()])
assert panL.query('ENSP!=cds').shape[0] == 0,'ensembl下载的ENSP和cds序列文件的start，end对应不上，请检查原始下载文件'

# 给cds注释基因组的信息
logger.info('给cds注释基因组的信息')
def getSeqOfCDS(enst):
	cs = cds[cds['transcript_id'] == enst][['seqname','start','end','strand']].apply(tuple,axis=1).to_list()
	return [h[i[0]][i[1]-1:i[2]] if i[3]=='+' else str(Seq(h[i[0]][i[1]-1:i[2]]).reverse_complement()) for i in cs]

def concatCds(x):#注释cds的序列
	x = sorted(x,key=lambda x:int(x[3]))
	s = []
	for i in x:
		if i[-1] == '-':
			s.append(str(Seq(h[i[0]][i[1]-1:i[2]]).reverse_complement()))
		else:
			s.append(h[i[0]][i[1]-1:i[2]])
	return '|'.join(s)

def concatCdsRange(x):#注释cds的range
	x = sorted(x,key=lambda x:int(x[3]))
	s = '|'
	for i in x:
		s += str(i[0])+":"+str(i[1])+"-"+str(i[2])+"|"
	return s
	
def pan(x):#判断cds的染色体是否再参考基因组的染色体中
	s = set([i[0] for i in x])
	return all([i in h for i in s])

#提取基因组上的cds区域
cds1 = cds[['seqname','start','end','strand','frame','gene_id','gene_name','gene_biotype','transcript_id','transcript_biotype','protein_id','exon_number']]
cdsInfo = {k:df.apply(lambda x:tuple([x['seqname'],x['start'],x['end'],x['exon_number'],x['strand']]),axis=1).to_list() for k,df in tqdm(cds1.groupby('transcript_id'))}
cdsSeq = {k:concatCds(v) for k,v in tqdm(cdsInfo.items()) if pan(v)}
cdsSeqRange = {k:concatCdsRange(v) for k,v in tqdm(cdsInfo.items()) if pan(v)}

df['cdsSeq'] = df['id'].apply(lambda x:cdsSeq[x.split('.')[0]])
df['cdsSeqRange'] = df['id'].apply(lambda x:cdsSeqRange[x.split('.')[0]])
df['cdsTranslated'] = df['seq'].apply(lambda x:str(Seq(x).translate()))
pan = (df['cdsTranslated'].str.strip('*')==df['ENSP'])|(df['cdsTranslated'] == df['ENSP'])
logger.info(f'{df[~pan].shape[0]} cds_ok is False')

# 对每个转录本的最后一个cds区域的end注释
def parsePos(x):#解析类似14:22438547-22438554这样的字符串，解析为chr,start,end
	c,p = x.split(':')
	return (c,*p.split('-'))
# 如果是正链，选择最后一个cds的end向后延长5000bp；如果是正链，选择最后一个cds的start向前延长5000bp；
df['last_cds_end'] = df.apply(lambda x:parsePos(x['cdsSeqRange'].split('|')[-2])[-1] if x['strand'] == '1' else parsePos(x['cdsSeqRange'].split('|')[-2])[1],axis=1)

# 注释mode信息
def giveMode(x):
	a,b = x['cdsSeq'].replace('|',''),x['seq']
	if a in b:
		return b.replace(a,'X')
	elif(len(a)==len(b)):
		return 'X'
	elif(len(a)==(len(b)-3)):
		return 'X'+b[-3:]
	else:
		return np.nan
df['mode'] = df.apply(giveMode,axis=1)
print(pd.DataFrame([{'mode':k,'num':v} for k,v in Counter(df['mode']).items()]).T)

#检查vcf中文件和基因组序列是否匹配上
logger.info('检查vcf中文件和基因组序列是否匹配上')
import vcf
vcfRead = vcf.Reader(filename = vcfFile)
muts = [(r.CHROM,r.POS,r.REF) for r in tqdm(vcfRead)]
#[t  for t in muts if h[t[0]][t[1]-1] != t[2][0]]
assert all([h[t[0]][t[1]-1:t[1]+len(t[2])-1] == t[2] for t in muts]),'some position not match'

logger.info('写出为bed文件和fa文件')
#写出为bed文件和fa文件
if not os.path.exists(os.path.join(refDir,'db')):
	os.mkdir(os.path.join(refDir,'db'))
def parsePos(x):
	c,p = x.split(':')
	return (c,*p.split('-'))
f = open(os.path.join(refDir,'db','cds.bed'),'w+')
for i,r in tqdm(df.iterrows()):
	csr = r['cdsSeqRange'].split('|')[1:-1]
	css = r['cdsSeq'].split('|')
	assert len(csr) == len(css),f'{r["id"]}'
	stopWithCodon = r['seq'][-3:] in ['TAG','TGA','TAA']
	addHead = r['mode'].split('X')[0] is isinstance(r['mode'],str) else ';'
	for j,v in enumerate(csr):
		c,s,e = parsePos(v)
		seq = css[j] if r['strand']=='1' else str(Seq(css[j]).reverse_complement())
		f.write(f'{c}\t{s}\t{e}\t{j+1}\t{r["id"]}\t{r["strand"]}\t{stopWithCodon}\t{addHead}\t{seq}\n')
	end = int(r["last_cds_end"])
	s,e = (end,end+EXPAND) if r['strand'] == '1' else (end-EXPAND,end)
	seq = h[r["chr"]][s-1:e]# if r['strand']=='1' else str(Seq(h[r["chr"]][s:e]).reverse_complement())
	f.write(f'{c}\t{s}\t{e}\t{len(csr)+1}\t{r["id"]}\t{r["strand"]}\t{stopWithCodon}\t{addHead}\t{seq}\n')
	f.flush()
f.close()
