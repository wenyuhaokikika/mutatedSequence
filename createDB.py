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
import re
import traceback
pattren = re.compile(r'[\W+\w+]*?get_variable_name\((\w+)\)')
__get_variable_name__ = []
def get_variable_name(x):#获得变量名为字符串
    global __get_variable_name__
    if not __get_variable_name__:
        __get_variable_name__ = pattren.findall(traceback.extract_stack(limit=2)[0][3])
    return __get_variable_name__.pop(0)

import argparse
parser = argparse.ArgumentParser(description='create database for ',prog='mutatedSequence')
parser.add_argument('-d','--DB',type=str,help='database path for cds and transcripts')
parser.add_argument('-r','--resource',type=str,help='ensembl database path for resource')
parser.add_argument('-e','--expand',default=5000,type=int,help='extend 5000bp for last cds region,make the 5kbp region as the last "cds" region for a transcript')
#parser.add_argument('-m','--mutations',type=str,help='vcf file of mutations')
parser.add_argument('-R','--release',type=str,help='release version')
args = parser.parse_args()

# DIR=/data/wenyuhao/55/ensembl/ensembl109
#progDir=/data/wenyuhao/55/ensembl/ensembl109/db/mutatedSequence
# python $progDir/createDB.py -d $DIR/db1 -r $DIR -R 109 -e 5000 -m /data/wenyuhao/tmp/runWGS/tmp/all.vcf

# def findResourceFiles(path,key,standard):
#   if path
release = args.release  #release in ensembl

#参数
refDir = args.DB #'./'
EXPAND = args.expand #5000 #最后一个终止子向后延伸的长度
resource = args.resource
#vcfFile = args.mutations #'/data/wenyuhao/55/ensembl/ensembl109/db/all.vcf.gz'

dnaFa = os.path.join(resource,'Homo_sapiens.GRCh38.dna.toplevel.fa.gz')
gtf = os.path.join(resource,f'Homo_sapiens.GRCh38.{release}.chr_patch_hapl_scaff.gtf.gz')
cdsFa = os.path.join(resource,'Homo_sapiens.GRCh38.cds.all.fa.gz')
enspFa = os.path.join(resource,'Homo_sapiens.GRCh38.pep.all.fa.gz')
for f in [dnaFa,gtf,cdsFa,enspFa]:
    assert os.path.exists(f),'resource or vcf files not exists'

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
assert panL.query('ENSP!=cds').shape[0] == 0,'ensembl下载的ENSP和cds序列文件的数目对应不上，请检查原始下载文件'

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
# 如果是正链，选择最后一个cds的end向后延长5000bp；如果负链，选择最后一个cds的start向前延长5000bp；
df['last_cds_end'] = df.apply(lambda x:parsePos(x['cdsSeqRange'].split('|')[-2])[-1] if x['strand'] == '1' else parsePos(x['cdsSeqRange'].split('|')[-2])[1],axis=1)

# 注释mode信息，mode中X代表核酸序列从人类的注释文件和基因组序列上切取的序列，一般ensembl提供的ENSP序列会在X的头部加'N'，尾部加入终止密码子，这个mode来说明这两个序列之间的关系。
# 如果这两个序列并没有包含关系的话，我们使用NaN来代替，这种情况并不多。在ensembl109版本中有21个转录本是这样的，如果突变出现在这样的转录本上的话
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

def tailCdsSeq(x):
    '''
    判断一下ENSP的最后三位是否和基因组对应的位置匹配的上。
    '''
    if isinstance(x['mode'],str):
        if x['mode'][-1]!='X':
            inChrom = h[x["chr"]][end:end+3] if x['strand'] == '1' else str(Seq(h[x["chr"]][end-3:end]).reverse_complement())
            return x['mode'][-3:] == inChrom
        else:return True
    else:return False

logger.info('写出为bed文件和fa文件')
#写出为bed文件和fa文件
if not os.path.exists(refDir):
	os.mkdir(refDir)
def parsePos(x):
	c,p = x.split(':')
	return (c,*p.split('-'))
df.to_csv(os.path.join(refDir,'transcipt'),index=False,sep='\t')
# f.write('chr,start,end,transcript,rank,strand,stopWithCodon,addHead,seq\n'.replace(',','\t'))
'''
输出文件解释:
1,base:chr,start,end,transcript
2,rank:cds的rank，转录本根据这些rank顺序将cds粘贴起来
3，strand:转录本的正负链
4，stopWithCodon:转录本是否是以终止密码子终止的，如果是的话，我们在转录本的最后一个cds的end位点后面加上5000bp的序列；如果不是的话，我们在转录本的最后一个cds的end位点后面不加序列，他是因为别的原因停止转录的。
5，addHead:如果某个转录本序列以;开头，则认为这个转录本在基因组上的序列和ensembl提供的cds序列不一致，这种情况下，我们在基因组上的序列前面加上;，如果为字符串，则在进行翻译的时候cds的head前需要加上这一段。
'''
def isStopWithCodon(r):
    end = int(r["last_cds_end"])
    if isinstance(r['mode'],str):
        if r['mode'][-1]=='X':
            return (len(r['seq'])%3==0)&(r['seq'][-3:] in ['TAG','TGA','TAA'])
        else:
            return (len(r['seq'])%3==0)&tailCdsSeq(r)&(r['seq'][-3:] in ['TAG','TGA','TAA'])
    # stopWithCodon = h[r["chr"]][end:end+3] in ['TAG','TGA','TAA'] if r['strand'] == '1' else h[r["chr"]][end-3:end] in ['CTA','TCA','TTA']
    # stopWithCodon = (r['seq'][-3:] in ['TAG','TGA','TAA'])|stopWithCodon
    # stopWithCodon = stopWithCodon&(len(r['seq'])%3==0)&tailCdsSeq(r)
f = open(os.path.join(refDir,'cds.bed'),'w+')
for i,r in tqdm(df.iterrows()):
	end = int(r["last_cds_end"])#+1
	csr = r['cdsSeqRange'].split('|')[1:-1]
	css = r['cdsSeq'].split('|')
	assert len(csr) == len(css),f'{r["id"]}'
	#stopWithCodon = r['seq'][-3:] in ['TAG','TGA','TAA']
	# 判断最后一个终止子是否是终止密码子的条件两种：
	# 1，最后一个cds的终止子是终止密码子；
	# 2，最后一个cds的结束位点后三位碱基是终止密码子, 比如cdsSeq为X，seq为XTAG,cdsSeq最后一个cds的end位点往后移码三位和TAG要对的上.
	# 3, 整个转录本的长度要被3整除。一些转录本的长度并不能被3整除，但是他的最后三位位点是标注终止密码子，实际上如果在其后面接上DNA序列，它并不会终止，这些只是看起来像的“伪终止子”,实际上他们的stopWithCodon依然为False。
	stopWithCodon =isStopWithCodon(r)
	addHead = r['mode'].split('X')[0] if isinstance(r['mode'],str) else '|'
	#如果某个转录本序列以;开头，则认为这个转录本在基因组上的序列和ensembl提供的cds序列不一致，这种情况下，我们在基因组上的序列前面加上|，
	# 这样在后续的处理中，我们就可以知道这个转录本的序列是不一致的。
	for j,v in enumerate(csr):
		c,s,e = parsePos(v)
		seq = css[j] if r['strand']=='1' else str(Seq(css[j]).reverse_complement())
		f.write(f'{c}\t{int(s)}\t{int(e)}\t{j+1}\t{r["id"]}\t{r["strand"]}\t{stopWithCodon}\t{addHead}\t{seq}\n')
	s,e = (end,end+EXPAND) if r['strand'] == '1' else (end-EXPAND,end)
	seq = h[r["chr"]][s-1:e] # if r['strand']=='1' else str(Seq(h[r["chr"]][s:e]).reverse_complement())
	f.write(f'{c}\t{s}\t{e}\t{len(csr)+1}\t{r["id"]}\t{r["strand"]}\t{stopWithCodon}\t{addHead}\t{seq}\n')
	f.flush()
f.close()

# 拼接的蛋白和ensembl的蛋白序列不一致
# 1，Seq比cdsSeq后面多了一个终止密码子，但是基因组上cdsSeq后三位并不是终止子，这种情况下，认为它停止并不是因为终止子停止的，不要在cds后面append额外的序列，将stopWithCodon变为False