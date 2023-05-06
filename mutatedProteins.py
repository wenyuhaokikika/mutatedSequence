import argparse
import pandas as pd
import os
parser = argparse.ArgumentParser(description='join all cds and translate to protein')
parser.add_argument('-b','--bed',type=str,help='bed file file')
parser.add_argument('-o','--out',type=str,help='bed file file')

args = parser.parse_args()
outPath = args.out
#mutatedBED = './db/mutated.bed'
mutatedBED = args.bed
mutatedSeq = pd.read_table(mutatedBED,header=None)
def parse(x):
	x = x.replace('addHead:;','addHead:|')
	p,i = x.split(' ')
	c,s,e = parsePos(p)
	trans,rank,strand,stopWithCodon,addHead = [i.split(':')[1] for i in i.split(';')]
	return c,s,e,trans,rank,strand,stopWithCodon,addHead
mutatedSeq['chr'],mutatedSeq['start'],mutatedSeq['end'],mutatedSeq['transcript'],mutatedSeq['rank'],mutatedSeq['strand'],mutatedSeq['stopWithCodon'],mutatedSeq['addHead'] = zip(*mutatedSeq[0].apply(parse))
mutatedSeq['stopWithCodon'] = mutatedSeq['stopWithCodon'].apply(eval)
mutatedSeq = mutatedSeq.rename({1:'seq'},axis=1)[['chr', 'start', 'end', 'transcript', 'rank', 'strand','stopWithCodon','seq','addHead']]
mutatedSeq.to_csv(mutatedBED,sep='\t',index=False)

def genateFa(h):
	h = h.sort_values('rank')
	if h['stopWithCodon'].iloc[0] == False:
		h = h.iloc[:-1,:]
	ss = h.apply(lambda x:x['seq'] if x['strand']=='1' else str(Seq(x['seq']).reverse_complement()),axis=1).to_list()
	ss = ''.join(ss)
	addHead = h['addHead'].iloc[0]
	ss = ss+addHead
	return ss
sequences = [SeqRecord(Seq(genateFa(df)),id=trans,description=f'strand:{df.iloc[0,:]["strand"]},start:{df["start"].min()},end:{df["start"].max()}') for trans,df in tqdm(mutatedSeq.groupby('transcript'))]
with open(os.path.join(outPath,'mutatedCDS.fasta'), "w") as output_handle:
    SeqIO.write(sequences, output_handle, "fasta")
def recSeq(x):
	if '|' in x:
		tmp = str(x.seq).replace('|','')
		return '|'+Seq(tmp).translate().split('*',1)[0]
	else:
		return x.seq.translate().split('*',1)[0]
sequences = [SeqRecord(recSeq(i),id=i.id,description = i.description) for i in sequences]
with open(os.path.join(outPath,'mutatedProteins.fasta'), "w") as output_handle:
    SeqIO.write(sequences, output_handle, "fasta")
