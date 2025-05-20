import gzip
import itertools
import os
import numpy as np
import pandas as pd
from weblogo import read_seq_data, LogoData, LogoOptions, LogoFormat , eps_formatter , png_formatter
from weblogo.colorscheme import ColorScheme, SymbolColor
def mkdir(path): 
    folder = os.path.exists(path) 
    if not folder:
        os.makedirs(path)

def get_sample(sample_info_file,PAM_length,PAM_orientation):
	bases=['A', 'T', 'C', 'G']
	barcode_dic={}
	PAM_dic={}
	with open(sample_info_file) as ff:
		for line in ff:
			sp=line.strip().split('\t')
			sample=sp[0]
			PAM_dic[sample]={}
			F_barcode=sp[1]
			R_barcode=sp[2]
			barcode_dic[F_barcode,R_barcode]=sample
			PAM_list=[''.join(p) for p in itertools.product(bases, repeat=PAM_length)]
			for PAM in PAM_list:
				# if PAM_orientation==5:
				# 	target=PAM
				# 	PAM_dic[sample][PAM]=0
				# elif PAM_orientation==3:
				# 	target=PAM
				PAM_dic[sample][PAM]=0
	return barcode_dic,PAM_dic

def reverse_complement(seq):
    nt_complement = dict({'A':'T','C':'G','G':'C','T':'A'})
    return "".join([nt_complement[c] for c in seq.upper()[-1::-1]])

def count_5(sgRNA,read_sequenceR1,read_sequenceR2,barcode_dic,PAM_dic,PAM_length=7):
	#PAM_length=
	F_seq=''
	R_seq=''

	rc_sgRNA=reverse_complement(sgRNA)

	find_in_R1=[read_sequenceR1.find(sgRNA),read_sequenceR2.find(rc_sgRNA)]
	find_in_R2=[read_sequenceR2.find(sgRNA),read_sequenceR1.find(rc_sgRNA)]

	if (find_in_R1[0]!=-1) and (find_in_R1[1]!=-1) :
		PAM=read_sequenceR1[find_in_R1[0]-PAM_length:find_in_R1[0]]
		F_seq=read_sequenceR1
		R_seq=read_sequenceR2
		# if len(PAM)!=PAM_length:
		# 	PAM=R_seq[find_in_R1[1]+len(sgRNA):find_in_R1[1]+len(sgRNA)+7]
		#print(F_seq)
		#print(R_seq)
	elif (find_in_R2[0]!=-1) and (find_in_R2[1]!=-1):
		PAM=read_sequenceR2[find_in_R2[0]-PAM_length:find_in_R2[0]]
		F_seq=read_sequenceR2
		R_seq=read_sequenceR1
		# if len(PAM)!=PAM_length:
		# 	PAM=R_seq[find_in_R1[1]+len(sgRNA):find_in_R1[1]+len(sgRNA)+PAM_length]
		#print(F_seq)
		#print(R_seq)
	if F_seq=='' or ('N' in PAM):
		pass
	else:
		for barcode,sample in barcode_dic.items():
			if (F_seq.find(barcode[0])!=-1) & (R_seq.find(barcode[1])!=-1):
				#print(sample,PAM)
				if len(PAM)==PAM_length:

					PAM_dic[sample][PAM]+=1

				break

	return(PAM_dic)

def count_3(sgRNA,read_sequenceR1,read_sequenceR2,barcode_dic,PAM_dic,PAM_length=7):
	#PAM_length=
	F_seq=''
	R_seq=''

	rc_sgRNA=reverse_complement(sgRNA)

	find_in_R1=[read_sequenceR1.find(sgRNA),read_sequenceR2.find(rc_sgRNA)]
	find_in_R2=[read_sequenceR2.find(sgRNA),read_sequenceR1.find(rc_sgRNA)]

	if (find_in_R1[0]!=-1) and (find_in_R1[1]!=-1) :
		PAM=read_sequenceR1[find_in_R1[0]+len(sgRNA):find_in_R1[0]+len(sgRNA)+PAM_length]
		F_seq=read_sequenceR1
		R_seq=read_sequenceR2
		# if len(PAM)!=PAM_length:
		# 	PAM=R_seq[find_in_R1[1]+len(sgRNA):find_in_R1[1]+len(sgRNA)+7]
		#print(F_seq)
		#print(R_seq)
	elif (find_in_R2[0]!=-1) and (find_in_R2[1]!=-1):
		PAM=read_sequenceR2[find_in_R2[0]+len(sgRNA):find_in_R2[0]+len(sgRNA)+PAM_length]
		F_seq=read_sequenceR2
		R_seq=read_sequenceR1
		# if len(PAM)!=PAM_length:
		# 	PAM=R_seq[find_in_R1[1]+len(sgRNA):find_in_R1[1]+len(sgRNA)+PAM_length]
		#print(F_seq)
		#print(R_seq)
	if F_seq=='' or ('N' in PAM):
		pass
	else:
		for barcode,sample in barcode_dic.items():
			if (F_seq.find(barcode[0])!=-1) & (R_seq.find(barcode[1])!=-1):
				#print(sample,PAM)
				if len(PAM)==PAM_length:

					PAM_dic[sample][PAM]+=1

				break

	return(PAM_dic)


def read_fastq(fastqR1,fastqR2,sgRNA,barcode_dic,PAM_dic,PAM_orientation,PAM_length):
	infileR1=gzip.open(fastqR1, 'rt')
	infileR2=gzip.open(fastqR2, 'rt')
	total_reads=0
	if PAM_orientation=='5':
		while infileR1.readline() and infileR2.readline():
			read_sequenceR1 = infileR1.readline().strip()
			infileR1.readline()
			infileR1.readline()
			read_sequenceR2 = infileR2.readline().strip()
			infileR2.readline()
			infileR2.readline()
			total_reads += 1
			PAM_dic=count_5(sgRNA,read_sequenceR1,read_sequenceR2,barcode_dic,PAM_dic)
	elif PAM_orientation=='3':
		while infileR1.readline() and infileR2.readline():
			read_sequenceR1 = infileR1.readline().strip()
			infileR1.readline()
			infileR1.readline()
			read_sequenceR2 = infileR2.readline().strip()
			infileR2.readline()
			infileR2.readline()
			total_reads += 1
			PAM_dic=count_3(sgRNA,read_sequenceR1,read_sequenceR2,barcode_dic,PAM_dic)


	infileR2.close()
	infileR1.close()

	return(PAM_dic)

def read_config(config_file):
	ff=open(config_file)
	f=ff.readlines()
	ff.close()
	f=[i.strip().split('\t') for i in f]
	config={i[0]:i[1] for i in f}
	return config

def make_diff(tsv_dir,control_sample):
	files=os.listdir(tsv_dir)
	files=[i for i in files if '.tsv' in i]
	files.remove(control_sample+'.tsv')
	for cas in files:
		name=cas
		cas=os.path.join(tsv_dir,cas)
		df2=pd.read_table(cas,header=None)
		total=df2[1].sum()
		df2[1]=df2[1]/total
		con=os.path.join(tsv_dir,control_sample+'.tsv')
		df1=pd.read_table(con,header=None)
		total=df1[1].sum()
		df1[1]=df1[1]/total
		df=pd.merge(df1,df2,on=0)
		df['log2']=(df['1_x']/df['1_y']).apply(np.log2)
		df=df.sort_values(by='log2',ascending=False)
		df.columns=['PAM','control','cas','log2']
		df.to_csv(cas.replace('.tsv','.csv'),index=0)

def draw_logo(tsv_dir,cutoff):
	files=os.listdir(tsv_dir)
	files=[i for i in files if '.csv' in i]
	for file in files:
		file=os.path.join(tsv_dir,file)
		tmp=open(f'./tmp.fa','w')
		with open(file) as ff:
			n=0
			for line in ff:
				break
			for line in ff:
				sp=line.split(',')
				if len(sp)==4:
					seq=sp[0]
				else:
					seq=sp[1]
				log2_fold_change=sp[-1].strip()
				try:
					if log2_fold_change=='inf' or float(log2_fold_change)>=cutoff: #3.5
						n+=1
						tmp.write(f">{n}\n{seq}\n")
				except:
					pass
		tmp.close()
		with open(f'./tmp.fa','r') as ff:
			seq_data = read_seq_data(ff)
		data=LogoData.from_seqs(seq_data)
		options = LogoOptions()
		options.fineprint = ""
		options.number_interval = 1
		options.resolution = 1024 
		#options.tic_length = 10
		groups = [
			SymbolColor("G", "orange", "Guanine"),
			SymbolColor("A", "green", "Adenine"),
			SymbolColor("T", "red", "Thymine"),
			SymbolColor("C", "blue", "Cytosine")
		]
		options.show_boxes = False
		dna_scheme = ColorScheme(groups, "DNA")
		options.color_scheme = dna_scheme
		logo_format = LogoFormat(data, options)
		formatter_pdf = eps_formatter(data,logo_format)
		formatter_png = png_formatter(data,logo_format)
		with open(f"{file.replace('.csv','.eps')}", "wb") as f:
			f.write(formatter_pdf)
		with open(f"{file.replace('.csv','.png')}", "wb") as f:
			f.write(formatter_png)
		os.remove('./tmp.fa')
