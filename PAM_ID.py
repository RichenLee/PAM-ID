import argparse
import os
from utils import *


def parse_arguments() -> argparse.Namespace:
    #"""解析命令行参数"""
    parser = argparse.ArgumentParser(description='Usage of PAM-ID')
    parser.add_argument('-c', '--config', required=True,
                        help='Parameters configuration file')
    parser.add_argument('-b', '--barcode', required=True,
                        help='Barcode sequence file')
    parser.add_argument('-o', '--output', required=True,
                        help='Output directory')
    return parser.parse_args()

def main():
    args = parse_arguments()
    config_file = args.config
    barcode_file = args.barcode
    output_dir = args.output
    flags=True
    error=''
    for file in [config_file,barcode_file]:
        if not os.path.exists(file):
            flags=False
            error+=file+' not found\n'
    if os.path.exists(output_dir):
        flags=False
        error+=output_dir+' already exists\n'
    else:
        os.makedirs(output_dir)
    if flags==False:
        print('Please check the input files and output directory:\n'+error)
        exit()
    config = read_config(config_file)
    sgRNA_seq=config['sgRNA_seq']
    fastqR1=config['fastqR1']
    fastqR2=config['fastqR2']
    PAM_length=int(config['PAM_length'])
    PAM_orientation=config['PAM_orientation']
    control_sample=config['control_sample']
    cutoff=float(config['log2fc_cutoff'])
    barcode_dic,PAM_dic=get_sample(barcode_file,PAM_length,PAM_orientation)
    print('PAM_ID is running...')
    PAM_dic=read_fastq(fastqR1,fastqR2,sgRNA_seq,barcode_dic,PAM_dic,PAM_orientation,PAM_length)
    for sample,PAM in PAM_dic.items():
        result_file=os.path.join(output_dir,sample+'.tsv')
        with open(result_file,'w') as ot:
            for target,count in PAM.items():
                ot.write(target+'\t'+str(count)+'\n')
    make_diff(output_dir,control_sample)
    draw_logo(output_dir,cutoff)
    print('PAM_ID is done!')
    
if __name__ == '__main__':
    main()