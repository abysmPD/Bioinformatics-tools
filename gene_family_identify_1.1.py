from optparse import OptionParser
import os
import pandas as pd
from Bio import SeqIO
import re

parser = OptionParser()
parser.add_option("--hmm", action='store', type="string", default='blank')
parser.add_option("--blastp", action='store', type="string", default='blank')
parser.add_option("--pfam", action='store', type="string",
                  default='blank')  # 输入pfam号，如PF00504. 如果有两个，输入PF00504-PF00504
parser.add_option("--gff", action='store', type="string", default='blank')
parser.add_option("--pep", action='store', type="string")
parser.add_option("--blastp_ts", action='store', type="string", default='40')

opts, args = parser.parse_args()

def biopython_fasta_to_df(_address : str):  #提取fasta文件生成一个数据框
    seq = []
    names = []
    for seq_record in SeqIO.parse(_address, "fasta"):
        seq.append(str(seq_record.seq))
        names.append(seq_record.id)
    df_out = pd.DataFrame({'name': names, 'seq': seq})
    return df_out

def gff2df(_gff):
    os.system('''awk -F '\t' '$1!~/#/ && $3=="mRNA" {print $0}' '''+_gff+'> mrna.gff')
    file_in = open('mrna.gff', 'r')
    mrna_list = []
    gene_list = []
    for line in file_in:
        mrna = re.findall(r'ID=(.*?)[;\n]', line)[0]
        gene = re.findall(r'Parent=(.*?)[;\n]', line)[0]
        mrna_list.append(mrna)
        gene_list.append(gene)
    file_in.close()
    df_out = pd.DataFrame({'gene':gene_list, 'mrna':mrna_list})
    return df_out

def extract_domain_seq(_hmm_out):
    os.system('''awk -F ' ' '$1!~/#/ {print $1"\t"$7"\t"$18"\t"$19}' '''+_hmm_out+'> first_hmm_filter.out')
    df_hmm = pd.read_table('first_hmm_filter.out', header=None, sep='\t')
    df_hmm = df_hmm.loc[df_hmm.loc[:, 1]<1e-5, :]
    df_pep = biopython_fasta_to_df(opts.pep)
    df_use = pd.merge(left=df_hmm, right=df_pep, left_on=0, right_on='name', how='left')
    file_out = open('first_hmm_domain_seqs.fa', 'w')
    for i, line in enumerate(df_use.loc[:, 'seq']):
        start = int(df_use.loc[i, 2])-1
        end = int(df_use.loc[i, 3])-1
        file_out.write('>'+df_use.loc[i, 'name']+'_'+str(i+1)+'\n')
        file_out.write(df_use.loc[i, 'seq'][start:end]+'\n')
    file_out.close()
    os.system('rm first_hmm_filter.out')

def hmm_search(_hmm, _pep):
    os.system('hmmsearch --domtblout first_hmm.out --cut_tc '+_hmm+' '+_pep)
    extract_domain_seq('first_hmm.out')
    os.system('clustalo -i first_hmm_domain_seqs.fa -o first_hmm_domain_seqs_align.fa')
    os.system('hmmbuild second.hmm first_hmm_domain_seqs_align.fa')
    os.system('hmmsearch --domtblout second_hmm.out second.hmm '+_pep)
    os.system('''awk -F ' ' '$1!~/#/ && $7<0.001 {print $1}' second_hmm.out | sort -u'''+'> second_hmm_gene.out')
    
def blastp_search(_blastp, _pep):
    os.system('makeblastdb -in '+_pep+' -input_type fasta -dbtype prot -title Base -parse_seqids -out Base')
    os.system('blastp -query '+_blastp+' -db Base -out blastp.out -evalue 1e-5 -num_threads 20 -outfmt 6')
    os.system(''' awk -F '\t' '$3>''' + opts.blastp_ts +
              ''' && $11<1e-5 {print $2}' blastp.out | sort -u > blastp_gene.out ''')

def extract_primary_mrna():
    #准备df文件
    df_mrna_gene = gff2df(opts.gff)
    df_pep = biopython_fasta_to_df(opts.pep)
    df_pep.loc[:, 'length'] = df_pep.loc[:, 'seq'].apply(len)
    df_merge = pd.merge(left=df_mrna_gene,right=df_pep,left_on='mrna',right_on='name',how='left')
    #df_merge_sort = df_merge.sort_values(by=['gene', 'length'], ascending=True).reset_index(drop=True)
    #df_merge_sort_duplicated = df_merge_sort.drop_duplicates(subset='gene', keep='last').reset_index(drop=True)
    #对final_gene去重
    df_gene = pd.read_table('final_genes_interpro_filter_genes.txt', header=None)
    #df_mrna = df_merge.loc[:, ['gene','mrna']]
    #df_primary = df_merge_sort_duplicated.loc[:, ['gene','mrna']]
    gene_merge_1 = pd.merge(left=df_gene, right=df_merge, left_on=0, right_on='mrna', how='left')
    gene_merge_2 = gene_merge_1.sort_values(by=['gene', 'length'], ascending=True).drop_duplicates(subset='gene', keep='last').reset_index(drop=True)
    #gene_final = pd.DataFrame(set(gene_merge_1.loc[:, 'gene'].tolist()))
    #gene_merge_2 = pd.merge(left=gene_final,right=df_primary,left_on=0,right_on='gene',how='left')
    out = pd.DataFrame(gene_merge_2.loc[:, 'mrna'].tolist())
    return out

def interpro_extract_domain():
    df = pd.read_table('interpro.tsv', header=None)
    final = pd.read_table('Final_primary_interpro.txt', header=None, sep='\t')
    df_final1 = pd.merge(left=final,right=df,on=0,how='left')
    df_final2 = df_final1.loc[df_final1.loc[:, 4] == opts.pfam, :]
    df_final2_group = df_final2.groupby(0)
    df_pep = biopython_fasta_to_df(opts.pep)
    pep = pd.merge(left=final,right=df_pep,left_on=0,right_on='name',how='left')
    file_out = open('result/Final_primary_interpro_filter_domain.fa', 'w')
    for name, group in df_final2_group:
        group.reset_index(drop=True, inplace=True)
        if group.shape[0] == 1:
            file_out.write('>'+name+'\n')
            start = group.loc[0, 6]
            end = group.loc[0, 7]
            file_out.write(pep.loc[pep.loc[:, 'name'] == name, 'seq'].iloc[0, ][start-1:end]+'\n')
        else:
            num = 0
            for i, line in enumerate(group.loc[:, 0]):
                num = num + 1
                file_out.write('>' + name + '-' + str(num) + '\n')
                start = group.loc[i, 6]
                end = group.loc[i, 7]
                file_out.write(pep.loc[pep.loc[:, 'name'] == name, 'seq'].iloc[0, ][start - 1:end] + '\n')
    file_out.close()

def main():
    if opts.hmm != 'blank':
        hmm_search(opts.hmm, opts.pep)
        blastp_search(opts.blastp, opts.pep)
        os.system('cat blastp_gene.out second_hmm_gene.out | sort -u > final_genes.txt')
    else:
        blastp_search(opts.blastp, opts.pep)
        os.system('cat blastp_gene.out | sort -u > final_genes.txt')

    os.system('seqkit grep -f final_genes.txt ' + opts.pep + ' > final_genes.fa')
    #interproscan对所有转录本筛选结构域
    df_pep = biopython_fasta_to_df(opts.pep)
    os.system(
        '/home/shuoliu/program/interproscan-5.39-77.0/interproscan.sh -i final_genes.fa -f tsv -appl Pfam -dp -b interpro')
    if opts.pfam != 'blank':
        if '-' in opts.pfam:
            pfam1 =opts.pfam.split('-')[0]
            pfam2 = opts.pfam.split('-')[1]
            inter_df = pd.read_table('interpro.tsv', header=None)
            inter_df_group = inter_df.groupby(0)
            inter_genes = []
            for name, group in inter_df_group:
                pfam_list = group.loc[:, 4].tolist()
                if (pfam1 in pfam_list) and (pfam2 in pfam_list):
                    inter_genes.append(name)
            inter_out = pd.DataFrame(inter_genes)
            inter_out.to_csv('Final_primary_interpro.txt', header=None, index=None)
            final_1 = pd.read_table('Final_primary_interpro.txt', header=None, sep='\t')
            final_2 = pd.merge(left=final_1, right=df_pep, left_on=0, right_on='name', how='left')
            #os.system('mkdir result')
            # file_out_2 = open('final_genes_interpro_filter.fa', 'w')
            # for i, line in enumerate(final_2.loc[:, 'name']):
            #     file_out_2.write('>' + final_2.loc[i, 'name'] + '\n')
            #     file_out_2.write(final_2.loc[i, 'seq'] + '\n')
            # file_out_2.close()
            file_out_3 = open('final_genes_interpro_filter_genes.txt', 'w')
            for i, line in enumerate(final_2.loc[:, 'name']):
                file_out_3.write(final_2.loc[i, 'name'] + '\n')
            file_out_3.close()
            # 只提取结构域
            #interpro_extract_domain()
        else:
            # os.system(
            #     ''' awk -F '\t' '$5=="''' + opts.pfam + '''" && $9 < 1e-10 {print $1}' interpro.tsv | sort -u > Final_primary_interpro.txt''')
            os.system(
                ''' awk -F '\t' '$5=="''' + opts.pfam + '''" {print $1}' interpro.tsv | sort -u > Final_primary_interpro.txt''')
            final_1 = pd.read_table('Final_primary_interpro.txt', header=None, sep='\t')
            final_2 = pd.merge(left=final_1, right=df_pep, left_on=0, right_on='name', how='left')
            #os.system('mkdir result')
            # file_out_2 = open('final_genes_interpro_filter.fa', 'w')
            # for i, line in enumerate(final_2.loc[:, 'name']):
            #     file_out_2.write('>' + final_2.loc[i, 'name'] + '\n')
            #     file_out_2.write(final_2.loc[i, 'seq'] + '\n')
            # file_out_2.close()
            file_out_3 = open('final_genes_interpro_filter_genes.txt', 'w')
            for i, line in enumerate(final_2.loc[:, 'name']):
                file_out_3.write(final_2.loc[i, 'name'] + '\n')
            file_out_3.close()
            # 只提取结构域
            #interpro_extract_domain()
    else:
        print('No pfam!')
    #筛选包含结构域结果中的最长转录本
    if opts.gff == 'blank':
        final_primary = pd.read_table('final_genes_interpro_filter_genes.txt', header=None)
    else:
        final_primary = extract_primary_mrna()  #对blastp和hmmsearch结果只保留最长转录本
    df_out = pd.merge(left=final_primary, right=df_pep, left_on=0, right_on='name', how='left')
    os.system('mkdir result')
    file_out = open('result/Final_genes.fa', 'w')
    for i, line in enumerate(df_out.loc[:, 'name']):
        file_out.write('>'+df_out.loc[i, 'name']+'\n')
        file_out.write(df_out.loc[i, 'seq']+'\n')
    file_out.close()
    file_out_5 = open('result/Final_genes.gene', 'w')
    for i, line in enumerate(df_out.loc[:, 'name']):
        file_out_5.write(df_out.loc[i, 'name']+'\n')
    file_out_5.close()

if __name__ == "__main__":
    main()