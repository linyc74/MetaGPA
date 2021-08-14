import pandas as pd 
import argparse

# read pfam.gtf into dictionary {'contig_1':['Thymidylat_synt','AAA',...]; 'contig_2':['PRA-PH']; ...}
def read():
    pfam_list = []
    dic = {}
    with open(input,'r') as f:
        for line in f:
            line = line.strip()
            contig_id = line.split('\t')[0]
            pfam = line.split('\t')[8].split('name "')[1].split(' ')[0]
            dic.setdefault(contig_id,[]).append(pfam)
            if not pfam in pfam_list:
                pfam_list.append(pfam)
    return pfam_list, dic

def subset_enriched_pfams(pfam_number, output):
    with open('pfam_a_statistics.csv','r') as f:
        enriched_pfam = []
        header = f.readline()
        for i in range(pfam_number):      # get the top 20 pfams by p-value
            line = f.readline()
            pfam_name = line.strip().split(',')[0].split(' ')[0]
            enriched_pfam.append(pfam_name)
    subset_df = pd.DataFrame(index=enriched_pfam, columns=list(dic))
    for i in enriched_pfam:
        for j in list(dic):
            if i in dic[j]:
                subset_df.at[i,j] = 1
            else:
                subset_df.at[i,j] = 0
    subset_df.to_csv(output, index=True, header=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--number', help='generate matrix for top # pfams', type=int, dest='pfam_number', default=20)
    args = parser.parse_args()
    
    input = 'modified_contigs_pfam_a.gtf'
    output = f'enriched_pfam_a_matrix.csv'

    pfam_list, dic = read()
    subset_enriched_pfams(args.pfam_number, output)