import pandas as pd 

input = 'selected_contigs_pfam_a_95.gtf'
output = 'selected_contigs_pfam_a_PA.csv'

# read pfam.gtf into dictionary {'contig_1':['Thymidylat_synt','AAA',...]; 'contig_2':['PRA-PH']; ...}
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

# generate a binary dataframe with pfams as rows and contigs as columns       
df = pd.DataFrame(index=pfam_list, columns=list(dic))
for i in pfam_list:
    for j in list(dic):
        if i in dic[j]:
            df.at[i,j]=1
        else:
            df.at[i,j]=0

def subset_enriched_pfams(pfam_number, output):
    with open('pfam_a_enrichment_corrected_95.csv','r') as f:
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

df.to_csv(output, index=True, header=True)

subset_enriched_pfams(20, 'top_20_enriched_pfam_PA.csv')