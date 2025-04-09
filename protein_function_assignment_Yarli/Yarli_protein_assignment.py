import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

### Define global variables here
# Define the color palette
BNP_colors = ['cornflowerblue' ,'#ffa510',plt.cm.Greens(np.linspace(0.2, 1, 5))[2]]

### Define the global functions here

def export_color_code(cstring):
    """
    Export color code in the form of f8f9fa-e9ecef-dee2e6-ced4da-adb5bd-6c757d-495057-343a40-212529
    Into usable color lists
    """
    return ['#'+ x for x in cstring.split('-')]

def find_replicates(lst):
    """
    Takes a list of strings and returns a list of sets where each set represents a set of replicates.
    """
    groups = defaultdict(set)
    for item in lst:
        # Extract the prefix and suffix from the string
        prefix, suffix = item.rsplit('_', 1)
        # Add the item to the corresponding group based on its prefix
        groups[prefix].add(item)

    # Identify sets of replicates and groups that contain all replicates
    all_replicates = set()
    replicates = []
    for group in groups.values():
        if len(group) > 0:
            replicates.append(group)
            if len(group) == len(lst):
                all_replicates = group
    
    # Return the replicates and the group that contains all replicates
    return replicates, all_replicates

# ----------- Import GEM model, Uniprot and proteomics data -----------
### GEM model from Nelson lab, 2021
GEM_rxns = pd.read_excel('iYali_nelson_GEM.xlsx', sheet_name='RXNS')
## Uniprot data
unprotdf = pd.read_csv('uniprotkb_proteome_UP000001300_2024_11_19.tsv',sep='\t').set_index('Entry')
unprotdf.insert(3, 'gene_name_yali', np.nan)
unprotdf['gene_name_yali'] = unprotdf.apply(
    lambda row: next((gene for gene in [row['Gene Names (ordered locus)'], row['Gene Names (ORF)'], row['Gene Names (primary)']] if isinstance(gene, str) and 'YALI0' in gene), np.nan), axis=1)
unprotdf['gene_name_yali'] = unprotdf['gene_name_yali'].apply(lambda x: x.replace('_', '') if isinstance(x, str) else x)
### Proteomics with absolute quantification
proteomics = pd.read_excel('absolute_yarli.xlsx')

# Mapped the gene name for each entry in the proteomics data, from Uniprot 0000001300
# Some missing value needed to be manually cured
# Saved as mappped_proteomics.xlsx. But again, will be manually edited
# Drop the rows where 'Protein Id' contains 'HUMAN' and print the dropped rows
dropped_rows = proteomics[proteomics['Protein Id'].str.contains('HUMAN')]
proteomics = proteomics[~proteomics['Protein Id'].str.contains('HUMAN')]
proteomics = proteomics[['Protein Id',"fitted Con nM",'Relative abundance']]
# Add a new column 'gene name' after 'Protein Id'
proteomics.insert(proteomics.columns.get_loc('Protein Id') + 1, 'entry', np.nan)
proteomics.insert(proteomics.columns.get_loc('Protein Id') + 2, 'gene name', np.nan)
proteomics.insert(proteomics.columns.get_loc('Protein Id') + 3, 'protein_name_uniprot', np.nan)
proteomics.insert(proteomics.columns.get_loc('Protein Id') + 4, 'protein_name_GEM', np.nan)
proteomics.insert(proteomics.columns.get_loc('Protein Id') + 5, 'Gene Ontology (GO)', np.nan)
proteomics.insert(proteomics.columns.get_loc('Protein Id') + 6, 'subsystem', np.nan)
proteomics.insert(proteomics.columns.get_loc('Protein Id') + 7, 'Length', np.nan)

proteomics['entry'] = proteomics['Protein Id'].apply(lambda x: x.split('_')[1] if '_' in x else np.nan)

## Mapping proteomics to uniprot, using entry, attaching the gene name
proteomics['gene name'] = proteomics['entry'].map(unprotdf['gene_name_yali'])
proteomics['Gene Ontology (GO)'] = proteomics['entry'].map(unprotdf['Gene Ontology (GO)'])
proteomics['Gene Ontology (GO)'] = proteomics['Gene Ontology (GO)'].fillna('')

proteomics['protein_name_uniprot'] = proteomics['entry'].map(unprotdf['Protein names'])
proteomics['Length'] = proteomics['entry'].map(unprotdf['Length'])
proteomics['Fitcon_multi_length'] = proteomics['Length']*proteomics['fitted Con nM']
proteomics['Relative abundance'] = proteomics['Fitcon_multi_length'].div(proteomics['Fitcon_multi_length'].sum())

proteomics.to_excel('20250408mapped_yl_abs_proteomics_wlength.xlsx', index=False) 

#----------- Mapping proteomics to GEM, using gene name, attahcing the GEM protein name and subsystem
#proteomics['protein_name_GEM'] = proteomics['gene name'].map(GEM_rxns.set_index('GENE ASSOCIATION')['NAME'])
# Mapping proteomics to GEM, using gene name, attaching the GEM protein name
proteomics['protein_name_GEM'] = proteomics['gene name'].apply(
    lambda gene: GEM_rxns[GEM_rxns['GENE ASSOCIATION'].str.contains(gene, na=False)]['NAME'].tolist() if pd.notna(gene) else np.nan
)
proteomics['subsystem'] = proteomics['gene name'].apply(
    lambda gene: GEM_rxns[GEM_rxns['GENE ASSOCIATION'].str.contains(gene, na=False)]['SUBSYSTEM'].tolist() if pd.notna(gene) else np.nan
)
proteomics['subsystem'] = proteomics['subsystem'].apply(lambda x: list(set(x)) if isinstance(x, list) else x)
proteomics['subsystem'] = proteomics['subsystem'].apply(lambda x: [] if x == [np.nan] else x)


#-----------Assign GEM's ubsystem to functional sectors ----------------------------
### Logic. Map from GEM model first, then protein names/Gene Go term 
proteomics['protein_name_GEM'] = proteomics['protein_name_GEM'].apply(lambda x: ';'.join([str(i) for i in x]) if isinstance(x, list) else ('' if isinstance(x, float) else x))
# Iterate over each row in the proteomics DataFrame
for index, row in proteomics.iterrows():
    if not row['subsystem'] or row['subsystem'] == [''] or row['subsystem'] is np.nan:
   
        if 'atp synthase' in row['protein_name_uniprot'].lower() or 'cytochrome c' in row['protein_name_uniprot'].lower() or 'cytochrome b-c' in row['protein_name_uniprot'].lower():
            proteomics.at[index, 'subsystem'] = ['oxidative phosphorylation']
        elif 'pyruvate dehydrogenase' in row['protein_name_uniprot'].lower():
            proteomics.at[index, 'subsystem'] = ['TCA']
        elif  any(keyword in row['protein_name_uniprot'].lower() for keyword in ['translation', 'ribosom', 'elongation factor','elongator','trna']):
            proteomics.at[index, 'subsystem'] = ['translation/ribosome']
        elif  any(keyword in row['protein_name_uniprot'].lower() for keyword in ['rna helicase','pre-rrna-processing']):
            proteomics.at[index, 'subsystem'] = ['translation/ribosome']
        elif  any(keyword in row['Gene Ontology (GO)'].lower() for keyword in ['trna synthase']):
            proteomics.at[index, 'subsystem'] = ['translation/ribosome']
        elif any(keyword in row['protein_name_uniprot'].lower() for keyword in ['lipid', 'fatty acid', 'fatty acyl','aceyl-coa','acylglycerol synthase']):
            proteomics.at[index, 'subsystem'] = ['lipid/fatty acid metabolism']
        elif any(keyword in row['protein_name_uniprot'].lower() for keyword in ['transcription', 'mrna splicing','mrna-splic','chromatin']):
            proteomics.at[index, 'subsystem'] = ['transcription']
        elif any(keyword in row['protein_name_uniprot'].lower() for keyword in ['ase','biosynthesis','biogenesis','sulfiredoxin']) or  any(keyword in row['protein_name_GEM'].lower() for keyword in ['ase']):#or any(keyword in row['Gene Ontology (GO)'].lower() for keyword in ['ase activity']):
            proteomics.at[index, 'subsystem'] = ['other enzymes']
        elif any(keyword in row['protein_name_uniprot'].lower() for keyword in ['transport', 'traffic','translocation']):
            proteomics.at[index, 'subsystem'] = ['transporting']    
        elif any(keyword in row['protein_name_uniprot'].lower() for keyword in ['histone', 'actin', 'tubulin', 'myosin','ubiquitin','signal']):
            proteomics.at[index, 'subsystem'] = ['other proteins']
        elif any(keyword in row['Gene Ontology (GO)'].lower() for keyword in ['fatty acid metabolic','lipid metabolic']):
            proteomics.at[index, 'subsystem'] = ['lipid/fatty acid metabolism']
        elif any(keyword in row['Gene Ontology (GO)'].lower() for keyword in ['lipid','fatty acid','fatty acyl']) and all(keyword not in row['Gene Ontology (GO)'].lower() for keyword in ['membrane']):
            proteomics.at[index, 'subsystem'] = ['lipid/fatty acid metabolism']
        elif any(keyword in row['Gene Ontology (GO)'].lower() for keyword in ['translation','ribosom','elongation factor','elongator']):
            proteomics.at[index, 'subsystem'] = ['translation/ribosome']
        elif any(keyword in row['Gene Ontology (GO)'].lower() for keyword in ['transcription','mrna splicing','mrna-splic']):
            proteomics.at[index, 'subsystem'] = ['transcription']
        elif any(keyword in row['protein_name_GEM'].lower() for keyword in ['transport','traffic','translocation','exchange']):
            proteomics.at[index, 'subsystem'] = ['transporting']
        elif any(keyword in row['protein_name_uniprot'].lower() for keyword in ['transport','traffic','translocation','exchange']):
            proteomics.at[index, 'subsystem'] = ['transporting']
        elif any(keyword in row['protein_name_uniprot'].lower() for keyword in ['autophagy','autophagic']):
            proteomics.at[index, 'subsystem'] = ['autophagy']
        elif any(keyword in row['protein_name_uniprot'].lower() for keyword in ['membrane']):
               proteomics.at[index, 'subsystem'] = ['membrane']
        elif 'YALI' in row['protein_name_uniprot'] and  row['Gene Ontology (GO)'] == 'membrane [GO:0016020]':
            proteomics.at[index, 'subsystem'] = ['membrane']
        elif any(keyword in row['protein_name_uniprot'].lower() for keyword in ['glycosylphosphatidylinositol-anchor biosynthesis']):
            proteomics.at[index, 'subsystem'] = ['Glycosylphosphatidylinositol (GPI)-anchor biosynthesis']
        elif any(keyword in row['protein_name_uniprot'].lower() for keyword in ['vesicular','vacuolar','endosomal','endoplasmic reticulum','response','proteasom','regulat','binding']):
            proteomics.at[index, 'subsystem'] = ['other proteins']
        elif 'YALI' not in row['protein_name_uniprot'] and row['protein_name_uniprot'] != '':
            proteomics.at[index, 'subsystem'] = ['other proteins']
        elif any(keyword in row['Gene Ontology (GO)'].lower() for keyword in ['translation']):
            proteomics.at[index, 'subsystem'] = ['translation/ribosome']
        elif any(keyword in row['Gene Ontology (GO)'].lower() for keyword in ['lipid','aceyl-coa','fatty acid']):
            proteomics.at[index, 'subsystem'] = ['lipid/fatty acid metabolism']
        elif any(keyword in row['Gene Ontology (GO)'].lower() for keyword in ['amino acid']):
            proteomics.at[index, 'subsystem'] = ['amino acid metabolism']
        elif any(keyword in row['Gene Ontology (GO)'].lower() for keyword in ['respiratory', 'cytochrome c']):
            proteomics.at[index, 'subsystem'] = ['respiration']
        elif any(keyword in row['Gene Ontology (GO)'].lower() for keyword in ['ase']) and any(keyword in row['Gene Ontology (GO)'].lower() for keyword in ['mitocho']):
            proteomics.at[index, 'subsystem'] = ['mitochondrial enzymes']
        elif any(keyword in row['Gene Ontology (GO)'].lower() for keyword in ['ase']):
            proteomics.at[index, 'subsystem'] = ['other enzymes']
        elif any(keyword in row['Gene Ontology (GO)'].lower() for keyword in ['mitochon']):
            proteomics.at[index, 'subsystem'] = ['other mitochondrial protein']
        elif any(keyword in row['Gene Ontology (GO)'].lower() for keyword in ['membrane']):
            proteomics.at[index, 'subsystem'] = ['membrane']
        elif row['Gene Ontology (GO)'] != '':
            proteomics.at[index, 'subsystem'] = ['other proteins']
        elif row['Gene Ontology (GO)'] == 'cytoplasm [GO:0005737]'  or row['Gene Ontology (GO)'] ==  'nucleus [GO:0005634]':
            proteomics.at[index, 'subsystem'] = ['other proteins']
        elif 'YALI' in row['protein_name_uniprot'] and row['Gene Ontology (GO)'] == '' and row['protein_name_GEM'] == '':
            proteomics.at[index, 'subsystem'] = 'unannotated'

proteomics['subsystem'] = proteomics['subsystem'].apply(lambda x: '|'.join([str(i) for i in x]) if isinstance(x, list) else ('' if isinstance(x, float) else x))

# -------------------------- Assign functional groups here -----------------------------


for index, row in proteomics.iterrows():
         # Assign 'ATP synthase' to the subsystem column
        if any(keyword in row['subsystem'].lower() for keyword in ['translation', 'ribosom', 'elongation factor','elongator','aminoacyl-trna biosynthesis']):
            proteomics.at[index,'functional group'] = 'translation'
        elif any(keyword in row['subsystem'].lower() for keyword in ['pentose phosphate pathway']):
            proteomics.at[index,'functional group'] = 'PPP'
        elif row['subsystem'] in ['fatty acid degradation', 'Fatty acid degradation','Triacylglycerol degradation']:
            proteomics.at[index, 'functional group'] = 'Lipid degradation'
       
        elif any(keyword in row['subsystem'].lower() for keyword in ['lipid', 'fatty acid', 'fatty acyl','glycosylphosphatidylinositol (gpi)-anchor biosynthesis','steroid biosynthesis','ergosterol biosynthesis','terpenoid backbone']):
            proteomics.at[index,'functional group'] = 'lipid pathway'

         #elif 'translation' in row['subsystem'].lower() or 'ribosom' in row['subsystem'].lower() or 'elongation factor' in row['subsystem'].lower():
        elif  any(keyword in row['subsystem'].lower() for keyword in ['glycolysis','tca','oxidative phospho','respirat','glyoxylate']):
            proteomics.at[index,'functional group'] = 'energy pathway'
        elif any(keyword in row['subsystem'].lower() for keyword in ['transcription', 'mrna splicing','mrna-splic']):
            proteomics.at[index,'functional group'] = 'transcription'


        elif any(keyword in row['subsystem'].lower() for keyword in ['transporting','other proteins','membrane','autophagy','mitochondrial protein']):
            proteomics.at[index,'functional group'] = 'other proteins'

        elif any(keyword in row['subsystem'].lower() for keyword in ['unannotated']):
            proteomics.at[index, 'functional group'] = 'unannotated'
        elif row['subsystem'] != '':
            proteomics.at[index, 'functional group'] = 'other enzymes'
      


### Do some manual assignment
##### Set ATP citrate lyase, the enzyme producing acetyl-CoA, as lipid pathway
proteomics.loc[ proteomics['protein_name_GEM'].notna() &  proteomics['protein_name_GEM'].str.contains('ATP-citrate lyase'), 'functional group'] = 'AcCoA biosynthesis'
proteomics.loc[ proteomics['protein_name_uniprot'].notna() &  proteomics['protein_name_uniprot'].str.contains('yruvate dehydrogenase'), 'functional group'] = 'AcCoA biosynthesis'
proteomics.loc[ proteomics['protein_name_uniprot'].notna() &  proteomics['protein_name_uniprot'].str.contains('Dihydrolipoyl dehydrogenase'), 'functional group'] = 'AcCoA biosynthesis'
proteomics.loc[ proteomics['protein_name_uniprot'].notna() &  proteomics['protein_name_uniprot'].str.contains('Pyruvate carboxylase'), 'functional group'] = 'AcCoA biosynthesis'
proteomics.loc[ proteomics['protein_name_uniprot'].notna() &  proteomics['protein_name_uniprot'].str.contains('Acetyl-CoA hydrolase'), 'functional group'] = 'AcCoA biosynthesis'
proteomics.loc[ proteomics['protein_name_uniprot'].notna() &  proteomics['protein_name_uniprot'].str.contains('Acetyl-coenzyme A synthetase'), 'functional group'] = 'AcCoA biosynthesis'
proteomics.loc[ proteomics['protein_name_uniprot'].notna() &  proteomics['protein_name_uniprot'].str.contains('hydroxymethylglutaryl-CoA lyase'), 'functional group'] = 'AcCoA biosynthesis'


# Export the rows, if funciontal group is lipid pathway, and fatty acid degradation is in the subsystem


# ---------------------- Differentiate fatty acid degradation and biosynthesis --------------------------
##### Screen the beta-oxidation pathway proteins here #####

temp1 = proteomics[(proteomics['functional group'] == 'lipid pathway') & (proteomics['subsystem'].str.contains('degradation')) & (~proteomics['subsystem'].str.contains('Biosynthesis of unsaturated fatty acids')& (~proteomics['subsystem'].str.contains('Fatty acid biosynthesis')))] # pick

temp2 = proteomics[(proteomics['functional group'] != 'lipid pathway') & (proteomics['protein_name_uniprot'].str.contains('CoA dehydrogenase')) ] # pick

temp3 = proteomics[(proteomics['functional group'] != 'lipid pathway') & (proteomics['Gene Ontology (GO)'].str.contains('acyl-CoA dehydrogenase')) ] # pick

temp4 = proteomics[(proteomics['functional group'] != 'lipid pathway') & (proteomics['protein_name_uniprot'].str.contains('enoyl-CoA'))] # pick

temp5  = proteomics[(proteomics['functional group'] != 'lipid pathway') & (proteomics['protein_name_uniprot'].str.contains('CoA thiolase'))] # pick

temp6 = proteomics[(proteomics['functional group'] == 'lipid pathway') & (proteomics['Gene Ontology (GO)'].str.contains('acyl-CoA dehydrogenase')) ] # pick

temp7 =  proteomics[( proteomics['protein_name_uniprot'].str.contains('Acyl-coenzyme A oxidase')) | ( proteomics['protein_name_GEM'].str.contains('Acyl-coenzyme A oxidase')) |(proteomics['protein_name_uniprot'].str.contains('3-ketoacyl-CoA thiolase')) |(proteomics['protein_name_uniprot'].str.contains('lipase'))] # pick # acylglycerol lipase'
# extract the index of the temp dataframes above, and then merge them together
temp1_index = temp1.index
temp2_index = temp2.index
temp3_index = temp3.index
temp4_index = temp4.index
temp5_index = temp5.index
temp6_index = temp6.index
temp7_index = temp7.index
temp_index = temp1_index.append(temp2_index)
temp_index = temp_index.append(temp3_index)
temp_index = temp_index.append(temp4_index)
temp_index = temp_index.append(temp5_index)
temp_index = temp_index.append(temp6_index)
temp_index = temp_index.append(temp7_index)
# delete the repeated index
temp_index = temp_index.drop_duplicates()
# extract the rows
temp = proteomics.loc[temp_index]
# assign the functional with the index in the proteomics
proteomics.loc[temp_index, 'functional group'] = 'Lipid degradation'
# Set GTPase activity like enzymes as other proteins --------------------------
temp = proteomics[(proteomics['protein_name_uniprot'].str.contains('YALI') )&( proteomics['Gene Ontology (GO)'].str.contains('GTPase'))& (proteomics['functional group'] == 'other enzymes')]
proteomics.loc[temp.index, 'functional group'] = 'other proteins'
#proteomics.head()


# ---------------------- Analyze relative proteomics here -----------------------------
# The limited culture is normalized to the bath culture, the protein abundance was refered from fold change
yl_rela_tsv = pd.read_csv('CanLipo_protein_quant_2049.tsv',sep = '\t')

prot_yarli_filtered = yl_rela_tsv.copy()

# Normalization starts here 
# Get the sample columns code
# Sample columns starts with rq_ and ends with ' sum'
sample_columns = [col for col in prot_yarli_filtered.columns if 'rq' in col and 'sum' in col]
prot_yarli_filtered['entry'] = prot_yarli_filtered['Protein Id'].apply(lambda x: x.split('|')[1] if '|' in x else np.nan)   
prot_yarli_filtered
prot_yarli_filtered = prot_yarli_filtered.merge(proteomics, on='entry', how='left', suffixes=('_filtered', '_proteomics'))
print('------- the relative measured accounts for the abundance of the proteomics data : -------')
print(prot_yarli_filtered['Relative abundance'].sum(axis = 0))
prot_yarli_filtered['B_scaled_abundance'] = prot_yarli_filtered['Relative abundance'].div( prot_yarli_filtered['Relative abundance'].sum(axis = 0),axis = 0)  
for sample in sample_columns:
    prot_yarli_filtered[sample + '_readnormtoB'] = [0]*len(prot_yarli_filtered)
    prot_yarli_filtered[sample + '_readnormtoB'] = prot_yarli_filtered[sample].div(prot_yarli_filtered[['rq_131n_sn sum','rq_131c_sn sum','rq_132n_sn sum']].mean(axis = 1))*(prot_yarli_filtered['B_scaled_abundance'])
for sample in sample_columns:
    prot_yarli_filtered[sample + '_frac'] = [0]*len(prot_yarli_filtered)
    prot_yarli_filtered[sample + '_frac'] = prot_yarli_filtered[sample + '_readnormtoB'].div(prot_yarli_filtered[sample + '_readnormtoB'].sum(axis = 0),axis = 0)

prot_yarli_filtered.to_excel('Yarli_result.xlsx', index=False,sheet_name='Yl_proteomics_normalized')    
