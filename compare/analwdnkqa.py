import pandas as pd

table1 = pd.read_excel('./other_23.xlsx')
table2 = pd.read_excel('./other_24.xlsx')

table1['Chromosome'] = table1['Chromosome'].astype(str)
table1['Position'] = table1['Position'].astype(str)
table1['Ref'] = table1['Ref'].astype(str)
table1['Alt'] = table1['Alt'].astype(str)
table2['Chromosome'] = table2['Chromosome'].astype(str)
table2['Position'] = table2['Position'].astype(str)
table2['Ref'] = table2['Ref'].astype(str)
table2['Alt'] = table2['Alt'].astype(str)

table = table1.merge(table2, on=['Gene', 'Chromosome', 'rsID', \
                                 'Position', 'Ref', 'Alt',\
                                 'EffType', 'IVS', 'OMIM',\
                                 'ClinVar', 'PROVEAN', 'SIFT',\
                                 'Polyphen2b', 'X1000G_AF',\
                                 'ExAC_AF', 'ESP6500_AF',\
                                 'our_AF', 'HGVS_name'],\
                     how = 'outer')
print(table)
writer = pd.ExcelWriter('./table.xlsx')
table.to_excel(writer, 'all')
writer.save()
