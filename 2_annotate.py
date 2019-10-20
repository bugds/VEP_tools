#!/usr/bin/python3

import pandas as pd
import sys
import datetime
import os

mutDF = dict()
patients = set()
SNPIDs = list()

suffix = 'TV'
VEPfile='VEP_results/VEP.' + suffix + '.vcf'
VCFfile='population_VCF/allVCF.' + suffix + '.vcf'
VAIfile='VAI_results/fromVAI.' + suffix + '.txt'

def fromVEP(VEPfile=VEPfile,
            output='tableVEP.txt'):
    '''Reformation of VEP output
    '''
    with open(VEPfile, 'r') as f, open(output, 'w') as o:
        line = '##'
        while line:
            if 'Format:' in line:
                # Parsing the INFO header, lacking tabs
                columns = line.count('|') + 1
                info = line.split('Format: ')[1][:-3].replace('|', '\t')
            if not '##' in line:
                if '#' in line:
                    # Writing new header with tabs
                    columns += line.count('\t')
                    o.write('IDENT' + '\t' + line[1:].replace('INFO', info))
                else:
                    toWrite = '_'.join(line.split('\t')[0:2]) + '_' \
                              + '/'.join(line.split('\t')[3:5]) + '\t' \
                              + line.replace('|', '\t')
                    o.write('\t'.join(toWrite.split('\t')[:columns]) + '\n')
            line = f.readline()

def ifNone(string):
    '''Unification of NA data
    '''
    if string == '':
        return 'NULL'
    else:
        return string

def parseInput(file, VAI=False):
    '''Parsing VCF, VEP, VAI files

    file - path to file
    VAI - indicator of VAI file to skip 17 rows

    Returns pandas dataframe
    '''
    if VAI:
        rowsToSkip = 17
    else:
        rowsToSkip = 0
    return pd.read_csv(file,
                       keep_default_na = False,
                       # index_col = 0, # r-g (reformatting)
                       sep = '\t',
                       skiprows = rowsToSkip,
                       dtype={'PATIENT': str})

def createDataBase(VCF, VEP, VAI):
    '''Merging dataframes from VCF and VEP and VAI outputs
    '''
    # patients = VCF['PATIENT'].unique() # r-g
    dataBase = pd.merge(VEP[['IDENT',
                             'AF', # 1KGenomes allele frequency
                             'EUR_AF', # European population of 1KGenomes
                             'gnomAD_AF', # GnomAD allele frequency
                             'CLIN_SIG', # Clinical significance
                             'Consequence', # Biological consequence
                             'SYMBOL', # Gene
                             'Existing_variation', # Known identificators
                             'MaxEntScan_diff', # Splicing predictor 2
                             'SIFT',
                             'PolyPhen'
                           ]],
                        VCF[['IDENT',
                             'PATIENT', # Probe ID
                             'FREQ', # Alt reads fraction
                             'DP', # Read depth
                             'FILTER' # VCF filter
                           ]],
                        on = 'IDENT', sort = 'False')
                        # right_index=True, left_index=True, sort='False') # r-g
    
    dataBase = \
        dataBase.rename(columns={'SIFT': 'SIFT_VEP', 'PolyPhen': 'PolyPhen_VEP'})
    
    infoVAI = ['VEST',
               'SIFT',
               'PP2HVAR',
               'PP2HDIV',
               'MUTTASTER',
               'MUTASSESSOR',
               'LRT']
    
    for info in infoVAI:
        dataBase[info] = ''

    for index, row in VAI.iterrows():
        for info in infoVAI:
            if info in row['Extra']:
                dataBase.loc[dataBase['IDENT'] == \
                             VAI.loc[index, 'Uploaded Variation'], info] = \
                    row['Extra'].split(info + '=')[1].split(';')[0]

    return dataBase

def reindexWithSamples(dataBase):
    for index, row in dataBase.iterrows():
        dataBase.loc[index, 'IDENT'] += '_' + str(row['PATIENT'])
    dataBase.set_index('IDENT', inplace=True, drop=True)
    return dataBase

def isPolymorphism(af, eur_af, gnomad):
    b = False
    percentage = 0.01
    try:
        if float(af) > percentage:
            b = True
    except ValueError:
        pass
    try:
        if float(eur_af) > percentage:
            b = True
    except ValueError:
        pass
    try:
        if float(gnomad) > percentage:
            b = True
    except ValueError:
        pass
    return b

def checkMinorAF():
    global SNPIDs
    global l2
    minorAFCandidates = list()
    Entrez.email = "bug.dmitrij@gmail.com"
    response = Entrez.efetch(db = "SNP",
                             id = ','.join(SNPIDs) \
                                     .replace('rs', ''))
    for SNPStr in response:
        if '[Homo sapiens]' in SNPStr:
            rs = SNPStr.split(' ')[1]
        elif '=' in SNPStr:
            if float(SNPStr.split('=')[1].split('/')[0]) > 0.05:
                minorAFCandidates.append(rs)
    return minorAFCandidates

fromVEP()
VCF = parseInput(VCFfile)
VEP = parseInput('tableVEP.txt')
#os.remove('tableVEP.txt')
VAI = parseInput(VAIfile, VAI=True)
dataBase = createDataBase(VCF,VEP,VAI)
dataBase['tier'] = ''
reindexWithSamples(dataBase)
for index, row in dataBase.iterrows():
    if isPolymorphism(row['AF'], row['EUR_AF'], row['gnomAD_AF']):
        dataBase.loc[index, 'tier'] = '5' 
    elif ('pathogen' in row['CLIN_SIG']) \
      or ('A' in row['MUTTASTER']):
        dataBase.loc[index, 'tier'] = '1'
    elif ('benign' in row['CLIN_SIG']) \
      or ('synon' in row['Consequence']) \
	  or ('P' in row['MUTTASTER']):
        dataBase.loc[index, 'tier'] = '4'
    elif ('D' in row['SIFT']) \
      or ('D' in row['PP2HVAR']) \
      or ('D' in row['PP2HDIV']) \
      or ('D' in row['MUTTASTER']) \
      or ('D' in row['LRT']) \
      or ('H' in row['MUTASSESSOR']) \
      or ('M' in row['MUTASSESSOR']) \
      or ('damag' in row['PolyPhen_VEP']) \
      or ('delet' in row['SIFT_VEP']):
        dataBase.loc[index, 'tier'] = '2'
    else:
        dataBase.loc[index, 'tier'] = '3'
writer = pd.ExcelWriter('./results/' \
                        + datetime.datetime.now().strftime("%Y_%m_%d") \
                        + '.' + suffix + '.xlsx')
dataBase.to_excel(writer, suffix)
writer.save()
