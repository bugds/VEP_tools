import requests, sys
import json
import pandas as pd
import re

# getting GRCh37 coordinates
'''
with open('all.txt', 'r') as variants, open('../clinDB/hg38.txt', 'w') as hg38:
    line = variants.readline()
    line = variants.readline()
    while line:
        ident = line.split('\t')[0]
        if not 'alt' in ident:
            hg38.write(ident.split('chr')[1].split('_')[0] + '\t' \
                       + ident.split('_')[1] + '\n')
        line = variants.readline()
'''

def readHg37():
    hg37Dict = {}
    line = hg.readline()
    while line:
        hg37Dict[line.split('\t')[0] + '_' + line.split('\t')[2]] = \
                                   line.split('\t')[1]
        line = hg.readline()
    return hg37Dict

def readHg38():
    hg38List = []
    with open('hg38.txt', 'r') as hg38:
        line = hg38.readline()
        while line:
            hg38List.append('_'.join(line.split('\t'))[:-1])
            line = hg38.readline()
    return hg38List
'''
def getHg37(chrom, posit):
    global hg37Dict
    key = chrom + '_' + posit
    if key in hgDict.keys():
        return hgDict[key]
    else:
        return 'NULL'
'''
def fromICGC(varICGC):
    global hg37Dict
    global clinData
    clinData.write('ICGC:\n')
    for hit in varICGC:
        chrom = str(hit['chromosome'])
        posit = str(hit['start'])
        if (chrom + '_' + posit) in hg37Dict.keys():
            clinData.write(chrom + '_' \
                           + hg37Dict[chrom + '_' + posit] \
                           + ': ' + hit['id'] + '\n')

def fromClinvar(varClinvar):
    global hg38List
    global clinData
    clinData.write('clinvar:\n')
    for index, row in varClinvar.iterrows():
        chrom = str(row['#CHROM'])
        posit = str(row['POS'])
        if (chrom + '_' + posit) in hg38List:
                clinData.write(chrom + '_' + posit \
                               + ': ' + row['ID'] + '\n')

def fromCGI(varCGI):
    global hg37Dict
    global clinData
    clinData.write('CGI:\n')
    for index, row in varCGI.iterrows():
        chrom = row['gdna'].split('chr')[1].split(':')[0]
        dirtyPos = row['gdna'].split(':g.')[1]
        posit = re.findall(r'\d+', dirtyPos)[0]
        if (chrom + '_' + posit) in hg37Dict.keys():
            clinData.write(chrom + '_' \
                           + hg37Dict[chrom + '_' + posit] \
                           + ': ' + row['reference'] + '\n')

def fromDoCM(varDoCM):
    global hg37Dict
    global clinData
    clinData.write('DoCM:\n')
    for index, row in varDoCM.iterrows():
        chrom = str(row['#CHROM'])
        posit = str(row['POS'])
        if (chrom + '_' + posit) in hg37Dict.keys():
            clinData.write(chrom + '_' \
                           + hg37Dict[chrom + '_' + posit] \
                           + ': ' + row['INFO'] + '\n')

with \
  open('hg37.txt', 'r') as hg, \
  open('CGI.tsv', 'r') as CGI, \
  open('DoCM.vcf', 'r') as DoCM, \
  open('ICGC.json', 'r') as ICGC, \
  open('OncoKB.txt', 'r') as OncoKB, \
  open('clinvar.vcf', 'r') as clinvar, \
  open('clinData.txt', 'w') as clinData, \
  open('../tableVEP.txt', 'r') as variants:
    hg37Dict = readHg37()
    hg38List = readHg38()
    print(hg38List)
    varICGC = json.loads(ICGC.read())['hits']
    varClinvar = pd.read_csv('clinvar.vcf', skiprows = 27, sep = '\t', \
                             dtype = 'str')
    varCGI = pd.read_csv('CGI.tsv', sep = '\t', dtype = 'str')
    varDoCM = pd.read_csv('DoCM.vcf', sep = '\t', skiprows = 6, dtype = 'str')
    fromICGC(varICGC)
    fromClinvar(varClinvar)
    fromCGI(varCGI)
    fromDoCM(varDoCM)
