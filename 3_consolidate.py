from Bio import Entrez
import pandas as pd
import requests, sys

mutDict = dict()
patients = set()
SNPIDs = list()

def ifNone(string):
    if string == '':
        return 'NULL'
    else:
        return string

def parseVCF():
    global vcf
    global mutDict
    global patients
    line = vcf.readline()
    line = vcf.readline()
    while line:
        if line.split('\t')[0] in mutDict:
            mutDict[line.split('\t')[0]]['patients']\
                                        .append(line.split('\t')[-1][:-1])
            patients.add(line.split('\t')[-1][:-1])
            mutDict[line.split('\t')[0]]['freq']\
                                        .append(line.split('\t')[14])
            mutDict[line.split('\t')[0]]['depth']\
                                        .append(line.split('\t')[11])
        else:
            mutDict[line.split('\t')[0]] = dict()
            mutDict[line.split('\t')[0]]['patients'] = \
                                        [line.split('\t')[-1][:-1]]
            patients.add(line.split('\t')[-1][:-1])
            mutDict[line.split('\t')[0]]['freq'] = [line.split('\t')[14]]
            mutDict[line.split('\t')[0]]['depth'] = [line.split('\t')[11]]
        line = vcf.readline()

def parseVEP():
    global vep
    global mutDict
    line = vep.readline()
    line = vep.readline()
    while line:
        mutDict[line.split('\t')[0]]['af'] = \
                                    ifNone(line.split('\t')[35])
        mutDict[line.split('\t')[0]]['eur_af'] = \
                                    ifNone(line.split('\t')[39])
        mutDict[line.split('\t')[0]]['gnomad'] = \
                                    ifNone(line.split('\t')[43])
        mutDict[line.split('\t')[0]]['clinsig'] = \
                                    ifNone(line.split('\t')[-14])
        mutDict[line.split('\t')[0]]['conseq'] = \
                                    ifNone(line.split('\t')[9])
        mutDict[line.split('\t')[0]]['symbol'] = \
                                    ifNone(line.split('\t')[11])
        mutDict[line.split('\t')[0]]['rs'] = \
                                    ifNone(line.split('\t')[25].split('&')[0])
        if mutDict[line.split('\t')[0]]['rs'].startswith('rs'):
            SNPIDs.append(mutDict[line.split('\t')[0]]['rs'])
        if 'COSM' in line.split('\t')[25]:
            mutDict[line.split('\t')[0]]['cosm'] = 'COSM' \
                  + line.split('\t')[25].split('COSM')[1].split('&')[0]
        else:
            mutDict[line.split('\t')[0]]['cosm'] = 'NULL'
        mutDict[line.split('\t')[0]]['sift'] = \
                                    ifNone(line.split('\t')[33])
        mutDict[line.split('\t')[0]]['polyphen'] = \
                                    ifNone(line.split('\t')[34])
        mutDict[line.split('\t')[0]]['APPRIS'] = \
                                    ifNone(line.split('\t')[32])
        mutDict[line.split('\t')[0]]['Loftool'] = \
                                    ifNone(line.split('\t')[-6])
        mutDict[line.split('\t')[0]]['MaxEnt'] = \
                                    ifNone(line.split('\t')[-2])
        line = vep.readline()

def parseVAI():
    global vai
    global mutDict
    line = '##'
    while '##' in line:
        line = vai.readline()
    line = vai.readline()
    while line:
        if 'MUTTASTER' in line:
            mutDict[line.split('\t')[0]]['muttaster'] = \
                                line.split('MUTTASTER=')[1].split(';')[0]
        if 'MUTASSESSOR' in line:
            mutDict[line.split('\t')[0]]['mutassessor'] = \
                                line.split('MUTASSESSOR=')[1][0].split(';')[0]
        if 'LRT' in line:
            mutDict[line.split('\t')[0]]['LRT'] = \
                                line.split('LRT=')[1][0].split(';')[0]
        line = vai.readline()

def createHeader(out):
    global patients
    line = 'key\trs\tcosm\tsymbol\tconseq\taf\taur_af\tgnomad\tclinsig\tsift\t'\
           + 'polyphen\tAPPRIS\tLoftool\tMaxEnt\tMutTaster\tMutAssessor\tLRT\t'
    patients = list(patients)
    patients.sort()
    for value in patients:
        line += value + 'freq' + '\t' + value + 'depth' + '\t'
    out.write(line[:-1] + '\n')

def isPolymorphism(af, eur_af, gnomad):
    b = False
    try:
        if float(af) > 0.05:
            b = True
    except ValueError:
        pass
    try:
        if float(eur_af) > 0.05:
            b = True
    except ValueError:
        pass
    try:
        if float(gnomad) > 0.05:
            b = True
    except ValueError:
        pass
    return b

def addRow(key, value):
    global patients
    toWrite = key + '\t' + value['rs'] \
                  + '\t' + value['cosm'] \
                  + '\t' + value['symbol'] \
                  + '\t' + value['conseq'] \
                  + '\t' + value['af'] \
                  + '\t' + value['eur_af'] \
                  + '\t' + value['gnomad'] \
                  + '\t' + value['clinsig'] \
                  + '\t' + value['sift'] \
                  + '\t' + value['polyphen'] \
                  + '\t' + value['APPRIS'] \
                  + '\t' + value['Loftool'] \
                  + '\t' + value['MaxEnt'] \
                  + '\t' + value.get('muttaster', 'NULL') \
                  + '\t' + value.get('mutassessor', 'NULL') \
                  + '\t' + value.get('LRT', 'NULL')
    for p in patients:
        if p in value['patients']:
            toWrite += '\t' + value['freq'][value['patients'].index(p)] \
                     + '\t' + value['depth'][value['patients'].index(p)]
        else:
            toWrite += '\t'*2
    return toWrite

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


def createExcel():
    everything = pd.read_csv('./results/all.txt', \
                         sep = '\t', \
                         # keep_default_na = False, \
                         index_col = 0)
    tier1 = pd.read_csv('./results/tier1.txt', \
                         sep = '\t', \
                         # keep_default_na = False, \
                         index_col = 0, \
                         # header = 0, \
                         names = ['key'] + list(everything.columns))
    tier2 = pd.read_csv('./results/tier2.txt', \
                         sep = '\t', \
                         # keep_default_na = False, \
                         index_col = 0, \
                         # header = 0, \
                         names = ['key'] + list(everything.columns))
    tier3 = pd.read_csv('./results/tier3.txt', \
                         sep = '\t', \
                         # keep_default_na = False, \
                         index_col = 0, \
                         # header = 0, \
                         names = ['key'] + list(everything.columns))
    tier4 = pd.read_csv('./results/tier4.txt', \
                         sep = '\t', \
                         # keep_default_na = False, \
                         index_col = 0, \
                         # header = 0, \
                         names = ['key'] + list(everything.columns))
    poly = pd.read_csv('./results/poly.txt', \
                         sep = '\t', \
                         # keep_default_na = False, \
                         index_col = 0, \
                         # header = 0, \
                         names = ['key'] + list(everything.columns))
    writer = pd.ExcelWriter('./results/tiers.xlsx')
    everything.to_excel(writer, 'all')
    tier1.to_excel(writer, 'I+II')
    tier2.to_excel(writer, 'II')
    tier3.to_excel(writer, 'III')
    tier4.to_excel(writer, 'IV')
    poly.to_excel(writer, 'Pol')
    writer.save()

def mutationAssessor(addition, mut):
    server = "http://mutationassessor.org/r3/?cm=var&var="
    ext = "hg38," + mut.replace('chr', '').replace('_', ',').replace('/', ',')
    r = requests.get(server+ext,
                     headers = {"Content-Type" : "application/json"},
                     params = {"frm" : "json", "fts" : "F_impact"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    decoded = r.json()
    print(mut + '\t' + decoded[0]['F_impact'])
    addition.write(mut + '\t' + repr(decoded[0]['F_impact']))

with \
  open('tableVEP.txt', 'r') as vep,\
  open('allVCF.vcf', 'r') as vcf,\
  open('fromVAI.txt', 'r') as vai,\
  open('./results/all.txt', 'w') as out,\
  open('./results/tier1.txt', 'w') as tier1,\
  open('./results/tier2.txt', 'w') as tier2,\
  open('./results/tier3.txt', 'w') as tier3,\
  open('./results/tier4.txt', 'w') as tier4,\
  open('./results/poly.txt', 'w') as poly:
    parseVCF()
    parseVEP()
    parseVAI()
    createHeader(out)
    for key, value in mutDict.items():
        toWrite = addRow(key, value) + '\n'
        out.write(toWrite)
        if isPolymorphism(value['af'], value['eur_af'], value['gnomad']):
            poly.write(toWrite)
        elif ('pathogen' in value['clinsig']) \
          or ('A' in value.get('muttaster', '0')):
            tier1.write(toWrite)
        elif ('benign' in value['clinsig']) \
          or ('synon' in value['conseq']) \
          or ('P' in value.get('muttaster', '0')):
            tier4.write(toWrite)
        elif ('delet' in value['sift']) \
          or ('damag' in value['polyphen']) \
          or ('D' in value.get('muttaster', '0')) \
          or ('D' in value.get('LRT', '0')) \
          or ('H' in value.get('mutassessor', '0')) \
          or ('M' in value.get('mutassessor', '0')):
            tier2.write(toWrite)
        elif (value['Loftool'] != 'NULL'):
          # or ('' in value['MaxEnt'])
          # or ('' in value['APPRIS'])
            if float(value['Loftool']) < 0.5:
                tier2.write(toWrite)
            else:
                tier3.write(toWrite)
        else:
            tier3.write(toWrite)
createExcel()
