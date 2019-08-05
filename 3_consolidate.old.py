from Bio import Entrez
# import xlsxwriter

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

def createHeader(out):
    global patients
    line = 'key\trs\tcosm\tsymbol\tconseq\taf\taur_af\tgnomad\tclinsig\tsift\t'\
           + 'polyphen\tAPPRIS\tLoftool\tMaxEnt\t'
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
                  + '\t' + value['MaxEnt']
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

'''
def createExcel():
    global mutDict
    global patients
    workbook = xlsxwriter.Workbook('./results/common.xlsx')
    sheets = {"All" : {}, \
              "Polymorphisms" : {}, \
              "Minor polymorphisms" : {}, \
              "Benign and synonymous" : {}, \
              "Unknown" : {}, \
              "Pathogenic" : {}, \
              "Predicted pathogenic" : {}}

    header = ["KEY", "RS", "COSM", "SYMBOL", "CONSEQ", "AF", "EUR_AF", \
        "GNOMAD", "CLINSIG", "SIFT", "POLYPHEN", "APPRIS", "LOFTOOL", "MAXENT"]

    for sheet in sheets.values():
        sheet['wsh'] = workbook.add_worksheet()
        sheet['row'] = 0
        for category in header:
            sheet['wsh'].write(counter, 0, category)
            counter += 1

    for key, value in mutDict.items():
        col = 0
        fillRow(sheets["All"])
        if isPolymorphism(value['af'], value['eur_af'], value['gnomad']):
            sheets[""]
        elif value['rs'] in minorAFCandidates:
            l2.write(toWrite)
        elif 'pathogen' in value['clinsig']:
            l5.write(toWrite)
        elif ('benign' in value['clinsig']) \
          or ('synon' in value['conseq']):
            l3.write(toWrite)
        elif ('delet' in value['sift']) \
          or ('damag' in value['polyphen']) \
          or (value['Loftool'] != 'NULL'):
          # or ('' in value['MaxEnt'])
          # or ('' in value['APPRIS'])
            try:
                if float(value['Loftool']) < 0.5:
                    l6.write(toWrite)
                else:
                    pass
            except ValueError:
                l6.write(toWrite)
        else:
            l4.write(toWrite)

        sheets["Tier I"]['wsh'].write(sheets["Tier I"]['row'], 0, key) 
        sheets["Tier I"]['wsh'].write(sheets["Tier I"]['row'], 1, value['rs'])
        sheets["Tier I"]['wsh'].write(sheets["Tier I"]['row'], 0, key) 
        sheets["Tier I"]['wsh'].write(sheets["Tier I"]['row'], 1, value['rs'])
        sheets["Tier I"]['row'] += 1

    workbook.close() 

def fillRow(sheet):
    header = ["KEY", "RS", "COSM", "SYMBOL", "CONSEQ", "AF", "EUR_AF", \
        "GNOMAD", "CLINSIG", "SIFT", "POLYPHEN", "APPRIS", "LOFTOOL", "MAXENT"]
    for 
'''

with \
  open('tableVEP.txt', 'r') as vep,\
  open('allVCF.vcf', 'r') as vcf,\
  open('./results/all.txt', 'w') as out,\
  open('./results/polymorphisms.txt', 'w') as l1,\
  open('./results/minor_AFdb.txt', 'w') as l2,\
  open('./results/ben_syn.txt', 'w') as l3,\
  open('./results/unknown.txt', 'w') as l4,\
  open('./results/pat.txt', 'w') as l5,\
  open('./results/pat_predict.txt', 'w') as l6:
    parseVCF()
    parseVEP()
    createHeader(out)
    minorAFCandidates = checkMinorAF()
    for key, value in mutDict.items():
        toWrite = addRow(key, value) + '\n'
        out.write(toWrite)
        if isPolymorphism(value['af'], value['eur_af'], value['gnomad']):
            l1.write(toWrite)
        elif value['rs'] in minorAFCandidates:
            l2.write(toWrite)
        elif 'pathogen' in value['clinsig']:
            l5.write(toWrite)
        elif ('benign' in value['clinsig']) \
          or ('synon' in value['conseq']):
            l3.write(toWrite)
        elif ('delet' in value['sift']) \
          or ('damag' in value['polyphen']) \
          or (value['Loftool'] != 'NULL'):
          # or ('' in value['MaxEnt'])
          # or ('' in value['APPRIS'])
            try:
                if float(value['Loftool']) < 0.5:
                    l6.write(toWrite)
                else:
                    pass
            except ValueError:
                l6.write(toWrite)
        else:
            l4.write(toWrite)
    # createExcel()
