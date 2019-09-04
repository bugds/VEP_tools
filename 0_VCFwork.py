#!/usr/bin/python3

import os

varSet = set()

def writeHeader(outputVCF, callerName, line):
    if callerName == 'Mutect2':
        specificHeaderPart = 'GT:AD:FREQ:DP'
    elif callerName == 'VarScan2':
        specificHeaderPart = \
            'GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR'
    elif callerName == 'Pisces':
        specificHeaderPart = 'GT:GQ:AD:DP:FREQ:NL:SB'
    outputVCF.write('IDENT\t'  \
                   + line[1:].replace(line.split('\t')[-1], 
                                      specificHeaderPart.replace(':', '\t') \
                                      + '\tPATIENT')
                             .replace('INFO\tFORMAT\t', 
                                      '') \
                   + '\n')

def getString(callerName, line, patient):
    if callerName == 'Mutect2':
        specificVariantPart = '\t'.join(line.split('\t')[9]\
                                            .split(':')[0:4])[:-1]
    else:
        specificVariantPart = '\t'.join(line.split('\t')[9]\
                                            .split(':'))[:-1]
    return '_'.join(line.split('\t')[0:2]) + '_' \
         + '/'.join(line.split('\t')[3:5]) + '\t' \
         + '\t'.join(line.split('\t')[0:7]) + '\t' \
         + specificVariantPart + '\t' + patient + '\n'

def main():
    firstHeaderFlag = True
    for f in os.listdir('./raw'):
        singleVCF = open('./raw/' + f, 'r')
        patient = f[3:5]
        line = '##'
        while line:
            if not ('##' in line):
                if ('#CHROM' in line) and firstHeaderFlag:
                    outputVCF = open('allVCF.' + suffix + '.vcf', 'w')
                    writeHeader(outputVCF, callerName, line)
                    firstHeaderFlag = False
                elif not ('#CHROM' in line):
                    string = getString(callerName, line, patient)
                    outputVCF.write(string)
                    varSet.add('\t'.join(string.split('\t')[1:6]) \
                             + '\t' + '.')
            elif ('##source=Mutect2' in line):
                callerName = 'Mutect2'
                suffix = 'M'
            elif ('##source=VarScan2' in line):
                callerName = 'VarScan2'
                suffix = 'V'
            elif ('##source=Pisces' in line):
                callerName = 'Pisces'
                suffix = 'TP'
            line = singleVCF.readline()
        singleVCF.close()

    outputVCF.close()
    fileForVEP = open('forVEP.' + suffix + '.vcf', 'w')
    for item in varSet:
        fileForVEP.write(item + '\n')
    fileForVEP.close()

main()