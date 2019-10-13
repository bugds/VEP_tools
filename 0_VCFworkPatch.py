#!/usr/bin/python3

import os

varSet = set()

def main():
    global suffix
    outputVCF = open('population_VCF/allVCF.' + suffix + '.vcf', 'r')
    newOutputVCF = open('population_VCF/allVCF.' + suffix + '.backspace.vcf', 'w')
    line = outputVCF.readline()
    while line:
        if ',' in line.split('\t')[0]:
            for altAllele in line.split('\t')[0].split('/')[1].split(','):
                toWrite = line.split('\t')[0].split('/')[0] + '/' + altAllele + '\t'
                toWrite += '\t'.join(line.split('\t')[1:5]) + '\t' + altAllele + '\t'
                toWrite += '\t'.join(line.split('\t')[6:8]) + '\t0/1\t'
                index = line.split('\t')[0].split('/')[1].split(',').index(altAllele)
                toWrite += line.split('\t')[9].split(',')[0] + ',' \
                           + line.split('\t')[9].split(',')[index] + '\t'
                toWrite += line.split('\t')[10].split(',')[index-1] + '\t'
                toWrite += line.split('\t')[11] + '\t' + line.split('\t')[12]
                newOutputVCF.write(toWrite)
        else:
            newOutputVCF.write(line)
        line = outputVCF.readline()
    outputVCF.close()
    newOutputVCF.close()
    
    newForVEP = open('VEP_results/forVEP.' + suffix + '.vcf', 'w')
    newOutputVCF = open('population_VCF/allVCF.' + suffix + '.backspace.vcf', 'r')
    line = newOutputVCF.readline()
    line = newOutputVCF.readline()
    while line:
        newForVEP.write('\t'.join(line.split('\t')[1:6]) + '\t' + '.' + '\n')
        line = newOutputVCF.readline()
    newOutputVCF.close()
    newForVEP.close()

# REPEATING LINES FOR VEP!!! MAKE A SET AND READ FROM IT!!!

s = ['V', 'TV', 'G', 'TG', 'P', 'TP']

for suffix in s:
    main()
main()
