import os

varSet = set()

output = open('allVCF.vcf', 'w')
bool = True

for f in os.listdir('./raw'):
        file = open('./raw/' + f, 'r')
        patient = f[-6:-4]
        partOfHeader = 'GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR'.replace(':', '\t')
        line = '##'
        while line:
                if not ('##' in line):
                        if '#CHROM' in line:
                                if bool:
                                        output.write('IDENT\t' + line[1:].replace('Sample1', partOfHeader + '\tPATIENT').replace('INFO\tFORMAT\t', ''))
                                        bool = False
                        else:
                                string = \
                                        '_'.join(line.split('\t')[0:2]) + '_' + '/'.join(line.split('\t')[3:5]) + \
                                        '\t' + '\t'.join(line.split('\t')[0:7]) + '\t' + '\t'.join(line.split('\t')[9].split(':'))[:-1] + '\t' + patient + '\n'
                                output.write(string)
                                varSet.add('\t'.join(string.split('\t')[1:7]))
                line = file.readline()
        file.close()

output.close()
file = open('forVEP.vcf', 'w')
for item in varSet:
        file.write(item + '\n')
file.close()
