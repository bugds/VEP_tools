with open('allVEP.vcf', 'r') as f, open('tableVEP.txt', 'w') as o:
    line = '##'
    while line:
        if 'Format:' in line:
            info = line.split('Format: ')[1][:-3].replace('|', '\t')
        if not '##' in line:
            if '#' in line:
                o.write('IDENT' + '\t' + line[1:].replace('INFO', info))
            else:
                o.write('_'.join(line.split('\t')[0:2]) + '_' + '/'.join(line.split('\t')[3:5]) + '\t' + line.replace('|', '\t'))
        line = f.readline()
