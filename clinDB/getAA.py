with open('clinData.txt', 'r') as clD, open('../results/all.txt', 'r') as res:
    clData = dict()
    line = clD.readline()
    while line:
        if '_' in line:
            clData[line.split(':')[0]] = line.split(':')[1][:-1]
        line = clD.readline()
    resData = dict()
    line = res.readline()
    line = res.readline()
    while line:
        if '_' in line:
            resData['_'.join(line.split('chr')[1].split('_')[:2])] = \
                        line.split('\t')[0]
        line = res.readline()
    for key in clData.keys():
        print(resData[key] + ':' + clData[key])
