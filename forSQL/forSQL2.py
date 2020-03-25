#!/usr/bin/python3

import os
import datetime

minVarCharLength = 50
naValue = '.'
ifUnknown = 'varchar({})'.format(minVarCharLength)

def valiDate(date_text):
    datetime.datetime.strptime(date_text, '%d.%m.%Y')
    return 'date'

def checkIf(strLine, valueType):
    exec(valueType + '(row[i])')
    return valueType                    

naList = ['NA', 'NULL', '-', '', '.']

for filename in os.listdir('./raw/'):
    file = open('./raw/' + filename, 'r')
    output = open('./results/sql_' + filename, 'w')
    line = 'insert into t' + filename.split('.')[0].capitalize() + '('
    line += '\t' + ',\t'.join(file.readline()[:-1].split('\t')) + ') values\n'
    output.write(line)
    line = file.readline()[:-1]
    typesArray = ['unknown']*len(line.split('\t'))
    while line:
        row = line.split('\t')
        # If len(row) < len(typesArray), all cells after it == 'NULL',
        # and need to be passed anyway
        for i in range(min(len(typesArray), len(row))):
            if row[i] in naList:
                pass
            elif typesArray[i] == 'unknown':
                try:
                    typesArray[i] = checkIf(row[i], 'int')
                except ValueError:
                    try:
                        typesArray[i] = checkIf(\
                            row[i].replace('%', '').replace(',', '.'), 'float')
                    except ValueError:
                        try:
                            typesArray[i] = valiDate(row[i])
                        except ValueError:
                                typesArray[i] = 'varchar({})'.format(len(row[i]))
            elif typesArray[i] == 'int':
                try:
                    typesArray[i] = checkIf(row[i], 'int')
                except ValueError:
                    try:
                        typesArray[i] = checkIf(\
                            row[i].replace('%', '').replace(',', '.'), 'float')
                    except ValueError:
                        typesArray[i] = 'varchar({})'.format(len(row[i]))
            elif typesArray[i] == 'float':
                try:
                    typesArray[i] = checkIf(\
                            row[i].replace('%', '').replace(',', '.'), 'float')
                except ValueError:
                    typesArray[i] = 'varchar({})'.format(len(row[i]))
            elif 'varchar' in typesArray[i]:
                length = int(typesArray[i].split('(')[1].split(')')[0])
                typesArray[i] = 'varchar({})'.format(max(length, len(row[i])))
        line = file.readline()[:-1]
    toWrite = '\t' + '\t'.join(['varchar(3)' if t == 'date' else t for t in typesArray]) + '\n'
    output.write(toWrite.replace('unknown', ifUnknown))
    file.close()
    file = open('./raw/' + filename, 'r')
    line = file.readline()
    line = file.readline()[:-1]
    while line:
        nextLine = file.readline()[:-1]
        c = ',\n'
        if nextLine == '':
            c = nextLine
        row = line.split('\t')
        toWrite = '(\t'
        if len(row) < len(typesArray):
            row.extend([naValue]*(len(typesArray) - len(row)))
        for i in range(len(typesArray)):
            if row[i] in naList:
                toWrite += naValue + ',\t'
            elif (typesArray[i] == 'int') or (typesArray[i] == 'float'):
                toWrite += row[i].replace('%', '').replace(',', '.') + ',\t'
            elif typesArray[i] == 'date':
                toWrite += 'convert(datetime,\'' + row[i] + '\',104),\t'
            else:
                toWrite += '\'' + row[i] + '\',\t'
        output.write(toWrite[:-2] + ')' + c)
        line = nextLine
    file.close()
    output.close()
