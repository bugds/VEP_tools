#!/usr/bin/python3

import os

varSet = set()

def main():
    global suffix
    linesSet = set()
    VEP = open('VEP.' + suffix + '.vcf', 'r')
    newVEP = open('newVEP.' + suffix + '.vcf', 'w')
    line = '##'
    while line:
        if line in linesSet:
            pass
        else:
            linesSet.add(line)
            newVEP.write(line)
        line = VEP.readline()
    VEP.close()
    newVEP.close()

s = ['V', 'TV', 'G', 'TG', 'P', 'TP']

for suffix in s:
    main()
main()
