'''
snp = open('snp_result.txt', 'r')
minorAFs = open('minorAF.txt', 'w')
line = snp.readline()
while line:
    if '[Homo sapiens]' in line:
        rs = line.split(' ')[1]
        mut = line = snp.readline()
    elif '=' in line:
        fr = float(line.split('=')[1].split('/')[0])
        if fr > 0.05:
            minorAFs.write(rs + '\t' + mut[:-1] + '\t' + line)
    line = snp.readline()
snp.close()
minorAFs.close()
'''

'''
import requests, sys

server = "https://rest.ensembl.org"
ext = "/variation/human/COSM13587?"
# "/overlap/id/ENSG00000157764?feature=somatic_variation"
# "/variation/human/rs34975585?"
 
r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})


if not r.ok:
  r.raise_for_status()
  sys.exit()
 
decoded = r.json()
with open('rar.txt', 'w') as file:
  file.write(repr(decoded))
'''

from Bio import Entrez
Entrez.email = "bug.dmitrij@gmail.com"
response = Entrez.efetch(db = "SNP", id = '201969638')
for i in response:
    print(i)
