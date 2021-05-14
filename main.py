from collections import OrderedDict


def parseSequencesToDict(lines):
  keys = []
  seqs = []
  for line in lines:
    if line[0] == '>':
      key = line[1:].strip()
      keys.append(key)
    else:
      sequence = line.strip()
      seqs.append(sequence)
  return dict(zip(keys, seqs))

def writeDictToFile(d, filename):
  f = open(filename, 'w')
  for key in d:
    f.write(key+'\n')
    f.write(d[key]+'\n')
  f.close()
  
filename = "pre_sorted_seqs.fasta"

with open(filename) as f:
  lines = f.readlines()

dict1 = parseSequencesToDict(lines)
dict1 = OrderedDict(sorted(dict1.items()))

for key in dict1:
  print(f'{key}: {dict1[key]}')

writeDictToFile(dict1, 'output.txt')