import sys
import os
from glob import glob

if not os.path.exists('splitted'):
    os.mkdir('splitted')

files = glob('outputs/*_best_trees.txt.gz')

for file in files:
    print(file)
    ph_id = file.rsplit('/', 1)[1].split('_', 1)[0]
    print(ph_id)
    fs = sys.getsizeof(open(file, 'r').read())
    opts = 500000
    if fs > opts:
        efno = fs / opts
        fmaxsize = efno / int(efno) * opts + 100

    olines = ''
    ofileno = 0
    for line in open(file, 'r'):
        print(ofileno)
        if sys.getsizeof(olines) < fmaxsize:
            olines += line
        else:
            ofile = open('splitted/%s_%s.txt' % (ph_id, ofileno), 'w')
            ofile.write(olines)
            ofile.close()
            olines = ''
            ofileno += 1

    ofile = open('splitted/%s_%s.txt' % (ph_id, ofileno), 'w')
    ofile.write(olines)
    ofile.close()
