#!/usr/bin/env python3

import sys
from utils import create_folder
from optparse import OptionParser


def main():
    parser = OptionParser()
    parser.add_option('-i', '--input', dest='ifile',
                      help=('File with multiple trees in newick, each line has'
                            ' the format: seed\tmodel\tlikelihood\tnewick'),
                      metavar='<file.nwk>')
    parser.add_option('-p', '--prefix', dest='prefix',
                      help='Output prefix.',
                      metavar='</path/to/dir/prefix> or <prefix>',
                      default='dist_output')
    (options, args) = parser.parse_args()

    # Creating the output directory
    odir = '../data/%s' % options.prefix
    create_folder(odir)

    file = options.ifile

    print('Processing: %s' % file)

    ph_id = file.rsplit('/', 1)[1].split('_', 1)[0]

    fs = sys.getsizeof(open(file, 'r').read())
    opts = 500000

    if fs > opts:
        efno = fs / opts
        fmaxsize = efno / int(efno) * opts + 100

        olines = ''
        ofileno = 0

        for line in open(file, 'r'):
            if sys.getsizeof(olines) < fmaxsize:
                olines += line  # .decode('utf-8')
            else:
                print(ofileno, 'written')
                ofile = open('%s/%s_%s.txt' % (odir, ph_id, ofileno), 'w')
                ofile.write(olines)
                ofile.close()
                olines = line
                ofileno += 1

        ofile = open('%s/%s_%s.txt' % (odir, ph_id, ofileno), 'w')
        ofile.write(olines)
        ofile.close()
    else:
        print('File is not greater than 500 KB.')


if __name__ == '__main__':
    main()
