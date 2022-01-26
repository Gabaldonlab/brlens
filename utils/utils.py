# -*- coding: utf-8 -*-

'''
utils.py -- Definition of useful basic functions

Written by Moisès Bernabeu <mail@mail.com>
January 2022
'''

import subprocess as sp
import sys
import os


# Running bash commands
def run_cmd(cmd, ommit=False):
    if ommit:
        try:
            process = sp.Popen(cmd, shell=True)
        except:
            pass
        process.communicate("Y\n")
        if process.wait() != 0:
            print("Error occurred, but you chose to ommit it")
    else:
        try:
            process = sp.Popen(cmd, shell=True)
        except OSError:
            sys.exit("Error: Execution cmd failed")
        process.communicate("Y\n")
        if process.wait() != 0:
            sys.exit("ERROR: Execution cmd failed")


def run_cmd_return(cmd):
    try:
        process = sp.Popen(cmd, shell=True, stdout=sp.PIPE).stdout
    except:
        sys.exit()
    return process.readlines()


# Creating folders
def create_folder(name):
    if not os.path.exists(name):
        try:
            os.mkdir(name)
        except Exception:
            print('Unable to create the directory: %s' % name)


def rmkey(dict, key):
    del dict[key]
    return dict


def rmbwbr(string, upper_del, lower_del):
    nsets = string.count(upper_del)
    ostr = ''
    if nsets == 1:
        ostr = string.split(upper_del)[0]
        ostr += string.split(lower_del)[1]
    else:
        parts = string.split(upper_del)
        ostrl = list()
        for i, item in enumerate(parts):
            if not i % 2:
                ostrl.append(item.replace(lower_del, ''))
        ostr = ' '.join(ostrl)
        ostr.replace('  ', ' ').replace('..', '.').replace(',,', ',')

    return ostr


def csv_to_dict(file, sep):
    '''
    It has to be a two columns csv.
    '''
    odict = dict()
    for line in open(file):
        if line != '' and '#' not in line:
            elements = line.replace('\n', '').split(sep)
            odict[elements[0]] = elements[1]

    return odict


def file_exists(filename):
    '''
    Check the existence and size of the file, returns True or False
    '''

    # Check file existence
    if os.path.isfile(filename):
        # Check the file is not empty
        if os.stat(filename).st_size != 0:
            return True
    else:
        return False


def hola():
    print('hola')


# def main():
#     cmd = 'flake8'
#     run_cmd(cmd)
#     miau = run_cmd_return(cmd)
#     print('OUT: %s' % miau)
#
#
# if __name__ == '__main__':
#     main()
