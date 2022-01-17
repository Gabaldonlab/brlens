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
