#!/usr/bin/env python3

import sys
import os
import argparse


def parseArguments():
    parser = argparse.ArgumentParser(description='transform file and header')
    parser.add_argument('--info_listfrs',type=str,help="input file",required=True)
    parser.add_argument('--sig_info',type=str,help="input file",required=True)
    args = parser.parse_args()
    return args


