#!/usr/bin/python3
import os

def add_home_command(subparsers):
    subparser = subparsers.add_parser("home",
        help="Display OGSync Homepage info and exit",
        usage='OGSync.py home',
        description="Open the OGSync Home Page and exit")
    subparser.set_defaults(callback=run_home_command)
    return subparser

def run_home_command(args):
    print('*************************************************')
    print('Please visit https://dvb.ac.cn/OGSync for detail.')
    print('*************************************************')