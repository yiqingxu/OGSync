#!/usr/bin/python3

import argparse, os, datetime, csv, json, io, time, hashlib, difflib

import oglib.og_def
from oglib.og_mongo import mongo_connector

from importlib import import_module
from progressbar import *

from Bio import SeqIO
from Bio import SeqFeature
from Bio.SeqRecord import SeqRecord

info_type_list = ['general', 'annotations', 'features', 'sequence']

def add_show_command(subparsers):
    subparser = subparsers.add_parser("show",
        help="Show detailed info of Organelle Genomes in the local database",
        usage='OGSync show refSeq infoType [options]',
        description='''
Show Organelle Genomes in the local database into tree view.
Currently, it has three type:
general: show the general information of the genome, such as name, id, description and etc.
sequence: show the nucleatide sequence of the genome.
annotations: show the annotations of the genome, such as molecule_type, taxonomy, source and etc.
features: show the features of the genome, including the qualifer info.
''')
    subparser.add_argument('refSeq', nargs=1,
        help = "the refSeq code of the Organcelle Genome in the database, NCNumber or INSDC."
    )
    subparser.add_argument('type', nargs=1,
        choices = info_type_list,
        help = "the type of info, currently it supports to show [general, annotations, features] info into tree view."
    )
    subparser.add_argument('-j', '--output-json',
        default=False, action='store_true',
        help="Output info in json format instead of tree view, False by default")
    subparser.set_defaults(callback=run_show_command)
    return subparser

class MyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, SeqFeature.FeatureLocation):
            return {"FeatureLocation":obj.__dict__}
        if isinstance(obj, SeqFeature.CompoundLocation):
            return {"CompoundLocation":obj.__dict__}
        if isinstance(obj, SeqFeature.Reference):
            return {"Reference":obj.__dict__}
        else:
            return super().default(obj)

def run_show_command(args):
    if args.debug:
        print("Connecting to local database.")
    ogsync_db = mongo_connector[oglib.og_def.OG_SYNC_NAME]
    data = ogsync_db[oglib.og_def.OG_SYNC_DATA_NAME]
    genbank = ogsync_db[oglib.og_def.OG_SYNC_GENBANK_NAME]

    # check the availability of the refseq
    status = oglib.og_def.check_one(args.refSeq[0], data)
    if status=="remote":
        print('Please use "OGSync add %s" to download, and try again.'%(args.refSeq[0]))
        return
    elif status==0 or status==2:
        return
    
    #show_dict = getattr( import_module("oglib.og_def"), "get_%s_from_gb"%(args.type[0]) )(gb_str)
    show_dict = oglib.og_def.show_genbank(args.refSeq[0], genbank, args.type[0])
    if args.debug:
        print(show_dict)

    if args.output_json:
        print(json.dumps(show_dict, indent=2, cls=MyEncoder))
    else:
        oglib.og_def.PyTreeView().display(show_dict)