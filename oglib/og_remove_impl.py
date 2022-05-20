#!/usr/bin/python3

import argparse, os, datetime, csv, json, io, time, hashlib, difflib
import oglib.og_def
from oglib.og_mongo import mongo_connector

from progressbar import *

def add_remove_command(subparsers):
    subparser = subparsers.add_parser("remove",
        help="Remove Organelle Genomes from local database",
        usage='OGSync remove refSeq1[,refSeq2,...] [options]',
        description="Remove Organelle Genomes from local database, but keep in the list")
    subparser.add_argument('refSeq', nargs='+',
        help = "the refSeq code of the Organcelle Genome in the database, NCNumber or INSDC."
    )
    subparser.set_defaults(callback=run_remove_command)
    return subparser

def run_remove_command(args):
    if args.debug:
        print("Connecting to local database.")
    ogsync_db = mongo_connector[oglib.og_def.OG_SYNC_NAME]
    data = ogsync_db[oglib.og_def.OG_SYNC_DATA_NAME]

    # check the availability of the list
    check_list = ",".join(args.refSeq).split(",")
    checked_list = oglib.og_def.check_list(check_list, data, False, args.debug)
    remove_list = checked_list["sync"]
    print("Warning: The following refseq was not added in local database.")
    print( ", ".join(checked_list["remote"]) )
    if args.debug:
        print("Removing list is: ",remove_list)

    # get ready to remove
    total = len(remove_list)
    widgets = ['Removing: (', Counter(), '/'+str(total)+') ', Percentage(), ' ', Bar('▒'),' ', Timer(), ' ', ETA()]
    pbar = ProgressBar(widgets=widgets, maxval=total)
    if not args.debug:
        pbar.start()
    
    # remove each from database
    for index,refseq in enumerate(remove_list):
        # construct the database operation
        query = {"og_type":"og_organism", "$or":[{"RefSeq":refseq}, {"INSDC":refseq}] }
        update = {"$set":{
            "status":"remote",
            "GI Number": "",
            "GeneBank": "",
            "Sequence": "",
            "CDS Nucleotid": "",
            "CDS Protein": ""
        }}

        # submit to database
        data.update_one(query, update)

        if not args.debug:
            pbar.update( index+1 )
        else:
            print(refseq," has been removed from database, (%d/%d)."%(index+1,total))
    
    if not args.debug:
        pbar.finish()