#!/usr/bin/python3

import argparse, os, datetime, csv, json, io, time, hashlib, difflib
import oglib.og_def
from oglib.og_mongo import mongo_connector
from prettytable import PrettyTable

from progressbar import *

def add_list_command(subparsers):
    subparser = subparsers.add_parser("list",
        help="List Organelle Genomes in local/remote database",
        usage='OGSync list [options]',
        description="List Organelle Genomes in local/remote database")
    subparser.add_argument('-r', '--remote',
        default=False, action='store_true',
        help="List Organelle Genomes in remote database, False by default.")
    subparser.add_argument('-l', '--list-RefSeq-only',
        default=False, action='store_true',
        help="List the RefSeq for Organelle Genomes, False by default")
    subparser.add_argument('-j', '--output-json',
        default=False, action='store_true',
        help="Output in json format, False by default")
    subparser.set_defaults(callback=run_list_command)
    return subparser

def run_list_command(args):
    if args.debug:
        print("Connecting to local database.")
    ogsync_db = mongo_connector[oglib.og_def.OG_SYNC_NAME]
    data = ogsync_db[oglib.og_def.OG_SYNC_DATA_NAME]

    query = {"og_type":"og_organism", "status":"sync"}
    if args.remote:
        query["status"]="remote"
    results = list( data.find(query) )

    if args.list_RefSeq_only:
        output_list = []
        for index,result in enumerate(results):
            output_list.append({
                "RefSeq":result["RefSeq"],
                "INSDC":result["INSDC"]
            })
        if args.output_json:
            print( json.dumps({"og_list":output_list}) )
        else:
            print( "\n".join( [ "(".join(list(x.values()))+")" for x in output_list ] ) )
    elif len(results)==0:
        print('There is no local Organelle Genomes, Try "OGSync add --help" to see detail.')
    else:
        tb = PrettyTable()
        tb.field_names = oglib.og_def.OG_SYNC_LIST_FIELD
        for index,result in enumerate(results):
            row = [ result[field] for field in tb.field_names[1:] ]
            tb.add_row( [index+1]+row )
        if args.output_json:
            print( tb.get_json_string() )
        else:
            print(tb)
    
    if args.debug:
        print("Listing is done.")