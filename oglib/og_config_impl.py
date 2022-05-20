#!/usr/bin/python3

import argparse
import json
import oglib.og_def
from oglib.og_mongo import mongo_connector

config_choices = oglib.og_def.OG_SYNC_CONFIGS

def add_config_command(subparsers):
    subparser = subparsers.add_parser("config",
        help="Configurate and customize OGSync",
        usage='OGSync config [options]',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''
Configurate and customize OGSync.
For example, use these commands to get/set NCBI_API_KEY
    OGSync config --get NCBI_API_KEY
    OGSync config --set '{"NCBI_API_KEY"}':"xxx"'
''')
    subparser.set_defaults(callback=run_config_command)
    subparser.add_argument('--get', choices=config_choices,
        help="get the config parameters")
    subparser.add_argument('--set', type=json.loads, 
        help="set the config parameters")
    return subparser

def run_config_command(args):
    ogsync_db = mongo_connector[oglib.og_def.OG_SYNC_NAME]
    config = ogsync_db[oglib.og_def.OG_SYNC_CONFIG_NAME]
    
    # get a property from db
    if args.get in config_choices:
        print( config.find_one({"_id":1})['config'][args.get]  )
    # set a property into db
    elif args.set:
        if config.find_one({"_id":1}):
            config.update_one({"_id": 1}, {"$set": {"config":args.set}} ) 
        else:
            config.insert_one({"_id":1, "config":args.set})