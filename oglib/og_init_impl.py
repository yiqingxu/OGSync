#!/usr/bin/python3

import oglib.og_def
from oglib.og_mongo import mongo_connector

def add_init_command(subparsers):
    subparser = subparsers.add_parser("init",
        help="Initialize the local OGSync database",
        usage='OGSync init [--reinstall]',
        description="Setup the local OGSync database in MongoDB")
    subparser.add_argument('-r', '--reinstall',
        default=False, action='store_true',
        help="force to reset the database. NOTE: THIS OPERATION WILL EMPTY THE DATABASE!!")
    subparser.set_defaults(callback=run_init_command)
    return subparser

def run_init_command(args):
    dbs = mongo_connector.list_database_names()
    if oglib.og_def.OG_SYNC_NAME not in dbs:
        print("OGSync is not initialized, seting up...")

    if args.reinstall:
        mongo_connector.drop_database(oglib.og_def.OG_SYNC_NAME)
        print("OGSync is erased.")

    # create database for OGSync if not exist
    ogsync_db = mongo_connector[oglib.og_def.OG_SYNC_NAME]
    if oglib.og_def.OG_SYNC_CONFIG_NAME not in ogsync_db.list_collection_names():
        print("OGSync_config is not configured, intializing...")
    
    # create config is not exist
    config = ogsync_db[oglib.og_def.OG_SYNC_CONFIG_NAME]
    if config.find_one({"_id":0}):
        config.update_one({"_id": 0}, {"$set": {"version":oglib.og_def.OG_SYNC_VERSION} } )
    else:
        config.insert_one({"_id": 0, "version":oglib.og_def.OG_SYNC_VERSION})
        config.insert_one({"_id": 1, "config":{"NCBI_API_KEY":""}})
    
    if oglib.og_def.OG_SYNC_DATA_NAME not in ogsync_db.list_collection_names():
        print("OGSync_data is empty, READY.")
    data = ogsync_db[oglib.og_def.OG_SYNC_DATA_NAME]
    if data.find_one({"_id":0}):
        data.update_one({"_id": 0}, {"$set": {"version":oglib.og_def.OG_SYNC_VERSION} } )
    else:
        data.insert_one({"_id": 0, "version":oglib.og_def.OG_SYNC_VERSION})
        data.create_index([('MD5', 1)])
    
    if oglib.og_def.OG_SYNC_GENEBANK_NAME not in ogsync_db.list_collection_names():
        print("OGSync_genebank is empty, READY.")
    genebank = ogsync_db[oglib.og_def.OG_SYNC_GENEBANK_NAME]
    if genebank.find_one({"_id":0}):
        genebank.update_one({"_id": 0}, {"$set": {"version":oglib.og_def.OG_SYNC_VERSION} } )
    else:
        genebank.insert_one({"_id": 0, "version":oglib.og_def.OG_SYNC_VERSION})
    
    if args.debug:
        print( str(config.find_one({"_id":0})) + " in "+oglib.og_def.OG_SYNC_CONFIG_NAME)
        print( str(data.find_one({"_id":0})) + " in "+oglib.og_def.OG_SYNC_DATA_NAME)

    print("OGSync is setup.")