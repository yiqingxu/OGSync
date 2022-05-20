#!/usr/bin/python3

import argparse, os, datetime, csv, json, io, time, hashlib, difflib
import oglib.og_def
from oglib.og_mongo import mongo_connector

from progressbar import *

__DEFAULT_NCBI_FULL_LIST_ADDRESS__ = 'https://www.ncbi.nlm.nih.gov/genomes/Genome2BE/genome2srv.cgi?action=download&orgn=&' \
    'report=organelles&group=--%20All%20Eukaryota%20--&subgroup=--%20All%20Eukaryota%20--&' \
    'host=All&format='

def add_update_command(subparsers):
    subparser = subparsers.add_parser("update",
        help="Sync the local OGSync list with NCBI Organelle Genome database",
        usage='OGSync update [options]',
        description="Sync the local OGSync list with NCBI Organelle Genome database")
    subparser.add_argument('-d', '--show-diff',
        default=False, action='store_true',
        help="if show the diffence of the genome list between current and latest, False by default")
    subparser.set_defaults(callback=run_update_command)
    return subparser

def run_update_command(args):
    #Download file from NCBI
    print("Syncing with NCBI database.")
    NCBI_CSV = oglib.og_def.load_from_url(__DEFAULT_NCBI_FULL_LIST_ADDRESS__)
    # NCBI_CSV = oglib.og_def.load_from_file("NCList.csv")
    ncbi_data = csv.DictReader(io.StringIO(NCBI_CSV), delimiter='\t')

    if args.debug:
        print("Connecting to local database.")
    ogsync_db = mongo_connector[oglib.og_def.OG_SYNC_NAME]
    data = ogsync_db[oglib.og_def.OG_SYNC_DATA_NAME]

    if args.debug:
        print("Re-formating dataset.")
    update_counter=0
    remain_counter=0
    data.update_many(
        {"latest":"latest"},
        {"$set":{"latest":"history"}}
    )

    # insert the csv_file as a record
    data_json =  json.dumps( [ row for row in ncbi_data ], ensure_ascii=False)
    og_list = {"og_type":"og_list", "raw_csv":NCBI_CSV, "raw_json":data_json}
    if args.debug:
        print("Updating the Organelle Genome List.")
    if data.find_one(og_list):
        remain_counter+=1
        data.update_one(
            og_list,
            {"$set":{"latest":"latest"}}
        )
    else:
        update_counter+=1
        og_list["latest"]="latest"
        data.insert_one(og_list)

    total = len(json.loads(data_json))
    # use the progressbar
    widgets = ['Updating: (', Counter(), '/'+str(total)+') ', Percentage(), ' ', Bar('▒'),' ', Timer(), ' ', ETA()]
    pbar = ProgressBar(widgets=widgets, maxval=total).start()

    # insert each row as a record
    for index,row in enumerate(json.loads(data_json)):
        og_organism = {"og_type":"og_organism"}
        og_organism.update(row)

        md5 = hashlib.md5(str(row).encode(encoding='utf-8')).hexdigest()
        og_organism["MD5"] = md5

        query = {"latest":"history", "MD5":md5}

        if data.find_one(query):
            remain_counter+=1
            data.update_one(
                {"MD5":md5},
                {"$set":{"latest":"latest"}}
            )
        else:
            update_counter+=1
            og_organism["latest"]="latest"
            og_organism["status"]="remote"
            data.insert_one(og_organism)

        pbar.update( index )
    pbar.finish()

    print("OGSync is updated: "+str(update_counter)+" Oganelle Genomes are updated.")
    if args.debug:
        print(remain_counter," were up-to-date and not modified.")

    if args.show_diff:
        ogsync_db = mongo_connector[oglib.og_def.OG_SYNC_NAME]
        data = ogsync_db[oglib.og_def.OG_SYNC_DATA_NAME]

        versions = data.find({"og_type":"og_list"}, {"raw_csv":1}).sort('_id',1).limit(2)
        versions = [x for x in versions]

        if len(versions)==2:
            change = difflib.unified_diff(versions[0]["raw_csv"].splitlines(), versions[1]["raw_csv"].splitlines(), lineterm='')
            print("\n".join(change))
        else:
            print("No history data is located to show the diff.")