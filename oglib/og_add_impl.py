#!/usr/bin/python3

import argparse, os, datetime, csv, json, io, time, hashlib, difflib, re
import oglib.og_def
from oglib.og_mongo import mongo_connector
from Bio import SeqIO
from Bio import SeqFeature

from progressbar import *

def add_add_command(subparsers):
    subparser = subparsers.add_parser("add",
        help="Add Organelle Genomes into local database and sync with NCBI",
        usage='OGSync add refSeq1[,refSeq2,...] [options]',
        description="Add Organelle Genomes into local database and sync with NCBI")
    subparser.add_argument('refSeq', nargs='+',
        help = "the refSeq code of the Organcelle Genome in the database, NCNumber or INSDC."
    )
    subparser.add_argument('-f', '--force-upgrade',
        default=False, action='store_true',
        help="if force to upgrade the local data, False by default")
    subparser.set_defaults(callback=run_add_command)
    return subparser

def run_add_command(args):
    if args.debug:
        print("Connecting to local database.")
    ogsync_db = mongo_connector[oglib.og_def.OG_SYNC_NAME]
    data = ogsync_db[oglib.og_def.OG_SYNC_DATA_NAME]
    genbank = ogsync_db[oglib.og_def.OG_SYNC_GENBANK_NAME]

    # check the availability of the list
    check_list = ",".join(args.refSeq).split(",")
    checked_list = oglib.og_def.check_list(check_list, data, args.force_upgrade, args.debug)
    add_list = checked_list["remote"]
    if len(checked_list["sync"])>0:
        print("Warning: The following refseq has already been added in local database, try -f/--force-upgrade to overwrite the data.")
        print( ", ".join(checked_list["sync"]) )
    if args.debug:
        print("Downloading list is: ",add_list)

    # get ready to download
    total = len(add_list)
    widgets = ['Adding: (', Counter(), '/'+str(total)+') ', Percentage(), ' ', Bar('▒'),' ', Timer(), ' ', ETA()]
    pbar = ProgressBar(widgets=widgets, maxval=total)
    if not args.debug:
        pbar.start()

    # insert each into database
    for index,NCNumber in enumerate(add_list):
        id = oglib.og_def.get_refseq_id(NCNumber, data)
        down_str = download_from_NCBI(NCNumber, data, args.debug)

        gb_str = down_str[0]
        add_genbank(id, gb_str, genbank, args.debug)
        cds_nucleotide_str = down_str[1]
        add_CDS(id, cds_nucleotide_str, "Nucleotide", genbank, args.debug )
        cds_protein_str = down_str[2]
        add_CDS(id, cds_protein_str, "Protein", genbank, args.debug )

        if not args.debug:
            pbar.update( index+1 )
        else:
            print("**       %6d / %-6d        **"%(index+1,total))

    if not args.debug:
        pbar.finish()
    else:
        print("**********************************")

def download_from_NCBI(refseq, db_connector, if_show_debug=False):
    if if_show_debug:
        print("**********************************")
        print("Downloading "+refseq+" from NCBI")

    # downlaod content from NCBI
    gi = oglib.og_def.get_gi_from_NC(refseq, if_show_debug)
    gb = oglib.og_def.get_genebank(refseq, if_show_debug)
    seq = oglib.og_def.get_sequence(refseq, if_show_debug)
    ft = oglib.og_def.get_feature_table(refseq, if_show_debug)
    cds_nucleotide = oglib.og_def.get_CDS_nucleotide(gi, if_show_debug)
    cds_protein = oglib.og_def.get_CDS_protein(gi, if_show_debug)

    # construct the database operation
    query = {"og_type":"og_organism", "$or":[{"RefSeq":refseq}, {"INSDC":refseq}] }
    update = {"$set":{
        "status":"sync",
        "GI Number": gi,
        "Features": ft,
        "GenBank": gb,
        "Sequence": seq,
        "CDS Nucleotid": cds_nucleotide,
        "CDS Protein": cds_protein
    }}

    # submit to database
    db_connector.update_one(query, update)
    return [gb,cds_nucleotide,cds_protein]

def add_genbank(id, gb_str, db_connector, if_show_debug=False):
    # extract info from gb_file
    seq = oglib.og_def.get_sequence_from_gb(gb_str)
    info = oglib.og_def.get_general_from_gb(gb_str)
    annotations = oglib.og_def.get_annotations_from_gb(gb_str)
    features = oglib.og_def.get_features_from_gb(gb_str)

    # construct the database operation
    query = {"og_type":"og_genbank", "$or":[{"RefSeq":id["RefSeq"]}, {"INSDC":id["INSDC"]}] }
    insert ={
        "og_type":"og_genbank",
        "RefSeq":id["RefSeq"],
        "INSDC":id["INSDC"],
        "GenBank": gb_str,
        "GenBank Info": json.dumps(info, cls=oglib.og_def.GBEncoder),
        "Sequence": json.dumps(seq, cls=oglib.og_def.GBEncoder),
        "Features": json.dumps(features, cls=oglib.og_def.GBEncoder),
        "Annotations": json.dumps(annotations, cls=oglib.og_def.GBEncoder)
    }
    update = {"$set": insert}

    # submit to database
    if db_connector.update_one(query, update).matched_count==0:
        db_connector.insert_one(insert)

    add_features(id, db_connector, features, if_show_debug)

    if if_show_debug:
        print(id["RefSeq"]+" has been added to database and genbank")
        print("**********************************")

def add_features(id, db_connector, features, if_show_debug=False):
    for index in features:
        # construct the database operation
        query = {"og_type":"og_feature", "$or":[{"RefSeq":id["RefSeq"]}, {"INSDC":id["INSDC"]}] }
        insert ={
            "og_type":"og_feature",
            "RefSeq":id["RefSeq"],
            "INSDC":id["INSDC"],
        }

        # submit to database
        for f in features[index]:
            insert.update( {
                "Index": index,
                f: json.dumps(features[index][f], cls=oglib.og_def.GBEncoder)
            } )
            query.update( {"Index": index} )
        update = {"$set": insert}

        if db_connector.update_one(query, update).matched_count==0:
            db_connector.insert_one(insert)
    if if_show_debug:
        print("%d features are added/updated."%(len(features)) )

def add_CDS(id, cds_str, cds_type, db_connector, if_show_debug=False):
    cds_record = SeqIO.parse(io.StringIO(cds_str), "fasta")
    for seq_record in cds_record:
        query = {"og_type":"og_cds", "Description": seq_record.description, "$or":[{"RefSeq":id["RefSeq"]}, {"INSDC":id["INSDC"]}] }
        insert ={
            "og_type":"og_cds",
            "RefSeq":id["RefSeq"],
            "INSDC":id["INSDC"],
            "Type": cds_type,
            "Name": seq_record.name,
            "CDS id": seq_record.id,
            "Description": seq_record.description,
            "Sequence": str(seq_record.seq)
        }

        ds = re.findall( r'\[(\S*)=(\S*)\]', seq_record.description )
        for d in ds:
            insert.update({d[0]: d[1]})

        update = {"$set": insert}

        if db_connector.update_one(query, update).matched_count==0:
            db_connector.insert_one(insert)

    if if_show_debug:
        print("%d %s CDSs are added/updated."%(len([r for r in SeqIO.parse(io.StringIO(cds_str), "fasta")]), cds_type) )
