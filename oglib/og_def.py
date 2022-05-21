#-*- coding: UTF-8 -*-

import requests, io, json
from oglib.og_mongo import mongo_connector

from Bio import SeqIO
from Bio import SeqFeature
from Bio.SeqRecord import SeqRecord

OG_SYNC_VERSION = "0.1.0"
OG_SYNC_NAME = 'OGSync'
OG_SYNC_CONFIG_NAME = OG_SYNC_NAME+"_config"
OG_SYNC_DATA_NAME = OG_SYNC_NAME+"_data"
OG_SYNC_GENBANK_NAME = OG_SYNC_NAME+"_genbank"
OG_SYNC_CITE = " OGSync, bioinformatics, 2022"

OG_SYNC_CONFIGS = ['NCBI_API_KEY']
OG_SYNC_LIST_FIELD = ["No", "RefSeq", "INSDC", "Group","SubGroup","Type","Size (Kb)", "#Organism/Name"]

__NCBI_EFETCH_URL__ = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi" #?api_key=f9f60bb9e0b84572e2674aeea0b9c81d8808&"

def load_from_url(url):
    try:
        _rawData = requests.get(url).text
    except Exception as e:
        print("Network Error, abort")
        print(e)
        exit()
    return _rawData

def load_from_file(fileName):
    with open(fileName, 'r') as f:
        return f.read()

def request_ncbi_data(ncbi_url, if_show_debug=False):
    if if_show_debug:
        print("Downloading NCBI content from: "+ncbi_url)
    ifOK = False
    content = ""
    while not ifOK:
        content = requests.get(ncbi_url).text
        if content.find("error\":\"API rate limit exceeded")!=-1:
            sleep( random()*3 )
            if if_show_debug:
                print(content)
                print(ncbi_url+"is re-downloading...")
            continue
        ifOK = True
    return content.strip()

def get_efetch_url(key=""):
    ogsync_db = mongo_connector[OG_SYNC_NAME]
    config = ogsync_db[OG_SYNC_CONFIG_NAME]
    result = config.find_one({"_id":1})['config'][ OG_SYNC_CONFIGS[0] ]
    if result:
        return __NCBI_EFETCH_URL__+"?api_key="+result+"&"
    else:
        return __NCBI_EFETCH_URL__+"?"

def get_gi_from_NC(NCNumber, if_show_debug):
    url = get_efetch_url() + "db=nucleotide&id=" + NCNumber + "&rettype=gi"
    return request_ncbi_data(url, if_show_debug)

def get_genebank(NCNumber, if_show_debug):
    url = get_efetch_url() + "db=nuccore&id=" + NCNumber + "&rettype=gb&retmode=text"
    return request_ncbi_data(url, if_show_debug)

def get_feature_table(NCNumber, if_show_debug):
    url = get_efetch_url() + "db=nuccore&id=" + NCNumber + "&rettype=ft&retmode=text"
    return request_ncbi_data(url, if_show_debug)

def get_sequence(NCNumber, if_show_debug):
    url = get_efetch_url() + "db=nuccore&id=" + NCNumber + "&rettype=fasta&retmode=text"
    return request_ncbi_data(url, if_show_debug)

def get_CDS_protein_by_NC(NCNumber, if_show_debug):
    giNumber = get_gi_from_NC(NCNumber, if_show_debug)
    return get_CDS_protein(giNumber, if_show_debug)

def get_CDS_nucleotide_by_NC(NCNumber, if_show_debug):
    giNumber = get_gi_from_NC(NCNumber, if_show_debug)
    return get_CDS_nucleotide(giNumber, if_show_debug)

def get_CDS_protein(GINumber, if_show_debug):
    url = get_efetch_url() + "db=nuccore&id=" + GINumber + "&rettype=fasta_cds_aa"
    return request_ncbi_data(url, if_show_debug)

def get_CDS_nucleotide(GINumber, if_show_debug):
    url = get_efetch_url() + "db=nuccore&id=" + GINumber + "&rettype=fasta_cds_na"
    return request_ncbi_data(url, if_show_debug)

def check_one(refseq, db_connector, force_upgrade=False):
    query = {"og_type":"og_organism", "$or":[{"RefSeq":refseq}, {"INSDC":refseq}] }
    result = list( db_connector.find(query) )
    if len(result)==0:
        print("Warning: RefSeq named "+refseq+" is not available, please double-check. It should be a refSeq 'NC_Number' or INSDC ID.")
        return 0
    if len(result)>=2:
        print("Warning: RefSeq named "+refseq+" matches multiple organism, please double-check.")
        return 2
    if force_upgrade:
        return "remote"
    else:
        return result[0]['status']

def check_list(check_list, db_connector, force_upgrade=False, if_show_debug=False):
    if if_show_debug:
        print("Checking RefSeq List: ",check_list)
    sync_list = set()
    remote_list = set()
    for refseq in check_list:
        result = check_one(refseq, db_connector, force_upgrade)
        if result==0 or result==2:
            continue
        elif result=="sync":
            sync_list.add(refseq)
        elif result=="remote":
            remote_list.add(refseq)
    re = {"sync":sync_list, "remote":remote_list}
    if if_show_debug:
        print("Checked RefSeq List: ",re)
    return re

# by github.com/basemax/PyTreeView
class PyTreeView():
    offset=0
    lines={}
    def parse(self, items):
        # global offset#, lines
        # if 'lines' in locals():
        #     self.lines={}
        start=0
        result=[]
        if isinstance(items, list):
            items={v: k for v, k in enumerate(items)}
        # print(items)
        dictPosition=1
        dictSize=len(items)
        for attr, value in items.items():
            result.append([self.offset, '', '', {}])
            if dictPosition == dictSize:
                result[-1][1]='└─'
                # if self.offset in self.lines:
                #     self.lines[self.offset]=self.lines[self.offset] + [len(result)]
                # else:
                #     self.lines.update({self.offset : [len(result)]})
                self.lines[self.offset]=len(result)
                # print("add", self.lines)
            else:
                result[-1][1]='├─'
                # self.lines.update({self.offset : [len(result)]})
                # print(self.lines, self.offset)
                if self.offset in self.lines:
                    del self.lines[self.offset]
            result[-1][2]=str(attr)
            result[-1][3]=self.lines.copy()
            if value == [] or value == {}:
                result[-1][2]+=': []'
            elif isinstance(value, list) or isinstance(value, dict):
                start=len(result)
                self.offset=self.offset+1
                result = result + self.parse(value)
                self.offset=self.offset-1
            else:
                result[-1][2]+=': '+str(value)
            dictPosition+=1
        return result

    def display(self, result):
        result=self.get(result)
        if result[-1] == '\n':
            result=result[:-1]
        print(result)

    def get(self, result):
        answer=''
        tree=self.parse(result)
        # print(tree)
        # for v in tree:
        #     print(v)
        # print(self.lines)
        i=0
        for v in tree:
            result=''
            if v[0] != 0:
                for i in range(v[0]):
                    # print(self.lines[i])
                    if i not in v[3]:
                        result+='│  '
                    else:
                        result+='   '
            else:
                result+=v[0]*'   '
            result+=v[1]+' '+v[2]
            # print(result)
            answer+=result+'\n'
            i+=1
        return answer

class GBEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, SeqFeature.FeatureLocation):
            return {"FeatureLocation":obj.__dict__}
        if isinstance(obj, SeqFeature.CompoundLocation):
            return {"CompoundLocation":obj.__dict__}
        if isinstance(obj, SeqFeature.Reference):
            return {"Reference":obj.__dict__}
        else:
            return super().default(obj)

def get_gb_str(refseq, db_connector):
    if check_one(refseq, db_connector)=="sync":
        gb = db_connector.find_one({"og_type":"og_organism", "$or":[{"RefSeq":refseq}, {"INSDC":refseq}]})
        return gb["GenBank"]
    else:
        return None

def get_refseq_id(refseq, db_connector):
    if check_one(refseq, db_connector)=="sync":
        record = db_connector.find_one({"og_type":"og_organism", "$or":[{"RefSeq":refseq}, {"INSDC":refseq}]}, ["RefSeq","INSDC"])
        return record
    else:
        return None

def get_sequence_from_gb(gb_str):
    gb_record = SeqIO.read(io.StringIO(gb_str), "genbank")
    return {gb_record.name:str(gb_record._seq)}

def show_genbank(refseq, db_connector, info_type):
    query = {"og_type":"og_genbank", "$or":[{"RefSeq":refseq}, {"INSDC":refseq}] }
    result = db_connector.find_one(query)

    print(info_type)
    trans={
        'general': "GenBank Info",
        'annotations': "Annotations",
        'features': "Features",
        'sequence': "Sequence"
    }
    return json.loads(result[trans[info_type]])

def get_general_from_gb(gb_str):
    gb_record = SeqIO.read(io.StringIO(gb_str), "genbank")
    show_gb = gb_record.__dict__
    show_mf = gb_record.features[0].__dict__
    pop_list=["_seq", "annotations","features","_per_letter_annotations"]
    for pop in pop_list:
        show_gb.pop(pop)
    show_gb.update({"main_feature":show_mf})
    return show_gb

def get_annotations_from_gb(gb_str):
    gb_record = SeqIO.read(io.StringIO(gb_str), "genbank")
    show_anno = gb_record.annotations

    jd = json.dumps(show_anno, indent=2, cls=GBEncoder)
    return json.loads(jd)

def get_features_from_gb(gb_str):
    gb_record = SeqIO.read(io.StringIO(gb_str), "genbank")
    ft_dict = {}
    for index,f in enumerate(gb_record.features[1:]) :
        ft_dict["Feature "+str(index+1)] = f.__dict__
    return ft_dict