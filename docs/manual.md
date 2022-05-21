# OGSync's Documentation

![GPL-3.0](https://img.shields.io/github/license/yiqingxu/OGSync) ![website](https://img.shields.io/website?down_color=lightgrey&down_message=offline&up_color=green&up_message=online&url=https%3A%2F%2Fyiqingxu.github.io%2FOGSync%2F) ![platform](https://img.shields.io/badge/platform-win--64%20%7C%20win--32%20%7C%20osx--arm64%20%7C%20osx--64%20%7C%20linux--64%20%20%7C%20linux--aarch64%20%7C%20linux--ppc64le-lightgrey)

![gnu](https://img.shields.io/badge/GNU-3-green?style=flat-square&logo=gnu) ![docker-compose](https://img.shields.io/badge/dockercompose-3-blue?style=flat-square&logo=docker) ![mongodb](https://img.shields.io/badge/mongo-4.4-green?style=flat-square&logo=mongodb)

![python](https://img.shields.io/badge/python-3-red?style=flat-square&logo=python) ![biopython](https://img.shields.io/pypi/status/biopython?label=biopython&style=flat-square) ![pymongo](https://img.shields.io/pypi/status/pymongo?label=pymongo&style=flat-square) ![progressbar](https://img.shields.io/pypi/status/progressbar?label=progressbar&style=flat-square) ![colorama](https://img.shields.io/pypi/status/colorama?label=colorama&style=flat-square)

![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/yiqingxu/OGSync) ![GitHub release (release name instead of tag name)](https://img.shields.io/github/v/release/yiqingxu/OGSync?display_name=release&include_prereleases) ![siz](https://img.shields.io/github/languages/code-size/yiqingxu/OGSync) ![repo](https://img.shields.io/github/repo-size/yiqingxu/OGSync) ![Docker Image Size (latest by date)](https://img.shields.io/docker/image-size/yiqingxu/ogsync)

## How to use

### Begin with `OGSync` CLI

The PATH for OGSync in the docker container is pre-configured, and user can just type `OGSync` to start.

### Get Help

You can type `-h/--help` after any command to get the help info:

``` sh
root@mst:/opt/OGSync# OGSync --help
usage: OGSync [-h] [-v] {home,init,add,remove,update,config,list,show} ...

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show OGSync version number and exit

avaiable commands:
  {home,init,add,remove,update,config,list,show}
    home                Display OGSync Homepage info and exit
    init                Initialize the local OGSync database
    add                 Add Organelle Genomes into local database and sync with NCBI
    remove              Remove Organelle Genomes from local database
    update              Sync the local OGSync list with NCBI Organelle Genome database
    config              Configurate and customize OGSync
    list                List Organelle Genomes in local/remote database
    show                Show detailed info of Organelle Genomes in the local database
```

### Check OGSync Version

Use `OGSync -v` or `OGSync --version` to see the version info.

``` sh
root@mst:/opt/OGSync# OGSync --version
OGSync 0.1.0
```

## Sub commands

Currently, OGSync provides 8 sub commands:

|  Sub command   | Description  |
|  ----  | ----  |
| [home](sub.md#ogsync-home) |Display OGSync Homepage info and exit |
| [init](sub.md#ogsync-init) | Initialize the local OGSync database |
| [add](sub.md#ogsync-add) | Add Organelle Genomes into local database and sync with NCBI |
| [remove](sub.md#ogsync-remove) | Remove Organelle Genomes from local database |
| [update](sub.md#ogsync-update) | Sync the local OGSync list with NCBI Organelle Genome database |
| [config](sub.md#ogsync-config) | Configurate and customize OGSync |
| [list](sub.md#ogsync-list) | List Organelle Genomes in local/remote database |
| [show](sub.md#ogsync-show) | Show detailed info of Organelle Genomes in the local database |

## For Developers

### Trouble Shooting

OGSync provide `-d`/`--debug` option as a global argument, and it turns on a verbose log for troubleshooting.

For example, when use `update` sub command, it will display a progress bar while updating, such as:

``` sh
root@mst:/opt/OGSync# OGSync update
Syncing with NCBI database.
Updating: (23343/23343) 100% |▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒| Elapsed Time: 0:00:10 Time: 0:00:10
OGSync is updated: 0 Oganelle Genomes are updated.
```

`--debug` will show a more detailed `update` progress, which normally consume more time, such as:

``` sh
root@mst:/opt/OGSync# OGSync update --debug
DEBUG MODE is set ON by options: {'show_diff': False, 'debug': True, 'callback': <function run_update_command at 0xffffa24ba040>}
Syncing with NCBI database.
Connecting to local database.
Re-formating dataset.
Updating the Organelle Genome List.
Updating: (23343/23343) 100% |▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒| Elapsed Time: 0:00:11 Time: 0:00:11
OGSync is updated: 0 Oganelle Genomes are updated.
23344  were up-to-date and not modified.
```

### OGSync API

Developers can use `OGSync` as a service by adding `-j`/`--json` argument.

For example, when use `show` sub command, it will display a tree view in terminal, such as:

``` sh
root@mst:/opt/OGSync# OGSync show NC_048988.1 general
general
├─ id: NC_048988.1
├─ name: NC_048988
├─ description: Abbottina binhi isolate AB-PEx-PiG01 mitochondrion, complete genome
├─ dbxrefs
│  └─ 0: BioProject:PRJNA642447
└─ main_feature
   ├─ location
   │  └─ FeatureLocation
   │     ├─ _start: 0
   │     ├─ _end: 16609
   │     ├─ _strand: 1
   │     ├─ ref: None
   │     └─ ref_db: None
   ├─ type: source
   ├─ id: <unknown id>
   └─ qualifiers
      ├─ organism
      │  └─ 0: Abbottina binhi
      ├─ organelle
      │  └─ 0: mitochondrion
      ├─ mol_type
      │  └─ 0: genomic DNA
      ├─ isolate
      │  └─ 0: AB-PEx-PiG01
      ├─ db_xref
      │  └─ 0: taxon:2745886
      └─ country
         └─ 0: China: Pingguo, Guangxi province
```

`-j` will parse the result as json string, which is more flexible for further development, such as:

``` sh
root@mst:/opt/OGSync# OGSync show NC_048988.1 general -j
general
{
  "id": "NC_048988.1",
  "name": "NC_048988",
  "description": "Abbottina binhi isolate AB-PEx-PiG01 mitochondrion, complete genome",
  "dbxrefs": [
    "BioProject:PRJNA642447"
  ],
  "main_feature": {
    "location": {
      "FeatureLocation": {
        "_start": 0,
        "_end": 16609,
        "_strand": 1,
        "ref": null,
        "ref_db": null
      }
    },
    "type": "source",
    "id": "<unknown id>",
    "qualifiers": {
      "organism": [
        "Abbottina binhi"
      ],
      "organelle": [
        "mitochondrion"
      ],
      "mol_type": [
        "genomic DNA"
      ],
      "isolate": [
        "AB-PEx-PiG01"
      ],
      "db_xref": [
        "taxon:2745886"
      ],
      "country": [
        "China: Pingguo, Guangxi province"
      ]
    }
  }
}
```
