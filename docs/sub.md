# OGSync sub commands

## OGSync home

Open the OGSync Home Page and exit

``` sh
OGSync.py home
```

### arguments

| argument | description |
| --- | --- |
| `-h`, `--help` | show this help message and exit |
| `--debug` | display debug info, for developers only |

### source

[oglib/og_home_impl.py](https://github.com/yiqingxu/OGSync/blob/main/oglib/og_home_impl.py)

---

## OGSync init

Setup the local OGSync database in MongoDB

``` sh
OGSync init [option]
```

### arguments

| argument | description |
| --- | --- |
| `-h`, `--help` | show this help message and exit |
| `-r`, `--reinstall` | force to reset the database. NOTE: THIS OPERATION WILL EMPTY THE DATABASE!! |
| `--debug` | display debug info, for developers only |

### source

[oglib/og_init_impl.py](https://github.com/yiqingxu/OGSync/blob/main/oglib/og_init_impl.py)

---

## OGSync update

Sync the local OGSync list with NCBI Organelle Genome database

``` sh
OGSync update [options]
```

### arguments

| argument | description |
| --- | --- |
| `-h`, `--help` | show this help message and exit |
| `-d`, `--show-diff` | if show the difference of the genome list between current and latest, False by default |
| `--debug` | display debug info, for developers only |

### source

[oglib/og_update_impl.py](https://github.com/yiqingxu/OGSync/blob/main/oglib/og__impl.py)

---

## OGSync config

Configurate and customize OGSync.

``` sh
OGSync config [options]
```

### arguments

| argument | description |
| --- | --- |
| -h, --help | show this help message and exit |
| --get {NCBI_API_KEY} | get the config parameters |
| --set SET | set the config parameters, in json format |
| --debug | display debug info, for developers only, in json format |

### optional config

| argument | description |
| --- | --- |
| NCBI_API_KEY | api key for NCBI E-utilities, see [this](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/) |

> For example, use these commands to get/set NCBI_API_KEY.
> Use `OGSync config --get NCBI_API_KEY` to get config value
> Use `OGSync config --set '{"NCBI_API_KEY"}':"xxx"'` to set config value

### source

[oglib/og_config_impl.py](https://github.com/yiqingxu/OGSync/blob/main/oglib/og_config_impl.py)

---

## OGSync add

Add Organelle Genomes into local database and sync with NCBI

``` sh
OGSync add refSeq1[,refSeq2,...] [options]
```

### arguments

| argument | description |
| --- | --- |
| `refSeq` | **REQUIRED**. the refSeq code of the Organelle Genome in the database, NCNumber or INSDC. |
| `-h`, `--help` | show this help message and exit |
| `-f`, `--force-upgrade` | if force to upgrade the local data, False by default |
| `--debug` | display debug info, for developers only |

### source

[oglib/og_add_impl.py](https://github.com/yiqingxu/OGSync/blob/main/oglib/og_add_impl.py)

---

## OGSync remove

Remove Organelle Genomes from local database

``` sh
OGSync remove refSeq1[,refSeq2,...] [options]
```

### arguments

| argument | description |
| --- | --- |
| `refSeq` | **REQUIRED**. the refSeq code of the Organelle Genome in the database, NCNumber or INSDC. |
| `-h`, `--help` | show this help message and exit |
| `--debug` | display debug info, for developers only |

### source
[oglib/og_remove_impl.py](https://github.com/yiqingxu/OGSync/blob/main/oglib/og_remove_impl.py)

## OGSync list

List Organelle Genomes in local/remote database

``` sh
OGSync list [options]
```

### arguments

| argument | description |
| --- | --- |
| `-h`, `--help` | show this help message and exit |
| `-r`, `--remote` | List Organelle Genomes in remote database, False by default. |
| `-l`, `--list-RefSeq-only` | List the RefSeq for Organelle Genomes, False by default |
| `-j`, `--output-json` | Output in json format, False by default|
| `--debug` | display debug info, for developers only |

### source
[oglib/og_list_impl.py](https://github.com/yiqingxu/OGSync/blob/main/oglib/og_list_impl.py)

---
## OGSync show

Show Organelle Genomes in the local database into tree/json view. 

``` sh
OGSync show refSeq infoType [options]
```

### info type

Currently, OGSync provides four types

| type | description |
| --- | --- |
| general | show the general information of the genome, such as name, id, description and etc. |
| sequence | show the nucleatide sequence of the genome. |
| annotations | show the annotations of the genome, such as molecule_type, taxonomy, source and etc. |
| features | show the features of the genome, including the qualifier info. |

### arguments

| argument | description |
| --- | --- |
| `refSeq` | **REQUIRED**. the refSeq code of the Organelle Genome in the database, NCNumber or INSDC. |
| `infoType` | **REQUIRED**. the type of info to show.
| `-h`, `--help` | show this help message and exit |
| `-j`, `--output-json` | Output in json format, False by default|
| `--debug` | display debug info, for developers only |

### source

[oglib/og_show_impl.py](https://github.com/yiqingxu/OGSync/blob/main/oglib/og_show_impl.py)
