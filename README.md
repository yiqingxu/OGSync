# OGSync

![gnu](https://img.shields.io/badge/GNU-3-green?style=flat-square&logo=gnu)
![docker-compose](https://img.shields.io/badge/dockercompose-3-blue?style=flat-square&logo=docker)
![python](https://img.shields.io/badge/python-3-red?style=flat-square&logo=python)

OGSync is an Organelle Genome Synchronizor to manage a local dataset for easy-access to NCBI Organelle Database.

## Quick Start

1. Download the repo and unzip
2. Start up the docker stack

``` sh
docker-compose up
```

3. Open termial and login OGSync container

``` sh
docker exec -it ogsync /bin/bash
```

4. Init the OGSync

``` sh
OGSync init
```

5. Update status and sync the list with NCBI database

``` sh
OGSync update
```

6. List all avaiable refSeq, use `-r` to display the `remote` genome code

``` sh
OGSync list -r
```

7. Add some refSeq by using `add` command. For example, if you want to sync the chloroplast genome of *Arabidopsis thaliana* to your local database.

``` sh
OGSync add NC_000932.1
```

8. Try more command to manage your local database and sync.

``` sh
OGSync --help
```