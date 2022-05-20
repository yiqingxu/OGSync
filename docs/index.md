# OGSync: A CLI Organelle Genome Database Synchronizer

![gnu](https://img.shields.io/badge/GNU-3-green?style=flat-square&logo=gnu)
![docker-compose](https://img.shields.io/badge/dockercompose-3-blue?style=flat-square&logo=docker)
![python](https://img.shields.io/badge/python-3-red?style=flat-square&logo=python)
![mongodb](https://img.shields.io/badge/mongo-4.4-green?style=flat-square&logo=mongodb)

OGSync is NOSQL-based, freely-available and user-friendly database which provides a command-line-interface platform for bioinformatics researchers to manage a local database to synchronize and analysis multiple complete genome sequences, gene sequences and feature annotations of species. Currently, this toolkit provides the function of managing local database, synchronizing data with NCBI Organelle Genome database, and a tree viewer of gene sequences and feature records in the database. High availability of distributed file system interface, extensive data analysis of feature records from GenBank files, visually human view and json-based data interface makes OGSync a valuable data management system for studies in organelle genomics.

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

## Motivation

After three-decade accumulation of high-throughput sequence data from various organisms, biology scientist have published over 23 thousands organelle genomes. These organelle DNAs are featured as a compact genome structure compared with nuclear genomes; thus, there are efficient molecular tools for the analysis of gene structure, genome structure, organelle function and evolution. However, an integrated organelle genome data management system, which could enable users to synchronize and analysis data in local computing system, has not previously been developed.

## Architecture

![architecture](https://raw.githubusercontent.com/yiqingxu/OGSync/main/img/OGSync.png)

The main module of the tool is implemented by the python and can be deployed by docker engine. Users can freely build their own database and easily manage organelle genome features such as nucleatide CDS, protein, sequence, annotations and etc.

## Reference

1. NCBI Organelle Genome Resources, <https://ncbi.nlm.nih.gov/genome/organelle/>