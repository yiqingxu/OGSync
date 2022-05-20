#!/usr/bin/python3
#-*- coding: UTF-8 -*-

import pymongo

MONGO_LINK = 'mongodb://OGSync:OGSyncDocker@mongo:27017/'
mongo_connector = pymongo.MongoClient(MONGO_LINK)


if __name__ == '__main__':
    print(mongo.list_database_names())