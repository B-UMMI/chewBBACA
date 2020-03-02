#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 00:17:16 2020

@author: rfm
"""


import pickle
import sqlite3
from sqlite3 import Error


def create_database(db_file):
    """
    """
    
    conn = None
    try:
        # creates db file if it does not exist
        conn = sqlite3.connect(db_file)
    except Error as e:
        print(e)
    finally:
        if conn:
            conn.close()


create_database('/home/rfm/Desktop/rfm/Lab_Software/CreateSchema_tests/new_create_schema_scripts/campy_schema/campy_test/.sqlite_profiles/profiles.db')


# define database tables

# date format YYYY-MM-DDTHH:MM:SS

# samples table
sql_create_samples_table = ('CREATE TABLE IF NOT EXISTS samples ('
                            'id integer PRIMARY KEY,'
                            'name text NOT NULL,'
                            'date text NOT NULL,'
                            'profile text NOT NULL,'
                            'FOREIGN KEY (profile) REFERENCES profiles (id)'
                            ');')

# loci table
sql_create_loci_table = ('CREATE TABLE IF NOT EXISTS loci ('
                         'id text PRIMARY KEY,'
                         'annotation text'
                         ');')

# profiles table
genes_list_file = '/home/rfm/Desktop/rfm/Lab_Software/Chewie_NS/NS_tests/saureus_config_schema/.genes_list'
with open(genes_list_file, 'rb') as glf:
    genes_list = pickle.load(glf)

genes_list = [locus.rstrip('.fasta') for locus in genes_list]

genes_sql_columns = ['{0} text'.format(locus) for locus in genes_list]
genes_sql_columns_text = ','.join(genes_sql_columns)

sql_create_profiles = ('CREATE TABLE IF NOT EXISTS profiles ('
                       'id text PRIMARY KEY,'
                       '{0}'
                       ');'.format(genes_sql_columns))


def create_connection(db_file):
    """ create a database connection to the SQLite database
        specified by db_file
    :param db_file: database file
    :return: Connection object or None
    """
    conn = None
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except Error as e:
        print(e)
 
    return conn


def create_table(conn, create_table_sql):
    """ create a table from the create_table_sql statement
    :param conn: Connection object
    :param create_table_sql: a CREATE TABLE statement
    :return:
    """
    try:
        c = conn.cursor()
        c.execute(create_table_sql)
    except Error as e:
        print(e)


conn = create_connection('/home/rfm/Desktop/rfm/Lab_Software/CreateSchema_tests/new_create_schema_scripts/campy_schema/campy_test/.sqlite_profiles/profiles.db')

if conn is not None:
    create_table(conn, sql_create_profiles)
    create_table(conn, sql_create_loci_table)
    create_table(conn, sql_create_samples_table)
else:
    print('No database connection.')


# select all tables
con = sqlite3.connect('/home/rfm/Desktop/rfm/Lab_Software/CreateSchema_tests/new_create_schema_scripts/campy_schema/campy_test/.sqlite_profiles/profiles.db')
cursor = con.cursor()
cursor.execute("SELECT * FROM sqlite_master WHERE type='table';")
print(cursor.fetchall())
cursor.close()
con.close()

# insert value into loci








