#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHOR

    Pedro Cerqueira
    github: @pedrorvc

    Rafael Mamede
    github: @rfm-targa

DESCRIPTION

"""


import pickle
import sqlite3
from sqlite3 import Error


def create_database_file(db_file):
    """
    """

    conn = None
    error = None
    try:
        # creates db file if it does not exist
        conn = sqlite3.connect(db_file)
    except Error as e:
        error = e
    finally:
        if conn:
            conn.close()

    return error


def create_connection(db_file):
    """ Creates a database connection to the SQLite database
        specified by db_file

        Args:
            param db_file: database file
        Returns:
            Connection object or None
    """
    conn = None
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except Error as e:
        print(e)

    return conn


def execute_statement(conn, create_table_sql, values):
    """ Creates a table from the create_table_sql statement

        Args:
            conn(): Connection object
            param create_table_sql(): a CREATE TABLE statement
        Returns:

    """
    try:
        c = conn.cursor()
        c.execute(create_table_sql, values)
    except Error as e:
        print(e)
        
    return c.lastrowid


# main
def create_database_structure(db_file, genes_list):
    """
    """

    message = create_database_file(db_file)
    print(message)

    # samples table
    # date format YYYY-MM-DDTHH:MM:SS
    sql_samples_table = ('CREATE TABLE IF NOT EXISTS samples ('
                             'id integer PRIMARY KEY,'
                             'name text NOT NULL,'
                             'date text NOT NULL,'
                             'profile text NOT NULL,'
                             'FOREIGN KEY (profile) REFERENCES profiles (id)'
                             ');')

    # loci table
    sql_loci_table = ('CREATE TABLE IF NOT EXISTS loci ('
                          'id text PRIMARY KEY,'
                          'annotation text,'
                          'custom_annotation text'
                          ');')
    
    # profiles table
    sql_profiles_table = ('CREATE TABLE IF NOT EXISTS profiles ('
                              'id text PRIMARY KEY,'
                              'type text'
                              ');')

    # profiles_alleles table
    sql_profiles_alleles_table = ('CREATE TABLE IF NOT EXISTS profiles_alleles ('
                                      'profile_id text,'
                                      'locus_id text,'
                                      'allele_id integer,'
                                      'FOREIGN KEY (profile_id) REFERENCES profiles (id),'
                                      'FOREIGN KEY (locus_id) REFERENCES loci (id),'
                                      'PRIMARY KEY (profile_id, locus_id)'
                                      ');')

    # create tables
    conn = create_connection(db_file)

    if conn is not None:
        row = execute_statement(conn, sql_samples_table)
        row = execute_statement(conn, sql_loci_table)
        row = execute_statement(conn, sql_profiles_table)
        row = execute_statement(conn, sql_profiles_alleles_table)
    else:
        print('No database connection.')

# database file
db_file = '/home/rfm/Desktop/rfm/Lab_Software/CreateSchema_tests/new_create_schema_scripts/campy_schema/campy_test/.profiles_sqlite/profiles.db'

# schema genes list
with open('/home/rfm/Desktop/rfm/Lab_Software/CreateSchema_tests/new_create_schema_scripts/campy_schema/campy_test/.genes_list', 'rb') as glf:
    genes_list = pickle.load(glf)

genes_list = [locus.rstrip('.fasta') for locus in genes_list]

# maximum number of columns is 2000
create_database_structure(db_file, genes_list)




def add_sample(connection, name, date, profile):
    """
    """

    insert_statement = ('INSERT INTO samples (name, date, profile)'
                        'VALUES(?, ?, ?);')

    row = execute_statement(conn, insert_statement, (name, date, profile))
    
    return row

conn = create_connection(db_file)
a = add_sample(conn, 'GIJO3', '20200131T042000', 'FAKEHASH')

def add_locus(locus_id, annotation=None, custom_annotation=None):
    """
    """

    pass


def add_profile(profile_id, profile_type):
    """
    """

    pass


def add_profile_allele(profile_id, locus_id, allele_id):
    """
    """








