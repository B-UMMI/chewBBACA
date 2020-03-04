#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHOR

    Pedro Cerqueira
    github: @pedrorvc

    Rafael Mamede
    github: @rfm-targa

DESCRIPTION


    Important comments:
        - the SQLite included with Python distribution (3.29.0)
        might be old and will not enforce the FOREIGN KEY constraint
        (this was only implemented in SQLite 3.6.19). We have to be
        sensible of that when altering the data in the database.
"""


import csv
import pickle
import hashlib
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


def create_table(conn, sql_statement):
    """
    """

    try:
        c = conn.cursor()
        c.execute(sql_statement)
    except Error as e:
        print(e)


def execute_statement(conn, sql_statement, values):
    """ Creates a table from the create_table_sql statement

        Args:
            conn(): Connection object
            param create_table_sql(): a CREATE TABLE statement
        Returns:

    """
    try:
        c = conn.cursor()
        c.execute(sql_statement, values)
    except Error as e:
        print(e)

    return c.lastrowid


# test inserting values into database
def select_all_rows(db_file, table):
    """
    Query all rows in the tasks table
    :param conn: the Connection object
    :return:
    """
    
    conn = create_connection(db_file)
    cur = conn.cursor()
    cur.execute('SELECT * FROM {0}'.format(table))

    rows = cur.fetchall()

    for row in rows:
        print(row)
        
    conn.close()


def create_insert_statement(table, columns):
    """
    """

    sql_insert = ('INSERT OR IGNORE INTO {0}({1}) '
                  'VALUES({2});'.format(table,
                                       ','.join(columns),
                                       ','.join('?'*len(columns))))

    return sql_insert


def insert_row(conn, sql_statement, row_data):
    """
    """

    with conn:

        sample_data = row_data
        sample_sql = sql_statement
        cur = conn.cursor()
        #print(row_data)
        cur.execute(sample_sql, sample_data)

    return cur.lastrowid


def insert_loci(db_file, loci_list):
    """
    """

    conn = create_connection(db_file)
    locus_sql = create_insert_statement('loci', ['locus_id'])
    print(locus_sql)
    loci = [(locus,) for locus in loci_list]
    print(loci)
    cur = conn.cursor()
    cur.executemany(locus_sql, loci)

    conn.commit()
    conn.close()


def insert_multiple(db_file, base_statement, data):
    """
    """

    conn = create_connection(db_file)
    cur = conn.cursor()
    cur.executemany(base_statement, data)

    conn.commit()
    conn.close()


def create_database_structure(db_file):
    """
    """

    message = create_database_file(db_file)
    print(message)

    # samples table
    # date format YYYY-MM-DDTHH:MM:SS
    sql_samples_table = ('CREATE TABLE IF NOT EXISTS samples ('
                             'id INTEGER PRIMARY KEY,'
                             'name TEXT NOT NULL,'
                             'date TEXT NOT NULL,'
                             'profile_id TEXT NOT NULL,'
                             'FOREIGN KEY (profile_id) REFERENCES profiles (profile_id)'
                             ');')

    # loci table
    sql_loci_table = ('CREATE TABLE IF NOT EXISTS loci ('
                          'locus_id text PRIMARY KEY,'
                          'annotation text,'
                          'custom_annotation text'
                          ');')

    # profiles table
    sql_profiles_table = ('CREATE TABLE IF NOT EXISTS profiles ('
                              'profile_id text PRIMARY KEY,'
                              'type text'
                              ');')

    # profiles_alleles table
    sql_profiles_alleles_table = ('CREATE TABLE IF NOT EXISTS profiles_alleles ('
                                      'profile_id text,'
                                      'locus_id text,'
                                      'allele_id integer,'
                                      'FOREIGN KEY (profile_id) REFERENCES profiles (profile_id),'
                                      'FOREIGN KEY (locus_id) REFERENCES loci (locus_id),'
                                      'PRIMARY KEY (profile_id, locus_id)'
                                      ');')

    # create tables
    conn = create_connection(db_file)

    if conn is not None:
        row = create_table(conn, sql_samples_table)
        row = create_table(conn, sql_loci_table)
        row = create_table(conn, sql_profiles_table)
        row = create_table(conn, sql_profiles_alleles_table)
    else:
        print('No database connection.')


# database file
db_file = '/home/rfm/Desktop/rfm/Lab_Software/CreateSchema_tests/new_create_schema_scripts/campy_schema/campy_test/.profiles_sqlite/profiles.db'

# schema genes list
with open('/home/rfm/Desktop/rfm/Lab_Software/CreateSchema_tests/new_create_schema_scripts/campy_schema/campy_test/.genes_list', 'rb') as glf:
    genes_list = pickle.load(glf)

genes_list = [locus.rstrip('.fasta') for locus in genes_list]

# maximum number of columns is 2000
create_database_structure(db_file)

# insert all loci identifiers into loci table
insert_loci(db_file, genes_list)



# test inserting whole AlleleCall matrix
# it is safer to import final matrix from file and insert into database

# read matrix
matrix = '/home/rfm/Desktop/rfm/Lab_Software/results_alleles.tsv'
with open(matrix, 'r') as m:
    matrix_lines = list(csv.reader(m, delimiter='\t'))

# get genes identifiers
loci_ids = [locus.rstrip('.fasta') for locus in matrix_lines[0][1:]]

# get sample identifiers
sample_ids = [l[0] for l in matrix_lines[1:]]

# get profiles
profiles = [l[1:] for l in matrix_lines[1:]]

# construct base statements

# profiles
profiles_concat = ['.'.join(p) for p in profiles]
profiles_hashes = [hashlib.sha256(p.encode('utf-8')).hexdigest() for p in profiles_concat]
profile_statement = create_insert_statement('profiles', ['profile_id', 'type'])

profiles_data = [(phash, 'wgMLST') for phash in set(profiles_hashes)]

# insert all profiles - it checks PK uniqueness condition
insert_multiple(db_file, profile_statement, profiles_data)

# samples
sample_statement = create_insert_statement('samples', ['name', 'date', 'profile_id'])
samples_data = [(sample_ids[i], 'FAKEDATE', profiles_hashes[i]) for i in range(len(sample_ids))]

# insert all samples - it checks PK uniqueness condition
insert_multiple(db_file, sample_statement, samples_data)

# profiles_alleles
profileid_statement = create_insert_statement('profiles_alleles', ['profile_id', 'locus_id', 'allele_id'])

# construct every data tuple with the allele for each profile!
profileid_data = []
for i in range(len(profiles)):
    profileid_tuples = [(profiles_hashes[i], loci_ids[j], profiles[i][j]) for j in range(len(loci_ids))]
    profileid_data += profileid_tuples

# insert all profileid - it checks PK uniqueness condition
insert_multiple(db_file, profileid_statement, profileid_data)


# we need to check which profiles from our matrix already are in the database
# remove those from data and only insert the new ones
# also insert samples and that's it!


# select all rows from tables
select_all_rows(db_file, 'loci')
select_all_rows(db_file, 'profiles')
select_all_rows(db_file, 'samples')
select_all_rows(db_file, 'profiles_alleles')




