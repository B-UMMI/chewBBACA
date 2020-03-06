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

        - json_extract is SQL function from JSON1 extension that is
        included in the SQL statements. It can be used to extract single
        values from a full JSON that has been inserted into a single table
        field/cell.
        - Changing special cases to shorter strings probably won't help
        because SQL has a fixed size for each table cell according to data
        type, and in the best scenario based on text length. Profiles are really
        big strings so it might be hard to lower space disk usage...
        
        TRY TO REMOVE LOCI THAT ARE LNF TO CONSIDERABLY REDUCE STRING LENGTH???
        
"""


import csv
import json
import pickle
import hashlib
import sqlite3
import sqlalchemy
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

    rows_list = []
    for row in rows:
        rows_list.append(row)

    conn.close()

    return rows_list


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
        cur.execute(sample_sql, sample_data)

    return cur.lastrowid


def insert_loci(db_file, matrix_file):
    """
    """

    matrix_lines = read_matrix(matrix_file)
    loci_list = [locus.rstrip('.fasta') for locus in matrix_lines[0][1:]]

    conn = create_connection(db_file)
    locus_sql = create_insert_statement('loci', ['name'])

    loci = [(locus,) for locus in loci_list]

    cur = conn.cursor()
    cur.executemany(locus_sql, loci)

    conn.commit()
    conn.close()

    print('Inserted {0} loci.'.format(len(loci_list)))


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
                          'id INTEGER PRIMARY KEY,'
                          'name TEXT,'
                          'annotation TEXT,'
                          'custom_annotation TEXT'
                          ');')

    # profiles table
    sql_profiles_table = ('CREATE TABLE IF NOT EXISTS profiles ('
                              'profile_id TEXT PRIMARY KEY,'
                              'type TEXT,'
                              'date TEXT,'
                              'profile_json JSON'
                              ');')

    # create tables
    conn = create_connection(db_file)

    if conn is not None:
        row = create_table(conn, sql_samples_table)
        row = create_table(conn, sql_loci_table)
        row = create_table(conn, sql_profiles_table)
    else:
        print('No database connection.')

    conn.commit()
    conn.close()


def read_matrix(matrix_file):
    """
    """

    with open(matrix_file, 'r') as m:
        matrix_lines = list(csv.reader(m, delimiter='\t'))

    return matrix_lines


def get_loci_ids(matrix_lines):
    """
    """

    loci_ids = [locus.rstrip('.fasta') for locus in matrix_lines[0][1:]]

    return loci_ids


def get_sample_ids(matrix_lines):
    """
    """

    sample_ids = [l[0] for l in matrix_lines[1:]]

    return sample_ids


def get_profiles(matrix_lines):
    """
    """

    profiles = [l[1:] for l in matrix_lines[1:]]

    return profiles


def profile_json(loci_ids, profile):
    """
    """

    json_profile = '{{"{0}":"{1}"'.format(loci_ids[0], profile[0])
    for l in loci_ids[1:]:
        json_profile += ', "{0}":"{1}"'.format(l, 1)
    json_profile += '}'

    return json_profile


def remove_inf(profile):
    """
    """

    new_profile = [a.lstrip('INF-') if 'INF-' in a else a for a in profile]

    return new_profile


def insert_allelecall_matrix(matrix_file, db_file):
    """
    """

    loci_list_db = select_all_rows(db_file, 'loci')
    loci_map = {t[1]: t[0] for t in loci_list_db}

    # read matrix
    matrix_lines = read_matrix(matrix_file)

    # get sample identifiers
    sample_ids = get_sample_ids(matrix_lines)

    # get profiles
    profiles = get_profiles(matrix_lines)

    profiles = [remove_inf(p) for p in profiles]

    # insert profiles
    # create JSON format of each profile
    conn = sqlite3.connect(db_file)
    # create a cursor
    c = conn.cursor()
    profiles_hashes = []
    for p in profiles:
        # add first entry to JSON only if locus value is not LNF
        i = 1
        not_lnf = False
        while not_lnf is False:
            if p[i-1] != 'LNF':
                json_profile = '{{"{0}":"{1}"'.format(i, p[i-1])
                not_lnf = True
            else:
                i += 1
        # start adding entries to JSON from index after non-LNF value
        for l in range(i, len(loci_map)):
            # only had if value is not LNF
            if p[l] != 'LNF':
                json_profile += ', "{0}":"{1}"'.format(l+1, p[l])
        json_profile += '}'

        profile_hash = hashlib.sha256(json_profile.encode('utf-8')).hexdigest()
        profiles_hashes.append(profile_hash)

        c.execute("INSERT OR IGNORE INTO profiles (profile_id, profile_json) VALUES (?, json(?));", (profile_hash, json_profile))

    conn.commit()
    conn.close()

    # insert samples
    sample_statement = create_insert_statement('samples', ['name', 'date', 'profile_id'])
    samples_data = [(sample_ids[i], 'FAKEDATE', profiles_hashes[i]) for i in range(len(sample_ids))]

    # insert all samples - it checks PK uniqueness condition
    insert_multiple(db_file, sample_statement, samples_data)

    return True


# database file
#db_file = '/home/rfm/Desktop/rfm/Lab_Software/CreateSchema_tests/new_create_schema_scripts/campy_schema/campy_test/.profiles_sqlite/profiles.db'
### maximum number of columns is 2000
#create_database_structure(db_file)
##
#matrix_file = '/home/rfm/Desktop/rfm/Lab_Analyses/Bacgentrack_data/7676_sagalactiae_allelecall_results/results_alleles_ns_version.tsv'
##
### insert all loci identifiers into loci table
#insert_loci(db_file, matrix_file)
#
## test inserting whole AlleleCall matrix
## it is safer to import final matrix from file and insert into database
#a = insert_allelecall_matrix(matrix_file, db_file)
#
#
#
#
#
## we need to check which profiles from our matrix already are in the database
## remove those from data and only insert the new ones
## also insert samples and that's it!
#
#
## select all rows from tables
#db_file = '/home/rfm/Desktop/ns_test/sagalactiae_test_schema/profiles_database/profiles.db'
#loci_list_db = select_all_rows(db_file, 'loci')
#profiles_list_db = select_all_rows(db_file, 'profiles')
#samples_list_db = select_all_rows(db_file, 'samples')
#
## select single field of profiles data
#conn = sqlite3.connect(db_file)
#
## create a cursor
#c = conn.cursor()
#
## json_extract is SQL function from JSON1 extension
#c.execute("select json_extract(profiles.profile_json, '$.1') from profiles;")
#
#rows = c.fetchall()
#
#rows_list = []
#for row in rows:
#    rows_list.append(row)
#
#conn.close()

