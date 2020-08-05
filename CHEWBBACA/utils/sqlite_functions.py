#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module

Notes
-----

The SQLite included with Python distribution (3.29.0)
might be old and will not enforce the FOREIGN KEY constraint
(this was only implemented in SQLite 3.6.19). We need to take
that into account when altering the data in the database.

"""


import os
import csv
import hashlib
import sqlite3


def create_database_file(db_file):
    """ Creates a SQLite database file.
        If the database file already exists,
        it will establish and close connection.

        Parameters
        ----------
        df_file : str
            Path to the SQLite database file.

        Returns
        -------
        error : None or sqlite3.OperationalError
            None if the SQLite database file was
            successfully created, OperationalError
            if it could not create/establish connection
    """

    conn = None
    error = None
    try:
        # creates db file if it does not exist
        conn = sqlite3.connect(db_file)
    except Exception as e:
        error = e
    finally:
        if conn:
            conn.close()

    return error


def create_connection(db_file):
    """ Creates a database connection to a SQLite
        database.

        Parameters
        ----------
        db_file: str
            Path to the SQLite database file.

        Returns
        -------
        conn : sqlite3.Connection or sqlite3.OperationalError
            SQLite Connection object if connection was
            successfull or error if it was not possible
            to connect to the database.
    """

    try:
        conn = sqlite3.connect(db_file)
    except Exception as e:
        conn = e

    return conn


def execute_statement(conn, statement):
    """ Executes a SQL statement.

        Parameters
        ----------
        conn : sqlite3.Connection
            SQLite Connection object.
        statement : str
            SQL statement to execute.

        Returns
        -------
        error : None or sqlite3.OperationalError
            None if the SQLite database file was
            successfully created, OperationalError
            if it could not create/establish connection
    """

    error = None
    try:
        c = conn.cursor()
        c.execute(statement)
        return c
    except Exception as e:
        error = e
        return error


def select_all_rows(db_file, table):
    """ Retrieves all rows in a table.

        Parameters
        ----------
        db_file : str
            Path to the SQLite database file.
        table : str
            Name of the table.

        Returns
        -------
        rows : list of tup
            List of all rows in the table. Each row
            is represented by a tuple with the values
            for all columns.
    """

    conn = create_connection(db_file)
    cur = conn.cursor()
    cur.execute('SELECT * FROM {0}'.format(table))

    result = cur.fetchall()

    rows = [r for r in result]

    conn.close()

    return rows


def create_insert_statement(table, columns, placeholders):
    """ Creates a base SQL insert statement.

        Parameters
        ----------
        table : str
            Name of the table.
        columns : list
            List with the names of the columns
            that values will be inserted into.

        Returns
        -------
        statement : str
            SQL insert statement that can be
            used to insert values into `columns`
            of a `table`.
    """

    statement = ('INSERT OR IGNORE INTO {0}({1}) '
                 'VALUES({2});'.format(table, ','.join(columns),
                                       ','.join(placeholders)))

    return statement


def insert_loci(db_file, matrix_file):
    """ Inserts loci into the loci table.

        Parameters
        ----------
        db_file : str
            Path to the SQLite database file.
        matrix_file : str
            Path to the TSV file with a matrix
            of allelic profiles.

        Returns
        -------
        The number of loci that were insert
        into the table.
    """

    matrix_lines = read_matrix(matrix_file)
    loci_list = [locus.rstrip('.fasta') for locus in matrix_lines[0][1:]]

    conn = create_connection(db_file)
    locus_sql = create_insert_statement('loci', ['name'],
                                        ['?'])

    loci = [(locus,) for locus in loci_list]

    cur = conn.cursor()
    cur.executemany(locus_sql, loci)

    conn.commit()
    conn.close()

    return len(loci_list)


def insert_multiple(db_file, base_statement, data):
    """ Executes several insert statements.

        Parameters
        ----------
        df_file : str
            Path to the SQLite database file.
        base_statement : str
            Base SQL insert statement to execute.
        data : list of tup
            A list with tuples that contain the
            column values to insert for each row.

        Returns
        -------
        error : None or sqlite3.OperationalError
            None if the SQL statement was successfully
            inserted, SQLite OperationalError otherwise.
    """

    error = None
    try:
        conn = create_connection(db_file)
        cur = conn.cursor()
        cur.executemany(base_statement, data)
        conn.commit()
        conn.close()
    except Exception as e:
        error = e

    return error


def create_database(db_file):
    """ Creates the database file and tables of a SQLite database
        that will store the allelic profiles determined with
        a schema.

        Parameters
        ----------
        db_file : str
            Path to the SQLite database file.

        Returns
        -------
        True if the SQLite database file and tables were
        successfully created, SQLite OperationalError otherwise.
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
                          'profile_json JSON,'
                          'subschema_id TEXT,'
                          'FOREIGN KEY (subschema_id) REFERENCES subschemas (subschema_id)'
                          ');')

    # subschemas table
    sql_subschemas_table = ('CREATE TABLE IF NOT EXISTS subschemas ('
                            'subschema_id TEXT PRIMARY KEY,'
                            'loci JSON'  # JSON with schema subset
                            ');')

    # create tables
    conn = create_connection(db_file)

    if isinstance(conn, str) is True:
        return conn
    else:
        row = execute_statement(conn, sql_samples_table)
        row = execute_statement(conn, sql_loci_table)
        row = execute_statement(conn, sql_profiles_table)
        row = execute_statement(conn, sql_subschemas_table)
        conn.commit()
        conn.close()

        return True


def read_matrix(matrix_file):
    """ Reads a TSV file that contains a matrix with
        allelic profiles.

        Parameters
        ----------
        matrix_file : str
            Path to the TSV file with the matrix.

        Returns
        -------
        matrix_lines : list of list
            A list with all the lines in the TSV file.
    """

    with open(matrix_file, 'r') as m:
        matrix_lines = list(csv.reader(m, delimiter='\t'))

    return matrix_lines


def get_loci_ids(matrix_lines):
    """ Extracts loci identifiers from a list
        with lines from a matrix of allelic profiles.

        Parameters
        ----------
        matrix_lines : list of list
            A list with all the lines read from a TSV
            file that contains a matrix of allelic profiles.

        Returns
        -------
        loci_ids : list
            List with the identifiers of all loci represented
            in the allelic profiles.
    """

    loci_ids = [locus.rstrip('.fasta') for locus in matrix_lines[0][1:]]

    return loci_ids


def get_sample_ids(matrix_lines):
    """ Extracts sample identifiers from a list
        with lines from a matrix of allelic profiles.

        Parameters
        ----------
        matrix_lines : list of list
            A list with all the lines read from a TSV
            file that contains a matrix of allelic profiles.

        Returns
        -------
        sample_ids : list
            List with the sample identifiers of all allelic
            profiles.
    """

    sample_ids = [l[0].rstrip('.fasta') for l in matrix_lines[1:]]

    return sample_ids


def get_profiles(matrix_lines):
    """ Extracts profiles from a list with lines from
        a matrix of allelic profiles.

        Parameters
        ----------
        matrix_lines : list of list
            A list with all the lines read from a TSV
            file that contains a matrix of allelic profiles.

        Returns
        -------
        profiles : list of dict
            List with one dictionary per allelic profile.
    """

    profiles = []
    loci_ids = matrix_lines[0][1:]
    for l in matrix_lines[1:]:
        lp = l[1:]
        profile = {loci_ids[i].rstrip('.fasta'): lp[i] for i in range(len(lp))}
        profiles.append(profile)

    return profiles


def remove_inf(profile):
    """ Remove 'INF-' prefix from inferred alleles.

        Parameters
        ----------
        profile : list
            List with the allele identifiers in an allelic
            profile.

        Returns
        -------
        clean_profile : list
            List with allele identifiers stripped of the
            'INF-' prefix.
    """

    clean_profile = [a.lstrip('INF-') if 'INF-' in a else a for a in profile]

    return clean_profile


def jsonify_profile(profile, loci):
    """
    """

    json_profile = ''
    for k, v in profile.items():
        # add first entry to JSON only if locus value is not LNF
        locus = k
        locus_id = loci[locus]
        allele_id = v
        if len(json_profile) == 0:
            if allele_id != 'LNF':
                json_profile += '{{"{0}":"{1}"'.format(locus_id, allele_id)
        else:
            if allele_id != 'LNF':
                json_profile += ', "{0}":"{1}"'.format(locus_id, allele_id)

    json_profile += '}'

    return json_profile


def insert_allelecall_matrix(matrix_file, db_file, insert_date):
    """ Inserts the data contained in a AlleleCall matrix into
        the SQLite database of the schema.

        Parameters
        ----------
        matrix_file : str
            Path to the TSV file with the matrix.
        db_file : str
            Path to the SQLite database file.
        insert_date : str
            Date in the name of the folder created
            by the AlleleCall process to save results
            (in the format Y-m-dTH:M:S).

        Returns
        -------
        A list with the following elements:

        - Number of inserted profiles.
        - Total number of profiles.
        - Number of unique profiles.
    """

    loci_list_db = select_all_rows(db_file, 'loci')
    loci_map = {t[1]: t[0] for t in loci_list_db}

    # read matrix
    matrix_lines = read_matrix(matrix_file)
    matrix_lines = [remove_inf(l) for l in matrix_lines]

    # get sample identifiers
    sample_ids = get_sample_ids(matrix_lines)

    # get profiles
    profiles = get_profiles(matrix_lines)

    # insert profiles
    # create JSON format of each profile
    profiles_hashes = []
    subschemas_hashes = []
    subschemas_loci = []
    inserted_profiles = 0
    profile_data = []
    for p in profiles:
        loci = [str(loci_map[locus]) for locus in p.keys()]
        loci_join = ','.join(loci)
        loci_hash = hashlib.sha256(loci_join.encode('utf-8')).hexdigest()
        if loci_hash not in subschemas_hashes:
            subschemas_loci.append(loci_join)
            subschemas_hashes.append(loci_hash)

        json_profile = jsonify_profile(p, loci_map)

        profile_hash = hashlib.sha256(json_profile.encode('utf-8')).hexdigest()
        profiles_hashes.append(profile_hash)

        profile_data.append((profile_hash, insert_date,
                             json_profile, loci_hash))

    # insert all profiles
    conn = create_connection(db_file)
    profile_statement = create_insert_statement('profiles',
                                                ['profile_id', 'date',
                                                 'profile_json', 'subschema_id'],
                                                ['?', '?', 'json(?)', '?'])
    previous_rcount = execute_statement(conn, "SELECT COUNT(*) "
                                              "FROM profiles;").fetchone()[0]
    profiles_res = insert_multiple(db_file, profile_statement, profile_data)

    posterior_rcount = execute_statement(conn, "SELECT COUNT(*) "
                                               "FROM profiles;").fetchone()[0]
    inserted_profiles = posterior_rcount - previous_rcount
    conn.close()

    # insert subschemas
    subschema_statement = create_insert_statement('subschemas',
                                                  ['subschema_id', 'loci'],
                                                  ['?', '?'])
    subschema_data = [(subschemas_hashes[i], subschemas_loci[i])
                      for i in range(len(subschemas_hashes))]

    # insert all subschemas
    subschema_res = insert_multiple(db_file, subschema_statement,
                                    subschema_data)

    # insert samples
    sample_statement = create_insert_statement('samples',
                                               ['name', 'date', 'profile_id'],
                                               ['?', '?', '?'])
    samples_data = [(sample_ids[i], insert_date, profiles_hashes[i])
                    for i in range(len(sample_ids))]

    # insert all samples - it checks PK uniqueness condition
    samples_res = insert_multiple(db_file, sample_statement, samples_data)

    return [inserted_profiles, len(profiles), len(set(profiles_hashes)),
            profiles_res, subschema_res, samples_res]


def select_outdated(loci, reassigned, cursor):
    """ Retrives the allelic profiles that have outdated
        allele identifiers.

        Parameters
        ----------
        loci : dict
            A dictionary with the numeric identifiers
            of loci as keys and the integer identifier
            of each locus in the local SQLite database.
        reassigned : dict
            A dictionary with loci identifiers as keys.
            Each key has a dictionary as value, with the
            current alleles identifiers as keys and the
            updated identifiers as values.
        cursor : sqlite3.Cursor
            SQLite cursor object.

        Returns
        -------
        profiles : dict of list
            A dictionary with profiles hashes as keys and
            lists as values. Each list contains the profile
            as str type and a variable number of tuples,
            one tuple per allele identifier that must be
            updated (tuples contain the locus identifier,
            the outdated allele identifier and the updated
            allele identifier).
    """

    profiles = {}
    for locus, alleles in reassigned.items():
        locus_id = loci[locus.split('-')[-1].rstrip('.fasta').lstrip('0')]
        for a1, a2 in alleles.items():
            # good idea to implement way to run many queries at once
            query = ("SELECT profiles.profile_id, profiles.profile_json "
                     "FROM profiles "
                     "WHERE json_extract(profiles.profile_json, '$.{0}') = '{1}';".format(locus_id, a1))
            cursor.execute(query)
            rows = [r for r in cursor.fetchall()]
            for r in rows:
                profile_hash = r[0]
                profile = r[1]
                if profile_hash in profiles:
                    profiles[profile_hash].append((locus_id, a1, a2))
                else:
                    profiles[profile_hash] = [profile, (locus_id, a1, a2)]

    return profiles


def alter_profiles(profiles, cursor):
    """ Alters allele identifiers in allelic profiles
        that are outdated.

        Parameters
        ----------
        profiles : dict
            A dictionary with profiles hashes as keys and
            lists as values. Each list contains the profile
            as str type and a variable number of tuples,
            one tuple per allele identifier that must be
            updated (tuples contain the locus identifier,
            the outdated allele identifier and the updated
            allele identifier).
        cursor : sqlite3.Cursor
            SQLite cursor object.

        Returns
        -------
        results : dict of str
            A dictionary with profiles hashes as keys and
            updated profiles as values.
    """

    results = {}
    for k, v in profiles.items():
        profile = v[0]
        to_replace = ["'$.{0}', '{1}'".format(r[0], r[2]) for r in v[1:]]
        # json_replace cannot receive a great number of arguments
        # change 50 identifiers at a time
        for i in range(0, len(to_replace), 50):
            to_replace_txt = ','.join(to_replace[i:i+50])
            query = ("SELECT json_replace('{0}', {1}) "
                     "FROM profiles;".format(profile, to_replace_txt))
            cursor.execute(query)
            res = cursor.fetchall()
            profile = res[0][0]
        results[k] = profile

    return results


def update_profiles(schema_directory, reassigned):
    """ Updates allele identifiers that have been changed.

        Parameters
        ----------
        schema_directory : str
            Path to the directory with the schema files.
        reassigned : dict of dict
            A dictionary with loci identifiers as keys.
            Each key has a dictionary as value, with the
            current alleles identifiers as keys and the
            updated identifiers as values.

        Returns
        -------
        The number of profiles that were altered.
    """

    db_file = os.path.join(schema_directory, 'profiles_database',
                           'profiles.db')

    if os.path.isfile(db_file) is False:
        return None

    # create connection to db
    conn = create_connection(db_file)
    cursor = conn.cursor()
    # get list of loci
    loci_list = select_all_rows(db_file, 'loci')
    loci = {l[1].split('-')[-1].lstrip('0'): l[0] for l in loci_list}
    # get profiles with identifiers that have to be changed
    profiles = select_outdated(loci, reassigned, cursor)

    # change identifiers in profiles
    updated_profiles = alter_profiles(profiles, cursor)

    conn.close()

    # change values in profiles table
    query = "update profiles set profile_json = ? where profile_id = ?;"
    params = [[p, h] for h, p in updated_profiles.items()]
    insert_multiple(db_file, query, params)

    return len(profiles)


# select all rows from tables
#db_file = ''
#loci_list_db = select_all_rows(db_file, 'loci')
#profiles_list_db = select_all_rows(db_file, 'profiles')
#samples_list_db = select_all_rows(db_file, 'samples')
#subschemas_list_db = select_all_rows(db_file, 'subschemas')
