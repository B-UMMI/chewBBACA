#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains functions to perform requests to
Chewie-NS (https://github.com/B-UMMI/Chewie-NS).

Code documentation
------------------
"""


import sys
import json
import requests
from getpass import getpass
from urllib.parse import urlparse, urlencode, urlsplit, parse_qs

try:
    from utils import constants as ct
except ModuleNotFoundError:
    from CHEWBBACA.utils import constants as ct


def simple_get_request(url, headers, endpoint_list=None,
                       parameters=None, verify=False, timeout=30):
    """Request data through a GET method sent to an endpoint.

    Construct an URI for an endpoint and request data through
    a GET method sent to that endpoint.

    Parameters
    ----------
    url : str
        The base URL for a Chewie-NS instance.
    headers : dict
        HTTP headers for GET requests.
    endpoint_list : list
        List with elements that will be concatenated to the
        base URI to create the URI for the API endpoint.
    parameters : dict or list or bytes
        Dictionary, list of tuples or bytes to send in the
        query string for the Request.
    verify : bool
        If the SSL certificates should be verified in HTTPS
        requests (False for no verification, True otherwise).
    timeout : int
        Maximum time in seconds to wait for response.

    Returns
    -------
    url : str
        Complete URI for the API endpoint.
    res : requests.models.Response
        The response returned through the GET method.
    """
    # unpack list of sequential endpoints and pass to create URI
    if endpoint_list is not None:
        url = make_url(url, *endpoint_list)

    res = requests.get(url, headers=headers, params=parameters,
                       timeout=timeout, verify=verify)

    # return queried endpoint URI and response
    return [url, res]


def simple_post_request(url, headers, endpoint_list=None,
                        data=None, files=None, verify=False):
    """Send data through a POST method.

    Construct an URI for an endpoint and send data through
    a POST method.

    Parameters
    ----------
    url : str
        The base URL for a Chewie-NS instance.
    headers : dict
        HTTP headers for POST requests.
    endpoint_list : list
        List with elements that will be concatenated to the
        base URI to create the URI for the API endpoint.
    data
        Data to send. Must be an object that can be JSON
        serialized.
    files : dict
        Dictionary with file data to send (see:
        https://requests.readthedocs.io/en/master/user
        /quickstart/#post-a-multipart-encoded-file).
    verify : bool
        If the SSL certificates should be verified in HTTPS
        requests (False for no verification, True otherwise).

    Returns
    -------
    url : str
        Complete URI for the API endpoint.
    res : requests.models.Response
        The response returned through the POST method.
    """
    # unpack list of sequential endpoints and pass to create URI
    if endpoint_list is not None:
        url = make_url(url, *endpoint_list)

    # if the content is defined as 'application/json'
    # the data has to be passed after being processed with
    # json.dumps()
    res = requests.post(url, headers=headers, data=data,
                        files=files, verify=verify)

    return [url, res]


def check_connection(url, headers=ct.HEADERS_GET_JSON):
    """Verify connection to a Chewie-NS instance.

    Perform a GET request to the stats/summary endpoint of a
    Chewie-NS instance to determine if it is possible to
    establish a connection.

    Parameters
    ----------
    url : str
       The base URL for a Chewie-NS instance.
    headers : dict
        HTTP headers for GET requests.

    Returns
    -------
    connection : bool
        True if the GET request to the stats/summary
        endpoint was successful (HTTP response status
        code 200 or 201), False otherwise.
    """
    try:
        url, res = simple_get_request(url, headers, ['stats', 'summary'])
        server_status = res.status_code
        if server_status in [200, 201]:
            connection = True
        else:
            connection = False
    except Exception:
        connection = False

    return connection


def user_info(base_url, headers_get):
    """Get the information about a registered user in Chewie-NS.

    Get the user identifier, role and permission level for a
    user registered in a Chewie-NS instance. The permission
    level is used to determine if the user can upload schemas
    to Chewie-NS.

    Parameters
    ----------
    base_url : str
        The base URL for a Chewie-NS instance.
    headers_get : dict
        HTTP headers for GET requests.

    Returns
    -------
    user_id : str
        User identifier in the Chewie-NS instance.
    user_role : str
        User role in the Chewie-NS instance.
    permission : bool
        True if the user has Admin or Contributor privileges,
        False otherwise.
    """
    # verify user role to check permission
    user_url, user_info = simple_get_request(base_url, headers_get,
                                             ['user', 'current_user'])
    user_info = user_info.json()

    user_id = str(user_info['id'])
    user_role = user_info['roles'].split(' ')[-1][:-1]
    permission = any(role in user_role for role in ['Admin', 'Contributor'])

    return [user_id, user_role, permission]


def species_ids(species_id, base_url, headers_get):
    """Get a species identifier or name in Chewie-NS.

    Parameters
    ----------
    species_id : str
        The identifier of the species in Chewie-NS or the
        name of the species.
    base_url : str
        The base URL for a Chewie-NS instance.
    headers_get : dict
        HTTP headers for GET requests.

    Returns
    -------
    A list with the species identifier and name if the species
    exists in Chewie-NS or a 404 HTTP status code if the species
    was not found.
    """
    try:
        # check if species unique integer identifier was provided
        int(species_id)
        species_url, species_info = simple_get_request(base_url, headers_get,
                                                       ['species', species_id])
        if species_info.status_code in [200, 201]:
            species_name = species_info.json()[0]['name']['value']
            return [species_id, species_name]
        else:
            return 404
    except ValueError:
        # check if the scientific name for species was provided
        species_name = species_id
        # get the list of species added to Chewie-NS
        ns_species = species_list(base_url, headers_get, ['species', 'list'])
        # try to get the species unique integer identifier
        species_id = ns_species.get(species_name, 'not_found')
        if species_id != 'not_found':
            return [species_id, species_name]
        else:
            return 404


def species_list(base_url, headers_get, endpoint_list):
    """Get the list of species in a Chewie-NS instance.

    Parameters
    ----------
    base_url : str
        The base URL for a Chewie Nomenclature Server
        instance.
    headers_get : dict
        HTTP headers for GET requests.
    endpoint_list : list
        List with elements that will be concatenated
        to the base URI to create the URI for the API endpoint.

    Returns
    -------
    species_lst : dict
        Dictionary with species names as keys and
        species identifiers as values.
    """
    species_url, res = simple_get_request(base_url, headers_get, endpoint_list)
    res = res.json()
    species_lst = {}
    for sp in res:
        # get species scientific name
        species = sp['name']['value']
        # get species unique integer identifier from URl
        species_url = sp['species']['value']
        species_id = species_url.split('/')[-1]
        species_lst[species] = species_id

    return species_lst


def is_url(url):
    """Check if a URL is valid.

    Parameters
    ----------
    url : str
        The URL to be validated.

    Returns
    -------
    True if the URL is valid, False otherwise.
    """
    try:
        result = urlparse(url)
        return all([result.scheme, result.netloc, result.path])
    except Exception:
        return False


def make_url(base_url, *res, **params):
    """Create a URL.

    Parameters
    ----------
    base_url : str
        The base URL.
    res : list
        List with elements that will be concatenated to the
        base URI to create the URI for the API endpoint.
    params : str
        Addtional parameters (WIP).

    Returns
    -------
    url : str
        URL constructed by adding endpoints and additional
        parameters. If the provided base URL is invalid
        returns 'An invalid URL was provided.'
    """
    url = base_url
    # Check if the url is valid
    if is_url(url):

        if url[-1] == '/':
            url = url[:-1]

        # Add the endpoints
        for r in res:
            url = '{0}/{1}'.format(url, r)

        # Add params if they are provided
        if params:
            url = '{0}?{1}'.format(url, urlencode(params))

        return url
    else:
        return 'An invalid URL was provided.'


def get_sequence_from_url(url):
    """Extract sequence from URL query component.

    Parameters
    ----------
    url : str
        URL that contains the sequence.

    Returns
    -------
    seq : str
        Sequence extracted from the URL.
    """
    seq = parse_qs(urlsplit(url).query)['sequence'][0]

    return seq


def login_user_to_NS(server_url, email, password):
    """Login user to a Chewie-NS instance.

    Parameters
    ----------
    server_url : str
        The base URL for a Chewie-NS instance.
    email : str
        Email linked to user's account in Chewie-NS.
    password : str
        Password for the user's account in Chewie-NS.

    Returns
    -------
    token : str
        Authorization token to perform requests to Chewie-NS
        instance if login is successful. Otherwise, returns
        False.
    """
    auth_params = {}
    auth_params['email'] = email
    auth_params['password'] = password

    auth_url, auth_r = simple_post_request(server_url,
                                           ct.HEADERS_POST_JSON,
                                           ['auth', 'login'],
                                           data=json.dumps(auth_params))

    auth_result = auth_r.json()
    if auth_result['status'] == 'success':
        token = auth_result['access_token']
    else:
        token = False

    return token


def capture_login_credentials(base_url):
    """Capture login credentials to login user into Chewie-NS.

    Parameters
    ----------
    base_url : str
        The base URL for a Chewie-NS instance.

    Returns
    -------
    token : str
        Authorization token to perform requests to Chewie-NS.
    """
    print('\nPlease provide login credentials:')
    user = input('USERNAME: ')
    password = getpass('PASSWORD: ')
    print()
    # get token
    token = login_user_to_NS(base_url, user, password)
    # if login was not successful, stop the program
    if token is False:
        sys.exit('Invalid credentials.')

    return token


def get_species_schemas(schema_id, species_id, base_url, headers_get):
    """Determine if there is a schema with specified identifier.

    Retrieves the list of schemas for a species in a Chewie-NS
    instance and determines if it has a schema with the
    specified identifier.

    Parameters
    ----------
    schema_id : str
        The identifier of the schema in Chewie-NS.
    species_id : str
        The identifier of the schema's species in Chewie-NS.
    base_url : str
        The base URL for the Chewie-NS instance.
    headers_get : dict
        HTTP headers for GET requests.

    Returns
    -------
    A list with the schema identifier, the schema URI, the
    schema name and a dictionary with all properties for
    the schema in Chewie-NS.

    Raises
    ------
    SystemExit
        - If the schema with the specified identifier does
          not exist.
        - If the process cannot retrieve the list of schemas
          for the species.
    """
    # get the list of schemas for the species
    schemas_url, schema_get = simple_get_request(base_url, headers_get,
                                                 ['species', species_id,
                                                  'schemas'])
    schema_get_status = schema_get.status_code
    if schema_get_status in [200, 201]:
        species_schemas = schema_get.json()

        # extract schemas identifiers, URIs and names from response
        schemas_info = []
        for s in species_schemas:
            schema_uri = s['schemas']['value']
            scid = schema_uri.split('/')[-1]
            schema_name = s['name']['value']
            schemas_info.append([scid, schema_uri, schema_name])

        # select schema with specified identifier
        schema = [s for s in schemas_info if schema_id in s]
        if len(schema) > 0:
            # get schema parameters
            schema = schema[0]
            url, schema_params = simple_get_request(schema[1], headers_get)
            schema_params = schema_params.json()[0]
            schema.append(schema_params)
            return schema
        else:
            sys.exit('Could not find a schema with provided identifier.')
    else:
        sys.exit('Could not retrieve schemas for current species.')


def upload_file(file, filename, url, headers, verify_ssl):
    """Upload a file to a Chewie-NS instance.

    Parameters
    ----------
    file : str
        Path to the file to upload.
    filename : str
        Name used to save the file in Chewie-NS.
    url : str
        Endpoint URL that receives the POST request.
    headers : dict
        HTTP POST request headers.
    verify_sll : bool
        If the SSL certificates should be verified in
        HTTPS requests (False for no verification, True
        otherwise).

    Returns
    -------
    response : requests.models.Response
        Response object from the 'requests' module.
    """
    file_handle = open(file, 'rb')
    files = {'file': (filename, file_handle)}
    up_url, response = simple_post_request(url, headers,
                                           files=files,
                                           verify=verify_ssl)

    file_handle.close()

    return response
