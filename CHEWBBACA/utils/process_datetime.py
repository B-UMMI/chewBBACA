#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains functions related to time measurement
and date formatting.

Code documentation
------------------
"""


import sys
import functools
import datetime as dt

try:
    from utils import parameters_validation as pv
except:
    from CHEWBBACA.utils import parameters_validation as pv


def get_datetime():
    """ Returns datetime module object with
        information about current date and hour.
    """

    current_datetime = dt.datetime.now()

    return current_datetime


def datetime_str(datetime_obj, date_format='%Y-%m-%dT%H:%M:%S'):
    """ Converts datetime module object to formatted string.

        Parameters
        ----------
        datetime_obj : datetime.datetime
            Datetime module object.
        date_format : str
            Format for string representation of the date
            object.

        Returns
        -------
        dt_str : str
            String representation of the date according
            to specified format.
    """

    dt_str = dt.datetime.strftime(datetime_obj, date_format)

    return dt_str


def datetime_obj(datetime_str, date_format='%Y-%m-%dT%H:%M:%S'):
    """ Converts formatted string to datetime module object.

        Parameters
        ----------
        datetime_str : str
            String to convert to datetime object.
        date_format : str
            Format of the string representation of the date
            object.

        Returns
        -------
        dt_obj : datetime.datetime
            Datetime module object.
    """

    dt_obj = dt.datetime.strptime(datetime_str, date_format)

    return dt_obj


def datetime_diff(sdate, edate):
    """ Returns the difference in minutes and the
        remaining in seconds between two dates.

        Parameters
        ----------
        sdate : datetime.datetime
            Datetime module object corresponding to
            the oldest date.
        edate : datetime.datetime
            Datetime module object corresponding to
            the most recent date.

        Returns
        -------
        A list with the following elements:
            minutes : float
                Time difference between the two dates
                in minutes.
            seconds : float
                The remaining time difference in seconds.
    """

    delta = edate - sdate
    minutes, seconds = divmod(delta.total_seconds(), 60)

    return [minutes, seconds]


def validate_date(date, date_format='%Y-%m-%dT%H:%M:%S.%f'):
    """ Checks if date is in specified format.

        Parameters
        ----------
        date : str
            String representing a date.
        date_format : str
            Date format.

        Returns
        -------
        valid : str
            The string representing the date if it
            is in the specified format.
    """

    valid = False
    try:
        date = dt.datetime.strptime(date, date_format)
        valid = date
    except ValueError:
        date = dt.datetime.strptime(date+'.0', date_format)
        valid = date

    return valid


# decorator to time main processes
def process_timer(func):
    # use functools to preserve info about wrapped function
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        # get process name and print header
        pv.process_header(sys.argv[1])

        start = get_datetime()
        start_str = datetime_str(start)
        print('Started at: {0}\n'.format(start_str))
        # run function
        func(*args, **kwargs)

        end = get_datetime()
        end_str = datetime_str(end)
        print('\nFinished at: {0}'.format(end_str))

        minutes, seconds = datetime_diff(start, end)
        print('Took {0: .0f}m{1: .0f}s.'.format(minutes, seconds))

    return wrapper
