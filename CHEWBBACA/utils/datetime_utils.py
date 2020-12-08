#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


DESCRIPTION

"""


import time
import datetime as dt


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
            String to conver to datetime object.
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


def validate_date(date):
    """ 
    """

    valid = False
    try:
        date = dt.datetime.strptime(date, '%Y-%m-%dT%H:%M:%S.%f')
        valid = date
    except ValueError:
        date = dt.datetime.strptime(date+'.0', '%Y-%m-%dT%H:%M:%S.%f')
        valid = date

    return valid
