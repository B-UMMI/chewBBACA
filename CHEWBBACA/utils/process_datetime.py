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


import datetime as dt


def datetime_str(datetime_obj, date_format='%Y-%m-%dT%H:%M:%S'):
	"""Convert datetime module object to formatted string.

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
	"""Convert formatted string to datetime module object.

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


def validate_date(date, date_format='%Y-%m-%dT%H:%M:%S.%f'):
	"""Check if date is in specified format.

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
