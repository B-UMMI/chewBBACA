#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module

Code documentation
------------------
"""


import os
import sys
import logging
import platform as pf

try:
    from utils import constants as ct
except ModuleNotFoundError:
    from CHEWBBACA.utils import constants as ct


# create logger
logger = logging.getLogger('LS')


def main(version, logfile):

    # handler for file logging
    FILE_HANDLER = logging.FileHandler(filename=logfile)
    FILE_HANDLER.setLevel(logging.DEBUG)
    FILE_FORMATTER = logging.Formatter(fmt=ct.FILE_HANDLER_FORMAT,
                                       datefmt=ct.FILE_HANLDER_DATEFMT)
    FILE_HANDLER.setFormatter(FILE_FORMATTER)

    # handler for stdout logging
    STDOUT_HANDLER = logging.StreamHandler()
    STDOUT_HANDLER.setLevel(logging.WARNING)
    STDOUT_FORMATTER = logging.Formatter(fmt=ct.STDOUT_HANDLER_FORMAT)
    STDOUT_HANDLER.setFormatter(STDOUT_FORMATTER)

    # config logging
    # this will apply config to root and affect all that is logged
    logging.basicConfig(level=logging.DEBUG,
                        handlers=[FILE_HANDLER, STDOUT_HANDLER])

    # log general info
    logger.info(f'chewBBACA version: {version}')
    logger.info(f'Authors: {ct.authors}')
    logger.info(f'Github: {ct.repository}')
    logger.info(f'Documentation: {ct.documentation}')
    logger.info(f'Contacts: {ct.contacts}')
    logger.info(f'Machine: {pf.processor()}, {os.cpu_count()}')
    logger.info(f'System: {pf.system()}, {pf.release()}')
    logger.info(f'Python version: {pf.python_version()}')
    logger.info(f'Command used: {" ".join(sys.argv)}')
