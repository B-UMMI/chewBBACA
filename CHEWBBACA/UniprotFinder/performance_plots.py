#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 13:35:31 2021

@author: rfm
"""


import os
import math
import random
from itertools import groupby

import plotly.graph_objs as go
from plotly.offline import plot
from plotly.subplots import make_subplots


species = ['S.agalactiae 32',
           'S.agalactiae 320',
           'L. pneumophila 32',
           'L. pneumophila 320',
           'C. difficile 32',
           'C. difficile 320']

legacy_time_seconds = [1010,
                      9033,
                      1284,
                      15880,
                      3512,
                      39002]

legacy_time_minutes = [v/60 for v in legacy_time_seconds]

legacy_memory_values = [261,
                        159,
                        116,
                        342,
                        171,
                        467]

latest_time_seconds = [56,
                      228,
                      81,
                      420,
                      158,
                      826]

latest_time_minutes = [v/60 for v in latest_time_seconds]

latest_memory_values = [237,
                        270,
                        187,
                        388,
                        221,
                        531]

# runtime bar plots
fig = go.Figure(data=[
    go.Bar(name='legacy',
           x=legacy_time_minutes,
           y=species,
           orientation='h',
           marker=dict(color='#4292c6')),
    go.Bar(name='latest',
           x=latest_time_minutes,
           y=species,
           orientation='h',
           marker=dict(color='#41ab5d'))])

fig.update_xaxes(type='log')
fig.update_layout(barmode='group')

plot(fig, filename='runtime.html')

# memory bar plots
fig = go.Figure(data=[
    go.Bar(name='legacy',
           x=legacy_memory_values,
           y=species,
           orientation='h',
           marker=dict(color='#4292c6')),
    go.Bar(name='latest',
           x=latest_memory_values,
           y=species,
           orientation='h',
           marker=dict(color='#41ab5d'))])

#fig.update_xaxes(type='log')
fig.update_layout(barmode='group')

plot(fig, filename='memory.html')



