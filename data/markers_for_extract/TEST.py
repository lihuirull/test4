# !/usr/bin/python 
# -*- coding: utf-8 -*-
# Author:lihuiru
# Created on 2023/12/10 17:44
import pandas as pd

import os

f = os.listdir("./")
lst = []
for file in f:
    if file.endswith(".xlsx"):
        df = pd.read_excel(file)
        lst.extend(df["Protein Type"].tolist())

print(set(lst))