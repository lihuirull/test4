# !/usr/bin/python 
# -*- coding: utf-8 -*-
# Author:lihuiru
# Created on 2023/12/10 17:44
import pandas as pd

import os


df = pd.read_csv(r"D:\user\data\fluphenotype\script\flupre\data\markers_for_extract\mammalian_virulence_formated.csv")
df = df[df.loc[:,"Protein"].isin([f"H{i}" for i in range(1,19)])]
df.loc[:,"Mutation"] = df.loc[:,"Mutation"].str.replace("HA\d-","")
print(df.Mutation)

# 使用更通用的正则表达式提取数字
df["Position"] = df["Mutation"].str.extract("(\d+)")
print(df)