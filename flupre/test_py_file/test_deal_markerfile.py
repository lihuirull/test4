# !/usr/bin/python 
# -*- coding: utf-8 -*-
# Author:lihuiru
# Created on 2023/11/29 9:49
from collections import defaultdict
from itertools import product

import pandas as pd
import re
import ast
import os
import numpy as np
pd.set_option('Max_columns', None)

HA_TYPES = [f"H{i}" for i in range(1, 17) if i != 3]
NA_TYPES = [f"N{i}" for i in range(1, 10) if i != 2]


# data = data[data["Protein Type"] != "combination"]
def split_residues(residues):
    """
    将残基字符串分割成位置和氨基酸部分。

    参数:
    - residues (str): 需要处理的残基字符串。

    返回:
    - List of tuples: 每个元组包含位置和氨基酸。
    """
    parts = []
    current_part = ''
    for char in residues:
        if char.isdigit() and current_part.isalpha():
            parts.append(current_part)
            current_part = char
        elif char.isalpha() and current_part.isdigit():
            parts.append(current_part)
            current_part = char
        else:
            current_part += char
    parts.append(current_part)  # 添加最后一部分
    return parts


def process_combinations_v2(residues, protein, index, *other_data):
    """
    处理残基中的单个蛋白的组合或者单个蛋白多种deletion构成的组合（用 '&' 分隔），并据此创建新的行。

    参数:
    - residues (str): 需要处理的残基字符串。
    - protein (str): 蛋白质类型。
    - index (int): 行索引。

    返回:
    - List of tuples: 每个元组包含处理后的 'Protein Type' 和 'Amino acid site'。
    """
    # 将残基分割成单独的部分
    split_residues_list = residues.split('&')
    new_rows = []
    # 遍历每个部分并创建新行
    # print(residues)
    for idx, residue in enumerate(split_residues_list):
        if "Deletion" in residue :
            deletion = residue.split("-", 1)[1]
            deletion = process_deletion(deletion)
            new_rows.append((f"{protein}-combination_{index}_{idx}", deletion, f"{protein}", *other_data))
            continue
        # elif "Deletion" in residue:
        #     deletion = process_deletion(residue)
        #     new_rows.append((f"{protein}-combination_{index}_{idx}", deletion, f"{protein}", *other_data))
        #     continue
        position = ''.join(filter(str.isdigit, residue)).strip()
        # 找到字符串中所有的数字并作为分割点
        split_points = [match.start() for match in re.finditer(r'\d', split_residues_list[-1])][-1] + 1

        amino_acids = split_residues_list[-1][split_points:]
        # 更新正则表达式，避免在匹配单个字符时产生空字符串
        updated_pattern = r"\[(.*?)\]|([^\[\]])"

        # 使用更新后的正则表达式进行拆分
        updated_split_result = [match for group in re.findall(updated_pattern, amino_acids)
                                for match in group if match]

        if '[' in amino_acids:
            for i in updated_split_result[idx]:
                for j in i:
                    # 打印每个变异
                    new_rows.append(
                        (f"{protein}-combination_{index}_{idx}", f"{position}{j}", f"{protein}", *other_data))

        else:
            # print('-'*50)

            # print(amino_acids)
            # print(idx)
            new_rows.append(
                (f"{protein}-combination_{index}_{idx}", f"{position}{amino_acids[idx]}", f"{protein}", *other_data))

    return new_rows


def process_deletion(residues):
    """
    处理删除类型的残基。

    参数:
    - residues (str): 需要处理的残基字符串。

    返回:
    - 经处理的残基数据。
    """
    positions = re.findall("\d+", residues)
    if "Complete" in residues or "Any" in residues:
        # 如果是完整删除或任意删除
        position_range = list(range(int(positions[0]), int(positions[1]) + 1))
        residues = [f"{i}-" for i in position_range]
        print(residues)
    else:
        # 如果不是完整删除或任意删除
        residues = f"{positions[0]}-"

    return residues


def process_human_residues_v3(row, index, other_columns):
    """
    根据更新的规则处理 'Amino acid site' 列，并保留其他列的内容。

    参数:
    - row (Series): 数据帧中的一行。
    - index (int): 行的索引。
    - other_columns (list): 需要保留的其他列的列名。

    返回:
    - List of tuples: 每个元组包含处理后的行信息。
    """
    protein = row['Protein Type']
    residues = row['Amino acid site']
    processed_rows = []

    # 获取其他列的数据
    other_data = {col: row[col] for col in other_columns}
    # 首先处理组合
    # print(residues)
    if '&' in residues:
        # 前面两种都是有“-”的
        if protein == "combination":
            processed_rows.extend(process_combinations_v4(residues, protein, index, *other_data.values()))
        elif protein in [f"H{i}" for i in range(1, 19)] and 'HA2' in residues:
            print(residues)
            processed_rows.extend(process_combinations_v4(residues, protein, index, *other_data.values()))
        # 此类为同一种蛋白的多个标志物构成，也不存在‘-’，唯一存在‘-’的可能是多个不同的deletion
        else:
            processed_rows.extend(process_combinations_v2(residues, protein, index, *other_data.values()))
    else:
        # 处理含有 '[]' 的残基
        if '[' in residues:
            prefix, rest = residues.split('[')
            variants = rest.split(']')[0]
            suffix = rest.split(']')[1] if ']' in rest else ''
            for variant in variants:
                new_residue = f"{prefix}{variant}{suffix}"
                processed_rows.append((protein, new_residue, protein, *other_data.values()))
        else:
            # 处理删除类型的残基
            if "Deletion" in residues:
                if "Any" in residues:
                    residues = process_deletion(residues)
                    for single_deletion in residues:
                        processed_rows.append((protein, single_deletion, protein, *other_data.values()))
                else:
                    residues = process_deletion(residues)
            processed_rows.append((protein, residues, protein, *other_data.values()))

    return processed_rows


def process_combinations_v4(residues, protein, index, *other_data):
    """
    处理残基中的combination组合和HA亚型的组合，因为只有这两种会存在-以表明是什么蛋白或者是HA2
    （用 '&' 分隔），并据此创建新的行。

    参数:
    - residues (str): 需要处理的残基字符串。
    - protein (str): 蛋白质类型。
    - index (int): 行索引。

    返回:
    - List of tuples: 每个元组包含处理后的 'Protein Type' 和 'Amino acid site'。
    """
    # 将残基分割成单独的部分
    split_residues_list = residues.split('&')

    new_rows = []
    # 遍历每个部分并创建新行
    for idx, residue in enumerate(split_residues_list):
        residue = residue.strip()
        protein_type = residue.split("-")[0]
        if "Deletion" in residue:
            deletion = residue.split("-", 1)[1]
            deletion = process_deletion(deletion)
            new_rows.append((f"{protein}-combination_{index}_{idx}", deletion, f"{protein_type}", *other_data))
            continue
        if "HA2" in residue or "HA1" in residue:
            new_rows.append(
                (f"{protein}-combination_{index}_{idx}", f"{residue.split('-',1)[1]}", f"{protein_type}", *other_data))
            continue

        position = ''.join(filter(str.isdigit, residue.split('-')[1])).strip()
        split_points = [match.start() for match in re.finditer(r'\d', residue)][-1] + 1
        amino_acids = residue[split_points:]

        if '[' in amino_acids:
            amino_acid = amino_acids.split("[")[1].split("]")[0]
            for i in amino_acid:
                new_rows.append(
                    (f"{protein}-combination_{index}_{idx}", f"{position}{i}", f"{protein_type}", *other_data))
        else:
            new_rows.append(
                (f"{protein}-combination_{index}_{idx}", f"{position}{amino_acids}", f"{protein_type}", *other_data))

    return new_rows


def merge_dicts_with_list(dict_list):
    """
    合并字典列表的函数，如果键重复，则将值合并为列表。

    参数:
    - dict_list (list): 包含字典的列表。

    返回:
    - 合并后的字典。
    """
    merged_dict = {}
    for d in dict_list:
        for key, value in d.items():
            if key in merged_dict:
                # 如果键已存在，则合并值为列表
                if not isinstance(merged_dict[key], list):
                    merged_dict[key] = [merged_dict[key]]
                merged_dict[key].append(value)
            else:
                # 如果键不存在，则直接添加
                merged_dict[key] = value
    return merged_dict


def generate_combinations(group):
    """
    根据特定类型对变异进行分组，并生成所有可能的组合。

    参数:
    - group (DataFrameGroupBy): 分组后的DataFrame。

    返回:
    - 所有可能的组合列表。
    """
    # 根据特定类型对变异进行分组
    spec_type_groups = group.groupby('Specific Type')
    # 创建字典以存储按特定类型映射到蛋白质的变异列表
    mutation_dict = defaultdict(list)
    for spec_type, g in spec_type_groups:
        for _, row in g.iterrows():
            if type(row["Mutation"]) == str:
                mutation_dict[spec_type].append({row['Protein']: row['Mutation'].strip()})
            else:
                mutation_dict[spec_type].append({row['Protein']: row['Mutation']})
    # 提取每个键对应的字典列表
    values_lists = [mutation_dict[key] for key in mutation_dict]
    # 生成所有可能的组合并合并字典
    combinations = [merge_dicts_with_list(comb) for comb in product(*values_lists)]
    return combinations


def read_and_prepare_data(file_name, column_names):
    """读取Excel文件并进行初步处理，包括调整列名"""
    data = pd.read_excel(file_name, engine = "openpyxl")
    # Replace empty strings with NaN
    data.replace('', np.nan, inplace = True)
    data = data.dropna(subset = ["Protein Type"],how = "all", axis = 0)
    data = data.dropna(how = "all", axis = 1)
    data.columns = column_names + data.columns[len(column_names):].tolist()
    data["Amino acid site"] = data["Amino acid site"].str.split('(', expand = True)[0]
    return data


def process_and_group_data(data, output_file):
    """处理数据并分组"""
    other_columns = data.columns.difference(['Protein Type', 'Amino acid site'])
    processed_data = []
    for index, row in data.iterrows():
        processed_rows = process_human_residues_v3(row, index, other_columns.tolist())
        processed_data.extend(processed_rows)

    column_names = ["Protein Type", "Mutation", "Protein"] + other_columns.tolist()
    processed_df = pd.DataFrame(processed_data, columns = column_names)

    processed_df.to_csv(output_file, index = False)

    processed_df["Specific Type"] = processed_df["Protein Type"].str.rsplit("_", n = 1).str[-1]
    processed_df['Protein Type'] = processed_df['Protein Type'].str.replace(r'_\d+$', '', regex = True)
    return processed_df.groupby('Protein Type')


def generate_protein_dict(grouped_data):
    """生成蛋白质字典"""
    new_protein_dict = defaultdict(list)
    for name, group in grouped_data:
        if 'combination' in name:
            new_protein_dict[name] = generate_combinations(group)
        else:
            new_protein_dict[name].extend(
                {row['Protein']: row['Mutation'].strip() if isinstance(row['Mutation'], str) else row['Mutation']}
                for _, row in group.iterrows()
            )
    return new_protein_dict


def main(input_file, output_file):
    """主函数"""
    column_names = ['Protein Type', 'Amino acid site']  # 列名调整
    data = read_and_prepare_data(input_file, column_names)
    grouped_data = process_and_group_data(data, output_file)
    new_protein_dict = generate_protein_dict(grouped_data)

    # 可以选择打印或返回结果
    print(len(new_protein_dict.keys()))
    print(new_protein_dict)
    return new_protein_dict


# 调用主函数
if __name__ == '__main__':
    phenotype = "mammalian_virulence"
    # phenotype = "human_adaptation"
    # phenotype = "transmissibility"
    # phenotype = "drug_resistance"
    # new_protein_dict = main(f"drug_resistance_std.xlsx",
    #                         f"../../data/markers_for_extract/{phenotype}_formated.csv")
    new_protein_dict = main(f"../../data/markers_for_extract/{phenotype}.xlsx",
        # f"../data/markers_for_extract/{phenotype}_formated.csv")
                            f"test_formated.csv")