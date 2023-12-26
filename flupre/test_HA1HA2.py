# !/usr/bin/python 
# -*- coding: utf-8 -*-
# Author:lihuiru
# Created on 2023/12/18 10:20
import os
import re

import pandas as pd
from Bio import SeqIO

pd.set_option('display.max_columns', None)

HA_TYPES = [f"H{i}" for i in range(1, 19) if i != 3]
NA_TYPES = [f"N{i}" for i in range(1, 10) if i != 2]

base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DB_PATH = os.path.join(base_dir, 'data', 'flu_db.dmnd')
STD_PATH = os.path.join(base_dir, 'data', 'std.fasta')
STRUCTURE_PATH = os.path.join(base_dir, 'data', 'HA_NA_mapdir')
STANDARD_PATH = os.path.join(base_dir, 'data', 'standard_seq_protein')
MODEL_PATH = os.path.join(base_dir, 'model')
DATA_PATH = os.path.join(base_dir, 'data')
MARKER_PATH = os.path.join(base_dir, 'data', 'markers_for_extract')


def load_markers(filepath):
    """
    Load and process markers from an input file.

    Parameters:
        filepath: Path to the Excel file containing virulence markers.

    Returns:
        Dictionary with protein types as keys and lists of virulence markers as values.
    """
    column_names = ['Protein Type', 'Amino acid site']
    data = pd.read_csv(filepath)
    data = data.dropna(how = "all", axis = 1)
    data.columns = column_names + data.columns[len(column_names):].tolist()
    data["Amino acid site"] = data["Amino acid site"].str.split('(', expand = True)[0]
    # data["Specific Type"] = data["Protein Type"].str.rsplit("_", n = 1).str[-1]
    # data['Protein Type'] = data['Protein Type'].str.replace(r'_\d+$', '', regex = True)
    return data.groupby('Protein')['Amino acid site'].apply(lambda x: list(set(x))).to_dict(), data


#
# def convert_HA_residues(marker_dict, structure_folder):
#     """
#     Converts HA/NA residues to H3/N2 numbering.
#
#     Parameters:
#         marker_dict: Dictionary with protein types as keys and marker lists as values.
#         structure_folder: Folder path where the structure mapping files are located.
#
#     Returns:
#         Updated marker_dict with HA/NA types converted to H3/N2 numbering.
#     """
#     updated_marker_dict = marker_dict.copy()  # Create copy
#     for protein in list(marker_dict.keys()):
#         if protein in HA_TYPES:
#             mapping_data = pd.read_csv(f"{structure_folder}/H3_{protein}.txt", sep = "\t", header = None,
#                                        names = ['H3', protein])
#             convert_to_h3_dict = dict(zip(mapping_data[protein], mapping_data['H3']))
#
#             residues = map_residues_to_h3(protein, marker_dict, convert_to_h3_dict)
#             # residues = [convert_to_H3_dict.get(re.search(r"\d+", i).group()) +
#             # re.search(r"[A-Z]", i).group() for i in
#             #             marker_dict[protein] if convert_to_H3_dict.get(re.search(r"\d+", i).group())]
#             if "H3" in updated_marker_dict:
#                 updated_marker_dict["H3"].extend(residues)
#             else:
#                 updated_marker_dict["H3"] = residues
#             del updated_marker_dict[protein]  # del key
#         elif protein in NA_TYPES:
#             if os.path.isfile(f"{structure_folder}/N2_{protein}.txt"):
#                 mapping_data = pd.read_csv(f"{structure_folder}/N2_{protein}.txt", sep = "\t", header = None,
#                                            names = ['N2', protein])
#             else:
#                 mapping_data = pd.read_csv(f"{structure_folder}/{protein}_N2.txt", sep = "\t", header = None,
#                                            names = [protein, 'N2'])
#             convert_to_n2_dict = dict(zip(mapping_data[protein], mapping_data['N2']))
#
#             residues = map_residues_to_h3(protein, marker_dict, convert_to_n2_dict)
#             if "N2" in updated_marker_dict:
#                 updated_marker_dict["N2"].extend(residues)
#             else:
#                 updated_marker_dict["N2"] = residues
#             del updated_marker_dict[protein]  # del key
#
#     return updated_marker_dict
#
# def map_residues_to_h3(protein, marker_dict, convert_to_h3_dict):
#     """
#     Maps the residue numbers for a given protein to the H3/N2 numbering system.
#
#     Parameters:
#         protein (str): The protein identifier.
#         marker_dict (dict): Dictionary containing markers for various proteins.
#         convert_to_h3_dict (dict): Dictionary that maps residue numbers to H3.
#
#     Returns:
#         list: A list of residues mapped to H3 numbering system.
#     """
#     markers = marker_dict[protein]
#
#     # 如果markers是字符串，将其转换为只含一个元素的列表
#     if isinstance(markers, str):
#         markers = [markers]
#
#     mapped_residues = []
#     for marker in markers:
#         # Ensure the marker is in the expected format (e.g., "12A")
#         marker_match = re.fullmatch(r"(\d+)([A-Z])", marker)
#         if not marker_match:
#             # print(f"Warning: marker '{marker}' is not in the correct format and will be skipped.")
#             continue
#
#         position, amino_acid = marker_match.groups()
#         h3_position = convert_to_h3_dict.get(position)
#         if h3_position is None:
#             # print(f"Warning: Position {position} does not have an H3 mapping "
#             #       f"in the structure comparison file and will be skipped.")
#             continue
#
#         mapped_residues.append(h3_position + amino_acid)
#
#     return mapped_residues

def get_h3_dict_and_hatype(protein, marker, convert_to_h3_dict):
    if protein in HA_TYPES:
        if "HA2" in marker:
            return convert_to_h3_dict["HA2"], "HA2"
        else:
            return convert_to_h3_dict["HA1"], "HA1" if "HA1" in marker else None
    elif protein in NA_TYPES:
        return convert_to_h3_dict, None
    return None, None


def adjust_position_and_get_h3_position(marker, hatype, H3_dict, protein):
    # marker_match = re.fullmatch(r"(\d+)([A-Z]|-)", marker)
    marker_match = re.search(r"(\d+)([A-Z]|-)", marker)

    if not marker_match:
        return None, None, hatype

    position, amino_acid = marker_match.groups()
    if not hatype and protein in HA_TYPES:
        minus = length_diffs[protein]
        position = str(int(position) - minus)
        hatype = "HA1"

    if H3_dict:
        # 处理除H3的情况
        return H3_dict.get(position), amino_acid, hatype
    else:
        # 处理H3（不需要位点转换）
        return f"HA1-{position}{amino_acid}"


def map_residues_to_h3(protein, marker_dict, convert_to_h3_dict, hatype = None):
    markers = [marker_dict[protein]] if isinstance(marker_dict[protein], str) else marker_dict[protein]
    map_dic = {"496S": "158S", "409P": "65P", "434G": "90G", "445G": "101G", "425G": "81G", "425M": "79M",
               "452T": "111T"}
    markers = [f"HA2-{map_dic[marker]}" if marker in map_dic else marker for marker in markers]

    mapped_residues = []

    for marker in markers:
        if hatype:
            H3_dict = convert_to_h3_dict[hatype]
        else:
            H3_dict, hatype = get_h3_dict_and_hatype(protein, marker, convert_to_h3_dict)

        if H3_dict is None:
            continue

        if ',' in marker:
            # 处理几种删除的情况
            for marker in marker.split(','):
                h3_position, amino_acid, updated_hatype = adjust_position_and_get_h3_position(marker, hatype, H3_dict,
                                                                                              protein)
                if h3_position is None:
                    continue

                hatype_prefix = f"{updated_hatype}-" if updated_hatype else ""
                mapped_residues.append(f"{hatype_prefix}{h3_position}{amino_acid}")
            continue
        elif not marker.endswith("-"):
            marker = marker.strip().split("-")[-1]

        h3_position, amino_acid, updated_hatype = adjust_position_and_get_h3_position(marker, hatype, H3_dict, protein)
        if h3_position is None:
            continue

        hatype_prefix = f"{updated_hatype}-" if updated_hatype else ""
        mapped_residues.append(f"{hatype_prefix}{h3_position}{amino_acid}")

    return mapped_residues


#
# def convert_HA_residues(marker_dict, structure_folder):
#     """
#     Converts HA/NA residues to H3/N2 numbering.
#
#     Parameters:
#         marker_dict: Dictionary with protein types as keys and marker lists as values.
#         structure_folder: Folder path where the structure mapping files are located.
#
#     Returns:
#         Updated marker_dict with HA/NA types converted to H3/N2 numbering.
#     """
#     updated_marker_dict = marker_dict.copy()  # Create copy
#     for protein in list(marker_dict.keys()):
#         if protein  == "H3":
#             res = []
#             for marker in marker_dict["H3"]:
#                 if "HA2" not in marker and "HA1" not in marker:
#                     marker = adjust_position_and_get_h3_position(marker, hatype = None, H3_dict = None,protein = "H3")
#                 res.append(marker)
#             updated_marker_dict["H3"] = res
#         # print(marker_dict["H3"])
#         if protein in HA_TYPES:
#
#             mapping_data_HA1 = pd.read_csv(
#                 rf"D:\user\data\fluphenotype\script\convert_site\HA_NA_mapdir\HA1/H3_{protein}.txt", sep = "\t",
#                 header = None,
#                 names = ['H3', protein])
#             convert_to_h3_dict_ha1 = dict(zip(mapping_data_HA1[protein], mapping_data_HA1['H3']))
#             # print(convert_to_h3_dict_ha1)
#             mapping_data_HA2 = pd.read_csv(
#                 rf"D:\user\data\fluphenotype\script\convert_site\HA_NA_mapdir\HA2/H3_{protein}.txt", sep = "\t",
#                 header = None,
#                 names = ['H3', protein])
#             convert_to_h3_dict_ha2 = dict(zip(mapping_data_HA2[protein], mapping_data_HA2['H3']))
#
#             combined_dict = {
#                 'HA1': convert_to_h3_dict_ha1,
#                 'HA2': convert_to_h3_dict_ha2
#             }
#             # print('-'*50)
#             # # print(convert_to_h3_dict_ha1)
#             # # print(convert_to_h3_dict_ha2)
#             # print(combined_dict)
#             # print(protein)
#             # print(combined_dict)
#             residues = map_residues_to_h3(protein, marker_dict, combined_dict)
#
#             if "H3" in updated_marker_dict:
#                 updated_marker_dict["H3"].extend(residues)
#             else:
#                 updated_marker_dict["H3"] = residues
#             del updated_marker_dict[protein]  # del key
#         elif protein in NA_TYPES:
#             if os.path.isfile(f"{structure_folder}/N2_{protein}.txt"):
#                 mapping_data = pd.read_csv(f"{structure_folder}/N2_{protein}.txt", sep = "\t", header = None,
#                                            names = ['N2', protein])
#             else:
#                 mapping_data = pd.read_csv(f"{structure_folder}/{protein}_N2.txt", sep = "\t", header = None,
#                                            names = [protein, 'N2'])
#             convert_to_n2_dict = dict(zip(mapping_data[protein], mapping_data['N2']))
#
#             residues = map_residues_to_h3(protein, marker_dict, convert_to_n2_dict)
#             if "N2" in updated_marker_dict:
#                 updated_marker_dict["N2"].extend(residues)
#             else:
#                 updated_marker_dict["N2"] = residues
#             del updated_marker_dict[protein]  # del key
#     updated_marker_dict = transform_marker_dict(updated_marker_dict)
#     return updated_marker_dict
import pandas as pd
import os


def load_mapping_data(filepath, column_names):
    if os.path.isfile(filepath):
        mapping_data = pd.read_csv(filepath, sep = "\t", header = None, names = column_names)
        return dict(zip(mapping_data[column_names[1]], mapping_data[column_names[0]]))
    return {}


def process_ha_type(protein, marker_dict, structure_folder, hatype):
    convert_to_h3_dict_ha1 = load_mapping_data(
        f"{STRUCTURE_PATH}/HA1/H3_{protein}.txt", ['H3', protein])
    convert_to_h3_dict_ha2 = load_mapping_data(
        f"{STRUCTURE_PATH}/HA2/H3_{protein}.txt", ['H3', protein])

    combined_dict = {'HA1': convert_to_h3_dict_ha1, 'HA2': convert_to_h3_dict_ha2}
    return map_residues_to_h3(protein, marker_dict, combined_dict, hatype)


def process_na_type(protein, marker_dict, structure_folder, hatype):
    convert_to_n2_dict = load_mapping_data(
        f"{STRUCTURE_PATH}/NA/N2_{protein}.txt", ['N2', protein])
    # print(convert_to_n2_dict)
    return map_residues_to_h3(protein, marker_dict, convert_to_n2_dict, hatype)


def convert_HA_residues(marker_dict, structure_folder, hatype):
    updated_marker_dict = marker_dict.copy()

    for protein in list(marker_dict.keys()):
        if protein == "H3":
            # 处理特殊情况
            res = []
            for marker in marker_dict["H3"]:
                if "HA2" not in marker and "HA1" not in marker:
                    # 假设 adjust_position_and_get_h3_position 函数适用于这种情况
                    marker = adjust_position_and_get_h3_position(marker, hatype = None, H3_dict = None, protein = "H3")
                res.append(marker)
            updated_marker_dict["H3"] = res

        if protein in HA_TYPES:
            residues = process_ha_type(protein, marker_dict, structure_folder, hatype)
            updated_marker_dict["H3"] = updated_marker_dict.get("H3", []) + residues
            del updated_marker_dict[protein]

        elif protein in NA_TYPES:
            residues = process_na_type(protein, marker_dict, structure_folder, hatype)
            updated_marker_dict["N2"] = updated_marker_dict.get("N2", []) + residues
            del updated_marker_dict[protein]

    return transform_marker_dict(updated_marker_dict)


def transform_marker_dict(marker_dict):
    transformed_data = {}
    for key, values in marker_dict.items():
        if key == 'H3':
            # 应用特定于 'H3' 的转换逻辑
            sub_dict = {}
            for value in values:
                prefix, suffix = value.split('-', 1)
                if prefix not in sub_dict:
                    sub_dict[prefix] = []
                sub_dict[prefix].append(suffix)
            transformed_data[key] = {k: v[0] if len(v) == 1 else v for k, v in sub_dict.items()}
        else:
            # 对其他键应用去重逻辑
            transformed_data[key] = list(set(values))
    return transformed_data


def read_fasta(file_path):
    """读取fasta文件并返回一个字典，其中键是描述（例如'H1'），值是序列"""
    sequences = {}
    for record in SeqIO.parse(file_path, "fasta"):
        description = record.description.split('_')[0]  # 假设描述的格式是 'H1_...'
        sequences[description] = str(record.seq)
    return sequences


def compare_sequences(seq_file1, seq_file2):
    """比较两个fasta文件中相同键的序列长度差异"""
    seq_dict1 = read_fasta(seq_file1)
    seq_dict2 = read_fasta(seq_file2)

    length_differences = {}
    for key in seq_dict1:
        if key in seq_dict2 and key in [f"H{i}" for i in range(1, 19)]:
            length_differences[key] = abs(len(seq_dict1[key]) - len(seq_dict2[key]))
    return length_differences


def annotate_markers(markers_path, STRUCTURE_PATH, hatype = None):
    """
    Annotate markers by loading virulence data, then converting HA types.

    Parameters:
        markers_path: Path to the Excel file with markers.

    Returns:
        A dictionary with annotated markers.
    """
    # Load markers from files
    marker_dict, data = load_markers(markers_path)

    # Duplicated
    marker_dict = {i: list(set(j)) for i, j in marker_dict.items()}
    print(f"之前的marker_dict{marker_dict}")
    # for i, j in marker_dict.items():
    #     if i in [f"H{i}" for i in range(1, 19)]:
    #         print(f"Protein:\n{i}")
    #         print(f"Markers:\n{j}")
    #         print("-" * 50)
    # Convert HA/NA residues to H3/N2 numbering and update marker_dict
    # Already Duplicated
    marker_dict = convert_HA_residues(marker_dict, STRUCTURE_PATH, hatype)

    # Duplicated
    # marker_dict = {i: list(set(j)) for i, j in marker_dict.items()}
    print(marker_dict)
    return marker_dict, data


# # 文件路径
# std_fasta_path = "test_py_file/std.fasta"
# complete_std_fasta_path = "test_py_file/complete_std.fasta"
#
# # 比较并获取长度差异
# length_diffs = compare_sequences(std_fasta_path, complete_std_fasta_path)
#
# s, data = annotate_markers(r"test_py_file/test_formated.csv",
#                            r"D:\user\data\fluphenotype\script\convert_site\HA_NA_mapdir")
# print(data)

# marker_dict = {'M2': '41C'}
# for protein in list(marker_dict.keys()):
#     values = marker_dict[protein]
#     if not isinstance(values, str):
#         print(list(set(values)))

# s =  {'H3': {'HA1': '218W', 'HA2': '156N'}}
# s1 =  {'H3': {'HA1': '218W', 'HA2': ['156N','222L']}}
# s = {'PB2': '158G', 'PA': '295P'}
# s1 = {'PB2': '158G', 'PA': ['295P', '225L']}
s = {'NP': '253I', 'PB2': '288Q'}
s1 = {'H3': {'HA1': ['184H', '98Y', '169P'], 'HA2': ['111T']},
      'NP': ['437T', '482S', '253I', '133L', '184K'],
      'PB2': ['504V', '191E', '288Q' ]}


# Revised function to correctly handle nested dictionaries
def is_subset_complex_revised(dict1, dict2):
    """
    Check if one dictionary is a complex subset of another, with revised logic for nested dictionaries.

    Parameters:
    - dict1, dict2: Dictionaries to be compared.

    Returns:
    - Boolean: True if dict1 is a subset of dict2, else False.
    """
    for key, value1 in dict1.items():
        if key not in dict2:
            return False

        value2 = dict2[key]

        # Check for nested dictionaries
        if isinstance(value1, dict) and isinstance(value2, dict):
            if not is_subset_complex_revised(value1, value2):
                return False
        # Check for list and string combinations
        elif isinstance(value1, list) and isinstance(value2, list):
            if not set(value1).issubset(set(value2)):
                return False
        elif isinstance(value1, str) and isinstance(value2, list):
            if value1 not in value2:
                return False
        elif isinstance(value1, str) and isinstance(value2, str):
            if value1 != value2:
                return False

    return True


# Test the revised function
print(is_subset_complex_revised(s, s1))

# def format_marker(marker, protein_prefix=''):
#     if '-' in marker:
#         amino_acid = marker.split('-')[0]
#         deletion_suffix = "Deletion"
#     else:
#         amino_acid = marker
#         deletion_suffix = ""
#
#     formatted_marker = f"{protein_prefix}-{amino_acid}{deletion_suffix}" \
#         if protein_prefix else f"{amino_acid}{deletion_suffix}"
#     return formatted_marker
#
# def format_marker_list(markers, protein_prefix=''):
#     if isinstance(markers, str):
#         return format_marker(markers, protein_prefix)
#
#     all_contain_dash = all('-' in marker for marker in markers)
#     if all_contain_dash:
#         start = markers[0].split('-')[0]
#         end = markers[-1].split('-')[0]
#         return f"{protein_prefix}-{start}-{end}CompleteDeletion"
#
#     return '&'.join(format_marker(marker, protein_prefix) for marker in markers)
#
# def process_dictionary(data_dict):
#     formatted_list = []
#     for protein, markers in data_dict.items():
#         if isinstance(markers, dict):
#             for sub_protein, sub_markers in markers.items():
#                 formatted_marker = format_marker_list(sub_markers, f"{protein}-{sub_protein}")
#                 formatted_list.append(formatted_marker)
#         else:
#             formatted_marker = format_marker_list(markers, protein)
#             formatted_list.append(formatted_marker)
#
#     return '&'.join(formatted_list)


def format_marker(marker, protein_prefix = ''):
    """
    Formats a single genetic marker. If the marker contains a hyphen ('-'),
    only the part before the hyphen is retained and appended with 'Deletion'.
    If a protein prefix is provided, it's added before the marker.

    Parameters:
        marker (str): The genetic marker to be formatted.
        protein_prefix (str): An optional prefix to be added before the marker.

    Returns:
        str: Formatted genetic marker.
    """
    # Check if the marker contains a hyphen and split accordingly.
    if '-' in marker:
        amino_acid = marker.split('-')[0]
        deletion_suffix = "Deletion"
    else:
        amino_acid = marker
        deletion_suffix = ""

    # Combine the protein prefix, amino acid, and deletion suffix.
    formatted_marker = f"{protein_prefix}-{amino_acid}{deletion_suffix}" \
        if protein_prefix else f"{amino_acid}{deletion_suffix}"
    return formatted_marker


def format_marker_list(markers, protein_prefix = ''):
    """
    Formats a list of markers or a single marker string.
    In case of a list where all elements contain a hyphen, a special formatted string is returned.
    Otherwise, each marker in the list is formatted individually.

    Parameters:
        markers (str or list): A string or list of strings representing genetic markers.
        protein_prefix (str): An optional prefix to be added before each marker.

    Returns:
        str: A single string representing the formatted markers, joined by '&'.
    """
    # Check if the input is a single string and format directly.
    if isinstance(markers, str):
        return format_marker(markers)

    # Determine if all markers in the list contain a hyphen.
    all_contain_dash = all('-' in marker for marker in markers)
    if all_contain_dash:
        # Create a special format string if all markers contain a hyphen.
        start = markers[0].split('-')[0]
        end = markers[-1].split('-')[0]
        return f"{start}-{end}CompleteDeletion"

    # Format each marker individually and join with '&'.
    return '&'.join(format_marker(marker, protein_prefix) for marker in markers)


# def process_dictionary(data_dict):
#     """
#     Processes a dictionary containing genetic markers.
#     If the dictionary has a single key-value pair, the value is formatted directly.
#     For multiple key-value pairs, each is formatted separately and joined by '&'.
#
#     Parameters:
#         data_dict (dict): A dictionary with protein names as keys and genetic markers as values.
#
#     Returns:
#         str: A single string representing the formatted contents of the dictionary.
#     """
#     # Process a single key-value pair directly.
#     if len(data_dict) == 1:
#         return format_marker_list(next(iter(data_dict.values())))
#
#     # Format each key-value pair separately if there are multiple.
#     return '&'.join(format_marker_list(markers, protein) for protein, markers in data_dict.items())

def process_dictionary(data_dict):
    formatted_list = []
    for protein, markers in data_dict.items():
        # 检查是否为嵌套字典
        if isinstance(markers, dict):
            for sub_protein, sub_markers in markers.items():
                # 将内部键作为额外的前缀
                formatted_marker = format_marker_list(sub_markers, f"{protein}-{sub_protein}")
                formatted_list.append(formatted_marker)
        else:
            # 非嵌套字典直接处理
            formatted_marker = format_marker_list(markers, protein)
            formatted_list.append(formatted_marker)

    return '&'.join(formatted_list)
# 测试用例
test_dicts = [
    {'H3': {'HA1': ['184R', '214L', '215N', '216S']}},
    {'H3': {'HA1': ['144E', '246K', '304T']}, 'PA': '615E'},
    {'PA': '653L', 'PB1': '229R'}
]
test_dicts = [{'H3': {'HA1': ['214L', '215N', '216S']}}]
# 测试结果
test_results = [process_dictionary(test_dict) for test_dict in test_dicts]
print(test_results)



# def map_residues_to_h3(protein, marker_dict, convert_to_h3_dict):
#     """
#     Maps the residue numbers for a given protein to the H3/N2 numbering system.
#
#     Parameters:
#         protein (str): The protein identifier.
#         marker_dict (dict): Dictionary containing markers for various proteins.
#         convert_to_h3_dict (dict): Dictionary that maps residue numbers to H3.
#
#     Returns:
#         list: A list of residues mapped to H3 numbering system.
#     """
#     markers = marker_dict[protein]
#
#     # 如果markers是字符串，将其转换为只含一个元素的列表
#     if isinstance(markers, str):
#         markers = [markers]
#
#     mapped_residues = []
#     for marker in markers:
#         # Ensure the marker is in the expected format (e.g., "12A")
#         marker_match = re.fullmatch(r"(\d+)([A-Z])", marker)
#         if not marker_match:
#             # print(f"Warning: marker '{marker}' is not in the correct format and will be skipped.")
#             continue
#
#         position, amino_acid = marker_match.groups()
#         h3_position = convert_to_h3_dict.get(position)
#         if h3_position is None:
#             # print(f"Warning: Position {position} does not have an H3 mapping "
#             #       f"in the structure comparison file and will be skipped.")
#             continue
#
#         mapped_residues.append(h3_position + amino_acid)
#
#     return mapped_residues
#
#
# def convert_HA_residues(marker_dict, structure_folder):
#     """
#     Converts HA/NA residues to H3/N2 numbering.
#
#     Parameters:
#         marker_dict: Dictionary with protein types as keys and marker lists as values.
#         structure_folder: Folder path where the structure mapping files are located.
#
#     Returns:
#         Updated marker_dict with HA/NA types converted to H3/N2 numbering.
#     """
#     updated_marker_dict = marker_dict.copy()  # Create copy
#     for protein in list(marker_dict.keys()):
#         if protein in HA_TYPES:
#             mapping_data = pd.read_csv(f"{structure_folder}/H3_{protein}.txt", sep = "\t", header = None,
#                                        names = ['H3', protein])
#             convert_to_h3_dict = dict(zip(mapping_data[protein], mapping_data['H3']))
#
#             residues = map_residues_to_h3(protein, marker_dict, convert_to_h3_dict)
#             # residues = [convert_to_H3_dict.get(re.search(r"\d+", i).group()) +
#             # re.search(r"[A-Z]", i).group() for i in
#             #             marker_dict[protein] if convert_to_H3_dict.get(re.search(r"\d+", i).group())]
#             if "H3" in updated_marker_dict:
#                 updated_marker_dict["H3"].extend(residues)
#             else:
#                 updated_marker_dict["H3"] = residues
#             del updated_marker_dict[protein]  # del key
#         elif protein in NA_TYPES:
#             if os.path.isfile(f"{structure_folder}/N2_{protein}.txt"):
#                 mapping_data = pd.read_csv(f"{structure_folder}/N2_{protein}.txt", sep = "\t", header = None,
#                                            names = ['N2', protein])
#             else:
#                 mapping_data = pd.read_csv(f"{structure_folder}/{protein}_N2.txt", sep = "\t", header = None,
#                                            names = [protein, 'N2'])
#             convert_to_n2_dict = dict(zip(mapping_data[protein], mapping_data['N2']))
#
#             residues = map_residues_to_h3(protein, marker_dict, convert_to_n2_dict)
#             if "N2" in updated_marker_dict:
#                 updated_marker_dict["N2"].extend(residues)
#             else:
#                 updated_marker_dict["N2"] = residues
#             del updated_marker_dict[protein]  # del key
#
#     return updated_marker_dict
#
#
# def annotate_markers(markers_path):
#     """
#     Annotate markers by loading virulence data, then converting HA types.
#
#     Parameters:
#         markers_path: Path to the Excel file with markers.
#
#     Returns:
#         A dictionary with annotated markers.
#     """
#     # Load markers from files
#     marker_dict, data = load_markers(markers_path)
#
#     # Duplicated
#     marker_dict = {i: list(set(j)) for i, j in marker_dict.items()}
#
#     # Convert HA/NA residues to H3/N2 numbering and update marker_dict
#     marker_dict = convert_HA_residues(marker_dict, STRUCTURE_PATH)
#
#     # Duplicated
#     marker_dict = {i: list(set(j)) for i, j in marker_dict.items()}
#
#     return marker_dict, data

# def adjust_position_and_get_h3_position(marker, hatype, H3_dict, protein):
#     # marker_match = re.fullmatch(r"(\d+)([A-Z]|-)", marker)
#     # length_diffs = compare_sequences(STD_PATH, COMPLETE_STD_PATH)
#     print(f"adjH2:\n{marker}")
#     marker_match = re.search(r"(\d+)([A-Z]|-)", marker)
#
#     if not marker_match:
#         return None, None, hatype
#
#     position, amino_acid = marker_match.groups()
#
#     if protein != "H3":
#         print(f"str\n{H3_dict.get(str(position))}")
#         print(f"int\n{H3_dict.get(int(position))}")
#         print(H3_dict)
#     if H3_dict:
#         # 处理除H3的情况
#         return H3_dict.get(int(position)), amino_acid, hatype
#     else:
#         if not hatype:
#             hatype = "HA1" #处理标志物字典的H3标志物
#         # 处理H3（不需要位点转换）
#         return f"{hatype}-{position}{amino_acid}"
# h3_position, amino_acid, updated_hatype  = adjust_position_and_get_h3_position("222R", "HA1", H3_dict, "H2")
# print(h3_position)
