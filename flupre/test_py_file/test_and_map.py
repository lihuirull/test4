# !/usr/bin/python 
# -*- coding: utf-8 -*-
# Author:lihuiru
# Created on 2023/12/18 10:49

import os
import re
from Bio import SeqIO
import pandas as pd

pd.set_option('Max_columns', None)

HA_TYPES = [f"H{i}" for i in range(1, 17) if i != 3]
NA_TYPES = [f"N{i}" for i in range(1, 10) if i != 2]

t = [{'N1': '106I'}, {'N1': '219Q'}, {'N1': '36-'}, {'N1': '44Q'}, {'N1': '49-'}, {'N1': '49-'}, {'N1': '50-'},
     {'N1': '51-'}, {'N1': '52-'}, {'N1': '53-'}, {'N1': '54-'}, {'N1': '55-'}, {'N1': '56-'}, {'N1': '57-'},
     {'N1': '58-'}, {'N1': '59-'}, {'N1': '60-'}, {'N1': '61-'}, {'N1': '62-'}, {'N1': '63-'}, {'N1': '64-'},
     {'N1': '65-'}, {'N1': '66-'}, {'N1': '67-'}, {'N1': '68-'}, {'N1': '69-'}, {'N1': '70-'}, {'N1': '71-'},
     {'N1': '72-'}, {
         'N1': ['49-', '50-', '51-', '52-', '53-', '54-', '55-', '56-', '57-', '58-', '59-', '60-', '61-', '62-', '63-',
                '64-', '65-', '66-', '67-', '68-', '69-', '70-', '71-', '72-']}, {'N1': '72Q'}]


# def format_markers(value, add_protein = ''):
#     if isinstance(value, str):
#         if add_protein:
#             return f'{add_protein}-{value}'
#         else:
#             return value
#     elif isinstance(value, list):
#         # 使用列表推导式来检查每个字符串中是否含有 '-'
#         contains_dash = ['-' in s for s in value]
#
#         # 检查是否所有的字符串都包含短横线 '-'
#         all_contain_dash = all(contains_dash)
#         if all_contain_dash:
#             return value[0].split('-')[0] + '-' + value[-1].split('-')[0] + "CompleteDeletion"
#         res = ''
#         for marker in value:
#             if '-' in marker:
#                 amino_acid = marker.split('-')[0]
#                 if add_protein:
#                     res += f"{add_protein}-{amino_acid}Deletion&"
#                 else:
#                     res += f"{amino_acid}Deletion&"
#             else:
#                 if add_protein:
#                     res += f"{add_protein}-{marker}&"
#                 else:
#                     res += f"{marker}&"
#         return res
# def deal_dict(i):
#
#     if len(i) == 1:
#         value = list(i.values())[0]
#         return format_markers(value).rsplit("&",1)[0]
#     else:
#         res = ""
#         for prot, value in i.items():
#             marker_symbol = format_markers(value,prot)
#             if marker_symbol.endswith("&"):
#                 res += f"{marker_symbol}"
#             else:
#                 res += f"{marker_symbol}&"
#         return res.rsplit("&",1)[0]

def format_marker(marker, protein_prefix = ''):
    """
    格式化单个遗传标记。如果标记中包含短横线（'-'），则仅保留短横线之前的部分，并附加'Deletion'。
    如果提供了蛋白质前缀，它将被添加到标记之前。

    参数:
        marker (str): 需要格式化的遗传标记。
        protein_prefix (str): 可选，添加到每个标记之前的前缀。

    返回:
        str: 格式化后的遗传标记。
    """
    # 检查标记是否包含短横线，并相应地分割。
    if '-' in marker:
        amino_acid = marker.split('-')[0]
        deletion_suffix = "Deletion"
    else:
        amino_acid = marker
        deletion_suffix = ""

    # 组合蛋白质前缀、氨基酸和删除后缀。
    formatted_marker = f"{protein_prefix}-{amino_acid}{deletion_suffix}" if protein_prefix else f"{amino_acid}{deletion_suffix}"
    return formatted_marker


def format_marker_list(markers, protein_prefix = ''):
    """
    格式化标记列表或单个标记字符串。
    如果是列表且所有元素都包含短横线，则返回特殊格式的字符串。
    否则，列表中的每个标记都将单独格式化。

    参数:
        markers (str 或 list): 表示遗传标记的字符串或字符串列表。
        protein_prefix (str): 可选，添加到每个标记之前的前缀。

    返回:
        str: 代表格式化后的标记的单个字符串，用'&'连接。
    """
    # 如果输入是单个字符串，直接格式化。
    if isinstance(markers, str):
        return format_marker(markers)

    # 确定列表中所有标记是否都包含短横线。
    all_contain_dash = all('-' in marker for marker in markers)
    if all_contain_dash:
        # 如果所有标记都包含短横线，创建特殊格式的字符串。
        start = markers[0].split('-')[0]
        end = markers[-1].split('-')[0]
        return f"{start}-{end}CompleteDeletion"

    # 单独格式化每个标记并用'&'连接。
    return '&'.join(format_marker(marker, protein_prefix) for marker in markers)


def process_dictionary(data_dict):
    """
    处理包含遗传标记的字典。
    如果字典只有一个键值对，直接格式化该值。
    对于多个键值对，分别格式化每个键值对并用'&'连接。

    参数:
        data_dict (dict): 以蛋白质名称为键，遗传标记为值的字典。

    返回:
        str: 代表字典内容格式化后的单个字符串。
    """
    # 如果只有一个键值对，直接处理这个值。
    if len(data_dict) == 1:
        return format_marker_list(next(iter(data_dict.values())))

    # 如果有多个键值对，分别处理每个键值对。
    return '&'.join(format_marker_list(markers, protein) for protein, markers in data_dict.items())


def compare_dicts_updated(dict1, dict2):
    for key, value1 in dict1.items():
        if key not in dict2:
            return False

        value2 = dict2[key]

        if isinstance(value1, list) and isinstance(value2, list):
            # 如果两个值都是列表，检查它们是否包含相同的元素（这里不考虑顺序）

            if not set(value1).issubset(set(value2)):
                return False
        elif isinstance(value1, str) and isinstance(value2, list):
            # 如果dict1中的值是字符串，而dict2中的值是列表，则检查字符串是否在列表中
            if value1 not in value2:
                return False
        elif value1 != value2:
            # 其他情况，直接比较值
            return False

    return True


dic2 = {'H3': ['216E', '223V', '146S', '263R', '225G', '229R'], 'M1': ['43M', '215A'],
        'M2': ['82S', '24D'], 'N1(N2 numbering)': [],
        'NP': ['482S', '184K', '437T', '105V', '253I', '373T', '133L', '286A'],
        'PA': ['383D', '224S', '190S', '550L', '237E', '321N', '149S', '295P', '409S', '394D', '330I', '100V'],
        'PA-X': [], 'PB1': ['298L', '652A', '115Q', '473V', '469T', '598L', '386R', '517I', '13P'],
        'PB1-F2': ['87E', '56V'],
        'PB2': ['627E', '715N', '191E', '661A', '504V', '559T', '495V', '283M', '339K', '309D', '66M', '89V', '133V',
                '389R', '598T', '288Q', '477G', '683T', '109V', '391E', '431M']}


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
    marker_match = re.fullmatch(r"(\d+)([A-Z])", marker)
    if not marker_match:
        return None, None, hatype

    position, amino_acid = marker_match.groups()
    if not hatype and protein in HA_TYPES:
        minus = length_diffs[protein]
        position = str(int(position) - minus)
        hatype = "HA1"

    return H3_dict.get(position), amino_acid, hatype


def map_residues_to_h3(protein, marker_dict, convert_to_h3_dict):
    markers = [marker_dict[protein]] if isinstance(marker_dict[protein], str) else marker_dict[protein]
    map_dic = {"496S": "158S", "409P": "65P", "434G": "90G", "445G": "101G", "425G": "81G", "425M": "79M", "452T": "111T"}
    markers = [f"HA2-{map_dic[marker]}" if marker in map_dic else marker for marker in markers]

    mapped_residues = []
    for marker in markers:
        H3_dict, hatype = get_h3_dict_and_hatype(protein, marker, convert_to_h3_dict)
        if H3_dict is None:
            continue

        h3_position, amino_acid, updated_hatype = adjust_position_and_get_h3_position(marker.split("-")[-1], hatype, H3_dict, protein)
        if h3_position is None:
            continue

        hatype_prefix = f"{updated_hatype}-" if updated_hatype else ""
        mapped_residues.append(f"{hatype_prefix}{h3_position}{amino_acid}")

    return mapped_residues




def convert_HA_residues(marker_dict, structure_folder):
    """
    Converts HA/NA residues to H3/N2 numbering.

    Parameters:
        marker_dict: Dictionary with protein types as keys and marker lists as values.
        structure_folder: Folder path where the structure mapping files are located.

    Returns:
        Updated marker_dict with HA/NA types converted to H3/N2 numbering.
    """
    updated_marker_dict = marker_dict.copy()  # Create copy
    for protein in list(marker_dict.keys()):
        if protein in HA_TYPES:

            mapping_data_HA1 = pd.read_csv(
                rf"D:\user\data\fluphenotype\script\convert_site\HA_NA_mapdir\HA1/H3_{protein}.txt", sep = "\t",
                header = None,
                names = ['H3', protein])
            convert_to_h3_dict_ha1 = dict(zip(mapping_data_HA1[protein], mapping_data_HA1['H3']))
            # print(convert_to_h3_dict_ha1)
            mapping_data_HA2 = pd.read_csv(
                rf"D:\user\data\fluphenotype\script\convert_site\HA_NA_mapdir\HA2/H3_{protein}.txt", sep = "\t",
                header = None,
                names = ['H3', protein])
            convert_to_h3_dict_ha2 = dict(zip(mapping_data_HA2[protein], mapping_data_HA2['H3']))

            combined_dict = {
                'HA1': convert_to_h3_dict_ha1,
                'HA2': convert_to_h3_dict_ha2
            }
            # print('-'*50)
            # # print(convert_to_h3_dict_ha1)
            # # print(convert_to_h3_dict_ha2)
            # print(combined_dict)
            # print(protein)
            # print(combined_dict)
            residues = map_residues_to_h3(protein, marker_dict, combined_dict)

            if "H3" in updated_marker_dict:
                updated_marker_dict["H3"].extend(residues)
            else:
                updated_marker_dict["H3"] = residues
            del updated_marker_dict[protein]  # del key
        elif protein in NA_TYPES:
            if os.path.isfile(f"D:/user/data/fluphenotype/script/convert_site/HA_NA_mapdir/NA/N2_{protein}.txt"):
                mapping_data = pd.read_csv(rf"D:\user\data\fluphenotype\script\convert_site\HA_NA_mapdir\NA/N2_{protein}.txt", sep = "\t", header = None,
                                           names = ['N2', protein])
            else:
                mapping_data = pd.read_csv(rf"D:\user\data\fluphenotype\script\convert_site\HA_NA_mapdir\/NA/{protein}_N2.txt", sep = "\t", header = None,
                                           names = [protein, 'N2'])
            convert_to_n2_dict = dict(zip(mapping_data[protein], mapping_data['N2']))

            residues = map_residues_to_h3(protein, marker_dict, convert_to_n2_dict)
            if "N2" in updated_marker_dict:
                updated_marker_dict["N2"].extend(residues)
            else:
                updated_marker_dict["N2"] = residues
            del updated_marker_dict[protein]  # del key
    return updated_marker_dict




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


# 用于测试和理解代码思路的
# total_markers = defaultdict(list)
# for pro, lst in new_protein_dict.items():
#     for dic in lst:
#         # 通过convert_HA_residues全部都会变成H3的，没有影响
#         total_markers[pro].append(convert_HA_residues(dic, "../data/structure"))
# print('-' * 50)
# print(total_markers)
# for i, j in total_markers.items():
#     for s in j:
#         # if (type(s) == str):
#         #     print(s)
#         # if compare_dicts_updated(s,dic2 ):
#         if s and all(s.values()):
#             markers_formated = process_dictionary(s)
#             print(markers_formated)
#             # print('-------------')
#             # print(s.values())
#             # print(s)

# res = process_dictionary(s)
# print(res)
# for s in t:
#     res = process_dictionary(s)
#     print(res)

"""把HA NA单个的标志物转变为H3 N2"""
# phenotype = "drug_resistance"
# data = pd.read_csv(f"../../data/markers_for_extract/{phenotype}_formated.csv")
phenotype = "mammalian_virulence"
data = pd.read_csv(f"test_formated.csv")
print(data)
print(data.loc[:, "Protein Type"].tolist())
# # 假设 HA_TYPES 和 NA_TYPES 是预定义的列表
# HA_TYPES.append('H3')
# NA_TYPES.append('N2')


data_with_hana = data[
    data.loc[:, "Protein Type"].isin(HA_TYPES) | data.loc[:, "Protein Type"].isin(NA_TYPES)]

data_without_hana = data[
    ~(data.loc[:, "Protein Type"].isin(HA_TYPES) | data.loc[:, "Protein Type"].isin(NA_TYPES))]
# data_with_combination.loc[:, "Residue"] = data_with_combination.loc[:, "Mutation"].apply(
#     lambda x: re.search('[A-Z-]', x).group())
# data_with_combination.loc[:, "Mutation"] = data_with_combination.loc[:, "Mutation"].apply(
#     lambda x: re.search('\d+', x).group())
# print(data_with_combination)


# 文件路径
std_fasta_path = "std.fasta"
complete_std_fasta_path = "complete_std.fasta"

# 比较并获取长度差异
length_diffs = compare_sequences(std_fasta_path, complete_std_fasta_path)

for idx, row in data_with_hana.iterrows():
    dic = {row["Protein Type"]: row['Mutation']}
    print(dic)
    newdic = convert_HA_residues(dic, r"D:/user/data/fluphenotype/script/flupre/data/HA_NA_mapdir")
    print(newdic)
    print('-' * 50)

    # 检查 newdic 是否为空
    if list(newdic.values())[0]:
        # print('-'*50)
        # print(dic)
        # print(newdic)
        data_with_hana.loc[idx, "Protein Type"] = list(newdic.keys())[0]
        data_with_hana.loc[idx, "Protein"] = list(newdic.keys())[0]
        data_with_hana.loc[idx, "Mutation"] = list(newdic.values())[0][0]
    else:
        # 如果 newdic 为空，可以选择跳过或赋予默认值
        # 例如：跳过当前迭代
        continue
        # 或者赋予默认值
        # data_with_hana.loc[idx, "Protein Type"] = 默认值
        # data_with_hana.loc[idx, "Mutation"] = 默认值

print(data_with_hana)

converted_data = pd.concat([data_with_hana,data_without_hana])
print(converted_data)
converted_data.to_csv(f"test2_formated.csv",index = False)
# converted_data.to_csv(f"../../data/markers_for_extract/{phenotype}_formated.csv",index = False)


# 以下为一些函数的探索过程
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
#     H3_dict = {}
#     # print(markers)
#     hatype = ""
#     for marker in markers:
#         # HA1也应该处理
#         if "HA2" in marker and protein in HA_TYPES:
#             H3_dict = convert_to_h3_dict["HA2"]
#             marker = marker.split("-")[-1]
#             hatype = "HA2"
#
#             # print('-' * 50)
#             # print(H3_dict)
#             # print(position)
#         elif protein in HA_TYPES:
#             H3_dict = convert_to_h3_dict["HA1"]
#             if "HA1" in marker:
#                 hatype = "HA1"
#             marker = marker.split("-")[-1]
#
#             # print('-' * 50)
#             # print(H3_dict)
#             # print(position)
#
#         elif protein in NA_TYPES:
#             H3_dict = convert_to_h3_dict
#         # 使用re.search查找匹配，并检查是否找到匹配项
#         search_result = re.search(r"\d+[A-Z]", marker)
#         if search_result:
#             # 获取匹配的字符串
#             marker = search_result.group()
#             map_dic = {"496S": "158S", "409P": "65P", "434G": "90G", "445G": "101G", "425G": "81G", "425M": "79M",
#                        "452T": "111T"}
#
#             # 检查marker是否在map_dic字典中
#             if marker in map_dic:
#                 marker = map_dic[marker]
#                 # 使用re.fullmatch来进一步验证格式
#                 marker_match = re.fullmatch(r"(\d+)([A-Z])", marker)
#
#
#         marker_match = re.fullmatch(r"(\d+)([A-Z])", marker)
#         if not marker_match:
#             # print(f"Warning: marker '{marker}' is not in the correct format and will be skipped.")
#             continue
#
#         position, amino_acid = marker_match.groups()
#         # position.map[]
#         # print(position)
#         # print(amino_acid)
#         if protein in HA_TYPES and not hatype:
#             minus = length_diffs[protein]
#             position = str(int(position) - minus)
#
#         # print(marker_match)
#         h3_position = H3_dict.get(position)
#         print(position)
#         print(h3_position)
#         # print(H3_dict)
#         if h3_position is None:
#             # print(f"Warning: Position {position} does not have an H3 mapping "
#             #       f"in the structure comparison file and will be skipped.")
#             continue
#         if hatype:
#             hatype = hatype+"-"
#         mapped_residues.append(f"{hatype}{h3_position}{amino_acid}")
#
#     return mapped_residues
# def map_residues_to_h3(protein, marker_dict, convert_to_h3_dict):
#     markers = marker_dict[protein]
#     if isinstance(markers, str):
#         markers = [markers]
#
#     map_dic = {"496S": "158S", "409P": "65P", "434G": "90G", "445G": "101G", "425G": "81G", "425M": "79M", "452T": "111T"}
#     markers = [map_dic.get(marker, marker) for marker in markers]
#
#     mapped_residues = []
#     for marker in markers:
#         hatype = None
#         if protein in HA_TYPES:
#             if "HA2" in marker:
#                 H3_dict = convert_to_h3_dict["HA2"]
#                 hatype = "HA2"
#             else:
#                 H3_dict = convert_to_h3_dict["HA1"]
#                 if "HA1" in marker:
#                     hatype = "HA1"
#         elif protein in NA_TYPES:
#             H3_dict = convert_to_h3_dict
#
#         marker = marker.split("-")[-1]
#         marker_match = re.fullmatch(r"(\d+)([A-Z])", marker)
#         if not marker_match:
#             continue
#
#         position, amino_acid = marker_match.groups()
#         if not hatype and protein in HA_TYPES:
#             minus = length_diffs[protein]
#             position = str(int(position) - minus)
#
#         h3_position = H3_dict.get(position)
#         if h3_position is None:
#             continue
#
#         hatype_prefix = f"{hatype}-" if hatype else ""
#         mapped_residues.append(f"{hatype_prefix}{h3_position}{amino_acid}")
#
#     return mapped_residues


# def renumber_proteins(fasta_path, acc_pro_dict, marker_dict):
#     """
#     Perform global alignment of protein sequences from a FASTA file against standard sequences and renumber them.
#     Specifically renumber the positions of proteins corresponding to H1-H18 (excluding H3)
#     based on their subtype standard sequences.
#
#     Parameters:
#         fasta_path (str): Path to the FASTA file containing protein sequences.
#         acc_pro_dict (dict): Dictionary mapping accession IDs to protein abbreviations.
#         marker_dict (dict): Dictionary with protein markers.
#
#     Returns:
#         dict: A dictionary with protein IDs and their renumbered positions.
#     """
#     fasta_sequences = SeqIO.parse(fasta_path, 'fasta')
#     renumbering_results = {}
#     proteins = ['HA1' for i in acc_pro_dict.keys() if i in [f"H{i}" for i in range(1, 19)]]
#     for record in fasta_sequences:
#         HA_results = defaultdict(dict)
#         protein_id = record.id
#         protein_abbr = acc_pro_dict.get(protein_id)
#         is_hana_type = protein_abbr in HA_TYPES or protein_abbr in NA_TYPES
#         if protein_abbr in marker_dict or is_hana_type:
#             try:
#                 if protein_abbr in [f"H{i}" for i in range(1, 19)]:
#                     # Construct the path to the standard sequence file
#                     standard_seq_path_HA1 = os.path.join(STANDARD_PATH, f"HA1/{protein_abbr}.fas")
#                     standard_seq_HA1 = next(SeqIO.parse(standard_seq_path_HA1, 'fasta')).seq
#
#                     # Perform global alignment
#                     alignments_HA1 = pairwise2.align.globalxx(standard_seq_HA1, record.seq)
#                     best_alignment_HA1 = max(alignments_HA1, key = lambda x: x.score)
#
#                     # Construct the path to the standard sequence file
#                     standard_seq_path_HA2 = os.path.join(STANDARD_PATH, f"HA2/{protein_abbr}.fas")
#                     standard_seq_HA2 = next(SeqIO.parse(standard_seq_path_HA2, 'fasta')).seq
#
#                     # Perform global alignment
#                     alignments_HA2 = pairwise2.align.globalxx(standard_seq_HA2, record.seq)
#                     best_alignment_HA2 = max(alignments_HA2, key = lambda x: x.score)
#                     # 返回的是列表
#                     HA_results["HA1"][protein_abbr] = renumber_sequence(best_alignment_HA1)
#                     HA_results["HA2"][protein_abbr] = renumber_sequence(best_alignment_HA2)
#
#                     continue
#                 path = f"{protein_abbr}.fas"
#                 if is_hana_type:
#                     path = f"NA/{protein_abbr}.fas"
#                 # Construct the path to the standard sequence file
#                 standard_seq_path = os.path.join(STANDARD_PATH, path)
#                 standard_seq = next(SeqIO.parse(standard_seq_path, 'fasta')).seq
#
#                 # Perform global alignment
#                 alignments = pairwise2.align.globalxx(standard_seq, record.seq)
#                 best_alignment = max(alignments, key = lambda x: x.score)
#
#                 if is_hana_type:
#                     # Store the renumbered sequence for HA/NA types
#                     HA_results[protein_abbr] = renumber_sequence(best_alignment)
#                 else:
#                     # Store the renumbered sequence for non-HA/NA types
#                     renumbering_results[protein_id] = renumber_sequence(best_alignment)
#
#             except Exception as e:
#                 print(f"An error occurred while processing {protein_id}: {str(e)}")
#         else:
#             print(f"No markers found for {protein_abbr} in the source data.")
#
#         # Convert other HA subtype numbering to H3
#         if is_hana_type:
#             renumbered_positions_HA1 = convert_HA_residues(HA_results["HA1"], STRUCTURE_PATH, hatype = "HA1")
#             renumbered_positions_HA2 = convert_HA_residues(HA_results["HA2"], STRUCTURE_PATH, hatype = "HA2")
#             renumbered_positions = merge_dictionaries(renumbered_positions_HA1, renumbered_positions_HA2)
#
#             pop_num = "H3" if protein_abbr in HA_TYPES else ("N2" if protein_abbr in NA_TYPES else None)
#             renumbered_positions[protein_id] = renumbered_positions.pop(pop_num)
#     renumbering_results.update(renumbered_positions)
#
#     return renumbering_results
