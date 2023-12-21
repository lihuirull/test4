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
    print(convert_to_n2_dict)
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


def annotate_markers(markers_path, STRUCTURE_PATH,hatype=None):
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
    for i, j in marker_dict.items():
        if i in [f"H{i}" for i in range(1, 19)]:
            print(f"Protein:\n{i}")
            print(f"Markers:\n{j}")
            print("-" * 50)
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

s = {'HA1': {'H2': ['1D', '2Q', '3I', '4C', '5I', '6G', '7Y', '8H', '9-', '10N', '11N', '12S', '13T', '14-', '15K', '16V', '17D', '18T', '19I', '20L', '21E', '22R', '23N', '24V', '25T', '26V', '27T', '28H', '29A', '30-', '31D', '32I', '33L', '34E', '35K', '36T', '37H', '38N', '39G', '40K', '41L', '42C', '43K', '44L', '45N', '46G', '47I', '48P', '49P', '50L', '51E', '52L', '53G', '54D', '55C', '56S', '57I', '58A', '59G', '60W', '61L', '62L', '63G', '64N', '65P', '66E', '67C', '68D', '69R', '70L', '71L', '72-', '73V', '74P', '75E', '76W', '77S', '78Y', '79I', '80M', '81E', '82-', '83E', '84N', '85P', '86R', '87-', '88-', '89L', '90C', '91Y', '92P', '93G', '94S', '95F', '96N', '97D', '98Y', '99E', '100E', '101L', '102K', '103H', '104L', '105L', '106S', '107S', '108V', '109-', '110H', '111F', '112E', '113K', '114V', '115K', '116I', '117L', '118P', '119-', '120D', '121-', '122W', '123T', '124Q', '125H', '126-', '127-', '128T', '129G', '130G', '131S', '132-', '133A', '134C', '135A', '136V', '137-', '138G', '139-', '140P', '141S', '142F', '143F', '144R', '145N', '146M', '147V', '148W', '149L', '150T', '151K', '152K', '153G', '154-', '155N', '156Y', '157P', '158V', '159A', '160K', '161-', '162S', '163Y', '164N', '165N', '166T', '167S', '168G', '169E', '170Q', '171M', '172L', '173-', '174I', '175W', '176G', '177V', '178H', '179H', '180P', '181N', '182D', '183E', '184A', '185E', '186Q', '187R', '188-', '189L', '190Y', '191Q', '192-', '193V', '194G', '195T', '196Y', '197V', '198S', '199-', '200-', '201T', '202S', '203-', '204L', '205N', '206K', '207R', '208S', '209-', '210P', '211E', '212I', '213A', '214-', '215R', '216P', '217-', '218V', '219-', '220-', '221-', '222G', '223G', '224R', '225M', '226E', '227F', '228S', '229W', '230T', '231-', '232L', '233-', '234-', '235L', '236D', '237T', '238I', '239-', '240F', '241E', '242S', '243T', '244G', '245N', '246L', '247V', '248A', '249P', '250E', '251Y', '252G', '253F', '254K', '255I', '256S', '257K', '258R', '259G', '260S', '261S', '262G', '263I', '264M', '265K', '266T', '267E', '268G', '269T', '270L', '271-', '272N', '273C', '274E', '275T', '276K', '277C', '278Q', '279T', '280P', '281L', '282G', '283A', '284I', '285N', '286T', '287T', '288L', '289P', '290F', '291H', '292N', '293-', '294H', '295P', '296L', '297T', '298I', '299G', '300E', '301C', '302P', '303K', '304Y', '305V', '306K', '307S', '308E', '309-', '310L', '311V', '312L', '313A', '314T', '315G', '316L', '317R', '318N', '319V', '320P', '321Q', '322I', '323E', '324S', '325R']}, 'HA2': {'H2': ['1G', '2L', '3F', '4G', '5A', '6I', '7A', '8G', '9F', '10I', '11E', '12G', '13G', '14W', '15Q', '16G', '17M', '18V', '19D', '20G', '21W', '22Y', '23G', '24Y', '25H', '26H', '27S', '28N', '29D', '30Q', '31G', '32S', '33G', '34Y', '35A', '36A', '37D', '38K', '39E', '40S', '41T', '42Q', '43K', '44A', '45-', '46D', '47-', '48I', '49T', '50N', '51K', '52V', '53N', '54S', '55V', '56I', '57E', '58K', '59M', '60N', '61T', '62Q', '63F', '64E', '65A', '66V', '67G', '68K', '69E', '70F', '71S', '72N', '73L', '74E', '75K', '76R', '77L', '78E', '79N', '80L', '81N', '82K', '83K', '84M', '85E', '86D', '87G', '88F', '89L', '90D', '91V', '92W', '93T', '94Y', '95N', '96A', '97E', '98L', '99L', '100V', '101L', '102M', '103E', '104N', '105E', '106R', '107T', '108L', '109D', '110F', '111H', '112D', '113S', '114N', '115V', '116K', '117N', '118L', '119Y', '120D', '121K', '122V', '123R', '124M', '125Q', '126L', '127R', '128D', '129N', '130-', '131K', '132E', '133L', '134G', '135N', '136G', '137C', '138F', '139E', '140F', '141Y', '142H', '143K', '144C', '145D', '146D', '147E', '148C', '149M', '150N', '151S', '152V', '153K', '154N', '155G', '156T', '157Y', '158D', '159Y', '160P', '161K', '162Y', '163E', '164E', '165E', '166S', '167K', '168L', '169-', '170R', '171N', '172E', '173I', '174K', '175G', '176V', '177K', '178L', '179-', '180-', '181M', '182G', '183V', '184Y', '185Q', '186I', '187L', '188A', '189I', '190Y', '191A', '192T', '193V', '194A', '195G', '196S', '197L', '198S', '199L', '200A', '201I', '202M', '203I', '204A', '205G', '206I', '207S', '208-', '209W', '210M', '211C', '212S', '213N', '214G', '215S', '216L', '217Q', '218C', '219R', '220I', '221C', '222I']}}
print(s["HA1"])
# marker_dict = {'M2': '41C'}
# for protein in list(marker_dict.keys()):
#     values = marker_dict[protein]
#     if not isinstance(values, str):
#         print(list(set(values)))

# s =  {'H3': {'HA1': '218W', 'HA2': '156N'}}
# s1 =  {'H3': {'HA1': '218W', 'HA2': ['156N','222L']}}
s = {'PB2': '158G', 'PA': '295P'}
s1 = {'PB2': '158G', 'PA': ['295P','225L']}
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
print(is_subset_complex_revised(s1, s))


def map_residues_to_h3(protein, marker_dict, convert_to_h3_dict):
    """
    Maps the residue numbers for a given protein to the H3/N2 numbering system.

    Parameters:
        protein (str): The protein identifier.
        marker_dict (dict): Dictionary containing markers for various proteins.
        convert_to_h3_dict (dict): Dictionary that maps residue numbers to H3.

    Returns:
        list: A list of residues mapped to H3 numbering system.
    """
    markers = marker_dict[protein]

    # 如果markers是字符串，将其转换为只含一个元素的列表
    if isinstance(markers, str):
        markers = [markers]

    mapped_residues = []
    for marker in markers:
        # Ensure the marker is in the expected format (e.g., "12A")
        marker_match = re.fullmatch(r"(\d+)([A-Z])", marker)
        if not marker_match:
            # print(f"Warning: marker '{marker}' is not in the correct format and will be skipped.")
            continue

        position, amino_acid = marker_match.groups()
        h3_position = convert_to_h3_dict.get(position)
        if h3_position is None:
            # print(f"Warning: Position {position} does not have an H3 mapping "
            #       f"in the structure comparison file and will be skipped.")
            continue

        mapped_residues.append(h3_position + amino_acid)

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
            mapping_data = pd.read_csv(f"{structure_folder}/H3_{protein}.txt", sep = "\t", header = None,
                                       names = ['H3', protein])
            convert_to_h3_dict = dict(zip(mapping_data[protein], mapping_data['H3']))

            residues = map_residues_to_h3(protein, marker_dict, convert_to_h3_dict)
            # residues = [convert_to_H3_dict.get(re.search(r"\d+", i).group()) +
            # re.search(r"[A-Z]", i).group() for i in
            #             marker_dict[protein] if convert_to_H3_dict.get(re.search(r"\d+", i).group())]
            if "H3" in updated_marker_dict:
                updated_marker_dict["H3"].extend(residues)
            else:
                updated_marker_dict["H3"] = residues
            del updated_marker_dict[protein]  # del key
        elif protein in NA_TYPES:
            if os.path.isfile(f"{structure_folder}/N2_{protein}.txt"):
                mapping_data = pd.read_csv(f"{structure_folder}/N2_{protein}.txt", sep = "\t", header = None,
                                           names = ['N2', protein])
            else:
                mapping_data = pd.read_csv(f"{structure_folder}/{protein}_N2.txt", sep = "\t", header = None,
                                           names = [protein, 'N2'])
            convert_to_n2_dict = dict(zip(mapping_data[protein], mapping_data['N2']))

            residues = map_residues_to_h3(protein, marker_dict, convert_to_n2_dict)
            if "N2" in updated_marker_dict:
                updated_marker_dict["N2"].extend(residues)
            else:
                updated_marker_dict["N2"] = residues
            del updated_marker_dict[protein]  # del key

    return updated_marker_dict


def annotate_markers(markers_path):
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

    # Convert HA/NA residues to H3/N2 numbering and update marker_dict
    marker_dict = convert_HA_residues(marker_dict, STRUCTURE_PATH)

    # Duplicated
    marker_dict = {i: list(set(j)) for i, j in marker_dict.items()}

    return marker_dict, data
