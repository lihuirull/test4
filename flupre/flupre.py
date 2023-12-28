#!/data/penglab3-20T/lihuiru/miniconda3/bin/python3
# -*coding: utf-8 -*-
# Author:lihuiru
# Created on 2023/11/7 10:58

import argparse
import os
import re
import subprocess
import sys
from collections import Counter
from collections import defaultdict
import datetime
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd

from .blastHASeq import blastHASeq
# 模块与函数重名时，最好用.模块的语法导入（这样不会显示不可callable）
from .blastSeq import blastSeq
from . import predict_host
from . import predict_virulence
from .changeStandardNameForMPNS import refreshStandardName
from .getBlastMostCommonHitProteinType import getMostCommonHitProtein, getMostCommonHitProteinLowLevelHost
from .getSeq import get
from .predictVirusHost import getVirusHost, getVirusHostLowLevel
from .standardize import standardize
from .translateDNA2Protein import translate, makeProteinFileForDownload

pd.set_option('display.max_columns', None)

HA_TYPES = [f"H{i}" for i in range(1, 19) if i != 3]
NA_TYPES = [f"N{i}" for i in range(1, 10) if i != 2]
length_diffs = {'H1': 17, 'H10': 17, 'H11': 16, 'H12': 17, 'H13': 18, 'H14': 17, 'H15': 18, 'H16': 19, 'H17': 18,
                'H18': 14, 'H2': 15, 'H3': 16, 'H4': 16, 'H5': 12, 'H6': 16, 'H7': 18, 'H8': 17, 'H9': 18}
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DB_PATH = os.path.join(base_dir, 'data', 'flu_db.dmnd')
STD_PATH = os.path.join(base_dir, 'data', 'std.fasta')
COMPLETE_STD_PATH = os.path.join(base_dir, 'data', 'complete_std.fasta')
STRUCTURE_PATH = os.path.join(base_dir, 'data', 'HA_NA_mapdir')
STANDARD_PATH = os.path.join(base_dir, 'data', 'standard_seq_protein')
MODEL_PATH = os.path.join(base_dir, 'model')
DATA_PATH = os.path.join(base_dir, 'data')
MARKER_PATH = os.path.join(base_dir, 'data', 'markers_for_extract')
RESULT_PATH = os.path.join(base_dir,'result')
TEMP_PATH = os.path.join(base_dir,'temp')

def process_antigen_text(file_path):
    """
    处理抗原信息文件的文本行。
    :param file_path: 抗原信息文件路径
    :return: 处理后的文本行列表
    """
    with open(file_path, 'r') as AntigenFile:
        text = AntigenFile.readlines()

    # 对text进行排序
    text.sort(key=lambda x: calculate_sort_key(x))

    # 第一次格式化text
    text = [format_line_stage1(line) for line in text]

    # 第二次格式化text
    text = [format_line_stage2(line) for line in text]

    return text

def calculate_sort_key(line):
    """
    根据提供的排序逻辑计算排序键。
    """
    parts = line.split('\t')
    try:
        value = parts[2]
        if "Max" in value:
            return 10000000000 + len(line)
        elif "Min" in value:
            return -100000 + len(line)
        else:
            return float('{:.10f}'.format(float(value)))
    except IndexError:
        return float('inf')

def format_line_stage1(line):
    """
    对文本行进行第一次格式化。
    """
    parts = line.split('\t')
    try:
        value = parts[2]
        if "Max" not in value and "Min" not in value:
            # 直接给parts[2]赋值，不额外添加'\t'
            parts[2] = '{0:.2e}'.format(float(value))
        # 使用'\t'.join(parts)来合并，自然会在每个部分之间插入'\t'
        return '\t'.join(parts)

    except IndexError:
        return line  # 如果格式不正确，返回原行

def format_line_stage2(line):
    """
    对文本行进行第二次格式化。
    """
    parts = line.split('\t')
    try:
        value = parts[2]
        if "Max" not in value and "Min" not in value:
            value_float = float(value)
            prefix = "Similar_" if value_float <= 1.0 else "Variant_"
            parts[2] = prefix + value
        elif "Min" in value:
            parts[2] = "similar_" + value
        elif "Max" in value:
            parts[2] = "variant_" + value
        return '\t'.join(parts)
    except IndexError:
        return line  # 如果格式不正确，返回原行


def ivew_task(resultFileDir, tempFileDir, inputFilePath):
    beginTime = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    DIR = base_dir
    print(DIR)
    workName = os.path.basename(inputFilePath)

    fileLog = open(resultFileDir + "/" + workName + ".result", "w", encoding = 'utf-8')
    fileLog.write("\n######################\t" + beginTime + "begin\n")
    ##########################################################################
    # 标准化序列并返回类型（DNA/蛋白质）
    inputFile, querySeqType, dic = standardize(inputFilePath, outFileDir = tempFileDir)

    fileLog.write("standard_name:original_name->" + str(dic) + "\n")
    ##########################################################################

    ##########################################################################
    dataBaseName = None
    if querySeqType == "nucleo":
        dataBaseName = "/protType_NA"  # blast
    if querySeqType == "protein":
        dataBaseName = "/protType_AA"  # blast
    querySeqFile = inputFile
    querySeqFileDir = tempFileDir
    DBDir = DIR + "/18Mid/standard_seq/allProteinTypeDB/"  # blast
    queryEValue = "1e-5"  # blast
    outfmt = "3"  # blast
    outName = "querySeqToProteinTypeDB"  # blast

    blastSeq(querySeqFile, querySeqType, DBDir, dataBaseName, queryEValue, outfmt, outName, querySeqFileDir,
             tempFileDir)  #
    predictedProteinType = getMostCommonHitProtein(outName, tempFileDir)
    print(predictedProteinType)  #
    fileLog.write("ProteinType->" + str(predictedProteinType) + "\n")
    for eachSeqType in predictedProteinType:
        if 'Unknown' in eachSeqType:
            fileLog.write('Attention! One or more sequence inputted may not belong to influenza viruses\n')
    times = Counter(num[1] for num in predictedProteinType).most_common()
    for time in times:
        if time[1] > 1:
            fileLog.write("Warning : You've entered " + str(time[1]) + " " + str(time[0]) + " proteins in\n")

    if querySeqType == "nucleo":
        dataBaseName = "/HostNuclDB"  # blast
    if querySeqType == "protein":
        dataBaseName = "/HostProtDB"  # blast
    querySeqFile = inputFile
    querySeqFileDir = tempFileDir
    DBDir = DIR + "/18Mid/standard_seq/allProteinTypeDB/"
    # print(DBDir)
    queryEValue = "1e-5"  # blast
    outfmt = "3"  # blas
    outName = "querySeqToHostDB"  # blast
    tempFileDir = tempFileDir  # blas
    blastSeq(querySeqFile, querySeqType, DBDir, dataBaseName, queryEValue, outfmt, outName, querySeqFileDir,
             tempFileDir)  #
    Host = getMostCommonHitProtein(outName, tempFileDir)
    HostLowLevel = getMostCommonHitProteinLowLevelHost(outName, tempFileDir, Host)
    virusHostLowLevel, HostNew2 = getVirusHostLowLevel(predictedProteinType, HostLowLevel)
    virusHost, HostNew = getVirusHost(predictedProteinType, Host)
    fileLog.write("virusHostLowLevel->" + str(virusHostLowLevel) + "\n")
    fileLog.write("VirusHost->" + str(virusHost) + "\n")
    ##########################################################################

    ##########################################################################
    # 获取对每条DNA序列翻译后注释的结果（包含蛋白类型和长度）
    # CDSList:[('querySeq1', 'HA(1..1698)'), ('querySeq2', 'NA(1..1419)'), ('querySeq3', 'NP(1..1494)'),
    # ('querySeq4', 'NS1(1..711),NS2(1..30,503..838)'), ('querySeq5', 'PA(1..2148),PA-X(1..570,572..760)'),
    # ('querySeq6', 'PB1(1..2274),PB1-F2(95..364)'), ('querySeq7', 'PB2(1..2277)'),
    # ('querySeq8', 'M1(1..816),M2(1..26,715..982)'), ('querySeq9', 'HA(1..1698)')]
    CDSList = []
    # proteinPath = ""
    if querySeqType == "nucleo":
        # 返回：outputName.replace(".fas", ".fas2"), outFileDir, dicCDS
        queryProteinSeqFile = translate(querySeqFile, querySeqFileDir, DIR, tempFileDir)
        querySeqFile = queryProteinSeqFile[0]

        # 新的蛋白文件路径
        # proteinPath = os.path.join(queryProteinSeqFile[1], querySeqFile)

        fileLog.write("CDS->" + str(queryProteinSeqFile[2]) + "\n")
        CDSList = queryProteinSeqFile[2]

    ###########################################################################
    proteinPath = makeProteinFileForDownload(tempFileDir, querySeqFile, resultFileDir, dic, predictedProteinType)
    print(f"proteinPath:\n{proteinPath}")
    ##########################################################################
    # 根据传入的dic，为query添加蛋白注释信息同时和原本序列id对应。
    # 拆解CDSList，这里主要是为了区分一条序列对应多个蛋白（NS, MP, PB1, PA)
    # 把('querySeq6', 'PB1(1..2274),PB1-F2(95..364)'),
    # 变成--》'>querySeq6_PB1': '>standard_AB434294', '>querySeq6_PB1-F2': '>standard_AB434294'
    # 不区分不同蛋白的无需为query添加蛋白注释信息

    # CDS->[('querySeq1', 'HA(1..1698)'), ('querySeq2', 'NA(1..1419)'), ('querySeq3', 'NP(1..1494)'), ('querySeq4', 'NS1(1..711),NS2(1..30,503..838)'), ('querySeq5', 'PA(1..2148),PA-X(1..570,572..760)'), ('querySeq6', 'PB1(1..2274),PB1-F2(95..364)'), ('querySeq7', 'PB2(1..2277)'), ('querySeq8', 'M1(1..816),M2(1..26,715..982)'), ('querySeq9', 'HA(1..1698)')]
    # standard_name_new:original_name->{'>querySeq1': '>standard_DQ992756', '>querySeq2': '>standard_CY078869_modified', '>querySeq3': '>standard_CY079270', '>querySeq7': '>standard_CY101570', '>querySeq9': '>standard_KJ200805', '>querySeq4_NS1': '>standard_CY004342', '>querySeq4_NS2': '>standard_CY004342', '>querySeq5_PA': '>standard_CY102193_modified', '>querySeq5_PA-X': '>standard_CY102193_modified', '>querySeq6_PB1': '>standard_AB434294', '>querySeq6_PB1-F2': '>standard_AB434294', '>querySeq8_M1': '>standard_CY005795', '>querySeq8_M2': '>standard_CY005795'}
    # ProteinType->[('querySeq1', 'H5'), ('querySeq2', 'N5'), ('querySeq3', 'NP'), ('querySeq4_NS1', 'NS1'), ('querySeq4_NS2', 'NS2'), ('querySeq5_PA', 'PA'), ('querySeq5_PA-X', 'PA-X'), ('querySeq6_PB1', 'PB1'), ('querySeq6_PB1-F2', 'PB1-F2'), ('querySeq7', 'PB2'), ('querySeq8_M1', 'M1'), ('querySeq8_M2', 'M2'), ('querySeq9', 'H6')]

    if querySeqType == "nucleo":
        predictedProteinType = refreshStandardName(predictedProteinType, dic)
        fileLog.write("standard_name_new:original_name->" + str(dic) + "\n")
        fileLog.write("ProteinType->" + str(predictedProteinType) + "\n")
    ##########################################################################
    CDSTypeList = []
    print(CDSList)
    for eachCDSType in CDSList:
        print(eachCDSType)
        for eachType in eachCDSType[1].split('),'):
            if eachType.split('(')[0] in ['NS1', 'NS2', 'M1', 'M2', 'PB1', 'PB1-F2', 'PA-X', 'PA']:
                CDSTypeList.append(eachCDSType[0] + '_' + eachType.split('(')[0])
            else:
                CDSTypeList.append(eachCDSType[0])
    print(CDSTypeList)
    # 完全一致
    # print([i[0] for i in predictedProteinType])
    ##########################################################################
    for eachQuery in predictedProteinType:
        if (querySeqType == "nucleo") and (eachQuery[0] not in CDSTypeList): continue
        print(dic)
        print(eachQuery)
        fileLog.write("\n" + eachQuery[0] + "\t" + dic[">" + eachQuery[0]] + "\t" + eachQuery[1].rstrip() + "\n")
        if str(eachQuery[1]).strip("\n") == "Unknown":
            fileLog.write("Virulent Site:\n" + "\tNo result" + "\n")
            continue
        querySeq, querySeqFileName = get(fileName = querySeqFile, dir = querySeqFileDir, seqName = eachQuery[0])
        mafftDir = DIR + "/app/mafft/mafft-7.158-without-extensions/scripts/"
        ##########################################################################
        # 计算和疫苗株的抗原距离
        if eachQuery[1] in ["H1", "H3", "H5", "H7", "H9"]:
            subType = eachQuery[1]
            querySeqFileDir = tempFileDir
            querySeqName = eachQuery[0]
            print(querySeqName)
            # 虽然此处的librarydir是/18Mid/antigen/modelData/，但其实predAV.pl只用了model和referseq目录
            # 而针对这五种亚型的HA，还需传入疫苗株的序列，去计算遗传距离
            os.system(
                "perl " + DIR + "/18Mid/antigen/predAV.pl" + " " + subType + " " + querySeqFileDir + querySeqFileName + " " + DIR + "/18Mid/antigen/modelData/vaccine/" + subType + " " + tempFileDir + querySeqFileName + ".antigenInfo" + " " + DIR + "/18Mid/antigen/modelData/" + " " + tempFileDir)
            try:
                AntigenFile = tempFileDir + querySeqFileName + ".antigenInfo"
                text = process_antigen_text(AntigenFile)
                fileLog.write("AntigenInfo\n")
                for each in text:
                    fileLog.write(each.replace('\tX', '\t-').replace('X\t', '-\t'))
            except Exception as e:
                print(e)
            # Add a module of antigen similarity calculation
        if eachQuery[1] in [f"H{i}" for i in range(1, 17)]:
            # 针对所有亚型，传入自己HA的序列在一个文件中），通过blastp比对确定五条最相似的序列并传入这五条序列（在一个文件中），
            # 计算两两遗传距离
            subType = eachQuery[1]
            querySeqFileDir = tempFileDir
            querySeqName = eachQuery[0]
            print(querySeqName)
            try:
                blastHASeq(prefixDir = DIR, blastQuerySeqHA = querySeqFileName, HAType = subType, tempDir = tempFileDir,
                           mafftDir = mafftDir)
                AntigenHAFile = tempFileDir + querySeqFileName + ".antigenInfo.blast"

                text = process_antigen_text(AntigenHAFile)

                fileLog.write("AntigenInfo_Blast\n")
                for each in text:
                    fileLog.write(each.replace('\tX', '\t-').replace('X\t', '-\t'))
            except Exception as e:
                print(e)
    return proteinPath, resultFileDir + workName + ".result"


def run_diamond_blast(input_fasta_file, output_path, threads, evalue = 1e-5, suppress_output = False):
    cmd = [
        "diamond", "blastp",
        "-d", DB_PATH,
        "-q", input_fasta_file,
        "-o", output_path,
        "-f", '6',
        "-p", str(threads),
        "-e", str(evalue),
        "--sallseqid",
        "--salltitles"
    ]
    if suppress_output:
        subprocess.call(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
    else:
        subprocess.call(cmd)
        print("\nDIAMOND BLAST completed.\n")


def read_annotation_results(output_path, threshold):
    # Read in the alignment file
    data = pd.read_csv(output_path, sep = "\t", header = None,
                       names = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                "qstart", "qend", "sstart", "send", "evalue", "bitscore"])

    # Filter based on the evalue threshold
    data = data[data["evalue"] <= threshold]
    # Extract the highest bitscore hits
    best_hits_idx = data.groupby("qseqid")["bitscore"].idxmax()
    best_hits = data.loc[best_hits_idx, ["qseqid", "sseqid"]]
    return best_hits


def map_accession_to_protein(best_hits):
    # Load the mapping dictionary
    protein_sequences = SeqIO.parse(STD_PATH, "fasta")
    id_pro_dic = {}
    for record in protein_sequences:
        protein_type = record.id.split("_")[0]
        id_pro_dic[record.id] = protein_type
    # Merge the best hits with the mapping
    best_hits = best_hits.merge(pd.DataFrame(list(id_pro_dic.items()), columns = ["sseqid", "Protein Abbreviation"]),
                                on = "sseqid", how = "left")
    return best_hits[['qseqid', 'Protein Abbreviation']]


def update_fasta_csv_with_annotations(input_fasta_path, annotations, output_directory, prefix, suppress_output):
    if not suppress_output:
        # Read the FASTA file
        records = list(SeqIO.parse(input_fasta_path, "fasta"))
        input_fasta_filename = os.path.split(input_fasta_path)[1]

        # Create a mapping dictionary from annotations DataFrame
        annotations_dict = annotations.set_index('qseqid')['Protein Abbreviation'].to_dict()

        # Update the record IDs with annotations
        for record in records:
            protein_info = annotations_dict.get(record.id)
            if protein_info:
                # Append the protein information only to the description field
                record.id = f"{record.id}_{protein_info} {record.description.split(' ', 1)[-1]}"
                record.description = ""

        # Write the annotated sequences to a new FASTA file
        output_fasta_filename = f"{prefix}{input_fasta_filename.split('.')[0]}_annotated.fasta"

        output_fasta_path = f"{output_directory}/{output_fasta_filename}"

        with open(output_fasta_path, "w") as output_handle:
            SeqIO.write(records, output_handle, "fasta")
        print("\nFASTA file updated with annotations.")


def annotate_fasta_file(input_fasta_path, output_directory = ".", prefix = "", evalue = 1e-5,
                        update_file = True, threads = 10, suppress_output = False):
    """
    Annotate a FASTA file using DIAMOND BLAST against a flu database.

    Parameters:
        input_fasta_path(str): Path to the input FASTA file or a directory.
        output_directory (str): Directory where the output files will be saved.
        prefix (str): Prefix to be added to the output filenames.
        evalue (float): E-value threshold for filtering BLAST hits.
        update_file (bool): If True, update the FASTA file with annotations.
        threads (int): Number of parallel threads.
        suppress_output (bool): If True, suppress output.
    """
    os.makedirs(output_directory, exist_ok = True)
    input_fasta_filename = os.path.split(input_fasta_path)[1]
    add_prefix = prefix + "_" if prefix else ""
    output_filename = f"{add_prefix}{input_fasta_filename.split('.')[0]}.aln"
    output_path = f"{output_directory}/{output_filename}"
    print(output_path)

    # Run DIAMOND BLAST
    run_diamond_blast(input_fasta_path, output_path, threads, evalue, suppress_output)

    # Read and process the BLAST results
    best_hits = read_annotation_results(output_path, evalue)
    annotations = map_accession_to_protein(best_hits)

    # Update the CSV file with annotations
    output_csv_filename = f"{add_prefix}{input_fasta_filename.split('.')[0]}_annotated.csv"
    output_csv_path = f"{output_directory}/{output_csv_filename}"
    annotations.to_csv(output_csv_path, index = False, encoding = 'utf-8-sig')
    print("CSV file updated with annotations.\n")

    if update_file:
        # Update the FASTA file with annotations
        update_fasta_csv_with_annotations(input_fasta_path, annotations, output_directory, add_prefix, suppress_output)

    return annotations


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


def get_h3_dict_and_hatype(protein, marker, convert_to_h3_dict):
    if protein in HA_TYPES:
        if "HA2" in marker:
            return convert_to_h3_dict["HA2"], "HA2"
        else:
            return convert_to_h3_dict["HA1"], "HA1" if "HA1" in marker else None
    elif protein in NA_TYPES:
        return convert_to_h3_dict, None
    return None, None


def adjust_position_type(position, H3_dict):
    if not H3_dict:
        return position  # 如果字典为空，则直接返回原始位置

    # 获取字典中第一个键的类型
    first_key_type = type(next(iter(H3_dict)))

    # 如果位置已经是正确的类型，则直接返回
    if isinstance(position, first_key_type):
        return position

    # 否则，尝试将位置转换为该类型
    try:
        return first_key_type(position)
    except ValueError:
        return None  # 或者处理转换失败的情况


def adjust_position_and_get_h3_position(marker, hatype, H3_dict, protein):
    # marker_match = re.fullmatch(r"(\d+)([A-Z]|-)", marker)
    # length_diffs = compare_sequences(STD_PATH, COMPLETE_STD_PATH)
    # print(f"adjH2:\n{marker}")
    marker_match = re.search(r"(\d+)([A-Z]|-)", marker)

    if not marker_match:
        return None, None, hatype

    position, amino_acid = marker_match.groups()

    # if not hatype and protein in HA_TYPES:
    if not hatype and protein in HA_TYPES:
        minus = length_diffs[protein]
        position = str(int(position) - minus)
        hatype = "HA1"

    if H3_dict:
        # 处理除H3的情况
        # return H3_dict.get(position), amino_acid, hatype
        adjusted_position = adjust_position_type(position, H3_dict)
        return H3_dict.get(adjusted_position), amino_acid, hatype
    else:
        if not hatype:
            hatype = "HA1"  # 处理标志物字典的H3标志物
        # 处理H3（不需要位点转换）
        return f"{hatype}-{position}{amino_acid}"


def map_residues_to_h3(protein, marker_dict, convert_to_h3_dict, hatype = None):
    markers = [marker_dict[protein]] if isinstance(marker_dict[protein], str) else marker_dict[protein]
    # print(f"H2protein:\n{markers}")
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
        # HA1 H2
        # print(f"con:\n{convert_to_h3_dict}")
        # print(f"H3_dict:\n{H3_dict}")
        h3_position, amino_acid, updated_hatype = adjust_position_and_get_h3_position(marker, hatype, H3_dict, protein)
        # print(f"H3_postion:\n{h3_position}")
        if h3_position is None:
            continue

        hatype_prefix = f"{updated_hatype}-" if updated_hatype else ""
        mapped_residues.append(f"{hatype_prefix}{h3_position}{amino_acid}")

    return mapped_residues


def load_mapping_data(filepath, column_names):
    if os.path.isfile(filepath):
        mapping_data = pd.read_csv(filepath, sep = "\t", header = None, names = column_names)
        return dict(zip(mapping_data[column_names[1]], mapping_data[column_names[0]]))
    return {}


def process_ha_type(protein, marker_dict, structure_folder, hatype):
    convert_to_h3_dict_ha1 = load_mapping_data(
        f"{structure_folder}/HA1/H3_{protein}.txt", ['H3', protein])
    convert_to_h3_dict_ha2 = load_mapping_data(
        f"{structure_folder}/HA2/H3_{protein}.txt", ['H3', protein])

    combined_dict = {'HA1': convert_to_h3_dict_ha1, 'HA2': convert_to_h3_dict_ha2}
    # print(f"combine\n{combined_dict}")
    return map_residues_to_h3(protein, marker_dict, combined_dict, hatype)


def process_na_type(protein, marker_dict, structure_folder, hatype):
    convert_to_n2_dict = load_mapping_data(
        f"{structure_folder}/NA/N2_{protein}.txt", ['N2', protein])
    # print(convert_to_n2_dict)
    return map_residues_to_h3(protein, marker_dict, convert_to_n2_dict, hatype)


def convert_HA_residues(marker_dict, structure_folder, hatype):
    # 传入都是{'H1': ['190R', 'HA2-64H', '186P', 445G', 'HA2-117S'], 'H10': ['220L'], 'PB1': ['250S', '170R']...}
    # 中间过程是：{'H3': ['HA1-190R', 'HA2-64H', 'HA1-186P', 'HA2-45G', 'HA2-117S'], 'PB1': ['220L']...}
    # 经过transform_marker_dict导出都是{'H3': {'HA1': ['138S', '225G', '304T'], 'HA2': ['156N', '107I']},'M1': ['192V']}
    # 此函数处理HA/NA，用于将基于标准序列的编号向H3/N2转换，其他蛋白重编号位点不变（即基于标注序列的编号）
    updated_marker_dict = marker_dict.copy()

    for protein in list(marker_dict.keys()):
        if protein == "H3":
            # 处理特殊情况
            res = []
            if not isinstance(marker_dict["H3"], list):
                marker_dict["H3"] = [marker_dict["H3"]]
            for marker in marker_dict["H3"]:
                if "HA2" not in marker and "HA1" not in marker:
                    # 假设 adjust_position_and_get_h3_position 函数适用于这种情况
                    # 用于识别标志物的序列的编号列表传入(此列表元素不会出现：HA1-XXX或者HA2-XXX，
                    # 因为此列表的整理就是HA1:[],HA2:[]，然后把对应的HA1的列表传入
                    # 不需要对此列表减去信号肽，所以编号时需指定hatype（hatype为HA1或HA）,不减信号肽
                    # 只是为了给这些marker字符串加上HA1/HA2的前缀

                    # 而传入标志物列表时会出现无需处理带有HA1和HA2标识的marker字符串，不带有的则需要减去信号肽，
                    # 给这些marker字符串加上HA1-的前缀
                    marker = adjust_position_and_get_h3_position(marker, hatype = hatype, H3_dict = None,
                                                                 protein = "H3")
                # print(marker)
                res.append(marker)
            updated_marker_dict["H3"] = res

        if protein in HA_TYPES:
            residues = process_ha_type(protein, marker_dict, structure_folder, hatype)
            updated_marker_dict["H3"] = updated_marker_dict.get("H3", []) + residues
            del updated_marker_dict[protein]

        elif protein in NA_TYPES:
            residues = process_na_type(protein, marker_dict, structure_folder, hatype)
            print(f"residues\n{residues}")
            if residues:
                updated_marker_dict["N2"] = updated_marker_dict.get("N2", []) + residues
            del updated_marker_dict[protein]

    # print(marker_dict)
    # print(updated_marker_dict)
    # print('-'*50)
    return transform_marker_dict(updated_marker_dict)


def transform_marker_dict(marker_dict):
    # 这一步只是为了把H3
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
            if not isinstance(values, str):
                transformed_data[key] = list(set(values))
            else:
                transformed_data[key] = values
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

    # Convert HA/NA residues to H3/N2 numbering and update marker_dict
    # Already Duplicated
    # 传入的marker_dict已经经过H3/N2的转换
    marker_dict = convert_HA_residues(marker_dict, STRUCTURE_PATH, hatype)

    # Duplicated
    # marker_dict = {i: list(set(j)) for i, j in marker_dict.items()}
    # print(marker_dict)
    return marker_dict, data


def renumber_sequence(best_alignment):
    """
    Renumber a protein sequence based on the best alignment result.

    Parameters:
        best_alignment (list of tuples): The alignment result between the standard and query sequences.

    Returns:
        list: A list of renumbered positions in the format f'{position}{amino_acid}'.
    """
    # Initialize the list for storing renumbered positions
    renumbered_positions = []
    count = 1  # Start counting positions from 1
    for std_char, query_char in zip(best_alignment[0], best_alignment[1]):
        if std_char != '-':  # Ignore gaps in the standard sequence
            renumbered_positions.append(f"{count}{query_char}")
            count += 1  # Increment the position counter for non-gap characters

    return renumbered_positions


def merge_dictionaries(dict1, dict2):
    merged_dict = {}
    for key in dict1:
        if key == 'H3':
            # 将 'H3' 键下的字典作为两个独立元素放入列表
            merged_h3_values = [dict1[key], dict2[key]]
            merged_dict[key] = merged_h3_values
        else:
            # 对于其他键，直接复制值
            merged_dict[key] = dict1[key]
    return merged_dict


from Bio import SeqIO, pairwise2


def perform_alignment_and_renumber(standard_seq_path, query_seq):
    standard_seq = next(SeqIO.parse(standard_seq_path, 'fasta')).seq
    alignments = pairwise2.align.globalxx(standard_seq, query_seq)
    best_alignment = max(alignments, key = lambda x: x.score)
    return renumber_sequence(best_alignment)


def process_ha_na(protein_abbr, sequence):
    ha_results = defaultdict(dict)
    # if protein_abbr != "H3":
    if protein_abbr in [f"H{i}" for i in range(1, 19)]:
        ha1_path = f"{STANDARD_PATH}/HA1/{protein_abbr}.fas"
        ha_results["HA1"][protein_abbr] = perform_alignment_and_renumber(ha1_path, sequence)

        ha2_path = f"{STANDARD_PATH}/HA2/{protein_abbr}.fas"
        ha_results["HA2"][protein_abbr] = perform_alignment_and_renumber(ha2_path, sequence)
    elif protein_abbr in [f"N{i}" for i in range(1, 10)]:
        na_path = f"{STANDARD_PATH}/{protein_abbr}.fas"
        ha_results[protein_abbr] = perform_alignment_and_renumber(na_path, sequence)
    return ha_results


def renumber_proteins(fasta_path, acc_pro_dict, marker_dict):
    # 针对待提取标志物的序列的重编号，
    fasta_sequences = SeqIO.parse(fasta_path, 'fasta')
    renumbering_results = {}

    for record in fasta_sequences:
        protein_id = record.id
        protein_abbr = acc_pro_dict.get(protein_id)
        is_hana_type = protein_abbr in HA_TYPES or protein_abbr in NA_TYPES
        # print("marker_dict")
        # print(marker_dict)
        if protein_abbr in marker_dict or is_hana_type:
            try:
                if protein_abbr in [f"H{i}" for i in range(1, 19)]:
                    # 给出在对应HA亚型的标准序列的编号
                    ha_results = process_ha_na(protein_abbr, record.seq)
                    # 对HA进行向H3映射的重编号
                    # print(f"ha_results:\n{ha_results}")
                    renumbered_positions_HA1 = convert_HA_residues(ha_results["HA1"], STRUCTURE_PATH, hatype = "HA1")
                    # print(f"RENUMBER:HA1\n{renumbered_positions_HA1}")
                    renumbered_positions_HA2 = convert_HA_residues(ha_results["HA2"], STRUCTURE_PATH, hatype = "HA2")
                    # print(f"RENUMBER:HA2\n{renumbered_positions_HA2}")
                    renumbered_positions = merge_dictionaries(renumbered_positions_HA1, renumbered_positions_HA2)
                    # pop_num = "H3" if protein_abbr in HA_TYPES else "N2"
                    renumbered_positions[protein_id] = renumbered_positions.pop("H3")
                    renumbering_results.update(renumbered_positions)
                elif protein_abbr in [f"N{i}" for i in range(1, 10)]:
                    ha_results = process_ha_na(protein_abbr, record.seq)
                    renumbered_positions = convert_HA_residues(ha_results, STRUCTURE_PATH, hatype = None)
                    renumbered_positions[protein_id] = renumbered_positions.pop("N2")
                    renumbering_results.update(renumbered_positions)
                else:
                    # 处理非 HA/NA 类型
                    standard_seq_path = os.path.join(STANDARD_PATH, f"{protein_abbr}.fas")
                    renumbering_results[protein_id] = perform_alignment_and_renumber(standard_seq_path, record.seq)
            except Exception as e:
                print(f"An error occurred while processing {protein_id}: {str(e)}")
        else:
            print(f"No markers found for {protein_abbr} in the source data.")
    # print(renumbering_results)
    # print('-' * 50)
    return renumbering_results


def merge_dicts_with_list(dict_list):
    """
    Function to merge a list of dictionaries. If keys are repeated, values are merged into a list.

    Parameters:
    - dict_list (list): List containing dictionaries.

    Returns:
    - Merged dictionary.
    """
    merged_dict = {}
    for d in dict_list:
        for key, value in d.items():
            if key in merged_dict:
                # Merge values into a list if the key already exists
                if not isinstance(merged_dict[key], list):
                    merged_dict[key] = [merged_dict[key]]
                merged_dict[key].append(value)
            else:
                # Directly add if the key does not exist
                merged_dict[key] = value
    return merged_dict


def generate_combinations(group):
    """
    Groups mutations by specific types and generates all possible combinations.

    Parameters:
    - group (DataFrameGroupBy): Grouped DataFrame.

    Returns:
    - List of all possible combinations.
    """
    # Group mutations by specific type
    spec_type_groups = group.groupby('Specific Type')
    # Create a dictionary to store lists of mutations mapped to proteins by specific type
    mutation_dict = defaultdict(list)
    for spec_type, g in spec_type_groups:
        for _, row in g.iterrows():
            if type(row["Amino acid site"]) == str:
                mutation_dict[spec_type].append({row['Protein']: row['Amino acid site'].strip()})
            else:
                mutation_dict[spec_type].append({row['Protein']: row['Amino acid site']})
    # Extract lists of dictionaries for each key
    values_lists = [mutation_dict[key] for key in mutation_dict]
    # Generate all possible combinations and merge dictionaries
    combinations = [merge_dicts_with_list(comb) for comb in product(*values_lists)]
    return combinations


def generate_protein_dict(grouped_data):
    """Generate protein dictionary"""
    new_protein_dict = defaultdict(list)
    for name, group in grouped_data:
        if 'combination' in name:
            new_protein_dict[name] = generate_combinations(group)
        else:
            new_protein_dict[name].extend(
                {row['Protein']: row['Amino acid site'].strip() if isinstance(row['Amino acid site'], str) else row[
                    'Amino acid site']}
                for _, row in group.iterrows()
            )
    return new_protein_dict


def load_total_markers(data):
    data["Specific Type"] = data["Protein Type"].str.rsplit("_", n = 1).str[-1]
    data['Protein Type'] = data['Protein Type'].str.replace(r'_\d+$', '', regex = True)
    return data.groupby('Protein Type')


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


def format_marker(marker, protein_prefix = ''):
    if '-' in marker:
        amino_acid = marker.split('-')[0]
        deletion_suffix = "Deletion"
    else:
        amino_acid = marker
        deletion_suffix = ""

    formatted_marker = f"{protein_prefix}-{amino_acid}{deletion_suffix}" \
        if protein_prefix else f"{amino_acid}{deletion_suffix}"
    return formatted_marker


def format_marker_list(markers, protein_prefix = ''):
    if isinstance(markers, str):
        return format_marker(markers, protein_prefix)

    all_contain_dash = all('-' in marker for marker in markers)
    if all_contain_dash:
        start = markers[0].split('-')[0]
        end = markers[-1].split('-')[0]
        return f"{protein_prefix}-{start}-{end}CompleteDeletion"

    return '&'.join(format_marker(marker, protein_prefix) for marker in markers)


def process_dictionary(data_dict):
    formatted_list = []
    for protein, markers in data_dict.items():
        if isinstance(markers, dict):
            for sub_protein, sub_markers in markers.items():
                formatted_marker = format_marker_list(sub_markers, f"{protein}-{sub_protein}")
                formatted_list.append(formatted_marker)
        else:
            formatted_marker = format_marker_list(markers, protein)
            formatted_list.append(formatted_marker)

    return '&'.join(formatted_list)


def process_protein_sequence(acc_id, renumbered_position, acc_pro_dic, marker_markers):
    # 检查输入序列的重编号位点残基有哪些是在marker_markers(标志物字典'H3':{'HA1':[],'HA2':[]})中存在的。
    protein_type = acc_pro_dic[acc_id]

    # Skip processing if protein type is unknown
    if protein_type == "Unknown":
        return None, None
    HA_TYPES_ALL = [f"H{i}" for i in range(1, 19)]
    NA_TYPES_ALL = [f"N{i}" for i in range(1, 10)]

    use_protein = "H3" if protein_type in HA_TYPES_ALL else protein_type

    expected_markers = marker_markers.get(use_protein, [])
    # print(f"expet:\n{expected_markers}")
    # print(f"renumber:\n{renumbered_position}")
    protein = f'H3' if protein_type in HA_TYPES_ALL else (
        f'N2' if protein_type in NA_TYPES_ALL else protein_type)
    markers = defaultdict(list)
    if use_protein == "H3":
        for hatype, ha_markers in expected_markers.items():
            for marker in ha_markers:
                match = re.match(r"(\d+)([A-Z]|-)", marker)
                index = 0 if hatype == "HA1" else 1
                if match and match.group() in renumbered_position[index][hatype]:
                    markers[hatype].append(match.group())
        return protein, markers
    # markers = defaultdict(list)
    # if use_protein == "H3":
    #     for hatype, ha_markers in expected_markers.items():
    #         for marker in ha_markers:
    #             match = re.match(r"(\d+)([A-Z])", marker)
    #             if match and match.group() in renumbered_position:
    #                 markers[hatype].append(match.group())
    #     return protein,markers
    markers_oth = []
    for marker in expected_markers:
        match = re.match(r"(\d+)([A-Z]|-)", marker)
        if match and match.group() in renumbered_position:
            markers_oth.append(match.group())
    return protein, markers_oth


def check_marker_combinations(total_markers, results_markers, markers_type, input_file_name, data, ha_type, na_type):
    results = []

    # Initialize results with empty/default values
    initial_data = {
        'Strain ID': '',  # or some default value
        'Amino acid site': '',  # or some default value
        'Protein Type': ''  # or some default value
    }
    results.append(initial_data)

    # Sequentially check if each type of marker for the phenotype is present in the identified marker dictionary.
    for marker_protein_type, marker_list in total_markers.items():
        # proba_comb is one of the multiple combinations of markers for each type, which is a dictionary.
        # 'combination-combination_449': [{'PB2': '158G', 'PA': '295P'}, {'PB2': '158A', 'PA': '295P'}]
        # {'H1-combination_409': [{'H3': {'HA1': ['149E', '208G']}}],
        # proba_comb is one of these dictionaries, if the dictionary is satisfied,
        # it is considered that there is a combination-combination_449 type of marker.
        for proba_comb in marker_list:
            print(f"子集比较：\n{proba_comb}\n{results_markers}")
            if proba_comb and is_subset_complex_revised(proba_comb, results_markers):
                # If the key-value pair in this dictionary exists in the identified marker dictionary,
                # return a more concise format.
                # if is_subset_complex_revised(proba_comb, results_markers):

                if proba_comb and all(proba_comb.values()):
                    markers_formated = process_dictionary(proba_comb)
                    print(f"proba_comb:\n{proba_comb}")
                    print(f"markers_formated:\n{markers_formated}")
                    results.append({
                        'Strain ID': input_file_name.split(".")[0],
                        'Amino acid site': markers_formated,
                        'Protein Type': marker_protein_type,
                    })

    results = pd.DataFrame(results)
    print(f"带合并：\n{results}")
    final_results = merge_dataframes(results, data, markers_type, ha_type, na_type)
    return final_results


def merge_dataframes(results, data, markers_type, ha_type, na_type):
    """
    Merge two DataFrames based on the 'Protein Type' column, handling rows with and without 'combination' differently.

    Args:
    results (pd.DataFrame): DataFrame containing results data.
    data (pd.DataFrame): DataFrame containing additional data to merge.
    markers_type (str): String indicating the type of markers to be used for renaming columns.

    Returns:
    pd.DataFrame: The merged DataFrame after processing.
    """
    # Pre-compile the regex pattern for performance
    combination_pattern = re.compile(r'combination')
    # Split the 'data' DataFrame into two based on 'combination' presence without using str.contains for each row
    data['HasCombination'] = data['Protein Type'].apply(lambda x: bool(combination_pattern.search(x)))
    data_with_combination = data[data['HasCombination']].drop_duplicates(subset = 'Protein Type')
    data_without_combination = data[~data['HasCombination']]

    # Drop the 'Amino acid site' column and the helper 'HasCombination' column to avoid merging them later
    data_with_combination = data_with_combination.drop(columns = ["Amino acid site", "HasCombination"])

    # Do the same split for the 'results' DataFrame
    results['HasCombination'] = results['Protein Type'].apply(lambda x: bool(combination_pattern.search(x)))
    results_with_combination = results[results['HasCombination']]
    results_without_combination = results[~results['HasCombination']]

    # Merge parts with and without 'combination' separately
    merged_with_combination = pd.merge(results_with_combination, data_with_combination, on = 'Protein Type',
                                       how = 'left')
    data_without_combination.loc[:, "Amino acid site"] = \
    data_without_combination.loc[:, "Amino acid site"].str.split("HA\d-").str[-1]
    merged_without_combination = pd.merge(results_without_combination, data_without_combination,
                                          on = ['Protein Type', 'Amino acid site'], how = 'left')

    # Concatenate the merged parts
    final_results = pd.concat([merged_with_combination, merged_without_combination], ignore_index = True)

    # Rename the 'Amino acid site' column and cleanup 'Protein Type'
    final_results.rename(columns = {'Amino acid site': f'{markers_type.title()} Markers'}, inplace = True)
    final_results['Protein Type'] = final_results['Protein Type'].str.replace('-combination.*', '', regex = True)

    final_results['Protein Type'] = final_results['Protein Type'].apply(lambda x: get_hana_string(x, ha_type, na_type))
    # Drop unnecessary columns and the helper 'HasCombination' column
    final_results.drop(columns = ["Specific Type", "Protein", "HasCombination_x", "HasCombination_y", "HasCombination"],
                       inplace = True)

    # Replace empty strings with NaN
    final_results.replace('', np.nan, inplace = True)

    # Drop rows where specific columns are NaN
    final_results.dropna(subset = ['Strain ID', f'{markers_type.title()} Markers', 'Protein Type'], how = "all",
                         inplace = True)

    final_results.drop_duplicates(inplace = True)
    return final_results


def get_hana_string(protein_type, ha_type, na_type):
    if protein_type in HA_TYPES or protein_type == "H3":
        return f'{ha_type}(H3 numbering)'
    elif protein_type in NA_TYPES or protein_type == "N2":
        return f'{na_type}(N2 numbering)'
    else:
        return protein_type


def identify_markers(input_file_path, renumbering_results, marker_markers, acc_pro_dic, markers_type, data,
                     output_directory = ".", prefix = ""):
    """
    Identifies virulence markers in protein sequences based on the provided marker markers
    and the renumbered sequences.

    Parameters:
    input_file_path (str): The path to the input file containing sequence data.
    renumbering_results (dict): Dictionary with accession ID as keys and renumbered positions as values.
    marker_markers (dict): Dictionary defining the virulence markers to be identified.
    acc_pro_dic (dict): Dictionary mapping accession IDs to protein types.
    markers_type (str): Type of markers to be identified (e.g., 'HA', 'NA').
    data (DataFrame): DataFrame containing additional data for marker identification.
    output_directory (str, optional): The directory where output files will be saved. Defaults to the current directory.
    prefix (str, optional): A prefix to be added to the output file names.

    Returns:
    DataFrame: A DataFrame with identified markers for each sequence.
    """
    # Ensure the output directory exists
    os.makedirs(output_directory, exist_ok = True)
    input_file_name = os.path.split(input_file_path)[1]
    results_markers = defaultdict(list)
    print(f"初始result——markers\n{renumbering_results}")
    # 有两个标志物列表，一个是用于识别具体单个的标志物的，一个是用于识别有哪些单个的或是组合的标志物的（是一个重新组织后的标志物列表形式）
    # 两个都经历了H3/NA转换
    # Process each accession ID and its renumbered position
    for acc_id, renumbered_position in renumbering_results.items():
        # renumbered_position是列表，除了HA外元素都是marker(226L)，HA是里面有两个字典，一个字典是HA1，对应一个列表（元素也是marker)
        # 另一个字典是HA2，对应一个列表（元素也是marker)
        # acc_id是accession号
        # 识别有哪些是存在于标志物列表的
        protein, markers = process_protein_sequence(acc_id, renumbered_position, acc_pro_dic, marker_markers)
        if protein:
            results_markers[protein] = markers

    # Initialize types
    ha_type = na_type = None

    # Identify HA and NA types based on acc_pro_dic
    print(acc_pro_dic)
    for acc, pro in acc_pro_dic.items():
        if pro in HA_TYPES or pro == "H3":
            ha_type = pro
        elif pro in NA_TYPES or pro == "N2":
            na_type = pro
    # S = {'H1-combination_409': [{'H1': ['163E', '222G']}],
    #  'M2': [{'M2': '41C'}, {'M2': '24D'}, {'M2': '82S'}],
    #  'combination-combination_460': [{'H3': ['299R', 'HA2-107I'], 'N2': '35R', 'M2': '41C'}]}

    # This is to handle each HA/NA, including those present in combinations,
    # the file only processed single HA/NA markers
    ori_markers = generate_protein_dict(load_total_markers(data))
    print(f"ori:\n{ori_markers}")
    total_markers = defaultdict(list)
    for pro, lst in ori_markers.items():
        for dic in lst:
            if dic and all(dic.values()):
                # After converting through convert_HA_residues, everything will become H3, so there's no impact
                total_markers[pro].append(convert_HA_residues(dic, STRUCTURE_PATH, hatype = None))
    print(f"total\n{total_markers}")
    print(f"result_markers\n{results_markers}")
    # Check marker combinations and merge results with data
    results_df = check_marker_combinations(total_markers, results_markers, markers_type,
                                           input_file_name, data, ha_type, na_type)

    # Add prefix to filename if provided
    add_prefix = prefix + "_" if prefix else ""
    filename = add_prefix + input_file_name.split(".")[0] + "_markers.csv"
    # Save the results to a CSV file
    results_df.to_csv(f"{output_directory}/{filename}", index = False, encoding = 'utf-8-sig')

    return results_df


def is_fasta_file(filename):
    # Check if input is a fasta file
    return re.search(r'\.(fasta|faa|fa)$', filename, re.IGNORECASE)


def find_files_with_string(directory, string):
    all_items = os.listdir(directory)
    files_with_string = [item for item in all_items
                         if string in item and os.path.isfile(os.path.join(directory, item))
                         and item.endswith("_annotated.csv")]

    return files_with_string[0]

def extract_protein_annotations(protein_path):
    """
    从FASTA文件提取蛋白类型和序列标识符，并将它们写入一个新的CSV文件。

    :param protein_path: FASTA文件的路径。
    :return: None
    """
    # 定义输出文件路径
    base_name = os.path.basename(protein_path).split('.')[0]
    result_dir_name = os.path.dirname(protein_path)
    output_path = os.path.join(result_dir_name, "{}_annotation.csv".format(base_name))

    # 解析FASTA文件
    sequences = SeqIO.parse(protein_path, "fasta")

    # 准备写入文件
    with open(output_path, 'w') as out_file:
        out_file.write("qseqid,Protein Abbreviation")
        for record in sequences:
            seq_id = record.id.split("|")[0]
            # 使用正则表达式从seq_id中提取蛋白类型
            protein_type_match = re.search("querySeq\d+_([^_]+)", seq_id)
            if protein_type_match:  # 确保找到了匹配
                protein_type = protein_type_match.group(1)
                # 将蛋白类型和记录ID作为一行写入文件
                out_file.write("{},{}\n".format(record.id, protein_type))

    # 输出文件操作完成
    print("Annotation written to {}".format(output_path))

def parse_args():
    parser = argparse.ArgumentParser(prog = 'flupre',
                                     description = 'flupre command line tool for flu markers '
                                                   'extraction, annotation, host and virulence level prediction.')
    subparsers = parser.add_subparsers(dest = 'subcommand', help = 'Sub-commands')

    # anno subcommand
    anno_parser = subparsers.add_parser('anno',
                                        help = 'Annotate a FASTA file or all FASTA files in a directory '
                                               'using DIAMOND BLAST against a flu database.')
    anno_parser.add_argument('-i', '--input', required = True,
                             help = 'Input FASTA file or directory containing FASTA files.')
    anno_parser.add_argument('-o', '--output_directory', type = str, default = RESULT_PATH,
                            help = 'Directory to save the output files. Defaults to the result directory.')
    anno_parser.add_argument('-temp', '--temp_directory', type = str, default = TEMP_PATH,
                             help = 'Directory to save the temp output files. Defaults to the temp directory.')
    # anno_parser.add_argument('-o', '--output_directory', type = str, default = '.',
    #                          help = 'Directory to save the output files. Defaults to the current directory.')
    # anno_parser.add_argument('-p', '--prefix', type = str, default = '', help = 'Prefix for the output filenames.')
    # anno_parser.add_argument('-e', '--evalue', type = float, default = 1e-5,
    #                          help = 'E-value threshold for DIAMOND BLAST hits. Defaults to 1e-5.')
    # anno_parser.add_argument('-u', '--update_file', action = 'store_true',
    #                          help = 'If set, updates the FASTA file with annotations.')
    # anno_parser.add_argument('-t', '--threads', type = int, default = 10,
    #                          help = 'Number of threads for DIAMOND BLAST. Defaults to 10.')

    # extract subcommand
    extract_parser = subparsers.add_parser('extract', help = 'Extract and process protein annotations.')
    extract_parser.add_argument('-i', '--input', required = True,
                                help = 'Input FASTA file or directory containing FASTA files.')
    extract_parser.add_argument('-a', '--anno_path', required = True,
                                help = 'Input annotation CSV file or directory containing annotation CSV files.')
    extract_parser.add_argument('-o', '--output_directory', type = str, default = '.',
                                help = 'Directory to save the output files. Defaults to the current directory.')
    extract_parser.add_argument('-p', '--prefix', type = str, default = '', help = 'Prefix for the output filenames.')
    extract_parser.add_argument('-type', '--markers_type', type = str, default = "virulence",
                                help = 'Type of markers. Defaults to virulence type.')
    # pred command
    pred_parser = subparsers.add_parser('predv', help = 'Predict new data labels using a trained model.')
    pred_parser.add_argument('-i', '--input', required = True, type = str,
                             help = 'Input CSV file with marker data or directory containing such files.')
    pred_parser.add_argument('-m', '--model_path', default = MODEL_PATH + '/random_forest_model.joblib', type = str,
                             help = 'Path to the trained model file.')
    pred_parser.add_argument('-th', '--threshold', default = 0.5, type = float,
                             help = 'Probability threshold for model prediction.')
    pred_parser.add_argument('-o', '--output_directory', type = str, default = '.',
                             help = 'Directory to save the prediction results. Defaults to the current directory.')
    pred_parser.add_argument('-p', '--prefix', type = str, default = '',
                             help = 'Prefix for the output filenames of the predictions.')

    # pred command
    pred_parser = subparsers.add_parser('predh', help = 'Predict new data labels using a trained model.')
    pred_parser.add_argument('-i', '--input', required = True, type = str,
                             help = 'Input CSV file with marker data or directory containing such files.')
    pred_parser.add_argument('-m', '--model_path', default = MODEL_PATH + '/gnb_model.joblib', type = str,
                             help = 'Path to the trained model file.')
    pred_parser.add_argument('-th', '--threshold_path', default = MODEL_PATH + '/optimal_threshold.joblib', type = str,
                             help = 'Path to the file containing the optimal threshold value.')
    pred_parser.add_argument('-f', '--top_features_path', default = MODEL_PATH + '/top_features.joblib', type = str,
                             help = 'Path to the file containing the top features.')
    pred_parser.add_argument('-o', '--output_directory', type = str, default = '.',
                             help = 'Directory to save the prediction results. Defaults to the current directory.')
    pred_parser.add_argument('-p', '--prefix', type = str, default = '',
                             help = 'Prefix for the output filenames of the predictions.')
    return parser.parse_args()


def process_anno_cmd(input_file, args):
    """
    Call the appropriate functions to process a single fasta file
    """
    # 在这里加ivew
    proteinPath, resultPath = ivew_task(args.output_directory,args.temp_directory, str(input_file))
    # 得到注释文件和包含抗原和blast预测的结果文件。
    extract_protein_annotations(proteinPath)
    # annotate_fasta_file(
    #     str(input_file),
    #     args.output_directory,
    #     args.prefix,
    #     args.evalue,
    #     args.update_file,
    #     args.threads
    # )


def process_extract_cmd(input_file, args, is_directory = True):
    """
        Call the appropriate functions to process a single fasta file
    """
    input_filename_pre = os.path.split(input_file)[1].split('.')[0]
    if is_directory:
        anno_filename = find_files_with_string(args.anno_path, input_filename_pre)
        annotations = pd.read_csv(f"{args.anno_path}/{anno_filename}")
    else:
        annotations = pd.read_csv(f"{args.anno_path}")
    # 在这里加上blast预测和抗原的必行代码
    acc_pro_dic = dict(zip(annotations.iloc[:, 0], annotations.iloc[:, 1]))
    for filename in os.listdir(MARKER_PATH):
        if filename.endswith("_formated.csv") and "lence" in filename:
            # marker_dict已经过H3/N2的位点转换，data未经过位点转换
            marker_dict, data = annotate_markers(MARKER_PATH + f"/{filename}", STRUCTURE_PATH)
            renumbering_results = renumber_proteins(
                fasta_path = str(input_file),
                acc_pro_dict = acc_pro_dic,
                marker_dict = marker_dict
            )
            markers_type = filename.split("_formated.csv")[0]
            markers_type = markers_type.split('_')[1] if "_" in markers_type else markers_type
            results_df = identify_markers(
                input_file_path = str(input_file),
                renumbering_results = renumbering_results,
                marker_markers = marker_dict,
                acc_pro_dic = acc_pro_dic,
                # Default directory
                output_directory = markers_type,
                prefix = args.prefix,
                markers_type = markers_type,
                data = data,
            )

            print("\nMarker extracted and saved to file.")


def process_directory(directory, args):
    for file in directory.iterdir():
        if is_fasta_file(str(file)):
            if args.subcommand == 'anno':
                process_anno_cmd(file, args)
            elif args.subcommand == 'extract':
                process_extract_cmd(file, args)


def process_single_file(file, args):
    if args.subcommand == 'anno':
        process_anno_cmd(file, args)
    elif args.subcommand == 'extract' and str(args.anno_path).endswith("_annotated.csv"):
        process_extract_cmd(file, args, is_directory = False)


def run_other_subcommand(args):
    input_path = Path(args.input)
    if input_path.is_dir():
        if (args.subcommand == "extract" and Path(args.anno_path).is_dir()) or (args.subcommand == "anno"):
            process_directory(input_path, args)
        else:
            print(f"Error: {args.anno_path} is not a valid directory")
    elif input_path.is_file():
        process_single_file(input_path, args)
    else:
        print(f"Error: {args.input} is not a valid file or directory", file = sys.stderr)


def main():
    args = parse_args()
    if args.subcommand == 'predv':
        predictions = predict_virulence.predict_new_data(
            str(Path(args.input)),
            args.model_path,
            args.threshold,
            DATA_PATH,
            args.output_directory,
            args.prefix,
        )
        print(predictions)
        print(f"Predictions completed.")
    elif args.subcommand == "predh":
        predictions = predict_host.predict_new_data(
            str(Path(args.input)),
            args.model_path,
            args.threshold_path,
            args.top_features_path,
            args.output_directory,
            args.prefix
        )
        print(predictions)
        print(f"Predictions completed.")

    else:
        run_other_subcommand(args)


if __name__ == '__main__':
    length_diffs = compare_sequences(STD_PATH, COMPLETE_STD_PATH)
    main()
