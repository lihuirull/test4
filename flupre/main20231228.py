# !/usr/bin/python 
# -*- coding: utf-8 -*-
# Author:lihuiru
# Created on 2023/12/28 10:24
import datetime
import os
from collections import Counter

from blastHASeq import blastHASeq
from blastSeq import blastSeq
from changeStandardNameForMPNS import refreshStandardName
from getBlastMostCommonHitProteinType import getMostCommonHitProtein, getMostCommonHitProteinLowLevelHost
from getSeq import get
from predictVirusHost import getVirusHost, getVirusHostLowLevel
from standardize import standardize
from translateDNA2Protein import translate, makeProteinFileForDownload

base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# ==D:/Learning/edegfile/think/platform/"

DB_PATH = os.path.join(base_dir, 'data', 'flu_db.dmnd')


def process_antigen_text(file_path):
    """
    处理抗原信息文件的文本行。
    :param file_path: 抗原信息文件路径
    :return: 处理后的文本行列表
    """
    with open(file_path, 'r') as AntigenFile:
        text = AntigenFile.readlines()

    # 对text进行排序
    text.sort(key = lambda x: calculate_sort_key(x))

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


def ivew_task(resultFileDir, inputFileDir, tempFileDir, inputFile):
    beginTime = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    DIR = base_dir
    print(DIR)
    workName = inputFile

    fileLog = open(resultFileDir + workName + ".result", "w", encoding = 'utf-8')
    fileLog.write("\n######################\t" + beginTime + "begin\n")
    ##########################################################################
    # 标准化序列并返回类型（DNA/蛋白质）
    inputFile, querySeqType, dic = standardize(inputFile, inputFileDir, outFileDir = tempFileDir)

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
    proteinPath = ""
    if querySeqType == "nucleo":
        # 返回：outputName.replace(".fas", ".fas2"), outFileDir, dicCDS
        queryProteinSeqFile = translate(querySeqFile, querySeqFileDir, DIR, tempFileDir)
        querySeqFile = queryProteinSeqFile[0]

        # 新的蛋白文件路径
        proteinPath = os.path.join(queryProteinSeqFile[1], querySeqFile)

        fileLog.write("CDS->" + str(queryProteinSeqFile[2]) + "\n")
        CDSList = queryProteinSeqFile[2]

    ###########################################################################
    makeProteinFileForDownload(tempFileDir, querySeqFile, resultFileDir, dic, predictedProteinType)

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


if __name__ == "__main__":
    # processed_text = process_antigen_text(r"D:\Learning\edegfile\think\platform\temp\querySeq1_HA.antigenInfo")
    # print(processed_text)
    dirUser = base_dir
    s,j = ivew_task(dirUser + '/result/', dirUser + '/querySeq/', dirUser + '/temp/', 'test_query_seq20.fas')
    print(s)
    print(j)
    print('done')
# AntigenHAFile = open(tempFileDir + querySeqFileName + ".antigenInfo.blast", 'r')
# text = AntigenHAFile.readlines()
# """
# 作用：对text列表根据第三列进行排序。
# 解释: 使用复杂的lambda函数作为排序关键字。它根据每行文本中第三列（由x.split('\t')[2]获取）的内容来决定排序顺序。
# 如果该列不包含"Max"且不包含"Min"，则按照它的浮点数值（格式化为10位小数）进行正常排序。
# 如果包含"Max"或"Min"，则进行特殊处理，以确保它们位于排序的特定位置。
# """
# text.sort(key = lambda x: ("Max" not in x.split('\t')[2] and ("Min" not in x.split('\t')[2]) and float(
#     '{:.10f}'.format(float(x.split('\t')[2])))) or (
#                                   ('Max' in x.split('\t')[2] and float(10000000000 + len(x))) or (
#                                   'Min' in x.split('\t')[2] and float(-100000 + len(x)))))
# """
# 作用: 格式化text列表中每行的第三列数据。
# 解释: 遍历text中的每一行，如果第三列不包含"Max"和"Min"，则将该列的数值转换为科学计数法表示。否则，保持原样。
# """
# text = list(map(
#     lambda x: ("Max" not in x.split('\t')[2] and "Min" not in x.split('\t')[2]) and x.replace(
#         '\t' + x.split('\t')[2] + '\t',
#         '\t' + str('{0:.2e}'.format(float(x.split('\t')[2]))) + '\t') or x, text))
#
# """
# 作用: 进一步格式化text列表中的数据，标注序列的相似度或变异。
# 解释: 对每行的第三列数据进行条件判断和替换：
# 如果不包含"Max"和"Min"：
# 如果值小于或等于1.0，添加前缀"Similar_"。
# 如果值大于1.0，添加前缀"variant_"。
# 如果包含"Min"，添加前缀"similar_"。
# 如果包含"Max"，添加前缀"variant_"。
# """
# text = list(map(lambda x: ("Max" not in x.split('\t')[2] and "Min" not in x.split('\t')[2]) and ((float(
#     x.split('\t')[2]) <= 1.0 and x.replace(x.split('\t')[2], "Similar_" + x.split('\t')[2])) or (float(
#     x.split('\t')[2]) > 1.0 and x.replace(x.split('\t')[2], "variant_" + x.split('\t')[2]))) or (
#                                   "Min" in x.split('\t')[2] and x.replace(x.split('\t')[2],
#                                                                           "similar_" +
#                                                                           x.split('\t')[2])) or (
#                                   "Max" in x.split('\t')[2] and x.replace(x.split('\t')[2],
#                                                                           "variant_" +
#                                                                           x.split('\t')[2])), text))
# # print(text)
