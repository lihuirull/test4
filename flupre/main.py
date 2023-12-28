import datetime
import os
from collections import Counter

from addClusterTransitionSubstitution import addClusterTransitionSubstitution
from analysisDBSite import analysisDBSite
from analysisDBSiteHostMarker import analysisDBSiteHostMarker
from blastHASeq import blastHASeq
from blastSeq import blastSeq
from changeStandardNameForMPNS import refreshStandardName
from getAlignedSite import getAlignedSite
from getBlastMostCommonHitProteinType import getMostCommonHitProtein, getMostCommonHitProteinLowLevelHost
from getDBSiteInfo import getDBSiteInfo
from getDifferentSite import getdifferentSite
from getSeq import get
from margeSeqAlignedByMafft import margeSeqAlignedSeq
from predictVirusHost import getVirusHost, getVirusHostLowLevel
from sortResult import sortResult
from standardize import standardize
from translateDNA2Protein import translate, makeProteinFileForDownload


# import codecs
# sys.stdout = codecs.getwriter("utf-8")(sys.stdout.detach())

def ivew_task(resultFileDir, inputFileDir, tempFileDir, inputFile):
    beginTime = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    # DIR = "/home/think/platform/"
    DIR = r"D:/Learning/edegfile/think/platform/"
    workName = inputFile

    fileLog = open(resultFileDir + workName + ".result", "w", encoding = 'utf-8')
    fileLog.write("\n######################\t" + beginTime + "begin\n")
    ##########################################################################
    # 标准化序列并返回类型（DNA/蛋白质）
    inputFile, querySeqType, dic = standardize(inputFile, inputFileDir, outFileDir = tempFileDir)

    fileLog.write("standard_name:original_name->" + str(dic) + "\n")
    ##########################################################################

    ##########################################################################
    if querySeqType == "nucleo":
        dataBaseName = "/protType_NA"  # blast
    if querySeqType == "protein":
        dataBaseName = "/protType_AA"  # blast
    querySeqFile = inputFile
    querySeqFileDir = tempFileDir
    # dirBlast = DIR+"app/blast+/bin/"
    dirBlast = DIR + "app/ncbi-blast-2.9.0+/bin/"
    DBDir = DIR + "18Mid/standard_seq/allProteinTypeDB/"  # blast
    queryEValue = "1e-5"  # blast
    outfmt = "3"  # blast
    outName = "querySeqToProteinTypeDB"  # blast
    blastSeq(dirBlast, querySeqFile, querySeqType, DBDir, dataBaseName, queryEValue, outfmt, outName, querySeqFileDir,
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
    DBDir = DIR + "18Mid/standard_seq/allProteinTypeDB/"
    # print(DBDir)
    queryEValue = "1e-5"  # blast
    outfmt = "3"  # blas
    outName = "querySeqToHostDB"  # blast
    tempFileDir = tempFileDir  # blas
    blastSeq(dirBlast, querySeqFile, querySeqType, DBDir, dataBaseName, queryEValue, outfmt, outName, querySeqFileDir,
             tempFileDir)  #
    Host = getMostCommonHitProtein(outName, tempFileDir)
    HostLowLevel = getMostCommonHitProteinLowLevelHost(outName, tempFileDir, Host)
    virusHostLowLevel, HostNew2 = getVirusHostLowLevel(predictedProteinType, HostLowLevel)
    virusHost, HostNew = getVirusHost(predictedProteinType, Host)
    fileLog.write("virusHostLowLevel->" + str(virusHostLowLevel) + "\n")
    fileLog.write("VirusHost->" + str(virusHost) + "\n")
    ##########################################################################

    ##########################################################################
    CDSList = []
    if querySeqType == "nucleo":
        # 返回：outputName.replace(".fas", ".fas2"), outFileDir, dicCDS
        queryProteinSeqFile = translate(querySeqFile, querySeqFileDir, DIR, tempFileDir)
        querySeqFile = queryProteinSeqFile[0]
        fileLog.write("CDS->" + str(queryProteinSeqFile[2]) + "\n")
        CDSList = queryProteinSeqFile[2]
    # 获取对每条DNA序列翻译后注释的结果（包含蛋白类型和长度）
    # CDSList:[('querySeq1', 'HA(1..1698)'), ('querySeq2', 'NA(1..1419)'), ('querySeq3', 'NP(1..1494)'),
    # ('querySeq4', 'NS1(1..711),NS2(1..30,503..838)'), ('querySeq5', 'PA(1..2148),PA-X(1..570,572..760)'),
    # ('querySeq6', 'PB1(1..2274),PB1-F2(95..364)'), ('querySeq7', 'PB2(1..2277)'),
    # ('querySeq8', 'M1(1..816),M2(1..26,715..982)'), ('querySeq9', 'HA(1..1698)')]
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
        # if querySeqType=="protein":
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
        mafftDir = DIR + "app/mafft/mafft-7.158-without-extensions/scripts/"
        # queryFileName = eachQuery[0]
        # standardDBDir = DIR + "18Mid/translatePerl/standard_seq_protein/"
        # outDir = tempFileDir
        # standardDB = eachQuery[1] + ".fas"
        #
        # outFile = margeSeqAlignedSeq(querySeq, queryFileName, standardDBDir, standardDB, mafftDir, outDir)
        # transposedFileName, transposedFileNameDir = getAlignedSite(file = outFile, dir = tempFileDir)
        # DBInfoVirulence = getDBSiteInfo(file = "DBSiteInfoVirulence", dir = "/home/think/platform/18Mid/",
        #                                 proteinType = eachQuery[1])
        # DBInfoVirulence = str(DBInfoVirulence)
        # differentSiteFileName, differentSiteInfo, DBSeqName, DBSite = getdifferentSite(transposedFileName,
        #                                                                                transposedFileNameDir,
        #                                                                                DBInfoVirulence)
        # # fileLog.write(str(DBSeqName)+str(DBSite).replace("{","\n{").replace("}","}\n"))
        # VirulentInfo = analysisDBSite(differentSiteFileName, outDir, str(DBInfoVirulence), 'Virulence')
        # fileLog.write("Virulent Site:\n" + VirulentInfo)
        # ##########################################################################
        # if eachQuery[1] in "N1 N2 N3 N4 N5 N6 N7 N8 N9 M2 PA":
        #     DBInfoResistance = getDBSiteInfo(file = "DBSiteInfoResistance", dir = "/home/think/platform/18Mid/",
        #                                      proteinType = eachQuery[1])
        #     DBInfoResistance = str(DBInfoResistance)
        #     differentSiteFileName, differentSiteInfo, DBSeqName, DBSite = getdifferentSite(transposedFileName,
        #                                                                                    transposedFileNameDir,
        #                                                                                    DBInfoResistance)
        #     ResistanceInfo = analysisDBSite(differentSiteFileName, outDir, str(DBInfoResistance), 'Resistance')
        #     fileLog.write("Resistance Site:\n" + ResistanceInfo)
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
                "perl " + DIR + "18Mid/antigen/predAV.pl" + " " + subType + " " + querySeqFileDir + querySeqFileName + " " + DIR + "/18Mid/antigen/modelData/vaccine/" + subType + " " + tempFileDir + querySeqFileName + ".antigenInfo" + " " + DIR + "/18Mid/antigen/modelData/" + " " + tempFileDir + " " + mafftDir)
            try:
                AntigenFile = open(tempFileDir + querySeqFileName + ".antigenInfo")
                text = AntigenFile.readlines()
                text.sort(key = lambda x: ("Max" not in x.split('\t')[2] and ("Min" not in x.split('\t')[2]) and float(
                    '{:.10f}'.format(float(x.split('\t')[2])))) or (
                                                  ('Max' in x.split('\t')[2] and float(10000000000 + len(x))) or (
                                                  'Min' in x.split('\t')[2] and float(-100000 + len(x)))))
                text = list(map(
                    lambda x: ("Max" not in x.split('\t')[2] and "Min" not in x.split('\t')[2]) and x.replace(
                        '\t' + x.split('\t')[2] + '\t',
                        '\t' + str('{0:.2e}'.format(float(x.split('\t')[2]))) + '\t') or x, text))
                text = list(map(lambda x: ("Max" not in x.split('\t')[2] and "Min" not in x.split('\t')[2]) and ((float(
                    x.split('\t')[2]) <= 1.0 and x.replace(x.split('\t')[2], "Similar_" + x.split('\t')[2])) or (float(
                    x.split('\t')[2]) > 1.0 and x.replace(x.split('\t')[2], "variant_" + x.split('\t')[2]))) or (
                                                  "Min" in x.split('\t')[2] and x.replace(x.split('\t')[2],
                                                                                          "similar_" +
                                                                                          x.split('\t')[2])) or (
                                                  "Max" in x.split('\t')[2] and x.replace(x.split('\t')[2],
                                                                                          "variant_" +
                                                                                          x.split('\t')[2])), text))
                fileLog.write("AntigenInfo\n")
                for each in text:
                    fileLog.write(each.replace('\tX', '\t-').replace('X\t', '-\t'))
                AntigenFile.close()
            except Exception as e:
                print(e)
            # Add a module of antigen similarity calculation
        if eachQuery[1] in ['H1', 'H2', 'H3', "H4", 'H5', "H6", 'H7', "H8", 'H9', "H10", "H11", "H12", "H13", "H14",
                            "H15", "H16"]:
            # 针对所有亚型，传入自己HA的序列在一个文件中），通过blastp比对确定五条最相似的序列并传入这五条序列（在一个文件中），
            # 计算两两遗传距离
            subType = eachQuery[1]
            querySeqFileDir = tempFileDir
            querySeqName = eachQuery[0]
            print(querySeqName)
            try:
                blastHASeq(blastQuerySeqHA = querySeqFileName, HAType = subType, tempDir = tempFileDir,
                           mafftDir = mafftDir)
                AntigenHAFile = open(tempFileDir + querySeqFileName + ".antigenInfo.blast", 'r')
                text = AntigenHAFile.readlines()
                """
                作用：对text列表根据第三列进行排序。
                解释: 使用复杂的lambda函数作为排序关键字。它根据每行文本中第三列（由x.split('\t')[2]获取）的内容来决定排序顺序。
                如果该列不包含"Max"且不包含"Min"，则按照它的浮点数值（格式化为10位小数）进行正常排序。
                如果包含"Max"或"Min"，则进行特殊处理，以确保它们位于排序的特定位置。
                """
                text.sort(key = lambda x: ("Max" not in x.split('\t')[2] and ("Min" not in x.split('\t')[2]) and float(
                    '{:.10f}'.format(float(x.split('\t')[2])))) or (
                                                  ('Max' in x.split('\t')[2] and float(10000000000 + len(x))) or (
                                                  'Min' in x.split('\t')[2] and float(-100000 + len(x)))))
                """
                作用: 格式化text列表中每行的第三列数据。
                解释: 遍历text中的每一行，如果第三列不包含"Max"和"Min"，则将该列的数值转换为科学计数法表示。否则，保持原样。
                """
                text = list(map(
                    lambda x: ("Max" not in x.split('\t')[2] and "Min" not in x.split('\t')[2]) and x.replace(
                        '\t' + x.split('\t')[2] + '\t',
                        '\t' + str('{0:.2e}'.format(float(x.split('\t')[2]))) + '\t') or x, text))

                """
                作用: 进一步格式化text列表中的数据，标注序列的相似度或变异。
                解释: 对每行的第三列数据进行条件判断和替换：
                如果不包含"Max"和"Min"：
                如果值小于或等于1.0，添加前缀"Similar_"。
                如果值大于1.0，添加前缀"variant_"。
                如果包含"Min"，添加前缀"similar_"。
                如果包含"Max"，添加前缀"variant_"。
                """
                text = list(map(lambda x: ("Max" not in x.split('\t')[2] and "Min" not in x.split('\t')[2]) and ((float(
                    x.split('\t')[2]) <= 1.0 and x.replace(x.split('\t')[2], "Similar_" + x.split('\t')[2])) or (float(
                    x.split('\t')[2]) > 1.0 and x.replace(x.split('\t')[2], "variant_" + x.split('\t')[2]))) or (
                                                  "Min" in x.split('\t')[2] and x.replace(x.split('\t')[2],
                                                                                          "similar_" +
                                                                                          x.split('\t')[2])) or (
                                                  "Max" in x.split('\t')[2] and x.replace(x.split('\t')[2],
                                                                                          "variant_" +
                                                                                          x.split('\t')[2])), text))
                # print(text)
                AntigenHAFile.close()
                fileLog.write("AntigenInfo_Blast\n")

                for each in text:
                    # print(each)
                    fileLog.write(each.replace('\tX', '\t-').replace('X\t', '-\t'))
            except Exception as e:
                print(e)

    ##########################################################################
    # downTime = (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    # fileLog.write("\n######################\t" + str(downTime) + "down\n")
    # # obtain the host marker
    # for eachQuery in predictedProteinType:
    #     if (querySeqType == "nucleo") and (eachQuery[0] not in CDSTypeList): continue
    #     if eachQuery[1] not in ['H1', 'M1', 'M2', 'NP', 'NS1', 'NS2', 'PA', 'PB1', 'PB1-F2', 'PB2']: continue
    #     print(eachQuery)
    #     if str(eachQuery[1]).strip("\n") == "Unknown":
    #         continue
    #     querySeq, querySeqFileName = get(fileName = querySeqFile, dir = querySeqFileDir, seqName = eachQuery[0])
    #     mafftDir = DIR + "app/mafft/mafft-7.158-without-extensions/scripts/"
    #     queryFileName = eachQuery[0]
    #     standardDBDir = DIR + "18Mid/translatePerl/standard_seq_protein/"
    #     outDir = tempFileDir
    #     standardDB = eachQuery[1] + ".fas"
    #
    #     outFile = margeSeqAlignedSeq(querySeq, queryFileName, standardDBDir, standardDB, mafftDir, outDir)
    #     transposedFileName, transposedFileNameDir = getAlignedSite(file = outFile, dir = tempFileDir)
    #     DBInfoVirulence = getDBSiteInfo(file = "DBSiteInfoHostMarker", dir = "/home/think/platform/18Mid/",
    #                                     proteinType = eachQuery[1])
    #     DBInfoVirulence = str(DBInfoVirulence)
    #     differentSiteFileName, differentSiteInfo, DBSeqName, DBSite = getdifferentSite(transposedFileName,
    #                                                                                    transposedFileNameDir,
    #                                                                                    DBInfoVirulence)
    #     VirulentInfo = analysisDBSiteHostMarker(differentSiteFileName, outDir, str(DBInfoVirulence), 'hostMarker')
    #     if VirulentInfo == '': continue
    #     fileLog.write("\n" + eachQuery[0] + "\t" + dic[">" + eachQuery[0]] + "\t" + eachQuery[1].rstrip() + "\n")
    #     fileLog.write(VirulentInfo)
    #
    # fileLog.close()
    # sortResult(resultFileDir, workName + ".result")
    # addClusterTransitionSubstitution(resultFileDir = resultFileDir, resultFileName = workName + ".result",
    #                                  DBFileName = "DBSiteInfoClusterTransition", dir = "/home/think/platform/18Mid/",
    #                                  proteinTypeList = predictedProteinType)
    # # finishedDir = os.path.join(resultFileDir, 'finished')
    # # os.makedirs(finishedDir)
    #
    # ##########################################################################
    # # os.system("cd "+tempFileDir+"\nrm -r *Seq* *.stdName")


if __name__ == "__main__":
    # ivew_task('/home/think/platform/result/', '/home/think/platform/querySeq/', '/home/think/platform/temp/', 'test_query_seq1.fas')
    # dirUser = '/home/work/IVEW/python_ivew/tempfiles/0338aaa1-0059-91af-0a88-6caefb1e091f/1602692159826'
    # ivew_task(dirUser+'/result/', dirUser+'/querySeq/', dirUser+'/temp/','Isolate1')
    # print('down2')

    # dirUser = '/home/cong/windows/FluPhenotype/temp/DNA'
    dirUser = '..'
    ivew_task(dirUser + '/result/', dirUser + '/querySeq/', dirUser + '/temp/', 'test_query_seq20.fas')
    # ivew_task(dirUser+'/result/', dirUser+'/querySeq/', dirUser+'/temp/','Isolate1')
    print('down2')
