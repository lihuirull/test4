
from blastSeq import blastSeq
from getBlastMostCommonHitProteinType import getMostCommonHitProtein
from getSeq import get
from getAlignedSite import getAlignedSite
from margeSeqAlignedByMafft import margeSeqAlignedSeq
from translateDNA2Protein import translate
from standardize import standardize
from getDifferentSite import getdifferentSite
from analysisDBSite import analysisDBSite
from getDBSiteInfo import getDBSiteInfo
from changeStandardNameForMPNS import refreshStandardName
from predictVirusHost import getVirusHost
from collections import Counter
from sortResult import sortResult
from blastHASeq import blastHASeq
from analysisDBSiteHostMarker import analysisDBSiteHostMarker
import datetime,os,sys,getopt
#import codecs
#sys.stdout = codecs.getwriter("utf-8")(sys.stdout.detach())

def ivew_task(resultFileDir, inputFileDir, tempFileDir, inputFile):
    beginTime = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    DIR = "/home/think/platform/"
    #resultFileDir = "/home/think/platform/result/"
    #tempFileDir = "/home/think/platform/temp/"
    #inputFile = "test_query_seq12.fas"
    workName = inputFile
    #inputFileDir = DIR+"/querySeq/"
    #cd /home/think/platform/PycharmProjects/virusDB
    # python main.py --folderDir /home/think/platform/ --resultFileDir  /home/think/platform/result/ --resultFileName user2virus3  --tempFileDir  /home/think/platform/temp/ --inputFileDir /home/think/platform/querySeq/ --inputFile test_query_seq11.fas
    #opts, args = getopt.getopt(sys.argv[1:], "hi:o:",["inputFileDir=","inputFile=","tempFileDir=","resultFileDir=","resultFileName=","folderDir="])
    #for op, value in opts:
    #    if op == "--inputFileDir":
    #        inputFileDir = value
    #    elif op == "--inputFile":
    #        inputFile = value
    #    elif op == "--tempFileDir":
    #        tempFileDir = value
    #    elif op == "--resultFileDir":
    #        resultFileDir = value
    #    elif op == "--resultFileName":
    #        workName = value
    #    elif op == "--folderDir":
    #        DIR = value

    fileLog = open(resultFileDir+workName+".result","w",encoding='utf-8')
    fileLog.write("\n######################\t"+beginTime+"begin\n")
    ##########################################################################
    inputFile,querySeqType,dic = standardize(inputFile,inputFileDir,outFileDir=tempFileDir)
    # print("standard_name:original_name->"+str(dic))
    #print('liuil---' + str(dic))
    fileLog.write("standard_name:original_name->"+str(dic)+"\n")
    ##########################################################################

    ##########################################################################
    if querySeqType=="nucleo":
        dataBaseName = "/protType_NA"                        #blast
    if querySeqType=="protein":
        dataBaseName = "/protType_AA"                        #blast
    querySeqFile = inputFile
    querySeqFileDir = tempFileDir
    dirBlast = DIR+"app/blast+/bin/"
    DBDir = DIR+"18Mid/standard_seq/allProteinTypeDB/"     #blast
    queryEValue = "1e-5"                                            #blast
    outfmt = "3"                                                    #blast
    outName = "querySeqToProteinTypeDB"                                        #blast
    # tempFileDir = DIR+"18Mid//standard_seq/out//outout//"      # blas
    # tempFileDir = tempFileDir
    blastSeq(dirBlast,querySeqFile,querySeqType,DBDir,dataBaseName,queryEValue,outfmt,outName,querySeqFileDir,tempFileDir)  #
    predictedProteinType = getMostCommonHitProtein(outName,tempFileDir)      #
    # print("ProteinType->"+str(predictedProteinType))
    fileLog.write("ProteinType->"+str(predictedProteinType)+"\n")
    for eachSeqType in predictedProteinType:
        if 'Unknown' in eachSeqType:
            fileLog.write('Attention! One or more sequence inputted may not belong to influenza viruses\n')
    times= Counter(num[1] for num in predictedProteinType).most_common()
    for time in times:
        if time[1]>1:
            fileLog.write("Warning : You've entered "+str(time[1])+" "+str(time[0])+" proteins in\n")

    if querySeqType=="nucleo":
        dataBaseName = "/HostNuclDB"                        #blast
    if querySeqType=="protein":
        dataBaseName = "/HostProtDB"                        #blast
    querySeqFile = inputFile
    querySeqFileDir = tempFileDir
    DBDir = DIR+"18Mid/standard_seq/allProteinTypeDB/"
    # print(DBDir)
    queryEValue = "1e-5"  # blast
    outfmt = "3"  # blas
    outName = "querySeqToHostDB"  # blast
    tempFileDir = tempFileDir  # blas
    blastSeq(dirBlast,querySeqFile,querySeqType,DBDir,dataBaseName,queryEValue,outfmt,outName,querySeqFileDir,tempFileDir)  #
    Host = getMostCommonHitProtein(outName,tempFileDir)
    virusHost,HostNew=getVirusHost(predictedProteinType,Host)
    # print("VirusHost->"+str(virusHost)+str(HostNew)+"\n")
    # print("ProteinHost->"+str(Host))
    # fileLog.write("ProteinHost->"+str(Host)+"\n")
    # fileLog.write("ProteinHostUsed->"+str(HostNew)+"\n")
    fileLog.write("VirusHost->"+str(virusHost)+"\n")
    ##########################################################################

    ##########################################################################
    if querySeqType=="nucleo":
        queryProteinSeqFile = translate(querySeqFile,querySeqFileDir,DIR,tempFileDir)
        querySeqFile = queryProteinSeqFile[0]
        fileLog.write("CDS->" + str(queryProteinSeqFile[2]) + "\n")
    ###########################################################################
    ##########################################################################
    if querySeqType=="nucleo":
        predictedProteinType = refreshStandardName(predictedProteinType,dic)
        fileLog.write("standard_name_new:original_name->"+str(dic)+"\n")
        # print("ProteinType->"+str(predictedProteinType))
        fileLog.write("ProteinType->"+str(predictedProteinType)+"\n")
    ##########################################################################

    ##########################################################################
    for eachQuery in predictedProteinType:
        #print("\nProcessing->"+eachQuery[0],dic[">"+eachQuery[0]],eachQuery[1])
        fileLog.write("\n"+eachQuery[0]+"\t"+dic[">"+eachQuery[0]]+"\t"+eachQuery[1].rstrip()+"\n")
        # queryProteinSeqFile = translate(querySeqFile,querySeqFileDir)
        if str(eachQuery[1]).strip("\n")=="Unknown":
            # print("No result->"+eachQuery[0])
            fileLog.write("Virulent Site:\n"+"\tNo result"+"\n")
            continue
        querySeq,querySeqFileName = get(fileName=querySeqFile, dir=querySeqFileDir, seqName=eachQuery[0])
        mafftDir = DIR+"app/mafft/mafft-7.158-without-extensions/scripts/"
        queryFileName = eachQuery[0]
        standardDBDir = DIR+"18Mid/translatePerl/standard_seq_protein/"
        outDir = tempFileDir
        standardDB = eachQuery[1]+".fas"

        outFile = margeSeqAlignedSeq(querySeq, queryFileName, standardDBDir, standardDB, mafftDir,outDir)
        transposedFileName, transposedFileNameDir= getAlignedSite(file=outFile, dir=tempFileDir)
        DBInfoVirulence=getDBSiteInfo(file="DBSiteInfoVirulence",dir="/home/think/platform/18Mid/",proteinType=eachQuery[1])
        DBInfoVirulence = str(DBInfoVirulence)
        differentSiteFileName,differentSiteInfo , DBSeqName ,DBSite= getdifferentSite(transposedFileName,transposedFileNameDir,DBInfoVirulence)
        # fileLog.write(str(DBSeqName)+str(DBSite).replace("{","\n{").replace("}","}\n"))
        VirulentInfo=analysisDBSite(differentSiteFileName,outDir,str(DBInfoVirulence),'Virulence')
        fileLog.write("Virulent Site:\n"+VirulentInfo)
        ##########################################################################
        if eachQuery[1] in "N1 N2 N3 N4 N5 N6 N7 N8 N9 M2 PA":
            DBInfoResistance=getDBSiteInfo(file="DBSiteInfoResistance",dir="/home/think/platform/18Mid/",proteinType=eachQuery[1])
            DBInfoResistance = str(DBInfoResistance)
            differentSiteFileName,differentSiteInfo , DBSeqName ,DBSite= getdifferentSite(transposedFileName,transposedFileNameDir,DBInfoResistance)
            # fileLog.write(str(DBSeqName)+str(DBSite).replace("{","\n{").replace("}","}\n"))
            ResistanceInfo=analysisDBSite(differentSiteFileName,outDir,str(DBInfoResistance),'Resistance')
            fileLog.write("Resistance Site:\n"+ResistanceInfo)
        ##########################################################################

        if eachQuery[1] in ["H1","H3","H5","H7","H9"]:
            subType = eachQuery[1]
            querySeqFileDir = tempFileDir
            querySeqName = eachQuery[0]
            print(querySeqName)
            os.system("perl "+DIR+"18Mid/antigen/predAV.pl"+" "+subType+" "+querySeqFileDir+querySeqFileName+" "+DIR+"/18Mid/antigen/modelData/vaccine/"+subType+" "+tempFileDir+querySeqFileName+".antigenInfo"+" "+DIR+"/18Mid/antigen/modelData/"+" "+tempFileDir+" "+mafftDir)
            try:
                AntigenFile = open(tempFileDir+querySeqFileName+".antigenInfo")
                text = AntigenFile.readlines()
                text.sort(key=lambda x: ("Max" not in x.split('\t')[2] and ("Min" not in x.split('\t')[2]) and float('{:.10f}'.format(float(x.split('\t')[2])))) or (('Max' in x.split('\t')[2] and float(10000000000 + len(x))) or ('Min' in x.split('\t')[2] and float(-100000 + len(x)))))
                text = list(map(lambda x: ("Max" not in x.split('\t')[2] and "Min" not in x.split('\t')[2]) and x.replace('\t'+x.split('\t')[2]+'\t','\t'+str('{0:.2e}'.format(float(x.split('\t')[2])))+'\t') or x,text))
                text = list(map(lambda x: ("Max" not in x.split('\t')[2] and "Min" not in x.split('\t')[2]) and ((float(x.split('\t')[2]) <= 1.0 and x.replace(x.split('\t')[2], "Similar_" + x.split('\t')[2])) or (float(x.split('\t')[2]) > 1.0 and x.replace(x.split('\t')[2], "varient_" + x.split('\t')[2]))) or ("Min" in x.split('\t')[2] and x.replace(x.split('\t')[2],"similar_" + x.split('\t')[2])) or ("Max" in x.split('\t')[2] and x.replace(x.split('\t')[2],"varient_" + x.split('\t')[2])), text))
                fileLog.write("AntigenInfo\n")
                for each in text:
                    fileLog.write(each)
                AntigenFile.close()
            except Exception as e:
                print(e)
            #增加一个抗原相似性计算的模块
        if eachQuery[1] in ['H1','H2','H3',"H4",'H5',"H6",'H7',"H8",'H9',"H10","H11","H12","H13","H14","H15","H16"]:
            subType = eachQuery[1]
            querySeqFileDir = tempFileDir
            querySeqName = eachQuery[0]
            print(querySeqName)
            try:
                blastHASeq(blastQuerySeqHA = querySeqFileName,HAType = subType,tempDir = tempFileDir,mafftDir = mafftDir)
                AntigenHAFile = open(tempFileDir + querySeqFileName + ".antigenInfo.blast",'r')
                text = AntigenHAFile.readlines()
                text.sort(key=lambda x: ("Max" not in x.split('\t')[2] and ("Min" not in x.split('\t')[2]) and float('{:.10f}'.format(float(x.split('\t')[2])))) or (('Max' in x.split('\t')[2] and float(10000000000 + len(x))) or ('Min' in x.split('\t')[2] and float(-100000 + len(x)))))
                text = list(map(lambda x: ("Max" not in x.split('\t')[2] and "Min" not in x.split('\t')[2]) and x.replace('\t'+x.split('\t')[2]+'\t','\t'+str('{0:.2e}'.format(float(x.split('\t')[2])))+'\t') or x,text))
                text = list(map(lambda x: ("Max" not in x.split('\t')[2] and "Min" not in x.split('\t')[2]) and ((float(x.split('\t')[2]) <= 1.0 and x.replace(x.split('\t')[2], "Similar_" + x.split('\t')[2])) or (float(x.split('\t')[2]) > 1.0 and x.replace(x.split('\t')[2], "varient_" + x.split('\t')[2]))) or ("Min" in x.split('\t')[2] and x.replace(x.split('\t')[2],"similar_" + x.split('\t')[2])) or ("Max" in x.split('\t')[2] and x.replace(x.split('\t')[2],"varient_" + x.split('\t')[2])), text))
                # print(text)
                AntigenHAFile.close()
                fileLog.write("AntigenInfo_Blast\n")

                for each in text:
                    # print(each)
                    fileLog.write(each)
            except Exception as e:
                print(e)


    ##########################################################################
    downTime = (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    fileLog.write("\n######################\t"+str(downTime)+"down\n")
    # 得到宿主标记物
    for eachQuery in predictedProteinType:
        if str(eachQuery[1]).strip("\n")=="Unknown":
            continue
        querySeq,querySeqFileName = get(fileName=querySeqFile, dir=querySeqFileDir, seqName=eachQuery[0])
        mafftDir = DIR+"app/mafft/mafft-7.158-without-extensions/scripts/"
        queryFileName = eachQuery[0]
        standardDBDir = DIR+"18Mid/translatePerl/standard_seq_protein/"
        outDir = tempFileDir
        standardDB = eachQuery[1]+".fas"

        outFile = margeSeqAlignedSeq(querySeq, queryFileName, standardDBDir, standardDB, mafftDir,outDir)
        transposedFileName, transposedFileNameDir= getAlignedSite(file=outFile, dir=tempFileDir)
        DBInfoVirulence=getDBSiteInfo(file="DBSiteInfoHostMarker",dir="/home/think/platform/18Mid/",proteinType=eachQuery[1])
        DBInfoVirulence = str(DBInfoVirulence)
        differentSiteFileName,differentSiteInfo , DBSeqName ,DBSite= getdifferentSite(transposedFileName,transposedFileNameDir,DBInfoVirulence)
        # fileLog.write(str(DBSeqName)+str(DBSite).replace("{","\n{").replace("}","}\n"))
        VirulentInfo=analysisDBSiteHostMarker(differentSiteFileName,outDir,str(DBInfoVirulence),'hostMarker')
        if VirulentInfo == '':continue
        fileLog.write("\n" + eachQuery[0] + "\t" + dic[">" + eachQuery[0]] + "\t" + eachQuery[1].rstrip() + "\n")
        fileLog.write(VirulentInfo)
    # fileLog.write(']\n\n')

    fileLog.close()
    sortResult(resultFileDir,workName+".result")
    # finishedDir = os.path.join(resultFileDir, 'finished')
    #os.makedirs(finishedDir)

    ##########################################################################
    #os.system("cd "+tempFileDir+"\nrm -r *Seq* *.stdName")



if __name__ == "__main__":
    # ivew_task('/home/think/platform/result/', '/home/think/platform/querySeq/', '/home/think/platform/temp/', 'test_query_seq1.fas')
    ivew_task('/home/work/IVEW/python_ivew/tempfiles/a9d6dba3-2048-c2b6-03df-a06f3349dd62/1556077782297/result/', '/home/work/IVEW/python_ivew/tempfiles/a9d6dba3-2048-c2b6-03df-a06f3349dd62/1556077782297/querySeq/', '/home/work/IVEW/python_ivew/tempfiles/a9d6dba3-2048-c2b6-03df-a06f3349dd62/1556077782297/temp/','Isolate1')
