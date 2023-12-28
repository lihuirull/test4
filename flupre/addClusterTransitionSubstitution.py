# -*- coding:UTF-8 -*-

def addClusterTransitionSubstitution(resultFileDir,resultFileName,DBFileName,dir,proteinTypeList):
    dic = {}
    dicFile = open(dir+DBFileName,'r',encoding='utf-8')
    for eachLine in dicFile:
        dic[eachLine.split('>')[0]] = eval(eachLine.split('>')[1].strip('\n'))
    dicFile.close()
    
    dicDBAllProteinType = {}
    fileInDic = open(dir+'DBSiteInfoClusterTransition_deleteionSite','r')
    for eachLine in fileInDic:
        dicDBAllProteinType[eachLine.split('>')[0]] = eachLine.split('>',1)[1].strip('\n')
    fileInDic.close()

    
    fileIn = open(resultFileDir+resultFileName)
    textIn = fileIn.readlines()
    fileIn.close()
    flag = 0
    blastHAFlag = 0
    fileOut = open(resultFileDir+resultFileName,'w')
    dicSite = {}
    proteinType = ''
    clusterTransitionInfo_deletionSite=''
    for eachLine in textIn:
        if eachLine.startswith('querySeq') and eachLine.split('\t')[2].strip('\n') in dic:
            proteinType = eachLine.split('\t')[2].strip('\n')
            dicSite = dic[proteinType]
        if eachLine.startswith('AntigenInfo') and proteinType!='':
            flag = 1
        if flag == 1 and eachLine.strip('\n')=='':
            flag = 0
        if flag == 1:
            if '\t>' in eachLine or eachLine.startswith('Virulent Site:') or eachLine.startswith('AntigenInfo'):
                fileOut.write(eachLine)
                if 'AntigenInfo_Blast' in eachLine:
                    blastHAFlag = 1
                continue
            if proteinType in dicDBAllProteinType:
                vaccineName = eachLine.split('\t')[1].split('____')[0]
                variantSiteList = eachLine.strip('\n').split('\t')[3:]
                clusterTransitionInfo_deletionSite = getDeletionSiteClusterTransitionSubstitution(proteinType,proteinTypeList,eval(dicDBAllProteinType[proteinType]),resultFileDir,vaccineName,blastHAFlag,variantSiteList)
            
            clusterTransitionInfo = analysesSite(eachLine,dicSite,proteinType)
            if clusterTransitionInfo_deletionSite != '':
                clusterTransitionInfo = clusterTransitionInfo+clusterTransitionInfo_deletionSite
            fileOut.write(eachLine.strip('\n')+clusterTransitionInfo+'\n')
        else:
            fileOut.write(eachLine)
            
    fileOut.close()

def analysesSite(eachLine,dicSite,proteinType):
    import re
    siteInfoDic = {}
    for eachSite in eachLine.strip('\n').split('\t')[3:]:
        site = ''
        AAChange = []
        for eachLetter in eachSite:
            if eachLetter in ['0','1','2','3','4','5','6','7','8','9']:
                site+=eachLetter
            else:
                AAChange.append(eachLetter)
        
        siteInfo = site+'__'+ AAChange[0]+'__'+AAChange[1]
        siteInfoDic[site] = siteInfo
        
    
    outInfo = ''
    for DBSite in dicSite:
        if '&' in DBSite:
            querySiteNum = ''
            AAChange = ['','']
            DBSiteSplited = DBSite.split('&')
            for eachDBSite in DBSiteSplited:
                if eachDBSite not in siteInfoDic:
                    break
                else:
                    AAChange[0] = AAChange[0]+siteInfoDic[eachDBSite].split('__')[1]
                    AAChange[1] = AAChange[1]+siteInfoDic[eachDBSite].split('__')[2]
            DBseq = dicSite[DBSite].split('~')[0]
            reference = '~'.join(dicSite[DBSite].split("~")[1:-1])
            describe = dicSite[DBSite].strip('\n').split("~")[-1]
            if (re.match(DBseq.split("|")[0],AAChange[0])  and re.match(DBseq.split("|")[1],AAChange[1])) or (re.match(DBseq.split("|")[0],AAChange[1])  and re.match(DBseq.split("|")[1],AAChange[0])):
                outInfo = outInfo+'|'+AAChange[0]+DBSite+AAChange[1]+'+'+reference+'+'+describe

        else:
            if DBSite in siteInfoDic:
                AAChange = ['','']
                AAChange[0] = siteInfoDic[DBSite].split('__')[1]
                AAChange[1] = siteInfoDic[DBSite].split('__')[2]
                DBseq = dicSite[DBSite].split('~')[0]
                reference = '~'.join(dicSite[DBSite].split("~")[1:-1])
                describe = dicSite[DBSite].strip('\n').split("~")[-1]
                if (re.match(DBseq.split("|")[0],AAChange[0])  and re.match(DBseq.split("|")[1],AAChange[1])) or (re.match(DBseq.split("|")[0],AAChange[1])  and re.match(DBseq.split("|")[1],AAChange[0])):
                    outInfo = outInfo+'|'+AAChange[0]+DBSite+AAChange[1]+'+'+reference+'+'+describe
            
            
            DicSiteH3= {'156': 'K|E~PMID:24264991~The single cluster-transition substitution for the HK68 to EN72 cluster transition was 155TY and 156KE alone was responsible for the antigenic difference between the TX77 and BK79 clusters.'} #a same site with DB
            if (proteinType == 'H3') and (DBSite in DicSiteH3) and (DBSite in siteInfoDic):
                AAChange = ['','']
                AAChange[0] = siteInfoDic[DBSite].split('__')[1]
                AAChange[1] = siteInfoDic[DBSite].split('__')[2]
                DBseq = DicSiteH3[DBSite].split('~')[0]
                reference = '~'.join(DicSiteH3[DBSite].split("~")[1:-1])
                describe = DicSiteH3[DBSite].strip('\n').split("~")[-1]
                if (re.match(DBseq.split("|")[0],AAChange[0])  and re.match(DBseq.split("|")[1],AAChange[1])) or (re.match(DBseq.split("|")[0],AAChange[1])  and re.match(DBseq.split("|")[1],AAChange[0])):
                    outInfo = outInfo+'|'+AAChange[0]+DBSite+AAChange[1]+'+'+reference+'+'+describe
    return outInfo

def getDeletionSiteClusterTransitionSubstitution(proteinType,proteinTypeList,dicDB,resultFileDir,vaccineName,blastHAFlag,variantSiteList):
    siteInfoDic = {}
    for eachSite in variantSiteList:
        site = ''
        AAChange = []
        for eachLetter in eachSite:
            if eachLetter in ['0','1','2','3','4','5','6','7','8','9']:
                site+=eachLetter
            else:
                AAChange.append(eachLetter)

        if site == '130' and (AAChange[1] not in ['-','X','x','K','k']):
            siteInfo =( AAChange[0] , '#' )
        else:
            siteInfo =( AAChange[0] , AAChange[1] )
        siteInfoDic[site] = siteInfo
    
    outInfo = ''
    for each in dicDB:
        flag = 0
        for eachLetter in each.split('&'):
            if (eachLetter in siteInfoDic) and (siteInfoDic[eachLetter][1] in ['#']):
                continue
            if eachLetter not in siteInfoDic:
                continue
            if siteInfoDic[eachLetter][0] in ['-','X','x']:
                if (eachLetter in siteInfoDic) and (siteInfoDic[eachLetter][1] not in ['-','X','x']):
                    flag = 1
            else:
                if (eachLetter in siteInfoDic) and (siteInfoDic[eachLetter][1] in ['-','X','x']):
                    flag = 1
            
        if flag == 1:
            DBseq = dicDB[each].split('~')[0]
            reference = '~'.join(dicDB[each].split("~")[1:-1])
            describe = dicDB[each].strip('\n').split("~")[-1]
            
            if each=='206&207&208&209&210&211&212&213&214&215':each='220-loop '
            outInfo = outInfo+'|'+DBseq.split('|')[0]+each+DBseq.split('|')[1]+'+'+reference+'+'+describe
            #print(outInfo)
    return outInfo
    

if __name__ == '__main__':
    dirUser = '../'
    addClusterTransitionSubstitution(resultFileDir=dirUser+'/result/',resultFileName = 'example'+'.result',DBFileName="DBSiteInfoClusterTransition",dir="/home/think/platform/18Mid/",proteinTypeList=eval("[('querySeq1', 'NS'), ('querySeq2', 'MP'), ('querySeq3', 'N6'), ('querySeq4', 'NP'), ('querySeq5', 'H1'), ('querySeq6', 'PA'), ('querySeq7', 'PB1'), ('querySeq8', 'PB2')]"))






