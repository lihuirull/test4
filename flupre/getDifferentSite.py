# -*- coding:UTF-8 -*-

def getdifferentSite(file,fileDir,DBSite):
    print("getdifferentSite ...")

    dicDB = eval(DBSite)
    DBSiteNum = ""
    for eachDBNum in dicDB:
        DBSiteNum =DBSiteNum+"&"+eachDBNum

    fileIn = open(fileDir+file)
    fileOut = open(fileDir+file+".siteInfo","w")
    text = fileIn.readlines()
    line = 0
    siteInfo = ""
    hyphenNumber = {}           #用来记录相对应列的空位数
    sifferentSiteDic = {}
    DBsiteDic = {}
    for each in text:
        line = line + 1
        if line == 1:
            # print(each.split('\t')[0])  #记录第一行，每一列的名称
            name = each.split()
            # hyphenNumber = dict().fromkeys([each for each in range(0, len(name))], 0)
            continue
        if line == 2:
            startSite = each.split("\t")
            for i in range(0,len(name)):
                # hyphenNumber[i]= -int(startSite[i])
                hyphenNumber[i]= 0
                sifferentSiteDic.setdefault(i,[])
                DBsiteDic.setdefault(i,[])
            # print(len(name), hyphenNumber)

            continue
        # print(each.split())
        colonm = 0
        for eachSite in each.strip("\n").replace("\t",''):          #对该行的每一列进行遍历
            if colonm == 0:                                         #得到查询序列在该行的氨基酸
                firstSite = eachSite
            # print(eachSite)
            if eachSite == "-":                                     #如果该位点是空位，相应的序列的hyphenNumber计数加一
                hyphenNumber[colonm] +=1
                # print(hyphenNumber[colonm],"--",colonm)
            # if eachSite != firstSite:
            #     # print(firstSite+"__"+eachSite+"__"+str(line)+"__"+str(line-hyphenNumber[colonm])+"__"+str(colonm)+"__"+str(name[colonm]))
            #     site = firstSite+"__"+eachSite+"__"+str(line)+"__"+str(line-hyphenNumber[colonm]-2)+"__"+str(colonm)+"__"+str(name[colonm])
            #     siteInfo = siteInfo+"\n"+site
            #     sifferentSiteDic[colonm].append(firstSite+"__"+eachSite+"__"+str(line-hyphenNumber[0]-2)+"__"+str(line-hyphenNumber[colonm]-2))
            if str(line-hyphenNumber[colonm]-2) in DBSiteNum.split("&"):
                DBsiteDic[colonm].append(firstSite + "__" + eachSite + "__" + str(line - hyphenNumber[0]-2) + "__" + str(line-hyphenNumber[colonm]-2))
                #(line - hyphenNumber[0]-2)行数减去空位数减2等于查询序列位点的序号，
                #(line-hyphenNumber[colonm]-2)行数减去相对应的列的空位数减2等于标准序列位点的序号

            colonm = colonm + 1

    fileIn.close()
    # print(hyphenNumber)
    # print(siteInfo)
    # fileOut.write(str(name)+str(DBsiteDic).replace("],","],\n").replace("{","\n{"))
    fileOut.write(str(name)+str(DBsiteDic).replace("{","\n{"))

    fileOut.close
    differentSiteFileName = file + ".siteInfo"
    return differentSiteFileName,sifferentSiteDic , name , DBsiteDic
# getdifferentSite("querySeq1.MargeQueryAndStandard.fas.mafft.out.fas.transposed","/home/think/18Mid/standard_seq/out/outout/","{'199': 'A|S','318': 'K|R','355':'K|RQ','508':'R|Q','627':'K|E','675':'L|I','683':'T|A','701':'N|D','627&701':'KN|ED'}")
