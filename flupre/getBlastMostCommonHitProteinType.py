# -*- coding:UTF-8 -*-
#
# def getMostCommonHitProtein(result,dir):
#     print("getMostCommonHitProtein ...")
#     NUM_OF_RESULT = 5     #用前NUM_OF_RESULT个结果中出现最多次数的蛋白进行判断query蛋白是属于哪一种蛋白
#     num = NUM_OF_RESULT
#     from collections import Counter
#     mostCommonClass = []
#     fileIn = open(dir+result,"r",encoding='utf-8')
#     text = fileIn.readlines()
#     flag = 0
#     count = 100
#     for each in text:
#         if each[0:6]=="Query=":
#             each = each.split("=")[-1].lstrip(" ")
#             queryName = each.rstrip('\n')
#             flag = 1                             #得到搜索序列
#             count = 0                            #判断结果序列位置
#             # print()
#             # print(each.rstrip('\n'))
#             proteinClassName = []
#             num = NUM_OF_RESULT                  #通过一个临时的变量num传递NUM_OF_RESULT的值，以便在之后的对结果个数判断时对NUM的更改，而不会改变NUM_OF_RESULT的值
#         if flag == 1 and count < 7+num:
#             count = count + 1
#             if count >= 7 and count < 7+num:
#                 # print(len(each.strip("\n")))
#                 if len(each.strip("\n"))==0:                            #判断是否为空行，判断查询结果的末尾，如果为空行，说明查询结果不够NUM_OF_RESULT个
#                     num = count - 7                           #如果查询结果不够NUM_OF_RESULT个时，NUM_OF_RESULT变为查询结果的个数，如果查询结果大于NUM_OF_RESULT，则不变
#                     # print(num)                                #若某一条查询序列一条结果都没有，那么在文本第七行为空，故num=0
#                 name = each.split("_")[0].strip('\n')
#                 if name=='':
#                     name='Unknown'
#                 # print(queryName,name,"annanannananana")
#                 # print(each.rstrip("\n"))
#                 proteinClassName.append(name)
#                 # print(proteinClassName)
#             if count== 7+num:
#                 names = Counter(proteinClassName).most_common(1)[0]   #得到出现次数最多的蛋白种类
#                 mostCommonClass.append((queryName,names[0]))
#                 # print(mostCommonClass)
#
#
#     fileIn.close()
#     return mostCommonClass
#
# def getMostCommonHitProteinLowLevelHost(result,dir,Host):
#     print("getMostCommonHitProtein ...")
#     dicHighHost = {}
#     for eachHighHost in Host:
#         dicHighHost[eachHighHost[0]] = eachHighHost[1]
#     NUM_OF_RESULT = 5     #用前NUM_OF_RESULT个结果中出现最多次数的蛋白进行判断query蛋白是属于哪一种蛋白
#     num = NUM_OF_RESULT
#     from collections import Counter
#     mostCommonClass = []
#     fileIn = open(dir+result,"r",encoding='utf-8')
#     text = fileIn.readlines()
#     flag = 0
#     count = 100
#     for each in text:
#         if each[0:6]=="Query=":
#             each = each.split("=")[-1].lstrip(" ")
#             queryName = each.rstrip('\n')
#             flag = 1                             #得到搜索序列
#             count = 0                            #判断结果序列位置
#             # print()
#             # print(each.rstrip('\n'))
#             proteinClassName = []
#             num = NUM_OF_RESULT                  #通过一个临时的变量num传递NUM_OF_RESULT的值，以便在之后的对结果个数判断时对NUM的更改，而不会改变NUM_OF_RESULT的值
#         if flag == 1 and count < 7+num:
#             count = count + 1
#             if count >= 7 and count < 7+num:
#                 # print(len(each.strip("\n")))
#                 if len(each.strip("\n"))==0:                            #判断是否为空行，判断查询结果的末尾，如果为空行，说明查询结果不够NUM_OF_RESULT个
#                     num = count - 7                           #如果查询结果不够NUM_OF_RESULT个时，NUM_OF_RESULT变为查询结果的个数，如果查询结果大于NUM_OF_RESULT，则不变
#                     # print(num)                                #若某一条查询序列一条结果都没有，那么在文本第七行为空，故num=0
#                 if len(each.split("_"))>1:
#                     name = each.split("_")[1].strip('\n')
#                     nameHighLevel = each.split("_")[0].strip('\n')
#                     if nameHighLevel != dicHighHost[queryName]:
#                         continue
#                 else:
#                     name=''
#                 if name=='':
#                     name='Unknown'
#                 # print(queryName,name,"annanannananana")
#                 # print(each.rstrip("\n"))
#                 proteinClassName.append(name)
#                 # print(proteinClassName)
#             if count== 7+num:
#                 names = Counter(proteinClassName).most_common(1)[0]   #得到出现次数最多的蛋白种类
#                 mostCommonClass.append((queryName,names[0]))
#                 # print(mostCommonClass)
#
#
#     fileIn.close()
#     return mostCommonClass
from collections import Counter
import os

def getMostCommonHitProtein(result, dir):
    """
    从BLAST结果中提取最常见的蛋白质类型。
    """
    print("getMostCommonHitProtein ...")
    NUM_OF_RESULT = 5
    mostCommonClass = []

    with open(os.path.join(dir, result), "r", encoding='utf-8') as fileIn:
        proteinClassName = []
        queryName = None

        for each in fileIn:
            if each.startswith("Query="):
                if queryName is not None and proteinClassName:
                    name_count = Counter(proteinClassName).most_common(1)
                    mostCommonClass.append((queryName, name_count[0][0] if name_count else 'Unknown'))
                queryName = each.split("=")[-1].strip()
                proteinClassName = []
            else:
                if len(proteinClassName) < NUM_OF_RESULT:
                    name = each.split("_")[0].strip()
                    if name:
                        proteinClassName.append(name)

        if queryName is not None and proteinClassName:
            name_count = Counter(proteinClassName).most_common(1)
            mostCommonClass.append((queryName, name_count[0][0] if name_count else 'Unknown'))

    return mostCommonClass
def getMostCommonHitProteinLowLevelHost(result, dir, Host):
    """
    根据低级宿主信息从BLAST结果中提取最常见的蛋白质类型。
    """
    print("getMostCommonHitProtein ...")
    dicHighHost = {h[0]: h[1] for h in Host}
    NUM_OF_RESULT = 5
    mostCommonClass = []

    with open(os.path.join(dir, result), "r", encoding='utf-8') as fileIn:
        proteinClassName = []
        queryName = None

        for each in fileIn:
            if each.startswith("Query="):
                if queryName is not None and proteinClassName:
                    name_count = Counter(proteinClassName).most_common(1)
                    mostCommonClass.append((queryName, name_count[0][0] if name_count else 'Unknown'))
                queryName = each.split("=")[-1].strip()
                proteinClassName = []
            else:
                if len(proteinClassName) < NUM_OF_RESULT:
                    parts = each.strip().split("_")
                    if len(parts) > 1:
                        nameHighLevel, name = parts[0], parts[1]
                        if nameHighLevel == dicHighHost.get(queryName, ''):
                            proteinClassName.append(name)

        if queryName is not None and proteinClassName:
            name_count = Counter(proteinClassName).most_common(1)
            mostCommonClass.append((queryName, name_count[0][0] if name_count else 'Unknown'))

    return mostCommonClass

# print("\n\n",getMostCommonHitProtein("新建文本.txt","/home/think/18Mid/standard_seq/"))
