# import subprocess
# c = subprocess.getoutput('ls\nls')
# print(c)
# print()
# a = subprocess.call('ls')
# print(a)



# a = [([""]*5)for i in range(2)]
# a[0][0] = "a"
# print(a)
#
# myList = [([] * 3) for i in range(4)]
# print(myList)
#
# myList[1][1] = "a"
# print(myList)

#
# aaa = ["asdfgh","zxcvbnm"]
# for each in zip(aaa):
#     print(each)

#
# name = ["aaaa","bbbb",'cccc','dddd']
# y = iter(name)
# print(next(y))
# print(name)
# print([next(iter(range(0,5)))])
# # dict1 = {name[for i in range(0,5)]}
# # print(dict1)
#
# dict2 = {y:""}
# print(dict2[y])


#
# a = "123456789\n"
# print(a[:-2])

#
# from Bio.Align.Applications import ClustalwCommandline
# in_file = "/home/think/18Mid/standard_seq/test_query_seq2.fas"
# clustalw_cline = ClustalwCommandline("clustalw", infile=in_file)
# print(clustalw_cline)

# file = open(r'./H6_standard_KJ200805','r')
# seq = file.readlines()
# print(seq)
# file.close()
#
# dic = {}
# dic[0] = 1
# dic[0] = dic[0] + 1
# print(dic[0])

# for each in range(0,32):
#     dic[each] = 1
# print(dic)
# print(list(each for each in range(0,32)))

# info = dict().fromkeys([each for each in range(0,32)], 0)
# info[2] += 1
# print(info)


# d1={}
# key=1
# value=2
# d1.setdefault(key,[]).append()
# print(d1)


# string = ""
# for i in range(49,73):
#     string = string +"&"+ str(i)
# print(string)

# import itertools
# lists = ['a','b','c','d']
# a  = list(itertools.permutations(lists,len(lists)))
# print(str(a))

# a = list((r"babababa"))
# print(a)
# c = a.index('a')
# print(c)
# import re
# m = re.match("aa[a-z]ca","aacca")
# print((m))#匹配不成功返回None,匹配成功返回位置详细信息
#
# DBSeq = "P[A-Z][AVLIM]P[A-Z][RK]"
# DBSeq = ""
# m = re.match(DBSeq,"PAAPAK")
# print((m!=None))
# if re.match(DBSeq,"PAAPAK"):
#     print(11)
# else:
#     print(22)




# e = "---"
# print(e.replace("-",'')=="")

# text = ""
# import os
# for root , dirs , files in os.walk("/home/think/18Mid/GisAidSeq"):
#     print(files)
#     print(root)
#
#     for each in files:
#         file = open(root+"/"+each,'r')
#         text = text + file.read()
#         file.close()
#     fileOut = open("marge.fa",'w')
#     fileOut.write(text)
#     fileOut.close()





# import datetime
# nowTime=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')#现在
# print(nowTime)



import sys
print((sys.argv))

import sys,getopt,os
def main(argv):
    print(sys.argv[0],sys.argv[1],len(argv))
    print(111)

if __name__ == "__main__":
    main(sys.argv)

