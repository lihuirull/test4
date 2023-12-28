# fileIn = open("/home/think/platform/18Mid/standard_seq/allProteinTypeDB/protType_AA.fas",'r')
# fileOut = open("/home/think/platform/18Mid/standard_seq/allProteinTypeDB/aaaaaaa.fas",'w')
# text = fileIn.readlines()
# for each in text:
#     if each[0] =='>' and "_PB1-F2(" in each or "_PA-X(" in each:
#         eachsplit = each.split("_")[-1].split("(")[0]
#
#         # each = each[3:]
#         each = each.replace(">PA",">PA-X").replace("PB1","PB1-F2(")
#         print(eachsplit,each)
#     fileOut.write(each)
# fileOut.close()
# fileIn.close()



#
# import sys
# print((sys.argv))
#
# import sys,getopt,os
# def main(argv):
#     print(sys.argv[0],sys.argv[1],len(argv))
#     print(111)
#
# if __name__ == "__main__":
#     main(sys.argv)


# fileIn =open("/home/think/1231.result",'r')
# text = fileIn.readlines()
# text.sort(key=lambda x:len(x))
# fileOut = open("/home/think/1231111111.result",'w')
# for each in text:
#
#     fileOut.write((each))
# fileOut.close()
# fileIn.close()


# def as_num(x):
#     y='{:.5f}'.format(x) # 5f表示保留5位小数点的float型
#     return(y)
# print(as_num(float('9.44462353318463e-05')))


# text.sort(key=lambda x: len(x))
# # text = (key=lambda x:"Max" not in x  and '{:.10f}'.format(float(x.split('\t')[2])),text)
# text.sort(key=lambda x: "Max" not in x ； '{:.10f}'.format(float(x.split('\t')[2



import sys, getopt
DIR = "/home/think/platform/"
resultFileDir = "/home/think/"
tempFileDir = "/home/think/platform/18Mid/standard_seq/temp/"
workName = "test_query_seq1000"
inputFile = "test_query_seq10.fas"
inputFileDir = DIR+"/18Mid/standard_seq/querySeq//"
opts, args = getopt.getopt(sys.argv[1:], "hi:o:",["inputFileDir=","inputFile=","tempFileDir=","resultFileDir=","resultFileName=","folderDir="])
for op, value in opts:
    if op == "--inputFileDir":
        inputFileDir = value
    elif op == "--inputFile":
        inputFile = value
    elif op == "--tempFileDir":
        tempFileDir = value
    elif op == "--resultFileDir":
        resultFileDir = value
    elif op == "--resultFileName":
        workName = value
    elif op == "--folderDir":
        DIR = value
        # sys.exit()
print(inputFileDir,inputFile,tempFileDir,resultFileDir,workName,DIR)