# -*- coding:UTF-8 -*-

def getAlignedSite(file="",dir=""):
    print("getAlignedSite ...")
    fileIn = open (dir+file,"r")
    text = fileIn.readlines()

    flag = 0
    name = []
    alignedSeqBeforeFilledStarKey = []      #初始化blast后提取的序列
    eachSeqStartAlignedPos = []             #初始化blast对齐的每一条序列的起点数组
    for each in text:
        if len(each)>2:
            flag = 1
        if flag == 1 and len(each)>2:
            name.append(each.split()[0])
            alignedSeqBeforeFilledStarKey.append(each.split("\t")[2])
            eachSeqStartAlignedPos.append(each.split()[1])

                                                             #blast中如果库中的序列在末尾部分与query序列不匹配，
                                                              # 在blast文件中在最后一行就直接省略掉了，不会用空格补充，
                                                                # 所以在此进行补全，避免转置矩阵时出现错误。
    alignedSeq = []
    flag1 = 0
    # fileseq = open(dir+file+".pureseq","w")
    for each in alignedSeqBeforeFilledStarKey:
        if flag1 == 0:
            flag1 = 1
            querySeqLength = len(each)
        each = str(each).rstrip("\n")+"*"*(querySeqLength-len(each))+"\n"
        alignedSeq.append(each)
        # fileseq.write(each)
    # fileseq.close()



    # seqNum = len(alignedSeq)
    seqNum = 2
    seqLenth = len(alignedSeq[0])-1           #因每行末尾有换行符，故减一
    # seqArrayTransposed = [([""]*seqNum)for i in range(seqLenth)]       #初始化转置序列数组
    # print(len(seqArrayTransposed),len(seqArrayTransposed[0]))
    seqArray = [([" "]*(seqLenth))for i in range(seqNum)]              #初始化拆分（blast后提取的序列）的数组
    for i in range(0,seqNum):
        for j in range(0,seqLenth):
            Column = alignedSeq[i][j]                  #在每一行得到每一列的元素
            # seqArrayTransposed[j][i] = Column       #转置原来的矩阵
            # seqArray[i].append(Column)            #将原来的字符串拆分成数组
            seqArray[i][j] = (Column)            #将原来的字符串拆分成数组
    # print(len(seqArray),seqArray)           #输出拆分后的字符串
    seqArrayTransposed = zip(*seqArray)         #转置拆分后的数组


    fileOut = open(dir +file+".transposed", "w")
    fileOut.write("\t".join(name)+"\n")
    fileOut.write("\t".join(eachSeqStartAlignedPos)+"\n")       #将对齐的每一条序列的起点 写入文件第一行中
    seqArrayTransposed = list(seqArrayTransposed)

    for l in seqArrayTransposed:
        fileOut.write("\t".join(l)+"\n")                #写入转置后的矩阵，每一行为对齐了的位点，由第一行可计算出每一个位点在每一条序列上的位置
    fileOut.close()
    fileIn.close()

    fileOutName = file+".transposed"
    return fileOutName , dir

# getAlignedSite("query___H6_standard_KJ200805___In___H6___DB.margeseq","/home/think/18Mid/standard_seq/out/")
# getAlignedSite("H6PlusQuerySeq.fas.fftns2.out.fas","/home/think/18Mid/standard_seq/out/")
