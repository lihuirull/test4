# def standardize(file,dir,outFileDir):
#     fileIn = open(dir+file,'r')
#     fileOut = open(outFileDir+file+".stdName",'w')
#     outFileName = file+".stdName"
#     text = fileIn.readlines()
#     n = 0
#     dic = {}
#
#     for each in text:
#         if each.find(">")!=-1:
#             n = n + 1
#             standardName = ">querySeq"+str(n)
#             dic[standardName]=str(each.strip())
#             # print(each.strip('\n'))
#             each = '\n'+standardName+"\n"
#
#             countATCG = 0
#             length = 0
#         else:
#             countATCG =countATCG + each.count("t")+each.count("a")+each.count("c")+each.count("g")+each.count("-")+each.count("T")+each.count("A")+each.count("C")+each.count("G")
#             length = length + len(each)-1
#             each = each.upper().strip('\n')
#             # print(each.count("t")+each.count("a")+each.count("c")+each.count("g")+each.count("-"),len(each)-1)
#         fileOut.write(each)
#     # dic = sorted(dic.items(), key=lambda d: d[0])
#     # print(dic)
#     fileOut.write('\n')
#     fileIn.close()
#     fileOut.close()
#     if (countATCG+1) / (length+1) > 0.9:
#         # print("DNA seq")
#         return outFileName,"nucleo",dic
#     else:
#         # print("Prot seq")
#         return outFileName,"protein",dic
from Bio import SeqIO
from Bio import SeqIO
import os
def standardize(filePath, outFileDir):
    """
    使用BioPython标准化基因序列文件，并判断序列是DNA序列还是蛋白质序列。
    :param file: 待处理的文件路径
    :param dir: 待处理文件的目录
    :param outFileDir: 输出文件的目录
    :return: 输出文件名、序列类型（'nucleo' 或 'protein'）、标准化后的序列字典
    """
    file = os.path.basename(filePath)
    outFileName = file + ".stdName"
    is_protein = False
    dic = {}
    n = 0
    # 读取并处理序列文件
    print(f"读入文件：{filePath}")
    print(f"输出文件：{outFileDir + outFileName}")

    with open(filePath, 'r') as inFile, open(outFileDir + "/" + outFileName, 'w') as outFile:
        for record in SeqIO.parse(inFile, "fasta"):
            n+=1
            standardName = ">querySeq" + str(n)
            sequence = str(record.seq).upper()

            dic[standardName] = f">{record.id}"  # 存储序列的ID和序列原始id
            outFile.write(f"{standardName}\n{sequence}\n")

            # 检查是否为蛋白质序列
            if not is_protein and any(char not in "ATCGN-" for char in sequence):
                is_protein = True

    sequence_type = "protein" if is_protein else "nucleo"
    return outFileName, sequence_type, dic


"""
打开输入和输出文件。
逐行读取输入文件。
如果行以“>”开始，表示这是一个新的序列的标记，生成一个标准化的序列名，并将其添加到字典中。
如果行不是序列标记，则将其转换为大写并去除两端空格，检查是否包含非ATCGN字符，如果有，则判定为蛋白质序列。
将处理后的行写入输出文件。
根据是否检测到蛋白质序列的特征，返回输出文件名、序列类型及序列字典。
"""
# standardize("test_query_seq2.fas","/home/think/18Mid/standard_seq/out/")