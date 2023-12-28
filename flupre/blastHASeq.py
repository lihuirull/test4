import os


def renameFastaSeq():
    # 把H1.fas-->H1
    for i in range(1, 17):
        fileIn = open('../18Mid/antigen/HA_blast/HA_protein/H' + str(i) + '.faa', 'r')
        textIn = fileIn.readlines()
        fileIn.close()

        fileOut = open('../18Mid/antigen/HA_blast/H' + str(i), 'w')
        for eachLine in textIn:
            # print(eachLine)
            eachLine = eachLine.strip('\n')
            if '>' in eachLine:
                print(eachLine.split(' | '))
                eachLine = '>' + eachLine.split(' | ')[1] + '____GISAID_' + eachLine.split(' | ')[0].strip(
                    '>') + '____PMID_0\n'
            else:
                eachLine = eachLine + '\n'
            fileOut.write(eachLine)
        fileOut.close()


def makeSeqDic():
    """
        为每个HA亚型的序列创建一个字典文件，方便之后的处理。
        id 和包含id的全序列
    """
    for i in range(1, 17):
        In = '../18Mid/antigen/HA_blast/HA_protein/H' + str(i)
        fileIn = open(In, 'r')
        textIn = fileIn.read()
        fileIn.close()
        dic = {}
        fileOut = open(In + '.dic', 'w')
        for eachSeq in textIn.split('>'):
            key = eachSeq.split('\n')[0].split('____')[0]
            value = '>' + eachSeq
            print(key)
            print(value)
            dic[key] = value
        fileOut.write(str(dic))
        fileOut.close()
        # break


def removeDup():
    """
    移除重复的序列，保留唯一的记录。
    """
    for i in range(1, 17):
        Out = '../18Mid/antigen/HA_blast/HA_protein/H' + str(i)
        fileIn = open(Out + '.dic', 'r')
        dic = eval(fileIn.read())
        fileIn.close()

        fileOut = open(Out, 'w')
        for each in dic:
            fileOut.write(dic[each])
        fileOut.close()


def makeBLastDB():
    """
        为每个HA亚型的序列构建BLAST数据库。得到HA_DB
    """
    import os
    for i in range(1, 17):
        In = '../18Mid/antigen/HA_blast/HA_protein/H' + str(i)
        Out = '../18Mid/antigen/HA_blast/HA_DB/H' + str(i)
        os.system('../app/blast+/bin/makeblastdb -in ' + In + ' -dbtype prot -parse_seqids -out ' + Out + '/H' + str(i))

def blastHASeq(prefixDir, blastQuerySeqHA, HAType, tempDir, mafftDir):
    """
    执行对HA亚型序列的BLAST搜索，并处理搜索结果。
    :param blastQuerySeqHA: BLAST查询序列文件
    :param HAType: HA亚型
    :param tempDir: 临时文件目录
    :param mafftDir: MAFFT工具目录
    """
    # prefixDir = "D:/Learning/edegfile/think/platform"
    with open(f'{prefixDir}/18Mid/antigen/HA_blast/HA_protein/{HAType}.dic', 'r') as fileDic:
        dic = eval(fileDic.read())

    # 执行BLAST搜索
    os.system(
        # D:/user/app/blast-BLAST_VERSION+/bin/blastp
        f'blastp -db {prefixDir}/18Mid/antigen/HA_blast/HA_DB/{HAType}/{HAType} -query {tempDir}{blastQuerySeqHA} -out {tempDir}{blastQuerySeqHA}.blastHAOut -evalue 1e-5 -num_threads 2 -outfmt 7 \n')

    # 处理BLAST结果
    with open(f'{tempDir}{blastQuerySeqHA}.blastHAOut', 'r') as fileIn:
        textIn = fileIn.read()
    blastResult = textIn.split(' hits found\n')[1].split('# BLAST processed')[0].strip('\n')

    with open(f'{tempDir}HAHitSeq', 'w') as fileHitSeq:
        lines = blastResult.split('\n')
        for eachLine in lines[:min(5, len(lines))]:
            key = eachLine.split('\t')[1].split('____')[0]
            # 通过前面得到的dic文件获取HA蛋白对应的Top5 hit的序列
            fileHitSeq.write(dic[key])

    # 使用Perl脚本进一步处理
    os.system(
        f"perl {prefixDir}/18Mid/antigen/predAV.pl {HAType} {tempDir}{blastQuerySeqHA}"
        f" {tempDir}HAHitSeq {tempDir}{blastQuerySeqHA}.antigenInfo.blast {prefixDir}/18Mid/antigen/modelData/ {tempDir} ")


# c = subprocess.getoutput(dirBlast+blastType+" -db "+DBDir+dataBaseName+" -query "+querySeqFileDir+querySeqFile+" -out "+outFileDir+outName+" -evalue "+eValue+" -num_threads 2 -outfmt "+outfmt+" -num_descriptions 5 -num_alignments 5\n"+"cd "+outFileDir+"\n")

if __name__ == "__main__":
    # 前四个函数之前已经运行过，得到了目前目录下的所有的结果
    # renameFastaSeq()
    # makeBLastDB()
    # removeDup()
    # makeSeqDic()
    # blastHASeq(blastQuerySeqHA = 'querySeq1_HA',HAType = 'H5',tempDir = '../temp/',mafftDir = '../app/mafft/mafft-7.158-without-extensions/scripts/')
    blastHASeq(blastQuerySeqHA = 'querySeq1', HAType = 'H1', tempDir = '/home/think/platform/temp/',
               mafftDir = '/home/think/platform/app/mafft/mafft-7.158-without-extensions/scripts/')
