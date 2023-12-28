import subprocess

def blastSeq(querySeqFile, seqType = "protein", DBDir = "", dataBaseName = "", eValue = "1e-5", outfmt = "3",
             outName = "", querySeqFileDir = "", outFileDir = ""):
    print("blastSeq ...")
    blastType = None
    if seqType == "protein":
        blastType = "blastp"
    elif seqType == "nucleo":
        blastType = "blastn"
    c = subprocess.getoutput(
        blastType + " -db " + DBDir + dataBaseName + " -query " + querySeqFileDir + querySeqFile + " -out " +
        outFileDir + outName + " -evalue " + eValue + " -num_threads 2 -outfmt " +
        outfmt + " -num_descriptions 5 -num_alignments 5\n" + "cd " + outFileDir + "\n")
    print(c)
