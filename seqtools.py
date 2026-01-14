def getCompBase(base):
    if (base == "A"):
        return "T"
    elif (base == "T"):
        return "A"
    elif (base == "G"):
        return "C"
    else:
        return "G"
    
def getCompSeq(seq):
    compSeq = ""
    for i in range(len(seq)):
        if seq[i] == "A":
            compSeq += "T"
        elif seq[i] == "T":
            compSeq += "A"
        elif seq[i] == "G":
            compSeq += "C"
        else:
            compSeq += "G"
    return compSeq

def getRevSeq(seq):
    return seq[::-1]
    
def checkPalindromicSeq(seq):
    return(getCompSeq(seq) == getRevSeq(seq))


            