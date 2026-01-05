def calcMeltTempPrimer(seq):
    # Count the number of different nucleotides
    numA = 0
    numT = 0
    numG = 0
    numC = 0
    for letter in seq:
        if (letter == "A"):
            numA += 1
        elif (letter == "T"):
            numT += 1
        elif (letter == "G"):
            numG += 1
        else:
            numC += 1
    # Calculate the melting temperature differently based on the size of the strand
    if (len(seq) <= 13):
        meltTemp = (numA + numT) * 2 + (numG + numC) * 4
    else:
        meltTemp = 64.9 + 41 * (numG + numC - 16.4) / (len(seq))

def calcGCContent(seq):
    numG = 0
    numC = 0
    for letter in seq:
        if (letter == "G"):
            numG += 1
        elif (letter == "C"):
            numC += 1
    return (numG + numC) / len(seq) * 100

def calcGCInClamp(seq):
    # Count the number of G's and C's in the last five nucleotides
    numGC = 0
    for i in range(-5, 0):
        if (seq[i] == "G" or seq[i] == "C"):
            numGC += 1
    return numGC

def getCompBase(base):
    if (base == "A"):
        return "T"
    elif (base == "T"):
        return "A"
    elif (base == "G"):
        return "C"
    else:
        return "G"

def calcScoreHairpin(hairpin):
    '''
    Checks a Hairpin object for basic viability
    
    hairpin: A DNA hairpin (Hairpin)
    Returns: A score representing the number of matches minus the number of mismatches in the stem (int)
    '''
    # Get the shortest stem length
    if (len(hairpin.stem1) < len(hairpin.stem2)):
        minLenStem = len(hairpin.stem1)
    else:
        minLenStem = len(hairpin.stem2)
    # Check the number of base pair matches and mismatches
    for i in range(minLenStem):
        if (getCompBase(hairpin.stem1[i]) == hairpin.stem2[i]):
            score += 1
        else:
            score -= 1
    return score



    



