import datastructures as ds

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
        # Reverse the order of one stem in order to align it with the other
        revStem1 = hairpin.stem1[::-1]
        if (getCompBase(revStem1[i]) == hairpin.stem2[i]):
            score += 1
        else:
            score -= 1
    return score

def createPossHairpins(seq):
    '''
    Finds possible hairpins that can form from a given sequence. A minimum stem length of two was chosen as well as a loop size 
    between 4 and 5 bases based on real world data. These hairpins are then found using a basic sliding window algorithm.

    seq: A DNA sequence (str)
    Returns: A list of Hairpin objects (Hairpin)
    '''
    minStemLength = 2
    minLoopSize = 4
    possHairpins = []
    # Create initial hairpin
    stem1 = seq[0:minStemLength]
    loop = seq[minStemLength: minStemLength + minLoopSize]
    stem2 = seq[minStemLength + minLoopSize:]
    possHairpins.append(ds.Hairpin(stem1, loop, stem2))
    while (len(stem2) >= minStemLength):
        if (len(loop) == minLoopSize):
            # Add a base to the loop
            loop += stem2[0]
            # Take off the first base of stem2
            stem2 = stem2[1:]
        else:
            # Add a base to stem1
            stem1 += loop[0]
            # Take off the first base of loop
            loop = loop[1:]
        # Add hairpin to list
        possHairpins.append(ds.Hairpin(stem1, loop, stem2))
    return possHairpins





    



