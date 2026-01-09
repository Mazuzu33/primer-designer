import datastructures as ds

def calcMeltTempSeq(seq):
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
    return meltTemp

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

def calcMaxVocab(alphaSize, wordLen):
    '''
    Helper function for calculating linguistic complexity

    alphaSize: The amount of letters in the alphabet (int)
    wordLen: The length of the word (int)
    Returns: The maximum vocabulary of words of length 1 to m that can be formed by taking substrings of a word of wordLen (int)
    ''' 
    maxVocab = 0
    # Calculate the maximum vocabulary of words of length i
    for i in range(1, wordLen + 1):
        # Represents the max number of unique words of length i with regards to the alphabet. 
        val1 = alphaSize ** i
        # Represents the max number of unique words of length i that can be obtained by taking substrings of a word of wordLen. 
        # Invalid when wordLen > alphaSize
        val2 = wordLen - i + 1
        # Choose the minimum of the two values
        maxVocab += min(val1, val2)
    return maxVocab

def calcLinComp(seq):
    '''
    Calculates the linguistic complexity (LC) of a sequence where LC is defined as the ratio of the number of substrings of any
    length in the given sequence to the maximum possible number of substrings obtainable from a sequence of the given sequence's length. 
    
    seq: A DNA sequence (string)
    '''
    trie = ds.SuffixTrie()
    trie.addSuffixes(seq)
    uniqueSubWords = trie.getEdges()
    maxVocab = calcMaxVocab(4, len(seq))
    return uniqueSubWords / maxVocab

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

def calcLoopGibbs(hairpin):
    '''
    Helper method for calculating the Gibbs free energy of a hairpin. Uses the nearest neighbor (NN) method.
    
    hairpin: A Hairpin object (Hairpin)
    Returns: The loop's Gibbs free energy contribution (int)
    '''
    G = 0
    # Dictionary describing the contribution of each possible nearest-neighbor pair in a loop towards the Gibbs free energy in kcal/mol
    nnLoopContribution = {"AA": 1.4, "TA": 1.2, "GA": 0.1, "CA": 0.4, "AT": 1.5, "TT": 1.0, "GT": 0.9, "CT": 0.2, "AG": 2.4, "TG": 1.7, 
                      "GG": 1.2, "CG": 0.0, "AC": 2.4, "TC": 1.9, "GC": 1.5, "CC": 1.0}
    # Combines the loop sequence with the base surrounding it on each side
    combinedSeq = hairpin.stem1[-1] + hairpin.loop + hairpin.stem2[0]
    for i in range(len(combinedSeq) - 1):
        G += nnLoopContribution[combinedSeq[i:i+2]]
    return G

def calcStemGibbs(hairpin):
    '''
    Helper method for calculating the Gibbs free energy of a hairpin. Uses the nearest neighbor (NN) method.
    
    hairpin: A Hairpin object (Hairpin)
    Returns: The stem's Gibbs free energy contribution (int)
    '''
    G = 0
    # Dictionary describing the contribution of each possible nearest-neighbor pair in a duplex towards the Gibbs free energy in kcal/mol
    nnStemContribution = {"AA/TT": -1.00, "TA/AT": -0.58, "GA/CT": -1.30, "CA/GT": -1.45, "AT/TA": -0.88, "TT/AA": -1.00,  "GT/CA": -1.44, "CT/GA": -1.28,
                        "AG/TC": -1.28, "TG/AC": -1.45, "GG/CC": -1.84, "CG/GC": -2.17, "AC/TG": -1.44,  "TC/AG": -1.30, "GC/CG": -2.24,  "CC/GG": -1.84}
    # Get the shortest stem length
    if (len(hairpin.stem1) < len(hairpin.stem2)):
        minLenStem = len(hairpin.stem1)
    else:
        minLenStem = len(hairpin.stem2)
    # Reverse the order of one stem in order to align it with the other
    revStem1 = hairpin.stem1[::-1]
    for i in range(minLenStem - 1):
        #Check if there is a mismatch anmd if so, assess a penalty of .438 kcal/mol
        if (getCompBase(revStem1[i]) != hairpin.stem2[i] or getCompBase(hairpin.stem2[i+1] != hairpin.stem2[i+1])):
            G += .438
        else:
            stemPair = (revStem1[i:i+1] + "/" + hairpin.stem2[i:i+1])
            G += nnStemContribution[stemPair]
    return G

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




    



