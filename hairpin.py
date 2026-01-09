import seqtools as st

# Dictionary describing the contribution of each possible nearest-neighbor pair in a loop towards the Gibbs free energy in kcal/mol
nnLoopContributions = {"AA": 1.4, "TA": 1.2, "GA": 0.1, "CA": 0.4, "AT": 1.5, "TT": 1.0, "GT": 0.9, "CT": 0.2, "AG": 2.4, "TG": 1.7, 
                      "GG": 1.2, "CG": 0.0, "AC": 2.4, "TC": 1.9, "GC": 1.5, "CC": 1.0}

# Dictionary describing the contribution of each possible nearest-neighbor pair in a duplex towards the Gibbs free energy in kcal/mol
nnDuplexContributions = {"AA/TT": -1.00, "TA/AT": -0.58, "GA/CT": -1.30, "CA/GT": -1.45, "AT/TA": -0.88, "TT/AA": -1.00,  "GT/CA": -1.44, 
                        "CT/GA": -1.28, "AG/TC": -1.28, "TG/AC": -1.45, "GG/CC": -1.84, "CG/GC": -2.17, "AC/TG": -1.44,  "TC/AG": -1.30,
                        "GC/CG": -2.24,  "CC/GG": -1.84}

# Initiaton paramters for Gibbs free energy in kcal/mol
termGC = -2.182
termAT = -1.653

# Mismatch penalty in kcal/mol
mismatchPenalty = 0.438

class Hairpin:
    """
    Describes a hairpin DNA structure

    stem1: A string representing the DNA bases from the 5' end of a hairpin up to the start of the loop (exclusive)
    loop: A string representing the DNA bases of the loop of a hairpin
    stem2: A string representing the DNA bases from the end of the loop (exclusive) to the 3' end of a hairpin
    """
    def __init__self(self, stem1, loop, stem2):
        self.stem1 = stem1
        self.loop = loop
        self.stem2 = stem2

    def getShortestStemLen(self):
        if (len(self.stem1) < len(self.stem2)):
            return len(self.stem1)
        else:
            return len(self.stem2)
        
    def __str__(self):
        return f"Stem1: {self.stem1}\nLoop: {self.loop}\nStem2: {self.stem2}"
    


def calcScoreHairpin(hairpin):
    '''
    Checks a Hairpin object for basic viability
    
    hairpin: A DNA hairpin (Hairpin)
    Returns: A score representing the number of matches minus the number of mismatches in the stem (int)
    '''
    # Get the shortest stem length
    minLen = hairpin.getShortestStemLen()
    # Check the number of base pair matches and mismatches
    for i in range(minLen):
        # Reverse the order of one stem in order to align it with the other
        revStem1 = st.getRevSeq(hairpin.stem1)
        if (st.getCompBase(revStem1[i]) == hairpin.stem2[i]):
            score += 1
        else:
            score -= 1
    return score

def calcLoopGibbs(hairpin):
    '''
    Helper method for calcHairpinGibbs. Uses the nearest neighbor (NN) method.
    
    hairpin: A Hairpin object (Hairpin)
    Returns: The loop's Gibbs free energy contribution (int)
    '''
    G = 0
    # Combines the loop sequence with adjacent bases
    combinedSeq = hairpin.stem1[-1] + hairpin.loop + hairpin.stem2[0]
    # Iterate through NN pairs
    for i in range(len(combinedSeq) - 1):
        G += nnLoopContributions[combinedSeq[i:i+2]]
    return G

def calcStemGibbs(hairpin):
    '''
    Helper method for calcHairpinGibbs. Uses the nearest neighbor (NN) method.
    
    hairpin: A Hairpin object (Hairpin)
    Returns: The stem's Gibbs free energy contribution (int)
    '''
    G = 0
    # Get the shortest stem length
    minLen = hairpin.getShortestStemLen()
    # Reverse the order of one stem in order to align it with the other
    revStem1 = hairpin.stem1[::-1]
    # Iterate through NN pairs
    for i in range(minLen - 1):
        #Assess a penalty of .438 kcal/mol for a mismatch
        if (st.getCompBase(revStem1[i]) != hairpin.stem2[i] or st.getCompBase(revStem1[i+1]) != hairpin.stem2[i+1]):
            G += mismatchPenalty
        # Otherwise add the NN pair's Gibbs free energy contribution
        else:
            duplexPair = (revStem1[i:i+1] + "/" + hairpin.stem2[i:i+1])
            G += nnDuplexContributions[duplexPair]
    return G

def calcHairGibbs(hairpin):
    '''
    Calculates the Gibbs free energy in kcal/mol of the hairpin
    
    hairpin: A Hairpin object (Hairpin)
    '''
    # Reverse the order of one stem in order to align it with the other
    revStem1 = hairpin.stem1[::-1]
    #Obtain the terminal pairs at both ends of the stem duplex
    minLen = hairpin.getShortestStemLen()
    termPair = revStem1[0] + "/" + hairpin.stem2[0]
    # Check for G/C or A/T pairs and add appropriate initiation paramater
    if (termPair == "G/C" or termPair == "C/G") {
        G += termGC
    }
    elif (termPair == "A/T" or termPair == "T/A"):
        G += termAT
    else:
        G += mismatchPenalty
    # Add in Gibbs free energy calculations for the stem and loop
    G += calcStemGibbs(hairpin) + calcLoopGibbs(hairpin)
    return G
    
def createPossHairpins(seq):
    '''
    Finds possible hairpins that can form from a given sequence. A minimum stem length of two was chosen as well as a loop size 
    between 4 and 5 bases based on real world data. These hairpins are then found using a basic sliding window algorithm.

    seq: A DNA sequence (str)
    Returns: A list of Hairpin objects (list)
    '''
    minStemLength = 2
    minLoopSize = 4
    possHairpins = []
    # Create initial hairpin
    stem1 = seq[0:minStemLength]
    loop = seq[minStemLength: minStemLength + minLoopSize]
    stem2 = seq[minStemLength + minLoopSize:]
    possHairpins.append(Hairpin(stem1, loop, stem2))
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
        possHairpins.append(Hairpin(stem1, loop, stem2))
    return possHairpins