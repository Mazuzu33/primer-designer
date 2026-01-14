import seqtools as st
import nndata as nn

class Hairpin:
    """
    Describes a hairpin DNA structure

    stem1: A string representing the DNA bases from the 5' end of a hairpin up to the start of the loop (exclusive)
    loop: A string representing the DNA bases of the loop of a hairpin
    stem2: A string representing the DNA bases from the end of the loop (exclusive) to the 3' end of a hairpin
    """
    def __init__(self, stem1, loop, stem2):
        self.stem1 = stem1
        self.loop = loop
        self.stem2 = stem2
        self.minStemLen = self.getMinStemLen()
        self.stemPairs = self.getStemPairs()
        self.revStem1 = stem1[::-1]
        self.score = self.calcScoreHairpin()
        self.loopGibbs = self.calcLoopGibbs()
        self.mismatchGibbs = self.calcMismatchGibbs()
        self.stemGibbs = self.calcStemGibbs() 
        self.hairGibbs = self.calcHairGibbs()
        self.stemEnthalpy = self.calcStemEnthalpy()
        self.stemEntropy = self.calcStemEntropy()
        self.loopEntropy = self.calcLoopEntropy()
        self.mismatchEnthalpy = self.calcMismatchEnthalpy()
        self.meltingTemperature = self.calcMeltingTemperature()

    def getMinStemLen(self):
        if (len(self.stem1) < len(self.stem2)):
            return len(self.stem1)
        else:
            return len(self.stem2)
        
    def getStemPairs(hairpin):
        '''
        Gets all valid nn pairs without mismatches for the stem of a hairpin
        
        hairpin: A Hairpin object (Hairpin)
        Returns: A list of nn pairs in the form of strings (list)
        '''
        nnPairs = []
        # Reverse the order of one stem in order to align it with the other
        # Iterate through possible NN pairs
        for i in range(hairpin.minStemLen - 1):
            #Check that there are no mismatches
            if (st.getCompBase(hairpin.revStem1[i]) == hairpin.stem2[i] or st.getCompBase(hairpin.revStem1[i+1]) == hairpin.stem2[i+1]):
                duplexPair = hairpin.revStem1[i:i+2] + "/" + hairpin.stem2[i:i+2]
                nnPairs.append(duplexPair)
        return nnPairs
        
    def calcScoreHairpin(hairpin):
        '''
        Checks a Hairpin object for basic viability
        
        hairpin: A DNA hairpin (Hairpin)
        Returns: A score representing the number of matches minus the number of mismatches in the stem (int)
        '''
        score = 0
        # Check the number of base pair matches and mismatches
        for i in range(hairpin.minStemLen):
            if (st.getCompBase(hairpin.revStem1[i]) == hairpin.stem2[i]):
                score += 1
            else:
                score -= 1
        return score
    
    def calcLoopGibbs(hairpin):
        '''
        Helper method for calcHairGibbs. Uses the nearest neighbor (NN) method.
        
        hairpin: A Hairpin object (Hairpin)
        Returns: The loop's Gibbs free energy contribution (int)
        '''
        G = 0
        # Combines the loop sequence with adjacent bases
        combinedSeq = hairpin.stem1[-1] + hairpin.loop + hairpin.stem2[0]
        # Iterate through NN pairs
        for i in range(len(combinedSeq) - 1):
            G += nn.loopPairToG[combinedSeq[i:i+2]]
        return G
    
    def calcMismatchGibbs(hairpin):
        '''
        Calculates the mismatch contribution to the Gibbs free energy of a hairpin
        
        hairpin: A Hairpin object (Hairpin)
        Returns: The mismatch contribution to the Gibbs free energy (float)
        '''
        G = 0
        for i in range(hairpin.minStemLen):
            # Reverse the order of one stem in order to align it with the other
            revStem1 = st.getRevSeq(hairpin.stem1)
            # Add a penalty for each mismatch
            if (st.getCompBase(revStem1[i]) != hairpin.stem2[i]):
                G += nn.mismatchPenalty
        return G
    
    def calcStemGibbs(hairpin):
        '''
        Helper method for calcHairGibbs. Uses the nearest neighbor (NN) method.
        
        hairpin: A Hairpin object (Hairpin)
        Returns: The stem's Gibbs free energy contribution (int)
        '''
        G = 0
        # Iterate through the valid NN pairs
        for nnPair in hairpin.stemPairs:
            G += nn.duplexPairToG[nnPair]
        # End correction based on the pair at the end of the stem and at the beginning of the loop
        # Checks that the base pair is matched correctly first
        if (st.getCompBase(hairpin.stem1[-1]) == hairpin.stem2[0]):
            # Add end correction based on the identity of pair
            if (hairpin.stem1 == "G" or hairpin.stem1 == "C"):
                G += nn.stemTermgcAmtG
            else:
                G += nn.stemTermatAmtG
        return G
    
    def calcHairGibbs(hairpin):
        '''
        Calculates the Gibbs free energy in kcal/mol of the hairpin
        
        hairpin: A Hairpin object (Hairpin)
        '''
        # Add in Gibbs free energy calculations for the stem and loop and mismatches
        G += hairpin.stemGibbs + hairpin.loopGibbs + hairpin.mismatchGibbs
        return G
    
    def calcStemEnthalpy(hairpin):
        '''
        Helper method for calcHairEnthalpy. Uses the nearest neighbor (NN) method.
        
        hairpin: A Hairpin object (Hairpin)
        Returns: The stem's enthalpy contribution (int)
        '''
        termgc = 0
        termat = 0
        H = 0
        # Iterate through the valid NN pairs
        for nnPair in hairpin.stemPairs:
            H += nn.duplexPairToH[nnPair]
        # End correction based on the pairs at the end of the stem
        # Calculate number of terminating GC and AT pairs
        if (st.getComp(hairpin.revStem1[0]) == hairpin.stem2[0]):
            if (hairpin.revStem1[0] == "G" or hairpin.revStem1[0] == "C"):
                termgc += 1
            else:
                termat += 1
        if (st.getCompBase(hairpin.revStem1[hairpin.minStemLen]) == hairpin.stem2[hairpin.minStemLen]):
            if (hairpin.revStem1[hairpin.minStemLen] == "G" or hairpin.revStem1[hairpin.minStemLen] == "C"):
                termgc += 1
            else:
                termat += 1
        # Add end correction based on the identity of pairs
        H += termgc * nn.duplexTermgcAmtH + termat * nn.duplexTermatAmtH
        return H
    
    def calcStemEntropy(hairpin):
        '''
        Helper method for calcHairEntropy. Uses the nearest neighbor (NN) method.
        
        hairpin: A Hairpin object (Hairpin)
        Returns: The stem's entropy contribution (int)
        '''
        S = 0
        # Iterate through the valid NN pairs
        for nnPair in hairpin.stemPairs:
            S += st.nnDuplexSContributions[nnPair]
        # End correction based on the pairs at the end of the stem
        # Calculate number of terminating GC and AT pairs
        if (st.getCompBase(hairpin.revStem1[0]) == hairpin.stem2[0]):
            if (hairpin.revStem1[0] == "G" or hairpin.revStem1[0] == "C"):
                termgc += 1
            else:
                termat += 1
        if (st.getComp(hairpin.revStem1[hairpin.minStemLen]) == hairpin.stem2[hairpin.minStemLen]):
            if (hairpin.revStem1[hairpin.minStemLen] == "G" or hairpin.revStem1[hairpin.minStemLen] == "C"):
                termgc += 1
            else:
                termat += 1
        # Add end correction based on the identity of pairs
        S += termgc * nn.duplexTermgcAmtS + termat * nn.duplexTermatAmtS

    def calcLoopEntropy(hairpin):
        '''
        Calculate the entropy contribution of the hairpin loop in kcal / mol K
        
        :param hairpin: Description
        '''
        # Add in the free energy contributed by the loop
        G = hairpin.stemGibbs + hairpin.loopGibbs
        H = hairpin.stemEnthalpy
        S = hairpin.stemEntropy
        T = 310
        # Adjust S appropriately to account for the change in free energy
        Snew = (H - G)/T
        # Get the contribution to entropy by the loop itself
        return Snew - S
    
    def calcMismatchEnthalpy(hairpin):
        '''
        Calculate the enthalpy contribution of the mismatches in kcal / mol
        
        :param hairpin: Description
        '''
        # Add in the free energy contributed by the mismatches
        G = hairpin.stemGibbs + hairpin.mismatchGibbs
        H = hairpin.stemEnthalpy
        S = hairpin.stemEntropy
        T = 310
        # Adjust H appropriately to account for the change in free energy
        Hnew = G + T*S
        # Get the contribution to enthalpy by the mismatcheds itself
        return Hnew - H
    
    def calcMeltingTemperature(hairpin):
        S = hairpin.stemEntropy + hairpin.loopEntropy
        H = hairpin.stemEnthalpy + hairpin.mismatchEnthalpy
        Tm = H/S
        return Tm

    def __str__(self):
        return f"Stem1: {self.stem1}\nLoop: {self.loop}\nStem2: {self.stem2}"

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
    # No hairpins possible if sequence is smaller than 8 bp
    if (len(seq) < 8):
        return possHairpins
    # Create initial hairpin
    stem1 = seq[0:minStemLength]
    loop = seq[minStemLength: minStemLength + minLoopSize]
    stem2 = seq[minStemLength + minLoopSize:]
    possHairpins.append(Hairpin(stem1, loop, stem2))
    while (len(stem2) > minStemLength or len(loop) > minLoopSize):
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

