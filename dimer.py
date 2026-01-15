import seqtools as st

# Dictionary describing the contribution of each possible nearest-neighbor pair in a duplex towards the Gibbs free energy in kcal/mol
nnDuplexGContributions = {"AA/TT": -1.00, "TA/AT": -0.58, "GA/CT": -1.30, "CA/GT": -1.45, "AT/TA": -0.88, "TT/AA": -1.00,  "GT/CA": -1.44, 
                        "CT/GA": -1.28, "AG/TC": -1.28, "TG/AC": -1.45, "GG/CC": -1.84, "CG/GC": -2.17, "AC/TG": -1.44,  "TC/AG": -1.30,
                        "GC/CG": -2.24,  "CC/GG": -1.84}

class Dimer:
    '''
    Describes a dimer DNA structure

    seq1: A string representing a DNA sequence
    seq2: A string representing another DNA sequence
    overlapStart1: An integer representing the index on seq1 of where the overlap starts
    overlapEnd1: An integer reporesenting the index on seq1 of where the overlap ends
    overlapStart2: An integer representing the index on seq2 of where the overlap starts
    overlapEnd2: An integer reporesenting the index on seq2 of where the overlap ends
    '''
    def __init__(self, seq1, seq2, overlapStart1, overlapEnd1, overlapStart2, overlapEnd2):
        self.seq1 = seq1
        self.seq2 = seq2
        self.overlapStart1 = overlapStart1
        self.overlapEnd1 = overlapEnd1
        self.overlapStart2 = overlapStart2
        self.overlapEnd2 = overlapEnd2

    def getOverlapLen(self):
        return self.overlapEnd1 - self.overlapStart1 + 1

def calcScoreDimer(dimer):
    '''
    Checks a Dimer object for basic viability
    
    dimer: A DNA dimer (Dimer)
    Returns: A score representing the number of matches minus the number of mismatches in the overlap region (int)
    '''
    score = 0
    overlapLen = dimer.getOverlapLen()
    # Get the start indexes of the overlap for both strands
    seq1Pos= dimer.overlapStart1
    seq2Pos = dimer.overlapStart2
    # Check the number of base pair matches and mismatches
    for i in range(overlapLen):
        if (st.getCompBase(dimer.seq1[seq1Pos]) == dimer.seq2[seq2Pos]):
            score += 1
        else:
            score -= 1
    return score

def calcDimerGibbs(dimer):
    G = 0
    overlapLen = dimer.getOverlapLen()
    # Get the start indexes of the overlap for both strands
    seq1Pos= dimer.overlapStart1
    seq2Pos = dimer.overlapStart2
    for i in range(overlapLen):
        if