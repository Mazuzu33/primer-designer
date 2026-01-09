class Primer:
    '''
    Describes a primer

    seq: A string representing the DNA sequence of the primer
    start: A integer representing the start location of the primer on the original inputted sequence
    end: A integer representing the end location of the primer on the original inputted sequence
    '''
    def __init__(self, seq, start, end):
        self.seq = seq
        self.start = start
        self.end = end
        self.meltTemp = self.calcMeltTemp(seq)
        self.gcContent = self.calcGCContent(seq)
        self.gcInClamp = self.calcGCInClamp(seq)

    def calcMeltTemp(seq):
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


    
