from dataclasses import dataclass

@dataclass
class Hairpin:
    """
    Describes a hairpin DNA structure

    stem1: A string representing the DNA bases from the 5' end of a hairpin up to the start of the loop (exclusive)
    loop: A string representing the DNA bases of the loop of a hairpin
    stem2: A string representing the DNA bases from the end of the loop (exclusive) to the 3' end of a hairpin
    """
    stem1: str
    loop: str
    stem2: str

@dataclass
class Primer:
    '''
    Describes a primer

    seq: A string representing the DNA bases of the primer
    start: A integer representing the start location of the primer on the original inputted sequence
    end: A integer representing the end location of the primer on the original inputted sequence
    meltTemp: A float representing the melting temperature of the primer
    gcContent: A float representing the gcContent of the primer
    gcInClamp: A integer representing the number of G's and C's in the GC clamp of the primer
    '''
    seq: str
    start: int
    end: int
    meltTemp: float
    gcContent: float
    gcInClamp: int



