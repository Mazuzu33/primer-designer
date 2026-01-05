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