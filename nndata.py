#Nearest neighbor data

# Dictionary describing the contribution of each possible nearest-neighbor pair in a loop towards the Gibbs free energy in kcal/mol
loopPairToG = {"AA": 1.4, "TA": 1.2, "GA": 0.1, "CA": 0.4, "AT": 1.5, "TT": 1.0, "GT": 0.9, "CT": 0.2, "AG": 2.4, "TG": 1.7, 
                      "GG": 1.2, "CG": 0.0, "AC": 2.4, "TC": 1.9, "GC": 1.5, "CC": 1.0}

# Dictionary describing the contribution of each possible nearest-neighbor pair in a duplex towards the Gibbs free energy in kcal/mol
duplexPairToG = {"AA/TT": -1.00, "TA/AT": -0.58, "GA/CT": -1.30, "CA/GT": -1.45, "AT/TA": -0.88, "TT/AA": -1.00,  "GT/CA": -1.44, 
                        "CT/GA": -1.28, "AG/TC": -1.28, "TG/AC": -1.45, "GG/CC": -1.84, "CG/GC": -2.17, "AC/TG": -1.44,  "TC/AG": -1.30,
                        "GC/CG": -2.24,  "CC/GG": -1.84}


# End correction for a hairpin stem for Gibbs free energy in kcal/mol
stemTermgcAmtG = -2.182
stemTermatAmtG = -1.653

#Initiation paramters for duplex for Gibbs free energy in kcal/mol
duplexTermgcAmtG = 0.98
duplexTermatAmtG = 1.03
symAmtG = 0.43

# Mismatch penalty for Gibbs free energy in kcal/mol
mismatchPenalty = 0.438