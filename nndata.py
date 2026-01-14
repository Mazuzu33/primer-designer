#Nearest neighbor data

# Dictionary describing the contribution of each possible nearest-neighbor pair in a loop towards the Gibbs free energy in kcal/mol
loopPairToG = {"AA": 1.4, "TA": 1.2, "GA": 0.1, "CA": 0.4, "AT": 1.5, "TT": 1.0, "GT": 0.9, "CT": 0.2, "AG": 2.4, "TG": 1.7, 
                      "GG": 1.2, "CG": 0.0, "AC": 2.4, "TC": 1.9, "GC": 1.5, "CC": 1.0}

# Dictionary describing the contribution of each possible nearest-neighbor pair in a duplex towards the Gibbs free energy in kcal/mol
duplexPairToG = {"AA/TT": -1.00, "TA/AT": -0.58, "GA/CT": -1.30, "CA/GT": -1.45, "AT/TA": -0.88, "TT/AA": -1.00,  "GT/CA": -1.44, 
                        "CT/GA": -1.28, "AG/TC": -1.28, "TG/AC": -1.45, "GG/CC": -1.84, "CG/GC": -2.17, "AC/TG": -1.44,  "TC/AG": -1.30,
                        "GC/CG": -2.24,  "CC/GG": -1.84}

# Dictionary describing the contribution of each possible nearest-neighbor pair in a duplex towards the enthalpy in kcal/mol
duplexPairToH = {"AA/TT": -7.9, "TA/AT": -7.2, "GA/CT": -8.2, "CA/GT": -8.5, "AT/TA": -7.2, "TT/AA": -7.9, "GT/CA": -8.4, 
                          "CT/GA": -7.8, "AG/TC": -7.8, "TG/AC": -8.5, "GG/CC": -8.0, "CG/GC": -10.6, "AC/TG": -8.4, "TC/AG": -8.2, 
                          "GC/CG": -9.8, "CC/GG": -8.0}

# Dictionary describing the contribution of each possible nearest-neighbor pair in a duplex towards the entropy in kcal/K mol
duplexPairToS = {"AA/TT": -0.0222, "TA/AT": -0.0213, "GA/CT": -0.0222, "CA/GT": -0.0227, "AT/TA": -0.0204, "TT/AA": -0.0222, 
                          "GT/CA": -0.0224, "CT/GA": -0.0210, "AG/TC": -0.0210, "TG/AC": -0.0227, "GG/CC": -0.0199, "CG/GC": -0.0272, 
                          "AC/TG": -0.0224, "TC/AG": -0.0222, "GC/CG": -0.0244, "CC/GG": -0.0199}

# End correction for a hairpin stem for Gibbs free energy in kcal/mol
stemTermgcAmtG = -2.182
stemTermatAmtG = -1.653

#Initiation paramters for duplex for Gibbs free energy in kcal/mol
duplexTermgcAmtG = 0.98
duplexTermatAmtG = 1.03
symAmtG = 0.43

#Initiation paramaters for duplex for enthalpy in kcal/mol
duplexTermgcAmtH = 0.1
duplexTermatAmtH = 2.3
symAmtH = 0

#Initiation paramaeters for duplex of entropy in kcal/K mol
duplexTermgcAmtS = -0.0028
duplexTermatAmtS = 0.0041
symAmtS = -0.0014

# Mismatch penalty for Gibbs free energy in kcal/mol
mismatchPenalty = 0.438