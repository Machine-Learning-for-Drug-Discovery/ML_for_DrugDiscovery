import os, csv
import numpy as np

def get3Let(residue):
    
    resDict = {
        "A" : "Ala",
        "C" : "Cys",
        "D" : "Asp",
        "E" : "Glu",
        "F" : "Phe",
        "G" : "Gly",
        "H" : "His",
        "I" : "Ile",
        "K" : "Lys",
        "L" : "Leu",
        "M" : "Met",
        "N" : "Asn",
        "P" : "Pro",
        "Q" : "Gln",
        "R" : "Arg",
        "S" : "Ser",
        "T" : "Thr",
        "V" : "Val",
        "W" : "Trp",
        "Y" : "Tyr"        
    }
    
    return resDict[residue]

def calcMeanHydroph(seq):
    
    # Hydrophobicity values (Adapted from Fauchère, J., and Pliska, V. 1983. Hydrophobic parameters {pi} of amino-acid side chains from the partitioning of N-acetyl-amino-acid amides. Eur. J. Med. Chem. 8: 369–375)
    hDict = {
        "ALA" : 0.310,
        "ASP" : -0.770,
        "GLU" : -0.640,
        "ILE" : 1.800,
        "MET" : 1.230,
        "SER" : -0.040,
        "TYR" : 0.960,
        "ARG" : -1.010,
        "CYS" : 1.540,
        "GLY" : 0.000,
        "LEU" : 1.700,
        "PHE" : 1.790,
        "THR" : 0.260,
        "VAL" : 1.220,
        "ASN" : -0.600,
        "GLN" : -0.220,
        "HIS" : 0.130,
        "LYS" : -0.990,
        "PRO" : 0.720,
        "TRP" : 2.250           
    }
    
    try:
        netHydrp = 0
        for res in seq:

            netHydrp += hDict[ get3Let(res).upper() ]
    except:
        print("Unkown residue(s) found in sequence:", seq)
    return netHydrp/len(seq)

def calcCharge(seq):
    netCharge = 0
    netCharge += (seq.count("E") + seq.count("D")) *-1
    netCharge += (seq.count("K") + seq.count("R")) * 1
    return netCharge

class PepEncode:
    """
    Peptide Encoding Class
    
    This class is responsible for creating a representation for 
    amino acid sequences using physicochemical properties.
    """
    
    def __init__(self, pathToScale="./"):
        
        self.scale          = "fpScales.csv"
        self.pathToScale   = pathToScale

        self.abbr3L   = []
        self.abbr1L   = []
        self._embed    = []
        
        self._loadScale()
    
    def _loadScale(self):
        
        # Loads amino acid scale
        with open(os.path.join(self.pathToScale, self.scale)) as csvfile:
            spamreader = csv.reader(csvfile, delimiter=',')
            for indx,row in enumerate(spamreader):
                if indx == 0:
                    # Skip headers
                    continue
                
                self.abbr3L.append(row[0])
                self.abbr1L.append(row[1])
                self._embed.append([float(x) for x in row[2:]])
        
        self.abbr3L = np.asarray(self.abbr3L)
        self.abbr1L = np.asarray(self.abbr1L)
        self._embed  = np.asarray(self._embed)
    
    def encode(self, seq):
        
        locSeq = seq.strip()
        retSeq = []
        
        for res in locSeq:
            try:
                indx = np.where(self.abbr1L == res)[0][0]
                retSeq.append( self._embed[indx] )
            except:
                print("Error! Residue {} not recognized.".format(res))
                print("Could not encode sequence: ->" + seq + "<- ")
                return
        
        return(np.asarray(retSeq).T)
