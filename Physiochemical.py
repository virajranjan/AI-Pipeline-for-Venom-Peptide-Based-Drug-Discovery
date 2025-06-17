hopp_woods = {
    'A': -0.5, 'R': 3.0, 'N': 0.2, 'D': 3.0, 'C': -1.0,
    'Q': 0.2, 'E': 3.0, 'G': 0.0, 'H': -0.5, 'I': -1.8,
    'L': -1.8, 'K': 3.0, 'M': -1.3, 'F': -2.5, 'P': 0.0,
    'S': 0.3, 'T': -0.4, 'W': -3.4, 'Y': -2.3, 'V': -1.5
}

#Hopp-Woods scales.
def Hydrophobicity(seq):
    score= [hopp_woods.get(aa,0.0) for aa in seq]
    return sum(score)/len(score) if score else 0.0

def Aliphatic_Index(seq):
    lenght  = len(seq)
    if lenght==0:
        return 0.00
    A = seq.count('A')/lenght*100
    V = seq.count('V')/lenght*100
    I = seq.count('I')/lenght*100
    L = seq.count('L')/lenght*100
    return A+2.9*V+3.9*(I+L)
