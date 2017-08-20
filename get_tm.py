import math


# Driver
def get_tm(primer, primer_conc, na_conc, mg_conc):
    
    enthalpy, entropy = get_dH_dS(primer, primer_conc, na_conc, mg_conc)
    uncorrected_tm = get_uncorrected_tm(enthalpy, entropy, primer_conc)
    tm = monovalent_correction(primer, na_conc, uncorrected_tm)

    return tm



def get_uncorrected_tm(enthalpy, entropy, primer_conc):

    term = 1.9865 * math.log(primer_conc/1000000.0)

    return (enthalpy / (entropy + term))
    


def get_dH_dS(primer, primer_conc, na_conc, mg_conc):


    enthalpy = 0
    entropy = 14
    

    enthalpy_vals = {"AA": 79, "AC": 84, "AG": 78, "AT": 72,
                     "CA": 85, "CC": 80, "CG":106, "CT": 78,
                     "GA": 82, "GC": 98, "GG": 80, "GT": 84,
                     "TA": 72, "TC": 82, "TG": 85, "TT": 79}


    entropy_vals = {"AA": 222, "AC": 224, "AG": 210, "AT": 204,
                    "CA": 227, "CC": 199, "CG": 272, "CT": 210,
                    "GA": 222, "GC": 244, "GG": 199, "GT": 224,
                    "TA": 213, "TC": 222, "TG": 227, "TT": 222}

        


    terminal = [primer[0], primer[-1]]

    for t in terminal:
        if t in "AT":
            entropy += -41
            enthalpy += -23
        elif t in "GC":
            entropy += 28
            enthalpy += -1
    


    for i in range(len(primer) - 1):
        enthalpy += enthalpy_vals[primer[i:i+2]]
        entropy += entropy_vals[primer[i:i+2]]


    enthalpy *= -100.0
    entropy *= -.1


    return enthalpy, entropy




def get_gc_content(primer):
    return (primer.count("G") + primer.count("C")) / len(primer)
                    


def monovalent_correction(primer, na_conc, temp):

    gc_content = get_gc_content(primer)
    
    term0 = 1/temp
    term1 = ((4.29 * gc_content) - 3.95) * .00001 * math.log(na_conc / 1000.0)
    term2 = 9.4 * .000001 * (math.log(na_conc / 1000.0)**2.0)

    tm_inverse = term0 + term1 + term2

    return 1/tm_inverse - 273.15





if __name__ == "__main__":
    print(get_tm("CTCTATCTAGCTCTCT", .25, 50, 0))
