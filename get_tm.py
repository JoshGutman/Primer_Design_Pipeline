import math


# Driver
def get_tm(primer, oligo_conc, na_conc, mg_conc):

    out = []
    avg = 0


    min_primer = min_primer_conversion(primer)
    max_primer = max_primer_conversion(primer)

    for p in [min_primer, max_primer]:
    
        enthalpy, entropy = get_dH_dS(p)
        uncorrected_tm = get_uncorrected_tm(enthalpy, entropy, oligo_conc)

        if na_conc == 0 and mg_conc == 0:
            out += uncorrected_tm - 273.15
            continue

        if na_conc == 0:
            ratio = 6.0
        else:
            ratio = math.sqrt(mg_conc) / (na_conc/1000.0)


        if ratio < .22:
            tm = monovalent_correction(primer, na_conc, uncorrected_tm)
        else:
            tm = divalent_correction(primer, na_conc, mg_conc, uncorrected_tm, ratio)

        out.append(round(tm, 2))
        avg += tm

    out.append(round(avg / 2, 2))
    return out




# Tm = deltaH / (deltaS - R*ln(Na+))
def get_uncorrected_tm(enthalpy, entropy, oligo_conc):

    term = 1.9865 * math.log(oligo_conc/1000000.0)

    return (enthalpy / (entropy + term))
    


# Get enthalpy and entropy according to nearest-neighbor
#parameters found in Santalucia, 1998
def get_dH_dS(primer):

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
                    


# Owczarzy et al., 2004
def monovalent_correction(primer, na_conc, tm):

    if na_conc == 0:
        na_conc = 150

    gc_content = get_gc_content(primer)
    
    term0 = 1/tm
    term1 = ((4.29 * gc_content) - 3.95) * .00001 * math.log(na_conc / 1000.0)
    term2 = 9.4 * .000001 * (math.log(na_conc / 1000.0)**2.0)

    tm_inverse = term0 + term1 + term2

    return 1/tm_inverse - 273.15



# Owczarzy et al., 2008
def divalent_correction(primer, na_conc, mg_conc, tm, ratio):

    if na_conc == 0:
        na_conc = 150

    if ratio > 6:
        a = 3.92
        d = 1.42
        g = 8.31
    else:
        a = 3.92 * (.843 - (.352 * math.sqrt(na_conc/1000.0) * math.log(na_conc/1000.0)))
        d = 1.42 * (1.279 - (.00403 * math.log(na_conc/1000.0)) - (.00803 * math.log(na_conc/1000.0)**2.0))
        g = 8.31 * (.486 - ((.258 * math.log(na_conc/1000.0)) + (.00525 * math.log(na_conc/1000.0)**3.0)))

    gc_content = get_gc_content(primer)

    term0 = 1/tm
    term1 = a - (.91*math.log(mg_conc/1000.0))
    term2 = gc_content * (6.26 + (d * math.log(mg_conc/1000.0)))
    term3 = (1/(2 * (len(primer)-1))) * (-48.2 + (52.5 * math.log(mg_conc/1000.0)) + (g * math.log(mg_conc/1000.0)**2.0))

    tm_inverse = term0 + ((term1 + term2 + term3) * .00001)

    return 1/tm_inverse - 273.15



# Convert degens to the base with the least energy
def min_primer_conversion(primer):

    # A, T, G, C
    min_degens = {"R":"A", "Y":"T", "S":"G", "W":"A",
                  "K":"T", "M":"A", "B":"T", "D":"A",
                  "H":"A", "V":"A", "N":"A"}

    for degen in min_degens:
        primer = primer.replace(degen, min_degens[degen])

    return primer


# Convert degens to the base with the most energy
def max_primer_conversion(primer):

    # C, G, T, A
    max_degens = {"R":"G", "Y":"C", "S":"C", "W":"T",
                  "K":"G", "M":"C", "B":"C", "D":"G",
                  "H":"C", "V":"C", "N":"C"}

    for degen in max_degens:
        primer = primer.replace(degen, max_degens[degen])

    return primer


if __name__ == "__main__":
    print(get_tm("GAAGCTTCRYTGGTCAGTTC", .2, 0, 3))
