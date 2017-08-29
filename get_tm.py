"""
Calculates the melting temperature of primers.

The driver function, get_tm(), estimates the minimum, maximum, and average
melting temperature in a PCR reaction of a primer.

The terms "oligonucleotide" and "primer" are used interchangeably throughout
this file.

"""
import math


def get_tm(primer, oligo_conc, na_conc, mg_conc):
    """Get the melting temperature of an oligonucleotide.

    Calculate the minimum, maximum, and average melting temperature for a
    given oligonucleotide. It uses the nearest-neighbor parameters described
    in Santalucia, 1998 [1]_, the monovalent salt (Na+) correction described
    in Owczarzy et al., 2004 [2]_, and the divalent salt (Mg++) correction
    described in Owczarzy et al., 2008 [3]_.

    Args:
        primer (str): A string of bases, including degeneracies.
        oligo_conc (float): Olgionucleotide concentration (μM) in PCR solution.
        na_conc (float): Na+ concentration (mM) in PCR solution.
        mg_conc (float): Mg++ concentration (mM) in PCR solution.

    Returns:
        :obj:`list` of :obj:`float`: List of minimum, maximum, and average
        melting temperatures.

    Notes:
        Driver function

        `oligo_conc`, `na_conc`, and `mg_conc` are inputted by the user in
        primer_design_pipeline.py.

    References:
        .. [1] SantaLucia, John. "A unified view of polymer, dumbbell, and
               oligonucleotide DNA nearest-neighbor thermodynamics." Proceedings
               of the National Academy of Sciences 95.4 (1998): 1460-1465.
        .. [2] Owczarzy, Richard, et al. "Effects of sodium ions on DNA duplex
               oligomers: improved predictions of melting temperatures."
               Biochemistry 43.12 (2004): 3537-3554.
        .. [3] Owczarzy, Richard, et al. "Predicting stability of DNA duplexes
               in solutions containing magnesium and monovalent cations."
               Biochemistry 47.19 (2008): 5336-5353.

    Examples:
        Default values for `oligo_conc`, `na_conc`, and `mg_conc`.

        >>> get_tm("GAAGCTTCRYTGGTCAGTTC", .25, 50, 0)
        [51.02, 56.56, 53.79]

    """
    out = []
    avg = 0

    min_primer = _min_primer_conversion(primer)
    max_primer = _max_primer_conversion(primer)

    for p in [min_primer, max_primer]:

        enthalpy, entropy = _get_dH_dS(p)
        uncorrected_tm = _get_uncorrected_tm(enthalpy, entropy, oligo_conc)

        # If no Na+ and Mg++, just return uncorrected_tm
        if na_conc == 0 and mg_conc == 0:
            out += uncorrected_tm - 273.15
            continue

        # Avoid ZeroDivisionError
        if na_conc == 0:
            ratio = 6.0
        else:
            ratio = math.sqrt(mg_conc) / (na_conc/1000.0)

        # Determine which correction function to use
        if ratio < .22:
            tm = _monovalent_correction(p,
                                        na_conc,
                                        uncorrected_tm)
        else:
            tm = _divalent_correction(p,
                                      na_conc,
                                      mg_conc,
                                      uncorrected_tm,
                                      ratio)

        out.append(round(tm, 2))
        avg += tm

    out.append(round(avg / 2, 2))
    return out


def _get_uncorrected_tm(enthalpy, entropy, oligo_conc):
    """Calculate the uncorrected temperature

    Calculate the uncorrected melting temperature of a oligonucleotide
    using nearest-neighbor parameters described in Santalucia, 1998 [1]_.

    Args:
        enthalpy (float): Enthalpy estimation
        entropy(float): Entropy estimation
        oligo_conc (float): Olgionucleotide concentration (μM) in PCR reaction.

    Returns:
        float: Uncorrected melting temperature in Kelvin.

    Notes:
        `enthalpy` and `entropy` are calculated in _get_dH_dS().

        `oligo_conc` is inputted by the user in primer_design_pipeline.py.

        The formula used is:

        .. math:: T_{M}(^o K) = \\frac{\\Delta H^o}{\\Delta S^o  - Rln[oligo]}

        where R is the ideal gas constant:

        .. math:: 1.9865\\;kcal\\;K^{-1}\\;mol^{-1}

    References:
        .. [1] SantaLucia, John. "A unified view of polymer, dumbbell, and
               oligonucleotide DNA nearest-neighbor thermodynamics." Proceedings
               of the National Academy of Sciences 95.4 (1998): 1460-1465.

    """
    term = 1.9865 * math.log(oligo_conc/1000000.0)

    return enthalpy / (entropy + term)


def _get_dH_dS(primer):
    """Calculate enthalpy and entropy.

    Calculate the estimated enthalpy and entropy of an oligonucleotide
    using the nearest-neighbor parameters described in Santalucia, 1998[1]_.

    Args:
        primer (str): A string of bases ("A", "C", "G", or "T").

    Returns:
        :obj:`tuple` of :obj:`float`: Tuple of the enthalpy and entropy.

    References:
        .. [1] SantaLucia, John. "A unified view of polymer, dumbbell, and
               oligonucleotide DNA nearest-neighbor thermodynamics." Proceedings
               of the National Academy of Sciences 95.4 (1998): 1460-1465.

    """
    enthalpy = 0
    entropy = 14   # Primers will always be symmetric

    # Avoid the use of floats for now
    enthalpy_vals = {"AA": 79, "AC": 84, "AG": 78, "AT": 72,
                     "CA": 85, "CC": 80, "CG": 106, "CT": 78,
                     "GA": 82, "GC": 98, "GG": 80, "GT": 84,
                     "TA": 72, "TC": 82, "TG": 85, "TT": 79}

    entropy_vals = {"AA": 222, "AC": 224, "AG": 210, "AT": 204,
                    "CA": 227, "CC": 199, "CG": 272, "CT": 210,
                    "GA": 222, "GC": 244, "GG": 199, "GT": 224,
                    "TA": 213, "TC": 222, "TG": 227, "TT": 222}

    # First and last bases in primer are terminal
    terminal = [primer[0], primer[-1]]
    for t in terminal:
        if t in "AT":
            entropy += -41
            enthalpy += -23
        elif t in "GC":
            entropy += 28
            enthalpy += -1

    # Get enthalpy and entropy for each consecutive pair of bases
    for i in range(len(primer) - 1):
        enthalpy += enthalpy_vals[primer[i:i+2]]
        entropy += entropy_vals[primer[i:i+2]]

    # Convert to cal
    enthalpy *= -100.0
    entropy *= -.1

    return enthalpy, entropy


def _get_gc_content(primer):
    """Get the ratio of Gs and Cs relative to the total amount of bases.

    Args:
        primer (str): A string of bases ("A", "C", "G", or "T")

    Returns:
        float

    """
    return (primer.count("G") + primer.count("C")) / len(primer)


def _monovalent_correction(primer, na_conc, uncorrected_tm):
    """Perform monovalent salt correction.

    Adjust the temperature using the monovalent salt (Na+) correction
    described in Owczarzy et al., 2004[1]_.

    Args:
        primer (str): A string of bases ("A", "C", "G", or "T").
        na_conc (float): Na+ concentration (mM) in PCR solution.
        uncorrected_tm (float): Temperature calculated using nearest-neighbor
            parameters described in Santalucia, 1998 [2]_.

    Returns:
        float: Corrected melting temperature in degrees Celcius.

    See Also:
        _divalent_correction : Correction for Mg++ concentration.

    Notes:
        `na_conc` is inputted by the user in primer_design_pipeline.py.

        `uncorrected_tm` is calculated in _get_uncorrected_tm().

        The formula is:

        .. math:: \\frac{1}{T_{M}(Na^+)} = \\frac{1}{T_{M}(1M\\;Na^+)} + \\Big[(4.29f_{GC}\\,-\\,3.95)\\;ln[Na^+]\\;+\\;0.940\\,ln^2[Na^+]\\Big]

    References:
        .. [1] Owczarzy, Richard, et al. "Effects of sodium ions on DNA duplex
               oligomers: improved predictions of melting temperatures."
               Biochemistry 43.12 (2004): 3537-3554.
        .. [2] SantaLucia, John. "A unified view of polymer, dumbbell, and
               oligonucleotide DNA nearest-neighbor thermodynamics." Proceedings
               of the National Academy of Sciences 95.4 (1998): 1460-1465.

    """
    # Avoid ZeroDivisionError
    if na_conc == 0:
        na_conc = 150

    gc_content = _get_gc_content(primer)

    term0 = 1/uncorrected_tm
    term1 = ((4.29 * gc_content) - 3.95) * .00001 * math.log(na_conc / 1000.0)
    term2 = 9.4 * .000001 * (math.log(na_conc / 1000.0)**2.0)

    tm_inverse = term0 + term1 + term2

    return 1/tm_inverse - 273.15


def _divalent_correction(primer, na_conc, mg_conc, uncorrected_tm, ratio):
    """Perform divalent salt correction.

    Adjust the temperature using the divalent salt (Mg++) correction
    described in Owczarzy et al., 2008[1]_.

    Parameters:
        primer (str): A string of bases ("A", "C", "G", or "T").
        na_conc (float): Na+ concentration (mM) in PCR solution.
        mg_conc (float): Mg++ concentration (mM) in PCR solution.
        uncorrected_tm (float): Temperature calculated using nearest-neighbor
            parameters described in Santalucia, 1998[2]_.
        ratio (float): Special ratio of Mg++ ions to Na+ ions.

    Returns:
        float: Corrected melting temperature in degrees Celcius.

    See Also:
        _monovalent_correction : Correction for Na+ concentration.

    Notes:
        `na_conc` and `mg_conc` are inputted by the user in
        primer_design_pipeline.py.

        `uncorrected_tm` is calculated in _get_uncorrected_tm().

        `ratio` uses the following formula:

        .. math:: \\frac{\\sqrt{[Mg^{++}]}}{[Na^+]}

        The formula is:

        .. math:: \\frac{1}{T_{M}(Mg^{++})} = \\frac{1}{T_{M}(1M\\;Na^+)} + \\Bigg[ a - 0.911\\,ln[Mg^{++} ] + f_{GC} \\times \\big(6.26 + d\\,ln[MG^{++}]\\big) + \\frac{1}{2(N_{bp}-1)} \\times \\big[-48.2 + 52.5\\,ln[Mg^{++}] + g\\,ln^2[Mg^{++}]\\big]\\Bigg] \\times 10^{-5}

        a, d, and g can either be constants or vary with the `na_conc` based
        on `ratio`.

    References:
        .. [1] Owczarzy, Richard, et al. "Predicting stability of DNA duplexes
               in solutions containing magnesium and monovalent cations."
               Biochemistry 47.19 (2008): 5336-5353.
        .. [2] SantaLucia, John. "A unified view of polymer, dumbbell, and
               oligonucleotide DNA nearest-neighbor thermodynamics." Proceedings of
               the National Academy of Sciences 95.4 (1998): 1460-1465.

    """
    # Avoid ZeroDivisionError
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

    gc_content = _get_gc_content(primer)

    term0 = 1/uncorrected_tm
    term1 = a - (.911*math.log(mg_conc/1000.0))
    term2 = gc_content * (6.26 + (d * math.log(mg_conc/1000.0)))
    term3 = (1/(2*(len(primer)-1))) * (-48.2 + (52.5 * math.log(mg_conc/1000.0)) + (g * math.log(mg_conc/1000.0)**2.0))

    tm_inverse = term0 + ((term1 + term2 + term3) * .00001)

    return 1/tm_inverse - 273.15


def _min_primer_conversion(primer):
    """Convert primer to sequence with least energy.

    Convert all degens within primer to the bases that have the least
    energy. This allows for the minimum possible temperature to be
    calculated.

    Args:
        primer (str): A string of bases, including degeneracies.

    Returns:
        str: A string of bases, with degeneracies converted to base with
        minimum energy.

    See Also:
        _max_primer_conversion : Convert primer to sequence with most energy.

    Notes:
        From least energy to most energy, the bases are A, T, G, C.

    """
    min_degens = {"R": "A", "Y": "T", "S": "G", "W": "A",
                  "K": "T", "M": "A", "B": "T", "D": "A",
                  "H": "A", "V": "A", "N": "A"}

    for degen in min_degens:
        primer = primer.replace(degen, min_degens[degen])

    return primer


def _max_primer_conversion(primer):
    """Convert primer to sequence with most energy.

    Convert all degens within primer to the bases that have the most
    energy. This allows for the maximum possible temperature to be
    calculated.

    Args:
        primer (str): A string of bases, including degeneracies.

    Returns:
        str: A string of bases, with degeneracies converted to base with
        maximum energy.

    See Also:
        _min_primer_conversion : Convert primer to sequence with most energy.

    Notes:
        From most energy to least energy, the bases are C, G, T, A.

    """
    max_degens = {"R": "G", "Y": "C", "S": "C", "W": "T",
                  "K": "G", "M": "C", "B": "C", "D": "G",
                  "H": "C", "V": "C", "N": "C"}

    for degen in max_degens:
        primer = primer.replace(degen, max_degens[degen])

    return primer


if __name__ == "__main__":
    print(get_tm("GAAGCTTCRYTGGTCAGTTC", .25, 50, 0))
