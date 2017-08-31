COMBINED_SEQS = None
TARGET_LIST = None
REFERENCE_FASTA = None
KEEP_FILES = None


def init(combined, target, reference, keep):
    
    global COMBINED_SEQS
    global TARGET_LIST
    global REFERENCE_FASTA
    global KEEP_FILES
    
    COMBINED_SEQS = combined
    TARGET_LIST = target
    REFERENCE_FASTA = reference
    KEEP_FILES = keep
    
