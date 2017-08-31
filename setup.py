def init(combined, target, reference, keep):
    Constants.assign_values(combined, target, reference, keep)


class Constants:

    combined_seqs = None
    target_list = None
    reference_fasta = None
    keep_files = None

    @classmethod
    def assign_values(cls, combined, target, reference, keep):
        cls.combined_seqs = combined
        cls.target_list = target
        cls.reference_fasta = reference
        cls.keep_files = keep
