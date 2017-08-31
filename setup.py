def init(combined, config, target, reference, keep):
    Constants.assign_values(combined, config, target, reference, keep)


class Constants:

    combined_seqs = None
    config_file = None
    target_list = None
    reference_fasta = None
    keep_files = None

    @classmethod
    def assign_values(cls, combined, config, target, reference, keep):
        cls.combined_seqs = combined
        cls.config_file = config
        cls.target_list = target
        cls.reference_fasta = reference
        cls.keep_files = keep
