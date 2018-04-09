import argparse
import pickle
import sys
import os


def get_probes(combo_pickle, combo_ids, outfile_name):
    ids = [int(i) for i in combo_ids.split(",")]
    all_combos = unpickle_combos(combo_pickle)

    combos = []
    for combo_id in ids:
        combos.append(get_combo_from_id(all_combos, combo_id))

    candidate_lists = []
    for combo in combos:

        # First, look for perfect probes of length 25
        candidates = []
        median_probes = get_probes_at_size(combo, 25)
        for probe in median_probes:
            if probe.combined_score == 0:
                candidates.append(probe)

        if len(candidates) > 0:
            candidate_lists.append(CandidateList(combo, candidates))
            continue

        # There are no perfect probes of length 25 -- look in range 20-30
        perfects = []
        imperfects = []
        for i in range(20, 31):
            if i == 25:
                continue
            
            probes = get_probes_at_size(combo, i)
            for probe in probes:
                if probe.combined_score == 0:
                    perfects.append(probe)
                else:
                    imperfects.append(probe)

        if len(perfects) > 0:
            candidate_lists.append(CandidateList(combo, perfects))
        else:
            score_probes(imperfects)
            candidate_lists.append(CandidateList(combo, imperfects))


    output_probes(candidate_lists, outfile_name)
                

        
        
    


def unpickle_combos(file):
    with open(file, "rb") as infile:
        combos = pickle.load(infile)
    return combos


def get_combo_from_id(combos, combo_id):
    for combo in combos:
        if combo.id == combo_id:
            return combo


def get_probes_at_size(combo, size):
    probes = []
    probe_size = min(size, len(combo.amplicon))

    for i in range(len(combo.amplicon)-probe_size):
        probe = Probe(combo.amplicon[i:i+probe_size], i, combo)
        
        if probe.seq[0] == "G":
            continue
        
        probe.add_distance_score()
        probe.add_gc_score()
        probe.add_tm_score()
        probe.add_combined_score()
        
        probes.append(probe)

    return probes


def score_probes(probes):
    
    def _create_reverse_dict(scores, instance_var):
        out = {}
        for score in scores:
            if score in out:
                for probe in probes:
                    if getattr(probe, instance_var) == score:
                        out[score].append(probe)
            else:
                out[score] = []
        return out

    def _add_score(scores):
        for idx, score in enumerate(sorted(list(scores.keys()), reverse=True)):
            for probe in scores[score]:
                probe.score += idx

    distance_scores = _create_reverse_dict(
        [probe.distance_score for probe in probes], "distance_score")

    gc_scores = _create_reverse_dict(
        [probe.gc_score for probe in probes], "gc_score")

    tm_scores = _create_reverse_dict(
        [probe.tm_score for probe in probes], "tm_score")

    _add_score(distance_scores)
    _add_score(gc_scores)
    _add_score(tm_scores)


def output_probes(candidate_lists, output_file):
    with open(output_file, "w") as outfile:

        fields = ["seq", "size", "start_idx", "distance_score", 
                  "gc_score", "tm_score", "score"]
        
        for field in fields:
            outfile.write(field + "\t")
        outfile.write("\n")

        for candidate_list in candidate_lists:
            
            outfile.write("\n\n{}.   {} - {}\n".format(
                candidate_list.combo.id,
                candidate_list.combo.forward.name,
                candidate_list.combo.reverse.name))
            outfile.write("-" *100 + "\n")
            
            for i in range(min(3, len(candidate_list.candidates))):
                for field in fields:
                    to_write = getattr(candidate_list.candidates[i], field)
                    if type(to_write) is float:
                        to_write = format(to_write, ".2f")
                    outfile.write("{}\t".format(to_write))
                if candidate_list.candidates[i].combined_score == 0:
                    outfile.write("[PERFECT]")
                outfile.write("\n")
    



class Probe:

    def __init__(self, seq, start_idx, combo):
        self.score = 0
        self.seq = seq
        self.size = len(seq)
        self.start_idx = start_idx
        self.combo = combo


    def add_distance_score(self):
        self.distance_score = min(
            self.start_idx,
            len(self.combo.amplicon) - (self.start_idx + self.size))

    def add_gc_score(self):
        self.gc_content = (self.seq.count("G") + self.seq.count("C"))/self.size
        if self.gc_content > .65:
            self.gc_score = (self.gc_content - .65) * 100
        elif self.gc_content < .35:
            self.gc_score = (.35 - self.gc_content) * 100
        else:
            self.gc_score = 0


    def add_tm_score(self):
        self.melting_tm = get_tm.get_tm(self.seq, self.combo.oligo_conc,
                                        self.combo.na_conc, self.combo.mg_conc)

        self.tm_diff = ((self.melting_tm[1] - self.combo.forward.tm[1]) +
                        (self.melting_tm[1] - self.combo.reverse.tm[1]) / 2)
        
        if 6 <= self.tm_diff <= 8:
            self.tm_score = 0
        elif self.tm_diff < 6:
            self.tm_score = (6 - abs(self.tm_diff))
        else:
            self.tm_score = (abs(self.tm_diff) - 8)


    def add_combined_score(self):
        self.combined_score = (self.distance_score +
                               self.gc_score + self.tm_score)


class CandidateList:

    def __init__(self, combo, candidates):
        self.combo = combo
        self.candidates = sorted(candidates, key=lambda probe: probe.score, reverse=True)
    


if __name__ == "__main__":
    # Python paths are annoying
    if __package__ is None:
        sys.path.append(os.path.abspath(os.path.join(
            os.path.dirname(__file__), '..')))
    from Primer_Design_Pipeline import combo_funcs
    from Primer_Design_Pipeline import get_tm

    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", help="Output file name (default probes.txt)")
    parser.add_argument("-p", "--primers", help="[REQUIRED] Comma-delimited IDs of primers to output", required=True)
    parser.add_argument("-i", "--input", help="Path to input .pickle file", default="combos.pickle")

    args = parser.parse_args()

    get_probes(args.input, args.primers, args.output)

    
    
