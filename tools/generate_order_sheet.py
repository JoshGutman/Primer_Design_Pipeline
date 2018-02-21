import argparse
import pickle
import csv
import sys
import io
import os


def combos_to_csv(outfile, infile, id_str, ut):
    ids = [int(i) for i in id_str.split(",")]
    combos = unpickle_combos(infile)
    for combo_id in ids:
        combo = get_combo_from_id(combos, combo_id)
        if combo is not None:
            order_info = combo.get_order_info(ut)
            amp_info = combo.get_amplicon_info()
            output_to_csv(outfile, order_info, amp_info, ut)


def unpickle_combos(file):
    with open(file, "rb") as infile:
        combos = pickle.load(infile)
    return combos


def get_combo_from_id(combos, combo_id):
    for combo in combos:
        if combo.id == combo_id:
            return combo


def output_to_csv(file, order_info, amplicon_info, ut):
    rows = read_csv(file)
    
    # If file is empty
    if not rows:
        if ut:
            columns = ["Target", "Primer", "Combined_name",
                       "Primer (5'-3')", "Final_name", "UT + Sequence", "To order",
                       "Tm","Amplicon", "Amplicon length", "Amplicon Length + UT"]
        else:
            columns = ["Target", "Primer", "Combined_name",
                       "Primer (5'-3')", "Final_name", "To order",
                       "Tm","Amplicon", "Amplicon length"]
            
        empty_row = ["" for item in columns]
        rows.append(columns)
        rows.append(empty_row)
    
    index = find_empty_row(rows)
    rows.insert(index, order_info[0])
    rows.insert(index+1, order_info[1])
    rows.append(amplicon_info)
    write_csv(file, rows)


def write_csv(file, data):
    with open(file, "w", newline="") as outfile:
        writer = csv.writer(outfile)
        for row in data:
            writer.writerow(row)


def read_csv(file):
    try:
        with open(file, newline="") as infile:
            reader = csv.reader(infile)
            return list(reader)
    except FileNotFoundError:
        return []


def find_empty_row(data):
    for index, row in enumerate(data):
        for item in row:
            if item != "":
                continue
            return index
    return len(data)


def csv_to_str(data):
    out = io.StringIO()
    writer = csv.writer(out)
    writer.writerow(data)
    return out.getvalue().strip().split(",")


if __name__ == "__main__":

    
    # Python paths are annoying
    if __package__ is None:
        sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
    from Primer_Design_Pipeline import combo_funcs
    
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", help="[REQUIRED] Path to output .csv file", required=True)
    parser.add_argument("-p", "--primers", help="[REQUIRED] Comma-delimited IDs of primers to output", required=True)
    parser.add_argument("-i", "--input", help="Path to input .pickle file", default="combos.pickle")
    parser.add_argument("-ut", "--universal_tail", help="Output Universal Tail information (default False)", const=True, nargs="?", default=False)

    args = parser.parse_args()
    combos_to_csv(args.output, args.input, args.primers, args.universal_tail)
