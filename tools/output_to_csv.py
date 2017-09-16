import pickle
import csv
import io
import argparse

from combo_funcs import Combo

COLUMNS = ["Type of Target", "Target", "Primer", "Combined_name",
           "Primer (5'-3')", "Final_name", "UT + Sequence", "To order", "Tm",
           "Amplicon length", "Amplicon Length + UT", "SNP position",
           "SNP call", "2nd SNP call", "OutGroup SNP call"]

EMPTY_ROW = ["", "", "", "", "", "", "", "", "", "", "", "", "", "", ""]

PICKLED_COMBOS = "combos.pickle"


def combos_to_csv(file, id_str):
    ids = [int(i) for i in id_str.split(",")]
    combos = unpickle_combos()
    for combo_id in ids:
        combo = get_combo_from_id(combos, combo_id)
        if combo is not None:
            output_to_csv(file, combo.get_order_info, combo.get_amplicon_info)


def unpickle_combos():
    global PICKLED_COMBOS
    with open(PICKLED_COMBOS, "rb") as infile:
        combos = pickle.load(infile)
    return combos


def get_combo_from_id(combos, combo_id):
    for combo in combos:
        if combo.id == combo_id:
            return combo


def output_to_csv(file, order_info, amplicon_info):
    rows = read_csv(file)

    # If file is empty
    if len(rows[0]) == 0:
        global COLUMNS
        global EMPTY_ROW
        write_csv(file, COLUMNS)
        write_csv(file, EMPTY_ROW)
    
    index = find_empty_row(rows)
    rows.insert(index, order_info)
    rows.append(amplicon_info)
    write_csv(file, rows)


def write_csv(file, data):
    with open(file, "w", newline="") as outfile:
        writer = csv.writer(outfile)
        for row in data:
            writer.writerow(row)


def read_csv(file):
    with open(file, newline="") as infile:
        reader = csv.reader(infile)
        return list(reader)


def find_empty_row(data):
    global EMPTY_ROW
    for index, row in enumerate(data):
        if row == EMPTY_ROW:
            break
    return index


def csv_to_str(data):
    out = io.StringIO()
    writer = csv.writer(out)
    writer.writerow(data)
    return out.getvalue().strip().split(",")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", help="[REQUIRED] Path to output .csv file", required=True)
    parser.add_argument("-i", "--ids", help="[REQUIRED] Comma-delimited IDs of primers to output", required=True)

    args = parser.parse_args()
    #combos_to_csv(
