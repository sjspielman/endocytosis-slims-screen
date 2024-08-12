#!/usr/bin/env python3

"""
Find all locations (possibly overlapping) for the `[LI]XXQXTG` motif in the
disordered regions of the human proteome
"""

import re
import sys

BOUND_BUFFER = 2 # consider +/-2 residues beyond the annotated disordered region as part of the disordered region
MIN_SIZE = 5     # only consider disordered regions (including the ^ buffer)that are >=5 residues

def main():

    input_json = "../reference/disorder_UP000005640.mjson"
    output_file = "temp.csv"

    # define regex for LI]\\w\\wQ\\wTG and allow overlapping matches
    compiled_motif = re.compile( "(?=([LI]\\w\\wQ\\wTG))" )

    null="NA" ## hack for eval() use below

    with open(input_json, "r") as f:
        input_lines = f.readlines()

    full_results = {}
    for line in input_lines:
        disorder_record = eval(line.strip())

        uniprot = disorder_record["acc"]
        sequence = disorder_record["sequence"]

        try:
            regions_raw = disorder_record["mobidb_consensus"]["disorder"]["predictors"][0]
        except:
            print(uniprot, "has no predictions. moving on.")
            continue

        if regions_raw["method"] != "simple":
            print(uniprot, "has no simple predictions. moving on.")
            continue
        regions = regions_raw["regions"] ## list of disordered regions, i.e. [[1, 24, 'D'], [40, 47, 'D'], [71, 85, 'D'], [87, 121, 'D'], [154, 157, 'D'], [167, 202, 'D'], [214, 221, 'D'], [251, 310, 'D']]


        detected_motifs = []
        for bounds in regions:
            start = bounds[0] - 1 + BOUND_BUFFER
            stop = bounds[1] - BOUND_BUFFER

            ### skip disordered regions that are too small
            if (stop < start) | (bounds[1] - bounds[0] < MIN_SIZE):
                continue
            for detected_motif in compiled_motif.finditer(sequence[start:stop]):
                ## save tuple per motif:
                # - starting index (from 0!) and the actual motif detected.
                # - the detected sequence
                this_slim = [str(detected_motif.start() + start +  1), detected_motif.group(1)] ## list of 2 items, starting index (from 0!) and the actual motif detected. the index start position gets bumped to HUMAN INTERPRETABLE (starts from 0) and as an index for the WHOLE sequence, not just for the substring we searched.
                detected_motifs.append(this_slim)
        if len(detected_motifs) > 0:
            full_results[uniprot] = detected_motifs

    # save results
    with open(output_file, "w") as f:
        f.write("uniprot_id,start_position,motif_sequence\n")

        for uniprot in full_results:
            for motif in full_results[uniprot]:
                row = ",".join([uniprot, motif[0], motif[1]]) + "\n"
                f.write(row)



if __name__ == "__main__":
    main()
