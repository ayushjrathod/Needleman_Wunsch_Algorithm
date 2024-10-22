from alignment import needleman_wunsch, trace_back
from file_manager import read_files

input_dir = './INPUTS/'

sequences = read_files(input_dir)

for seq1, seq2, filename in sequences:
    S, Ix, Iy, final_score = needleman_wunsch(seq1, seq2, match_score=7, mismatch_score=-5, gap_open=-3, gap_extend=-1)

    alignment_a, alignment_b = trace_back(S, seq1, seq2, match_score=7, mismatch_score=-5, gap_extend=-1)

    print(f"File: {filename}")
    print("Alignment Score:", final_score/len(alignment_a))

    print("Sequence 1:", alignment_a)
    print("Sequence 2:", alignment_b)
    print("")
