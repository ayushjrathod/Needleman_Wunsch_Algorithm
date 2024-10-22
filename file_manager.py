import os


def read_files(input_dir):
    sequences = []
    for filename in os.listdir(input_dir):
        if filename.endswith('.seq'):
            with open(os.path.join(input_dir, filename), 'r') as file:
                lines = file.read().splitlines()
                if len(lines) >= 2:
                    sequences.append((lines[0], lines[1], filename))
    return sequences
