import pandas as pd
import os
import sys


REMOVE_COLS = ('unshared', 'voteA', 'voteB', 'congD')


def remove_cols_from_file(fn, out_fn=None, cols=REMOVE_COLS):
    if not out_fn:
        out_fn, out_ext = os.path.splitext(fn)
        out_fn = out_fn + '_stripped' + out_ext
    with open(fn) as f:
        header_lines = f.readlines()[:2]
    df = pd.read_csv(fn, skiprows=2, sep='\t')
    df = df.set_index('Unnamed: 0')
    df.index.name = ' '
    for col in cols:
        if col in df:
            del df[col]
    with open(out_fn, 'w') as f:
        for line in header_lines:
            f.write(line)
        df.to_csv(f, sep='\t')


if __name__ == "__main__":
    fn = sys.argv[1]
    out_fn = None
    if len(sys.argv) > 2:
        out_fn = sys.argv[2]
    remove_cols_from_file(fn, out_fn, REMOVE_COLS)
