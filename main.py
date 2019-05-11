import numpy as np
import pandas as pd
from sig_modules.deconstructsigs_old import DeconstructSigs_old


def main():
    vector_list = [28, 12, 12, 21, 10, 4, 1, 13, 39, 27, 128, 23, 14, 11, 21,
                   16, 54, 20, 38, 40, 6, 4, 7, 7, 21, 22
        , 5, 18, 4, 8, 4, 13, 37, 28, 66, 50, 5, 5, 13, 9, 16, 8, 17, 13, 6,
                   11, 8, 7, 26, 11, 5, 15
        , 5, 7, 2, 6, 31, 35, 62, 32, 6, 1, 9, 7, 21, 4, 14, 16, 3, 3, 2, 7,
                   27, 18, 1, 30, 13, 11
        , 2, 22, 25, 39, 49, 33, 15, 3, 7, 18, 16, 18, 13, 19, 6, 8, 6, 21]
    colon1_list = [3, 0, 2, 0, 0, 1, 0, 0, 5, 4, 22, 4, 1, 1, 0, 1, 6, 0, 4, 6,
                   1, 0, 1, 2, 1, 2, 1, 1, 1, 0, 0, 1, 7, 2, 10, 2, 1, 0, 0, 0,
                   0, 1, 3, 2, 1
        , 1, 2, 1, 4, 2, 1, 2, 0, 0, 0, 0, 6, 6, 18, 4, 1, 0, 1, 1, 0, 0, 1, 2,
                   1, 0, 0, 0, 3, 3, 0, 3, 1, 0, 0, 1, 4, 6, 9, 2, 1, 1, 2, 2,
                   2, 2
        , 1, 1, 0, 0, 1, 1]
    liver2_list = [10, 8, 6, 12, 10, 3, 1, 9, 19, 17, 16, 9, 8, 7, 19, 6, 34,
                   14, 29, 26, 5, 4, 4, 5, 15, 16, 4, 13
        , 3, 7, 3, 10, 19, 19, 17, 36, 4, 3, 8, 5, 11, 5, 12, 6, 2, 10, 3, 5,
                   12, 7, 3, 7, 3, 6, 1, 4
        , 16, 15, 8, 19, 4, 1, 7, 4, 12, 4, 8, 8, 1, 2, 0, 3, 14, 10, 0, 9, 7,
                   10, 2, 14, 18, 19, 7, 22
        , 11, 1, 4, 11, 11, 9, 9, 11, 5, 7, 3, 14]
    vector = np.array(liver2_list)
    vector = vector / sum(vector)
    sig_df = pd.read_csv(
        '/home/egors/deconstruct_sigs/data/signatures_probabilities.txt',
        sep='\t', usecols=np.arange(3, 40))

    # ds = DeconstructSigs2(m96vector = vector,output_folder=
    # '/home/egors/Deconstruct_digs/deconstructSigs/outputs',analysis_handle = 'liver2',use_cosmic=True,signatures_df=sig_df)

    ds = DeconstructSigs_old(
        maf_file_path='/home/egors/deconstruct_sigs/sample_data/test.maf',
        skip_rows=5)

    weights = ds.which_signatures(verbose=True)
    # ds.plot_signatures(weights, explanations=True)
    ds.figures(weights, explanations=True)
    print()


if __name__ == '__main__':
    main()
