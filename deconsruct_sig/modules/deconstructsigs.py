import datetime as dt
import os
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.style.use('ggplot')
from collections import defaultdict
import math
from scipy.optimize import minimize_scalar
from matplotlib.font_manager import FontProperties

courier_font = FontProperties(family='courier new', weight='bold')


class DeconstructSigs:
    """A Python implementation of the DeconstructSigs2 algorithm described in
    https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0893-4. Modeled after the R implementation
    coded by Rachel Rosenthal which can be found at https://github.com/raerose01/deconstructSigs.

    From the GenomeBiology description:
        The deconstructSigs approach determines the linear combination of pre-defined signatures that most accurately
        reconstructs the mutational profile of a single tumor sample. It uses a multiple linear regression model with
        the caveat that any coefficient must be greater than 0, as negative contributions make no biological sense. """

    # base pairs complement mapping dict
    pair = {
        'A': 'T',
        'C': 'G',
        'T': 'A',
        'G': 'C'
    }

    # pyrimidine bases
    pyrimidines = ['C', 'T']

    # purine bases
    purines = ['G', 'A']

    def __init__(self, m96vector, use_cosmic=True, signatures_df=None,
                 annotation_df=None, cutoff=0.06, analysis_handle='testing',
                 verbose=False, threshold=1e-3):
        """
        Initialize a DeconstructSigs2 object.

        """
        assert  type(m96vector) ==  type(np.zeros(0))
        assert len(m96vector) == 96
        self.m96_vector = m96vector


        self.threshold = threshold

        self.signature_cutoff = cutoff
        self.analysis_handle = analysis_handle

        package_path = os.path.join(
            os.path.dirname(os.path.realpath(__file__)), os.path.pardir)
        self.cosmic_signatures_filepath = os.path.join(package_path,
                                                       'data/signatures_probabilities.txt')
        self.cosmic_signature_explanations_filepath = os.path.join(
            package_path, 'data/about_cosmic_sigs.txt')

        self.timestamp = dt.datetime.now().strftime("%Y_%m_%d-%H_%M_%S")
        self.verbose = verbose

        self.__setup_subs_dict()

        if use_cosmic:
            self.__load_cosmic_signatures()
            self.__load_cosmic_signature_explanations()
        else:
            if isinstance(signatures_df, pd.DataFrame):
                if signatures_df.isna().any().any():
                    raise ValueError('bad sig dataframe ')
            self.pre_defined_signatures = signatures_df
            if annotation_df:
                self.pre_defined_signature_explanations = annotation_df
            else:
                temdf = pd.DataFrame()
                signames = [col for col in signatures_df.columns.values
                            if 'Signature' in col]
                temdf['Signature'] = signames
                temdf['Association'] = 'no info'
                self.pre_defined_signature_explanations = temdf
        self.final_result_weights = None

        # Remove unnecessary columns from the cosmic signatures data and make the S matrix. Note: the substitution
        # contexts are in alphabetical order (A[C>A]A, A[C>A]C, A[C>A]G, A[C>A]T, A[C>G]A, A[C>G]C... etc.)
        self.S = np.array(self.pre_defined_signatures.select(
            lambda x: not re.search(
                "(Substitution Type)|"
                "(Trinucleotide)|"
                "(Somatic Mutation Type)"
                "|(Context)|(SUBS)|(Unnamed)", x), axis=1))
        self.signature_names = [signame for signame in
                                self.pre_defined_signature_explanations[
                                    'Signature'].values]
        # self.signature_names = ['Signature 1', 'Signature 2', 'Signature 3',
        #                         'Signature 4', 'Signature 5',
        #                         'Signature 6', 'Signature 7', 'Signature 8',
        #                         'Signature 9', 'Signature 10',
        #                         'Signature 11', 'Signature 12', 'Signature 13',
        #                         'Signature 14', 'Signature 15',
        #                         'Signature 16', 'Signature 17', 'Signature 18',
        #                         'Signature 19', 'Signature 20',
        #                         'Signature 21', 'Signature 22', 'Signature 23',
        #                         'Signature 24', 'Signature 25',
        #                         'Signature 26', 'Signature 27', 'Signature 28',
        #                         'Signature 29', 'Signature 30']

    def which_signatures(self, signatures_limit=None, associated=None,
                         verbose=False):
        """Wrapper on __which_signatures function. Return weights """
        w = self.__which_signatures(signatures_limit=signatures_limit,
                                    associated=associated)
        self.final_result_weights = w
        return w

    def __which_signatures(self, signatures_limit=None, associated=None):
        """Get the weights transformation vector. If a vector
        of associated indices is provided, only consider the weights at the indicated indices."""
        # If no signature limit is provided, simply set it to the number of signatures
        if signatures_limit is None:
            signatures_limit = len(self.S[0])

        # Зашлушка под искулючаемые индексы на будущее
        T = self.m96_vector
        ignorable_indices = []
        ignorable_signatures = []
        iteration = 0

        w = self.__seed_weights(T, self.S, ignorable_indices=ignorable_indices)

        error_diff = math.inf
        error_threshold = self.threshold
        while error_diff > error_threshold:
            iteration = iteration + 1
            error_pre = self.__get_error(T, self.S, w)
            if error_pre == 0:
                break
            self.__status(
                "Iter {}:\n\tPre error: {}".format(iteration, error_pre))
            w = self.__update_weights(T, self.S, w,
                                      signatures_limit=signatures_limit,
                                      ignorable_signature_indices=ignorable_indices)
            error_post = self.__get_error(T, self.S, w)
            self.__status("\tPost error: {}".format(error_post))
            error_diff = (error_pre - error_post) / error_pre
            self.__status("\tNew normalized weights: ")
            if self.verbose:
                self.__print_normalized_weights(w)

        normalized_weights = w / sum(w)

        # Filter out any weights less than 0.6
        np.place(normalized_weights,
                 normalized_weights < self.signature_cutoff, 0)
        return normalized_weights

    def __print_normalized_weights(self, w):
        """A standard way to print normalized weights given a vector of potentially not yet normalized weights"""
        normalized_weights = w / sum(w)
        for i, weight in enumerate(normalized_weights):
            if weight != 0:
                self.__status(
                    "\t\t{}: {}".format(self.signature_names[i], weight))

    def __status(self, text):
        if self.verbose:
            sys.stdout.write('{}\n'.format(text))

    def __setup_subs_dict(self):
        """
        A dictionary to keep track of the SNVs and the trinucleotide context in which SNVs occurred in order to
        build a mutational signature for the samples
        """
        self.subs_dict = defaultdict(lambda: defaultdict(int))
        for sub in [self.__standardize_subs('C', 'A'),
                    self.__standardize_subs('C', 'G'),
                    self.__standardize_subs('C', 'T'),
                    self.__standardize_subs('T', 'A'),
                    self.__standardize_subs('T', 'C'),
                    self.__standardize_subs('T', 'G')]:
            ref = sub[0]
            for left_bp in ['A', 'C', 'T', 'G']:
                for right_bp in ['A', 'C', 'T', 'G']:
                    trinuc_context = '{}{}{}'.format(left_bp, ref, right_bp)
                    self.subs_dict[sub][trinuc_context] = 0

    def __load_cosmic_signatures(self, ):
        """Load cosmic signatures file. Note that the mutation contexts are listed in alphabetical order:
        (A[C>A]A, A[C>A]C, A[C>A]G, A[C>A]T, A[C>G]A, A[C>G]C... etc.) """
        self.pre_defined_signatures = pd.read_csv(
            '{}'.format(self.cosmic_signatures_filepath), sep='\t',
            engine='python')

    def __load_cosmic_signature_explanations(self):
        """Load about_cosmic_sigs.txt file, which contains correlations and proposed etiologies for the cosmic
        signatures."""
        self.pre_defined_signature_explanations = pd.read_csv(
            '{}'.format(self.cosmic_signature_explanations_filepath),
            sep='\t', engine='python')

    @staticmethod
    def __standardize_subs(ref, alt):
        """
        A function that converts substitutions into their pyrimidine-based notation. Only C and T ref alleles.
        :param ref: The reference allele
        :param alt: The alternate allele
        :return If reference allele is pyrimidine, returns string in format 'ref>alt.' Otherwise, returns string in
        format 'ref_complement_base>alt_complement>base' such that the ref is always a pyrimidine in the return value.
        """
        if ref in DeconstructSigs.purines:
            if alt == "-":
                print()
            return '{}>{}'.format(DeconstructSigs.pair[ref],
                                  DeconstructSigs.pair[alt])
        else:
            return '{}>{}'.format(ref, alt)

    @staticmethod
    def __standardize_trinuc(trinuc):
        """
        A function that ensures trinucleotide contexts are centered around a pyrimidine, using reverse complementary
        sequence to achieve this if necessary.
        :param trinuc: A string representing a trinucleotide context, e.g. 'ACT' or 'GAT'
        :return: An uppercase representation of the given trinucleotide if the center base pair is a pyrimidine,
        otherwise an uppercase representation of the reverse complementary sequence to the given trinucleotide.
        """
        trinuc = trinuc.upper()
        if trinuc[1] in DeconstructSigs.purines:
            return '{}{}{}'.format(DeconstructSigs.pair[trinuc[2]],
                                   DeconstructSigs.pair[trinuc[1]],
                                   DeconstructSigs.pair[trinuc[0]])
        else:
            return trinuc

    @staticmethod
    def __get_reconstructed_tumor_profile(signatures, w):
        """Reconstruct a tumor profile given a set of signatures and a vector of signature weights"""
        w_norm = w / sum(w)
        return w_norm.dot(np.transpose(signatures))

    def __get_error(self, tumor, signatures, w):
        """
        Calculate the SSE between the true tumor signature and the calculated linear combination of different signatures
        :param tumor: normalized array of shape (1, 96) where each entry is a mutation context fraction for the tumor
        :param signatures: array of shape (96, num_signatures) where each row represents a mutation context and each
        column is a signature
        :param w: array of shape (num_signatures, 1) representing weight of each signature
        :return: sum of squares error between reconstructed tumor context fractions and actual tumor profile
        """
        tumor = tumor / sum(tumor)
        reconstructed_tumor_profile = self.__get_reconstructed_tumor_profile(
            signatures, w)
        error = tumor - reconstructed_tumor_profile
        squared_error_sum = np.sum(error.dot(np.transpose(error)))
        return squared_error_sum

    def get_difference_vector(self):
        return self.m96_vector - self.__get_reconstructed_tumor_profile(
            self.S, self.final_result_weights)

    def __update_weights(self, tumor, signatures, w, signatures_limit,
                         ignorable_signature_indices=None):
        """
        Given a set of initial weights, update the weights array with new values that shrink the sum of squares
        error metric.
        :param tumor: normalized array of shape (1, 96) where each entry is a mutation context fraction for the tumor
        :param signatures: signatures: array of shape (96, num_signatures) where each row represents a mutation context
        and each column is a signature
        :param w: array of shape (num_signatures, 1) representing weight of each signature
        :param signatures_limit: How many of the total signatures to consider when assigning weights
        :param ignorable_signature_indices: an array of indices into the signatures array indicating which to ignore
        :return: a new weights array, w.
        """
        if ignorable_signature_indices is None:
            ignorable_signature_indices = []

        # The number of signatures already being used in the current linear combination of signatures
        num_sigs_present = len([weight for weight in w if weight != 0])

        # The total number of signatures to choose from
        num_sigs = np.shape(signatures)[1]

        # The current sum of squares error given the present weights assigned for each signature
        error_old = self.__get_error(tumor, signatures, w)

        # Which weight indices to allow changes for; if we haven't reached the limit all weights are fair game
        if num_sigs_present < signatures_limit:
            changeable_indices = range(num_sigs)
        else:
            # Work with signatures already present if we have reached maximum number of contributing signatures allowed
            changeable_indices = np.nonzero(w)[0]
        changeable_indices = [i for i in changeable_indices if
                              i not in ignorable_signature_indices]

        # zero square matrix of num signatures dimensions
        v = np.zeros((num_sigs, num_sigs))

        # 1 * num signatures vector with values preset to infinity
        new_squared_errors = np.empty(num_sigs, )
        new_squared_errors.fill(math.inf)

        # Only consider adjusting the weights which are allowed to change
        for i in changeable_indices:
            # Find the weight x for the ith signature that minimizes the sum of squared error
            def to_minimize(x):
                # Initialize a temporary zero vector of length number of signatures
                tmp = np.zeros((1, num_sigs))
                tmp[0, i] = x
                return self.__get_error(tumor, signatures, w + tmp[0,])

            error_minimizer = minimize_scalar(to_minimize, bounds=(-w[i], 1),
                                              method="bounded").x
            v[i, i] = error_minimizer
            w_new = w + v[i]
            new_squared_errors[i] = self.__get_error(tumor, signatures, w_new)

        # Find which signature can be added to the weights vector to best reduce the error
        min_new_squared_error = min(new_squared_errors)

        # Update that signature within the weights vector with the new value that best reduces the overall error
        if min_new_squared_error < error_old:
            index_of_min = np.argmin(new_squared_errors, axis=0)
            w[index_of_min] = w[index_of_min] + v[index_of_min, index_of_min]

        return w

    def __seed_weights(self, tumor, signatures, ignorable_indices=None):
        """
        Find which of the cosmic signatures best approximates the tumor signature, and seed the weights such that that
        signature is assigned weight 1 and all other signatures are assigned weight zero. These are the seed weights
        upon which the algorithm will build as it tries to further reduce sum of squared error.
        :param tumor: normalized array of shape (1, 96) where each entry is a mutation context fraction
        :param signatures: array of shape (96, num_signatures) where each row represents a mutation context and each
        column is a signature
        :return: normalized array of shape (num_signatures, 1) representing weight of each signature
        """
        if ignorable_indices is None:
            ignorable_indices = []

        num_sigs = len(signatures[0])
        ss_errors = np.empty(num_sigs, )
        ss_errors.fill(math.inf)
        for i in range(num_sigs):
            if i not in ignorable_indices:
                tmp_weights = np.zeros((num_sigs,))
                tmp_weights[i] = 1
                error = self.__get_error(tumor, signatures, tmp_weights)
                ss_errors[i] = error
        # Seed index that minimizes sum of squared error metric
        seed_index = np.argmin(ss_errors, axis=0)
        final_weights = np.zeros(num_sigs)
        final_weights[seed_index] = 1
        return final_weights
