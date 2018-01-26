from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os, sys, h5py
import numpy as np

import wrangler

#-----------------------------------------------------------------------------------------------

# fixed_length window to extract about peaks
window = 100
tf_name = 'CTCF'
file_name = 'ENCFF252PLM.bed.gz'

# path to bed file
data_path = '.'
bed_path = os.path.join(data_path, file_name)

# set genome path
genome_path = os.path.join(data_path, 'hg38.fa')
if not os.path.isfile(genome_path):
    wrangler.download.reference_genome('hg38', data_path)

#-----------------------------------------------------------------------------------------------

# create new bed file with window enforced
file_path = os.path.join(data_path, tf_name + '.bed')
wrangler.bedtools.enforce_constant_size(bed_path, file_path, window, compression='gzip')

# extract sequences from bed file and save to fasta file
fasta_path = os.path.join(data_path, tf_name + '.fa')
wrangler.bedtools.to_fasta(file_path, fasta_path, genome_path)

# parse sequence and chromosome from fasta file
sequences = wrangler.fasta.parse_sequences(fasta_path)

# filter sequences with absent nucleotides
pos_sequences, _ = wrangler.munge.filter_nonsense_sequences(sequences)

# convert filtered sequences to one-hot representation
pos_one_hot = wrangler.munge.convert_one_hot(pos_sequences, max_length=window)

# create negative sequences by performing a di-nucleotide shuffle of positive sequences
neg_fasta_path = os.path.join(data_path, tf_name+'_background.fa')
wrangler.meme.shuffle(fasta_path, neg_fasta_path, kmer=2, verbose=1)

# parse sequence and chromosome from fasta file
neg_sequences = wrangler.fasta.parse_sequences(neg_fasta_path)

# convert filtered sequences to one-hot representation
neg_one_hot = wrangler.munge.convert_one_hot(neg_sequences, max_length=window)

# merge postive and negative sequences
one_hot = np.vstack([pos_one_hot, neg_one_hot])
labels = np.vstack([np.ones((len(pos_one_hot), 1)), np.zeros((len(neg_one_hot), 1))])

# split dataset into training set, cross-validation set, and test set
train, valid, test, _ = wrangler.munge.split_dataset(one_hot, labels, valid_frac=0.1, test_frac=0.2)

# save dataset to hdf5 file
print('saving dataset')
with h5py.File(os.path.join(data_path, tf_name+'_'+str(window)+'.h5'), 'w') as fout:
    X_train = fout.create_dataset('X_train', data=train[0], dtype='float32', compression="gzip")
    Y_train = fout.create_dataset('Y_train', data=train[1], dtype='int8', compression="gzip")
    X_valid = fout.create_dataset('X_valid', data=valid[0], dtype='float32', compression="gzip")
    Y_valid = fout.create_dataset('Y_valid', data=valid[1], dtype='int8', compression="gzip")
    X_test = fout.create_dataset('X_test', data=test[0], dtype='float32', compression="gzip")
    Y_test = fout.create_dataset('Y_test', data=test[1], dtype='int8', compression="gzip")
