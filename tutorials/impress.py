
import matplotlib.image as mpimg
from scipy.misc import imresize
import os, sys
import numpy as np

def normalize_pwm(pwm, factor=None):

    MAX = np.max(np.abs(pwm))
    pwm = pwm/MAX
    if factor:
        pwm = np.exp(pwm*factor)
    norm = np.outer(np.ones(pwm.shape[0]), np.sum(np.abs(pwm), axis=0))
    return pwm/norm


def load_alphabet(char_path, alphabet, colormap='standard'):

    def load_char(char_path, char, color):
        colors = {}
        colors['green'] = [10, 151, 21]
        colors['red'] = [204, 0, 0]
        colors['orange'] = [255, 153, 51]
        colors['blue'] = [0, 0, 204]
        colors['cyan'] = [153, 204, 255]
        colors['purple'] = [178, 102, 255]
        colors['grey'] = [160, 160, 160]
        colors['black'] = [0, 0, 0]

        img = mpimg.imread(os.path.join(char_path, char+'.eps'))
        img = np.mean(img, axis=2)
        x_index, y_index = np.where(img != 255)
        y = np.ones((img.shape[0], img.shape[1], 3))*255
        for i in range(3):
            y[x_index, y_index, i] = colors[color][i]
        return y.astype(np.uint8)


    colors = ['green', 'blue', 'orange', 'red']
    if alphabet == 'dna':
        letters = 'ACGT'
        if colormap == 'standard':
            colors = ['green', 'blue', 'orange', 'red']
        chars = []
        for i, char in enumerate(letters):
            chars.append(load_char(char_path, char, colors[i]))

    elif alphabet == 'rna':
        letters = 'ACGU'
        if colormap == 'standard':
            colors = ['green', 'blue', 'orange', 'red']
        chars = []
        for i, char in enumerate(letters):
            chars.append(load_char(char_path, char, colors[i]))


    elif alphabet == 'structure': # structural profile

        letters = 'PHIME'
        if colormap == 'standard':
            colors = ['blue', 'green', 'orange', 'red', 'cyan']
        chars = []
        for i, char in enumerate(letters):
            chars.append(load_char(char_path, char, colors[i]))

    elif alphabet == 'pu': # structural profile

        letters = 'PU'
        if colormap == 'standard':
            colors = ['cyan', 'purple']
        elif colormap == 'bw':
            colors = ['black', 'grey']
        chars = []
        for i, char in enumerate(letters):
            chars.append(load_char(char_path, char, colors[i]))

    return chars



def seq_logo(pwm, height=30, nt_width=10, norm=0, alphabet='dna', colormap='standard'):

    def get_nt_height(pwm, height, norm):

        def entropy(p):
            s = 0
            for i in range(len(p)):
                if p[i] > 0:
                    s -= p[i]*np.log2(p[i])
            return s

        num_nt, num_seq = pwm.shape
        heights = np.zeros((num_nt,num_seq));
        for i in range(num_seq):
            if norm == 1:
                total_height = height
            else:
                total_height = (np.log2(num_nt) - entropy(pwm[:, i]))*height;
            if alphabet == 'pu':
                heights[:,i] = np.floor(pwm[:,i]*np.minimum(total_height, height));
            else:
                heights[:,i] = np.floor(pwm[:,i]*np.minimum(total_height, height*2));

        return heights.astype(int)


    # get the alphabet images of each nucleotide
    package_directory = '.'
    char_path = os.path.join(package_directory,'chars')
    chars = load_alphabet(char_path, alphabet, colormap)

    # get the heights of each nucleotide
    heights = get_nt_height(pwm, height, norm)

    # resize nucleotide images for each base of sequence and stack
    num_nt, num_seq = pwm.shape
    width = np.ceil(nt_width*num_seq).astype(int)

    if alphabet == 'pu':
        max_height = height
    else:
        max_height = height*2
    #total_height = np.sum(heights,axis=0) # np.minimum(np.sum(heights,axis=0), max_height)
    logo = np.ones((max_height, width, 3)).astype(int)*255;
    for i in range(num_seq):
        nt_height = np.sort(heights[:,i]);
        index = np.argsort(heights[:,i])
        remaining_height = np.sum(heights[:,i]);
        offset = max_height-remaining_height

        for j in range(num_nt):
            if nt_height[j] > 0:
                # resized dimensions of image
                nt_img = imresize(chars[index[j]], (nt_height[j], nt_width))

                # determine location of image
                height_range = range(remaining_height-nt_height[j], remaining_height)
                width_range = range(i*nt_width, i*nt_width+nt_width)

                # 'annoying' way to broadcast resized nucleotide image
                if height_range:
                    for k in range(3):
                        for m in range(len(width_range)):
                            logo[height_range+offset, width_range[m],k] = nt_img[:,m,k];

                remaining_height -= nt_height[j]

    return logo.astype(np.uint8)