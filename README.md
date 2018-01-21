# Tutorials tensorflow

Contains 3 notebook-based tutorials:

	* Notebook Tensorflow basics: overview of graphs, variables, ops, placeholders, sessions

	* Notebook Tensorflow neural networks: overview of neural network components

	* Notebook CNN TF analysis: walk through of building, training, testing, and evaluating a CNN for a supervised classificaiton task

	* Notebook CNN RNAcompete analysis: walk through of building, training, testing, and evaluating a CNN for a supervised regression task


# Tutorials deepomics

Contains 3 notebook-based tutorials:

	* Notebook Deepomics CNN TF analysis: example of how to use deepomics to train, test, and evaluate a CNN for a supervised classification task

	* Notebook Deepomics CNN RNAcompete analysis: example of how to use deepomics to train, test, and evaluate a CNN for a supervised regression task
	
	* Notebook Deepomics VAE example Frey Faces: example of training a variational autoencoder to fit the distribution of the data.  



## Dependencies

#### python dependencies:
	- tensorflow (release > 1.0, preferrably r1.4)
	- numpy
	- scipy
	- matplotlib
	- jupyter-notebook
	- PILLOW
	- sklearn
	- h5py
	- six
	- pandas  


#### dependencies for preprocessing:
	- bedtools 
	- wget    (if requiring automatic downloads of files)
	- meme    (if requiring di-nucleotide shuffle)
	- Piranha (if calling peaks)
	- samtools
	- RNAplfold (if predicting RNA secondary structure, apart of Vienna package)


#