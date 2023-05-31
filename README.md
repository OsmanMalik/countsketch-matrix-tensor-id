# Fast randomized matrix and tensor interpolative decomposition using CountSketch
This repository provides the code we use in our paper 

O.A. Malik, S. Becker. *Fast randomized matrix and tensor interpolative decomposition using CountSketch*. **Adv Comput Math** 46, article number: 76, 2020. https://doi.org/10.1007/s10444-020-09816-9

You can also view the paper for free on SharedIt via this link: https://rdcu.be/b9eFU

## Some further details
In our paper, we propose a new fast randomized algorithm for interpolative decomposition (ID) of matrices utilizing the CountSketch technique. We then extend our method to the recently proposed tensor ID problem. This repository provides all the code we use in our experiments in the paper, including implementations of our proposed CountSketch matrix and tensor ID algorithms which are provided in the functions **CS_matrix_ID.m** and **CS_tensor_ID.m**, respectively.

## Requirements
Parts of our code relies on the following software:
* Tensor Toolbox by Bader, Kolda and others (available at https://www.tensortoolbox.org). 
We used version 2.6 of Tensor Toolbox in our work, but newer versions should also work fine.
* RSVDPACK by Voronin and Martinsson (available at https://github.com/sergeyvoronin/LowRankMatrixDecompositionCodes).
* Interpolative Decomposition based on Strong RRQR by Xin Xing (available at https://www.mathworks.com/matlabcentral/fileexchange/69352-interpolative-decomposition-based-on-strong-rrqr).

## Installation
The mex files can be compiled by using the standard Matlab mex commands. For the other software listed above, please see the instructions provided with each of them. 

## Experiments
The four script files **run_experiment_n.m** (where n = 1, 2, 3, 4), can be used to reproduce our experiments. Each experiment is described below.
* **Experiment 1:** This experiment compares standard matrix ID, Gaussian matrix ID, SRFT matrix ID, and two versions of our proposed CountSketch matrix ID (with strongly rank-revealing QR, or with column pivoted QR). It can be set up so that the test matrices are either dense with a specific structure, or just random sparse matrices with no particular structure aside from sparsity. Our paper does not include any results from this experiment.
* **Experiment 2:** This experiment compares standard matrix ID, Gaussian matrix ID, SRFT matrix ID, and our proposed CountSketch matrix ID (with column pivoted QR). The test matrices in this experiment are sparse, with a specific structure. The figure below shows results from this experiment that are included in our paper. In the experiment, we decompose sparse matrices with density 0.5%, a varying number of rows and 10,000 columns. The target rank is 1,000.
<p align="center">
	<img src="plot-experiment-2b-error.png" width="60%">
	<img src="plot-experiment-2b-time.png" width="60%">
</p>

* **Experiment 3:** This experiment compares Gaussian matrix ID, SRFT matrix ID, and our proposed CountSketch matrix ID (with column pivoted QR). The decomposed matrix is a real-world matrix. In our paper, we use the matrix *specular*, which is available at https://sparse.tamu.edu/Brogan/specular. It has 477,976 rows and 1,600 columns, contains 7,647,040 nonzero elements, and has a rank of 1,442. In the experiment, we set the target rank to 1,442. The table below shows results from this experiment that are included in our paper.

| Algorithm | Error | Run time (s) |
| :--- | :--- | :--- |
| Gaussian | 1.505e-15 | 20.38 |
| SRFT | 1.507e-15 | 18.40 |
| CountSketch (proposal) | 1.504e-15 | 0.59 |
* **Experiment 4:** This experiment compares tensor ID using the Gram matrix approach, Gaussian tensor ID, and CountSketch tensor ID. The test tensors in this experiment are sparse CP tensors with a specific structure. The figure below shows results from this experiment that are included in our paper. In the experiment, we consider 5-way CP tensors with each side of equal and varying size. The factor matrix density is 1%, the number of rank-1 terms used for representing each tensor is 10,000, and we use a target rank of 1,000.
<p align="center">
	<img src="plot-experiment-4d-error.png" width="60%">
	<img src="plot-experiment-4d-time.png" width="60%">
</p>

## Referencing this code
If you use our code in any of your own work, please reference our paper:
```
@article{malik-2020-fast-randomized,
  author    = {Osman Asif Malik and Stephen Becker},
  title     = {Fast randomized matrix and tensor interpolative decomposition using {C}ount{S}ketch},
  journal   = {Advances in Computational Mathematics},
  volume    = {46},
  year      = {2020},
  doi	    = {10.1007/s10444-020-09816-9},
}
```

Most of the code in this repository is implementations of algorithms invented by other researchers. We have done our best to include relevant references in the comments of our code.

## Author contact information
Please feel free to contact me at any time if you have any questions or would like to provide feedback on this code or on our paper. I can be reached at osman.malik@colorado.edu.

## Licenses
The code in the folder AlgorithmCUR is downloaded from Christos Boutsidis's website (http://www.boutsidis.org/software.html) and is included with this repository for convenience.

All other code in this project falls under the license in the root of this project.
