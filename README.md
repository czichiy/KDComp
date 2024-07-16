# KDComp

KDComp is a MATLAB project designed for computing kernel dictionaries and related matrices used in signal processing. The project includes functions for calculating kernel dictionaries (KD), RkkD matrices, TD matrices, and a testing script to evaluate the performance of these computations.


## Files Overview

### file 1. CalculateKD

This function computes the KD matrix.  
```sh
KD = CalculateKD(Dic, Nd, sigma, Ruu)

Dic: Dictionary matrix.
Nd: Dimension of the dictionary.
sigma: Standard deviation.
Ruu: Autocorrelation matrix.
```
### file 2. CalculateRkkD

This function computes the RkkD matrix.  
```sh
RkkD = CalculateRkkD(Dic, Nd, sigma, Ruu)
Dic: Dictionary matrix.
Nd: Dimension of the dictionary.
sigma: Standard deviation.
Ruu: Autocorrelation matrix.
```
### file 3. CalculateTD

This function computes the TD matrix.
```sh
TD = CalculateTD(KDvec, Cv, Nd)
KDvec: Vectorized KD matrix.
Cv: Covariance matrix.
Nd: Dimension of the dictionary.
```
Explanation:

	1.	Function Definition:
	•	function TD = CalculateTD(KDvec, CvD, Nd): Defines the function CalculateTD that takes KDvec, CvD, and Nd as inputs and computes the transformation matrix TD.
	2.	Vectorized Inner Product Method:
	•	TD = reshape(KDvec' * reshape(CvD', [], 1), Nd, Nd);
	•	reshape(CvD', [], 1): Converts the matrix CvD into a column vector.
	•	KDvec' * reshape(CvD', [], 1): Computes a vectorized inner product.
	•	reshape(..., Nd, Nd): Reshapes the resulting vector into an  Nd \times Nd  matrix TD, representing the transformation matrix.
	3.	Benefits of Vectorization:
	•	Efficiency: Avoids the inefficiencies of nested loops, using MATLAB’s optimized matrix operations.
	•	Scalability: Particularly advantageous for larger values of Nd, where traditional methods would be slower.
	4.	Uncommented Traditional Approach (for Reference):
	•	Provided as commented-out code, it demonstrates a nested loop method for computing TD. This approach is typically slower and less efficient compared to the vectorized method shown.

...

### File 4. Fix_dic_test  

The MATLAB file Fix_dic_test.m implements a Kernel Least Mean Square (KLMS) algorithm using a Gaussian kernel with a fixed uniform dictionary. The code involves several key components: parameter initialization, input signal generation, filter operation, theoretical model computation, and result plotting.


## Installation

To use KDComp, clone the repository to your local machine:

```sh
git clone https://github.com/czichiy/KDComp.git
```

Make sure you have MATLAB installed on your system.

