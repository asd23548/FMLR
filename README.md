# FMLR
Prepare your data into the format of a matlab struct:
- data.X : the covariate matrix with size n by p, where each row is a flatten image.
- data.Y : the reponse vector with size n by 1, where each element is reponse value in the regression model.
- data.mask : a binary matrix with the same size of the input image, e.g. 32 by 32 reprent a 32 by 32 pixel image. The mask image indicate which pixel should be included in the model, i.e. the number of 1's should be equal to p in data.X. The mask will be used to generate the Gaussian kernel and restore the full image estimation.

Prepare your algorithm setup into the format of a matlab struct:
- options.q_max : the number of maximum mixing conponents.
- options.band_width_s : the vector of band width in the Gassian Kernel function.
- options.lambda_s : the vector of penalty values in the regularization term.
- options.iter_max : the maximum of iterations for each parameter combination.

To run your code, please refer to the example in FMEM_SIM.m
