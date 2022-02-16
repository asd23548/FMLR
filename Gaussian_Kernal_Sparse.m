function K = Gaussian_Kernal_Sparse(mask, r, sigma)
% Generate the sparse Gaussian_Kernal Matrix under threshold r.

[Mask_Idx, ~, Nbr_Dist_Sqr] = Head_File_For_Mask(mask, r*sigma);
q = size(Mask_Idx, 1);
K_value = exp(-Nbr_Dist_Sqr./2/sigma^2);
Nbr_Size = length(Nbr_Dist_Sqr);
v1 = repmat(1:q, Nbr_Size, 1);
v1 = v1(:);
v2 = Mask_Idx(:, 2:end);
v2 = v2';
v2 = v2(:);
v2(v2 == 0) = q + 1;
K_value = repmat(K_value, 1, q);
K = sparse(v1, v2, K_value, q, q+1);
K(:, end) = [];
K = K + speye(q);