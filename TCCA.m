function [H, Z] = TCCA(X, var_mats, cov_ten, set, para)
% -------------------------------------------------------------------------
% Tensor CCA
% -------------------------------------------------------------------------

use_inv2 = 1;
rDim = para.rDim;

% -------------------------------------------------------------------------
% Compute the projection matrices for different views
% -------------------------------------------------------------------------
var_mats_inv2 = cell(set.nbV, 1);
for v = 1:set.nbV
    var_mats_inv2{v} = (var_mats{v} + para.epsilon*eye(size(var_mats{v})))^(-0.5);
end
M_ten = ttm(cov_ten, var_mats_inv2);

% -------------------------------------------------------------------------
% if para.rDim == 1, the results of different runs are almost the same,
% otherwise, we should set a random seed.
% -------------------------------------------------------------------------
rand('seed', 123);
P_kten = cp_als(M_ten, rDim);

% -------------------------------------------------------------------------
% Map the data into the low-dimensional space
% -------------------------------------------------------------------------
Z = cell(set.nbV, 1);
H = cell(set.nbV, 1);
for v = 1:set.nbV
    if use_inv2
        H{v} = var_mats_inv2{v} * P_kten.U{v};
    else
        H{v} = P_kten.U{v};
    end
    Z{v} = X{v} * H{v};
end

end

