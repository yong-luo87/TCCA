function [H, Spect] = KTCCA(K, var_gram_mats, cov_gram_ten, set, para)
% -------------------------------------------------------------------------
% Kernel TCCA
% -------------------------------------------------------------------------

epsilon = para.epsilon;
rDim = para.rDim;

K_ = cell(set.nbV, 1);
for v = 1:set.nbV
    K_{v} = K{v}(set.idxTrnTCCA, set.idxTrnTCCA);
end
clear K

L_mat_inv = cell(set.nbV, 1);
for v = 1:set.nbV
    % termK = (1 - epsilon) * var_gram_mats{v} + epsilon * K_{v};
    termK = var_gram_mats{v} + epsilon * K_{v}; % Here, var_gram_mats{v} = K_{v} * K_{v} / (nbSample-1);
    % termK = termK + 1e-10*eye(size(termK));
    L_mat = chol(termK);
    L_mat_inv{v} = L_mat^-1;
    L_mat_inv{v} = L_mat_inv{v}';
    clear L_mat termK
end

% S_ten = ttm(cov_gram_ten, L_mat_inv{1}', 1);
% S_ten = ttm(S_ten, L_mat_inv{2}', 2);
% cov_gram_mat = double(cov_gram_ten);

S_ten = ttm(cov_gram_ten, L_mat_inv);
% S_mat = double(S_ten);

% -------------------------------------------------------------------------
% if para.rDim == 1, the results of different runs are almost the same,
% otherwise, we should set a random seed.
% -------------------------------------------------------------------------
rand('seed', 123);
% P_kten = parafac_als(S_ten, rDim);
P_kten = cp_als(S_ten, rDim);
% [U, S, V] = svds(double(S_ten), rDim); P_kten.U{1} = U; P_kten.U{2} = V;
clear S_ten

H = cell(set.nbV, 1);
for v = 1:set.nbV
    H{v} = L_mat_inv{v}' * P_kten.U{v};
end
Spect = P_kten.lambda;

end

