function [var_gram_mats, cov_gram_ten] = var_cov_gram_ten_calculation(K, set)
% -------------------------------------------------------------------------
% Compute the self-variance gram matrices and covariance gram tensor
% -------------------------------------------------------------------------

nbView = length(K);
nbSample = set.nbTrnTCCA;

K_ = cell(nbView, 1);
for v = 1:nbView
    K_{v} = K{v}(set.idxTrnTCCA, set.idxTrnTCCA);
end
clear K

var_gram_mats = cell(nbView, 1);
for v = 1:nbView
    var_gram_mats{v} = K_{v} * K_{v};
    var_gram_mats{v} = var_gram_mats{v} / (nbSample-1);
end

fprintf('\n n: ');
for n = 1:nbSample
    u = cell(nbView,1);
    for v = 1:nbView
        u{v} = K_{v}(:,n);
    end
    cov_x = ktensor(u);
    if n == 1
        cov_gram_ten = full(cov_x);
    else
        cov_gram_ten = cov_gram_ten + full(cov_x);
    end
    clear u cov_x
    
    if rem(n, 100) == 0
        fprintf('%d ', n);
    end
    if rem(n, 1000) == 0
        fprintf('\n n: ');
    end
end
cov_gram_ten = cov_gram_ten / (nbSample - 1);

disp('finished!');

end

