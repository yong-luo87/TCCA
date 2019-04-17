function [var_mats, cov_ten] = var_cov_ten_calculation(X)
% -------------------------------------------------------------------------
% Compute the self variance matrices and covariance matrix
% -------------------------------------------------------------------------

nbView = length(X);
nbSample = size(X{1}, 1);

var_mats = cell(nbView,1);
for v = 1:nbView
    var_mats{v} = X{v}' * X{v};
    var_mats{v} = var_mats{v} / (nbSample-1);
    var_mats{v} = var_mats{v} + eps*eye(size(var_mats{v}));
end

fprintf('\n n: ');
for n = 1:nbSample
    u = cell(nbView,1);
    for v = 1:nbView
        u{v} = X{v}(n,:)';
    end
    cov_x = ktensor(u);
    if n == 1
        cov_ten = full(cov_x);
    else
        cov_ten = cov_ten + full(cov_x);
    end
    clear u cov_x
    
    if rem(n, 100) == 0
        fprintf('%d ', n);
    end
    if rem(n, 1000) == 0
        fprintf('\n n: ');
    end
end
cov_ten = cov_ten / (nbSample-1);

disp('finished!');

end

