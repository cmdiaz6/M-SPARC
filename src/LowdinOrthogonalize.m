% for Lowdin Orthogonalization Phi' = S^(-1/2) * Phi
% returns S^(-1/2) matrix
% final transformation to be done outside of function
function [F, U, evals] = LowdinOrthogonalize(S)

%diagonalize S
[U,D] = eig(S);

% check if any evals < 1e-8 and print warning
evals = diag(D);
%fprintf('lowdin sizes %d %d \n',size(S))
fprintf('LOWDEN OVERLAP EIGENVALUES ')
fprintf(' %f  ',evals);
fprintf('\nsum = %f\n',sum(evals))
%fprintf('\n');
if (any(evals<1e-8)) 
    fprintf('small eigenvalues found\n')
    fprintf('%d \n',evals(evals<1e-8))
    % correct to 1e-8?
end

% S^(-1/2) diagonal elements
D=diag(1./sqrt(evals));

% F = U * S^(-1/2)_diag * conj(U')
F=U*D*U';

end


