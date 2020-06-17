function dpsi = calc_gradients(S, psi)
% @brief    Calculate gradients of psi

%fprintf('calc_gradients %d %d %d %d\n',S.Nev,S.Nx,S.N,S.dx)

dpsi=zeros(size(psi,1),S.Nev,3);

%TODO: possible to use blochGradient (e.g. atomicForce.m)?
%   kpt=1; % only called for kpt=1 currently (flosicSCF.m).
%    kpt_vec = S.kptgrid(kpt,:);
%    dpsi(:,:,1) = blochGradient(S,kpt_vec,1)*psi(:,:);
%    dpsi(:,:,2) = blochGradient(S,kpt_vec,2)*psi(:,:);
%    dpsi(:,:,3) = blochGradient(S,kpt_vec,3)*psi(:,:);

for n = 1:S.Nev
    % reshape to 3d
    psi_n = reshape(psi(:,n),S.Nx,S.Ny,S.Nz);
    % calculate gradients
    [Dpsi_x,Dpsi_y,Dpsi_z] = gradient(psi_n,S.dx,S.dy,S.dz);
    % reshape back to 1d

%% make sure this orders dpsi in the same way
% swapping X and Y
    dpsi(:,n,2) = reshape(Dpsi_x,S.N,1);
    dpsi(:,n,1) = reshape(Dpsi_y,S.N,1);
    dpsi(:,n,3) = reshape(Dpsi_z,S.N,1);
end
fprintf('---calc_gradients is permuting X and Y gradients---\n')

end
