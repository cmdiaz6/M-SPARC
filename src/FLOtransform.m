
function [psi, FILO, Q] = FLOtransform(S,nfrm,psi,FOD,occ)
% @brief    FLOtransform transforms psi to Fermi-Lowdin orbitals
%
% @param S				A struct that contains the following (or more) fields:
% @param FOD          coordinates for Fermi-Orbital descriptors ai
% @param nfrm           Number of states
% @param occ            occupancies

% @param FILO           Matrix from Lowdin orthogonalization
% @param Q              Eigenvalues from Lowdin Orthogonalization
%
% @authors  Carlos M. Diaz <cmdiaz6@miners.utep.edu>

% create Fermi-Lowdin orbitals
% create T matrix
T=zeros(nfrm,nfrm);
for iwf=1:nfrm
    drho=0.0;
    %fprintf('interpolate FOD %d: %f %f %f \n',iwf,S.FOD(:,iwf))
    % evaluate rho at ai
    for jwf=1:nfrm
        psi_ai = interpolatePoint(S, psi(:,jwf), FOD(:,iwf));
        drho = drho + psi_ai * psi_ai;
%        drho = drho + psi_ai * psi_ai * occ(jwf); % testing partial occupation
    end
%    fprintf('rho @ ai %d is %f\n',iwf,drho);
    for jwf=1:nfrm
        % evaluate psi_j at ai
        psi_ai = interpolatePoint(S, psi(:,jwf), FOD(:,iwf));
%        fprintf('psi %d @ ai %d is %f\n', jwf, iwf, psi_ai);
        T(iwf,jwf) = psi_ai / sqrt(drho);
%        T(iwf,jwf) = sqrt(occ(jwf)) * psi_ai / sqrt(drho); % testing partial occupation
    end
end
%T' %print

% Lowdin-orthogonalize wavefunctions psi (on occupied space only)
%  closic: O=T'*T, psi=psi*T*F
%O = T'*T;
O = T*T';
%fprintf('O before LowdinOrthogonalize')
%O
[F, FILO, Q] = LowdinOrthogonalize(O);
%F' %print
%fprintf('transformation matrix T^T * F')
%T'*F

% apply transformation on psi (done as transpose)
%psi(:,1:nfrm)=psi(:,1:nfrm)*T*F; %use with O=T'T ?
psi(:,1:nfrm)=psi(:,1:nfrm)*T'*F; %use with O=T*T'
end

