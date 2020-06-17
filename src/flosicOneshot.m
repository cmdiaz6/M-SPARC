function [Esic, sic_results, fod_forces] = flosicOneshot(S, nfrm, wkpt, psi, occ, FOD)

fprintf('occupation\n');
occ(:,1)'
if S.nspin == 2
    occ(:,2)'
end
fprintf('eigenvalues\n');
S.EigVal

Esic=0;
sic_results=zeros(4,sum(nfrm));
fod_forces =zeros(3,sum(nfrm));

ks=0;
ibeg=1; iend=nfrm(1);
for ispin = 1:S.nspin
    fprintf('\nstarting SIC loop for spin %d\n',ispin);
    if nfrm(ispin) > 0
        % for spin: pass in FODs. Spin up stored after spin dn in S.FOD - S.FOD(1:3,nfrm(1)+1:end)
        for kpt = 1:S.tnkpt
            ks=ks+1;
% Calculate SIC energy
            %[Esic_ks, sic_results(:,ibeg:iend), SIC, AvgSICP(:,ispin)] = flosicEnergy(S, nfrm(ispin), wkpt(kpt), psi(:,:,ks), occ(:,ispin), FOD(:,ibeg:iend));
            [Esic_ks, sic_results(:,ibeg:iend), SIC, ~] = flosicEnergy(S, nfrm(ispin), wkpt(kpt), psi(:,:,ks), occ(:,ispin), FOD(:,ibeg:iend));
            Esic = Esic + Esic_ks;
% Calculate forces on FODs
            % skip for higher kpts. I don't know what FOD forces at other kpts would even mean yet.
            if kpt == 1
                [fod_forces(:,ibeg:iend)] = flosicForces(S, nfrm(ispin), psi(:,:,ks), SIC);
                %store spin-dependent SIC matrix for SCF and atomic Forces
            end
        end
        % send in next FODs for spin=2
        ibeg=ibeg+nfrm(1); iend=iend+nfrm(2);
    end
end

