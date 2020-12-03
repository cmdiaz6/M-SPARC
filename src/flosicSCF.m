function [Esic, sic_results, Veff] = flosicSCF(S, Veff, nfrm, wkpt, psi, rho, occ, FOD)

% Add Average SIC potential to S.Veff for psuedo-self-consistent FLOSIC
%bsav=S.b; %save S.b
Esic=0;
AvgSICP=zeros(size(S.psi,1),S.nspin);

ks=0;
ibeg=1; iend=nfrm(1);
for ispin = 1:S.nspin
    fprintf('\nstarting SIC loop for spin %d\n',ispin);
    if nfrm(ispin) > 0
        % for spin: pass in FODs. Spin up stored after spin dn in S.FOD - S.FOD(1:3,nfrm(1)+1:end)
        for kpt = 1:S.tnkpt
            ks=ks+1;
% Calculate SIC energy
% SIC matrix only needed for flosicForces post-SCF
            [Esic_ks, sic_results(:,ibeg:iend), ~, AvgSICP(:,ispin)] = flosicEnergy(S, nfrm(ispin), wkpt(kpt), psi(:,:,ks), occ(:,ispin), FOD(:,ibeg:iend));

% divide average SIC potential by rho. SICP.*rho_i done in flosicEnergy
            if S.nspin == 1
                % multiply by 2, to account for spin
                AvgSICP(:,ispin) = 2.0*AvgSICP(:,ispin)./rho(:);
            else
                AvgSICP(:,ispin) = AvgSICP(:,ispin)./rho(:,ispin+1);
            end

            % for testing SIC scaling factor
            %alpha = 0.5;
            %fprintf('scaling AvgSIC by %.2f\n',alpha);
            %AvgSICP(:,ispin) = alpha*AvgSICP(:,ispin);
            %fprintf('sum AvgSICP %f\n',sum(AvgSICP(:,ispin).*S.W));

            % minus sign gives right answer
            Veff(:,ispin) = Veff(:,ispin) - AvgSICP(:,ispin);

            Esic = Esic + Esic_ks;
        end
        % send in next FODs for spin=2
        ibeg=ibeg+nfrm(1); iend=iend+nfrm(2);
    end
end
fprintf(' Etot+Esic = %.8f\n', S.Etotal+Esic);

% transform back for SCF
%S.b=bsav; %necessary?

        %plot AvgSICP along Z-axis
        %vplot = reshape(AvgSICP,S.Nx,S.Ny,S.Nz);
        %loglog(1:round(S.Nz/2),squeeze(vplot(round(end/2),round(end/2),round(end/2):end)))
end
