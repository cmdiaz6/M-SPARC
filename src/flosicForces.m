%[fod_forces(:,ibeg:iend)] = flosicForces(S, nfrm(ispin), psi(:,:,ks), SIC);
%function [force_out] = flosicForces(S, nfrm, psi, FILO, Q, SIC)
function [force_out] = flosicForces(S, nfrm, psi, SIC)
% @brief    flosicForces calculates the forces on the flosic FODs
%
% @param S				A struct that contains the following (or more) fields:
% @param S.FOD          coordinates for Fermi-Orbital descriptors ai
% @param psi			The FLO-transformed psi
% @param FILO			Lowdin-Orthogonalized F matrix
% @param Q              Lowdin-Orthogonalization eigenvalues
% @param SIC            SIC matrix <psi_l | VSIC_k | psi_k>

fprintf('flosicForces %d\n',nfrm);
debdax = zeros(3,nfrm);

% construct FMAT and derivatives to get dSnm/dam = <Fn | dFm/dam>
    % Calculate gradient of psi  (grep for this)
Dpsi(:,:,1) = psi;
Dpsi(:,:,2:4) = calc_gradients(S, psi);

% create Fermi-Lowdin orbitals
% create T matrix
FMAT=zeros(nfrm,nfrm,4);
for iwf=1:nfrm
    drho(1:4)=0.0;
    % evaluate drho at ai
    for jwf=1:nfrm
        dpsi(1:4)=0.0;
        for kx=1:4
            dpsi(kx) = interpolatePoint(S, Dpsi(:,jwf,kx), S.FOD(:,iwf));
            drho(kx) = drho(kx) + dpsi(1) * dpsi(kx);
        end
    end
    fprintf('drho %f %f %f %f\n',drho(1),drho(2),drho(3),drho(4))
    for kx=1:4
        for jwf=1:nfrm
            % evaluate psi_j at ai
            psi_ai = interpolatePoint(S, Dpsi(:,jwf,kx), S.FOD(:,iwf));
            FMAT(jwf,iwf,kx) = psi_ai / sqrt(drho(1));
            if kx>=2
                %> \f$ \frac {d} {dx} \frac {\Psi_i} {\sqrt{\rho}} \f$
                FMAT(jwf,iwf,kx) = FMAT(jwf,iwf,kx) - FMAT(jwf,iwf,1)*drho(kx)/drho(1);
            else
                fprintf('%d %d fmat %.4f %.4f\n',iwf,jwf,FMAT(jwf,iwf,kx),psi_ai);
            end
        end
    end
end

for kx=1:4
fprintf('FMAT\n');
for iwf=1:nfrm
    fprintf(' %.4f  \t',FMAT(1:nfrm,iwf,kx))
    fprintf('\n')
end
end

% testing old FILO input against redoing it
O = FMAT(:,:,1)'*FMAT(:,:,1);
[~, FILO, Q] = LowdinOrthogonalize(O);
fprintf('new FILO=\n');
for iwf=1:nfrm
    fprintf(' %.6f ',FILO(:,iwf)) 
    fprintf('\n')
end
% end testing

% --- end FMAT part ---
% check if FMAT(:,:,1) = FILO

%!>  Equation: (22) 
%!>  \f$ ( T_ak T_bl - Tal Tbk ) \f$
%!>  \f$  \langle \phi_l | V_k | \phi_k \rangle \f$
HAM=zeros(nfrm,nfrm);
for ia=1:nfrm
    for ib=1:nfrm
        for kf=1:nfrm
            for lf=1:nfrm
                HAM(ib,ia)=HAM(ib,ia) + ...
                (FILO(kf,ia)*FILO(lf,ib)-FILO(lf,ia)*FILO(kf,ib)) * SIC(kf,lf);
            end
        end
%!>  \f$ \sqrt{Q_a*Q_b} \f$
        HAM(ib,ia)=HAM(ib,ia)/(sqrt(Q(ib))*sqrt(Q(ia)));
    end
end

%!>  \f$ HAM* T_{b,n}\f$
%! SC=FILO*HAM  = T * (T*T-T*T)/(sqrt(Q*Q)  * <phi_l | V_k | phi_k>
%! <Fn | dFm/dAmx>

SC = FILO*HAM;

for ia=1:nfrm
    for kx=1:3
%!>  \f$ dS_{n,m}/dam  \f$
%! OVER = this has to be <Fn|dFm/dAmx>
        OVER = FMAT(:,:,1)' * FMAT(:,:,kx+1);
        for mf=1:nfrm
            for nf=1:nfrm
                debdax(kx,mf)=debdax(kx,mf) - ...
                FILO(mf,ia)*SC(nf,ia)*OVER(nf,mf);
            end
        end
    end % kx=1:3
end % ia=1:nfrm


%!> \f$ (\sqrt{Qb}-\sqrt{Qa})/(\sqrt{Qb}+\sqrt{Qa})  \f$
for ia=1:nfrm
    for ib=1:nfrm
        HAM(ib,ia)=HAM(ib,ia)*( sqrt(Q(ib))-sqrt(Q(ia)) );
        HAM(ib,ia)=HAM(ib,ia)/( sqrt(Q(ib))+sqrt(Q(ia)) );
    end
end

for kx=1:3
    OVER = FMAT(:,:,kx+1)' * FMAT(:,:,1);
    %OVER = 2.0*OVER; % why doubled? was halved later anyways for SC2
    for ia=1:nfrm
        for ib=1:nfrm
            SC1=zeros(nfrm,1);
            for mf=1:nfrm
                for jf=1:nfrm
                    SC1(mf)=SC1(mf) - (FILO(jf,ib)*FILO(mf,ia) + FILO(mf,ib)*FILO(jf,ia)) * OVER(mf,jf);
                end
            end
            for mf=1:nfrm
                debdax(kx,mf)=debdax(kx,mf) - 0.5*SC1(mf)*HAM(ib,ia);
            end
        end
    end
end %kx=1:3

fprintf('FOD forces\n');
for mf=1:nfrm
   fprintf(' %d  %14.6e  %14.6e  %14.6e\n',mf,debdax(1:3,mf));
end

force_out=debdax;

%end function
end
