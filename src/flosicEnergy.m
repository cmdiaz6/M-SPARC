function [sic_energy, sic_results, SIC, AvgSICP] = flosicEnergy(S, nfrm, wkpt, psi, occ, FOD)
% @brief    flosicEnergy implements FLOSIC method energy correction
%
% @param S				A struct that contains the following (or more) fields:
% @param S.b 			Pseudocharge. should be zero'd for SIC loop, restored afterwards.
% @param S.Lap			Discretized Laplacian operator (matrix)
% @param S.LapPreconL	Left preconditioner for discretized Laplacian operator (matrix)
% @param S.LapPreconU	Rigtht preconditioner for discretized Laplacian operator (matrix)
% @param S.N			Total number of grid nodes
% @param S.Nev 			Number of states
% @param S.poisson_tol	Tolerance for solving poisson equation
% @param S.Relax_iter	The current relaxation iteration number
% @param S.W			Weights for spatial integration for Cuboidal domain
% @param FOD            coordinates for Fermi-Orbital descriptors ai
% @param nfrm           number of Fermi Orbital descriptors for each spin
% @param occ            occupancies

% @param F              Matrix from Lowdin orthogonalization (FILO)
% @param Q              Eigenvalues from Lowdin Orthogonalization
% @param SIC            SIC matrix : <psi_k | VSIC_l | psi_l>
%
% @authors  Carlos M. Diaz <cmdiaz6@miners.utep.edu>
%

%zero pseudo-charge density, b, for Poisson_RHS
S.b=0;

fprintf('inside sic_energy: nfrm = %d, Nelectron = %d \n',nfrm, S.Nelectron);
fprintf('printing FODs\n');
for iwf=1:nfrm
    fprintf(' %.6f ',FOD(:,iwf));
    fprintf('\n')
end

SIC=zeros(nfrm,nfrm);
AvgSICP=zeros(size(S.psi,1),1);

% transform psi to Fermi-Lowdin orbitals
[psi, ~, ~] = FLOtransform(S, nfrm, psi,FOD, occ);
%print FILO, Q later for debugging flosicForces
%[psi, FILO, Q] = FLOtransform(S, nfrm, psi,FOD, occ);

fprintf(2,' ------------------\n');
% loop through orbitals for SIC energies
sic_energy=0;
fprintf('SIC loop %d\n',nfrm);
for iwf=1:nfrm %(ispin)
    fprintf('starting sic loop %d occ %f\n',iwf,occ(iwf));

    % Electron density in spin-polarized form
    S.rho = zeros(S.N,3);
    %S.rho(:,1) = S.wkpt(kpt)* sum( bsxfun(@times,psi(:,iwf).*conj(psi(:,iwf)),occ(iwf)'),2);
    S.rho(:,1) = wkpt* sum( bsxfun(@times,psi(:,iwf).*conj(psi(:,iwf)),occ(iwf)'),2);
    S.rho(:,1) = real(S.rho(:,1));
    % store all orbital density in spin-up rho
    S.rho(:,2) = S.rho(:,1);
    S.rho(:,3) = 0;
    fprintf('charge for orbital %d: %f\n',iwf,sum(S.rho(:,1) .* S.W))

    % Electrostatic potential
    qsav=S.NetCharge;
    %S.NetCharge = -1; %try for PBC
    S = poissonSolve(S, S.poisson_tol, 10);
    S.NetCharge=qsav; %restore NetCharge
    fprintf('asy %f %f\n',S.phi(1)*sqrt(3*(S.L1/2)^2),S.phi(end)*sqrt(3*(S.L1/2)^2));

    % spin polarized forms always needed for Vxc/Exc
    spinsav=S.nspin;
    S.nspin=2;
    % Exchange-correlation potential
    S = exchangeCorrelationPotential(S);
    S.nspin=spinsav; %restore spin

    % Calculate energy
    %[Etotal,Eband,Exc,Exc_dc,Eelec_dc,Eent] = EvaluateTotalEnergy(S);
    % Exchange-correlation energy
    Exc   =      sum(S.e_xc.*S.rho(:,1).*S.W);
    % Electrostatic energy
    %S.phi = S.phi*(S.nspin/2); % why is this in NRLMOL?
    Eelec = -0.5*sum((S.rho(:,1)).*S.phi.*S.W);

    Etotal_i = Eelec - Exc;
    Etotal_i = Etotal_i*(2.0/S.nspin); % *S.occfac

    % for testing against own old implementation
    % exchange-only
    C2 = -(3/4)*(6/pi)^(1/3); % constant for spin-polarized exchange energy
    Ex = C2*sum((S.rho(:,1).^(4/3)).*S.W) ;
    %Vx = (4/3)*C2*S.rho(:,1).^(1/3) ;
    %correlation-only
    Ec = 0.0;
    %[Ec, ~] = pw91cor(S.rho(:,1),S); 
    %[Ec, Vc] = pw91cor(S.rho(:,1),S); 
    %fprintf('pw91 exchange %f  correlation %f\n',Ex,Ec);
    %fprintf('Exc 1: %f  pw91: %f  diff: %f\n',Exc,Ex+Ec,Exc-Ex-Ec);
    %fprintf('sum diff Vxc %f\n',sum(S.Vxc(:,1)-Vx-Vc));

    sic_results(:,iwf) = [Etotal_i, Eelec, Exc, Ex, Ec]; %, xdel];

    sic_energy = sic_energy + Etotal_i;

    fprintf('sic results: Esic_i, Eelec, Exc, LDA-Ex, LDA-Ec\n');
    fprintf(' %.8f  ',sic_results(:,iwf));
    fprintf(2,'\n ------------------\n');

%%% build SIC matrix column for FOD FORCES %%%
    %SICP = real(S.Vxc + S.phi);
    % only use spin-up part ov Vxc
    SICP = real(bsxfun(@plus,S.phi,S.Vxc(:,1)));
% IN NRLMOL this is over Nbas and includes unoccupied orbitals
% unoccupied orbitals are needed in the SCF section, but not for FOD forces
    for l=1:nfrm
        SIC(l,iwf) = -sum(psi(:,l).*SICP.*S.W.*psi(:,iwf));
    end
%%% Average SIC potential for non-variational SCF test
%%% return into spin dimension in SCF.m
    AvgSICP = AvgSICP + SICP.*S.rho(:,1); %./S.rho(:,1);
    % divide by total density outside
end %nfrm

fprintf('SIC results: (Esic, Ecoul, Exc, Ex, Ec): \n');
for ii = 1:nfrm
    fprintf(' %d   %.6f  %.6f  %.6f  %.6f  %.6f\n',ii,sic_results(1:5,ii));
end
fprintf('SIC energy:                    %.9f (Ha)\n',sic_energy);
%fprintf('Total energy:                  %.9f (Ha)\n',S.Etotal+sic_energy);
fprintf(' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n');

end
