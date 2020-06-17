function [occ] = fixOccupation(S,occ,nfrm)
% @brief    fixed occupation for SIC

fprintf('old occupation\n')
occ

%reset to integer occupation
%nocc=sum(occ);
%nocc=round(nocc);

% or use FOD input occupation
nocc=nfrm;

occ(:,:)=0.0;
occ(1:nocc(1),1)=1.0;
if S.nspin == 2
    occ(1:nocc(2),2)=1.0; 
end

fprintf('________________________________________________________________\n');
fprintf('                           Eigenvalues\n');
fprintf('________________________________________________________________\n');
fprintf(' Fermi energy = %f\n',S.lambda_f);
ks = 1;
for spin = 1:S.nspin
    if S.nspin == 2
        fprintf('spin %d:   ',spin);
    end
	for kpt = 1:S.tnkpt
        fprintf('%f  ', S.EigVal(:,ks) );
		ks = ks + 1;
	end
    fprintf('\n');
end
fprintf('________________________________________________________________\n');
fprintf('new occupation\n')
occ
end
