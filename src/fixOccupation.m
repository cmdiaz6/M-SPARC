function [occ] = scf(S,occ,nfrm)

fprintf('old occupation\n')
occ

%reset to integer occupation
nocc=sum(occ);
nocc=round(nocc);
% or use FOD input occupation
nocc=nfrm;

occ(:,:)=0.0;
occ(1:nocc(1),1)=1.0;
if S.nspin == 2
    occ(1:nocc(2),2)=1.0; 
end

fprintf('eigenvalues \n')
S.EigVal
fprintf(' Fermi energy = %f\n',S.lambda_f);

end
