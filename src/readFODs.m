function S = readFODs(S, filename)
% @brief	readFODs(filename) reads data from FOD positions from file 'filename'
%
% @param filename	The FOD filename, without suffix.
%
% @authors
%

% set flag for FLOSIC run if file is present
S.FlosicFlag = 0;

% open .fod file
fid1=fopen(strcat(filename,'.fod'),'r');
%fid1=fopen(filename,'r'); 
if (fid1 == -1) 
    return
%    error('\n Cannot open file "%s.fod"\n',filename);
%    fprintf('\n Cannot open file "%s.fod"\n',filename);
end 
S.FlosicFlag = 1;
fprintf(' Reading FOD Input file: %s\n',filename);

% # of FOD positions
textscan(fid1,'%s',1,'delimiter','\n');
C_nfrm = textscan(fid1,'%f %f',1,'delimiter','\n');
nfrm = [C_nfrm{1,1},C_nfrm{1,2}];

fprintf(' Number of FODs %d %d\n',nfrm(1),nfrm(2));

% FOD positions
textscan(fid1,'%s',2,'delimiter','\n');
C_FODs = textscan(fid1,'%f %f %f',nfrm(1)+nfrm(2),'delimiter','\n');
coords = [C_FODs{:,1},C_FODs{:,2},C_FODs{:,3}];

% translation vector
textscan(fid1,'%s',2,'delimiter','\n');
C_tvec = textscan(fid1,'%f %f %f',1,'delimiter','\n');
tvec = [C_tvec{1,1},C_tvec{1,2},C_tvec{1,3}];

fclose(fid1);

fprintf('translate vector %f %f %f\n',tvec);
FOD = coords + tvec;

% Transpose to get correct access
FOD = FOD';

% attach to S structure
S.nfrm = nfrm;
S.FOD = FOD;

fprintf('printing FODs: nfrm= %d  %d\n',nfrm);
for ifrm = 1:size(FOD,2)
    fprintf('%d     %.8f %.8f %.8f\n',ifrm,FOD(:,ifrm));
end

end
