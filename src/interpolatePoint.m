function psi_ai = interpolatePoint(S,psi,ai)
% @brief    interpolatePoint(psi,ai) finds psi at position ai by interpolation (cubic)
%
% @param S				A struct that contains the following (or more) fields:
% @param S.dx			The grid size in x-dimesion
% @param S.dy			The grid size in y-dimesion
% @param S.dz			The grid size in z-dimesion
% @param S.Nx			The number of nodes in [0,L1)
% @param S.Ny			The number of nodes in [0,L2)
% @param S.Nz			The number of nodes in [0,L3)
% @param psi            function values in linear format
% @param ai             position coordinates to calculate interpolation
% @param psi_ai         interpolated function value at point ai
%
% @authors	Carlos Diaz <cmdiaz6@miners.utep.edu>

% order of interpolation
% uses 2*np+1 points in each dimension
%np = 2;
np = 6;

% Modified from Calculate_b_guessRho_Eself.m: line 88
% Indices of closest grid point to atom
pos_ii = round(ai(1) / S.dx) + 1;
pos_jj = round(ai(2) / S.dy) + 1;
pos_kk = round(ai(3) / S.dz) + 1;

%fprintf('nearest grid point %d %d %d\n',pos_ii,pos_jj,pos_kk);

% get surrounding points
% Starting and ending indices of interpolation region
ii_s = pos_ii - np; 
ii_e = pos_ii + np;
jj_s = pos_jj - np;
jj_e = pos_jj + np;
kk_s = pos_kk - np;
kk_e = pos_kk + np;	
% Check if the region is inside the domain (for clusters)
isInside = (ii_s>1) && (ii_e<S.Nx) && (jj_s>1) && (jj_e<S.Ny) && (kk_s>1) && (kk_e<S.Nz);
assert(isInside,'Error: FOD too close to boundary for interpolation'); 

X = ii_s:ii_e;
Y = jj_s:jj_e;
Z = kk_s:kk_e;

% get psi at these points on a 3D grid
%reshape and get slices of psi
psi3d = reshape(psi,[S.Nx S.Ny S.Nz]);

% for testing: psi = square of distance to mid-point
%mid=26; %S.Nx/2;
%for i=1:S.Nx
% for j=1:S.Ny
%  for k=1:S.Nz
%   psi3d(i,j,k)=((i-mid)^2+(j-mid)^2+(k-mid)^2);
%  end
% end
%end

psi_sample = psi3d(X,Y,Z);

%interpolate points to get psi(ai)

%express ai location in terms of grid points
ax = ai(1) / S.dx + 1;
ay = ai(2) / S.dy + 1;
az = ai(3) / S.dz + 1;
%fprintf('ai as grid point at %f %f %f\n',ax, ay, az);

% This is the wrong way to do it.
% The X and Y dimensions seem to be permuted somehow
%[XX,YY,ZZ]=meshgrid(X,Y,Z);
%psi_ai = interp3(XX,YY,ZZ, psi_sample, ax, ay, az, 'cubic');

% This is the correct way to do it
[XX,YY,ZZ]=ndgrid(X,Y,Z);
%psi_ai = interpn(XX,YY,ZZ, psi_sample, ax, ay, az, 'linear');
psi_ai = interpn(XX,YY,ZZ, psi_sample, ax, ay, az, 'cubic');
%psi_ai = interpn(XX,YY,ZZ, psi_sample, ax, ay, az, 'spline');
end
