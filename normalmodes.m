% [kd,Gvec,Fvec] = NORMALMODES(dz,rho,F)
% 
% Calculate eigenvalues (deformation wavenumbers) and eigenvectors
% (normal modes) of the dual Sturm-Liouville problems
% 
%     (s*Fm')' = - km^2*Fm,     Fm'(-H) = Fm'(0) = 0
%     
%     s*Gm'' = - km^2*Gm,       Gm(-H) = Gm(0) = 0
% 
% where s = f^2/N^2 and m = 0, 1, 2, ... These are equivalent eigenvalue
% problems, as can be seen from the relations
% 
%     Gm' = - km*Fm,      s*Fm' = km*Gm. 
% 
% Specifically, Gm represents w-like modes (vertical velocity) and Gm
% represents psi-like modes (horizontal streamfunction). Note that the
% barotropic G-mode is identically zero (G1 = 0) for a rigid lid because 
% the depth-averaged vertical velocity must be zero in that case. By 
% contrast, the barotropic F-mode is simply constant (F1 = 1) for a rigid
% lid because there can still be a depth-averaged horizontal flow. In
% either case, the barotropic deformation wavenumber k1 = 0 for a rigid
% lid. 
% 
% The problem is solved by discretizing the vertical second derivative in
% the Gm equation using the density profile rho(n) and the corresponding 
% layers of thickness dz(n). The parameter F = f0^2*(L/2*pi)^2/(g'*H0). The
% vectors Fm are then calculated using Fm = - Gm'/km (plus F1 = 1). 
% 
% SEE ALSO: vmodes.m 

function [kd,Gvec,Fvec] = normalmodes(dz,rho,F)

nz = length(dz);

% Normalize layer thicknesses and depths by total depth
H = sum(abs(dz)); dz = abs(dz)/H;

% Density difference between layers normalized by average difference
drho = rho(2:nz)-rho(1:nz-1); drho0 = mean(drho); drho = drho/drho0;      
if any(drho<=0); error('drho must all be positive'); end

% components of tri-diagonal matrix for G-modes
a = F./drho(1:nz-1)./dz(1:nz-1); 
c = F./drho(1:nz-1)./dz(2:nz); 
b = - a - c;
A = diag(b) + diag(a(2:nz-1),- 1) + diag(c(1:nz-2),+1);
[V,D] = eig(A); [kd,ri] = sort(sqrt(-diag(D))); 
Gvec = V(:,ri); 

% Normalize G-modes so that <Gm*Gn/s(z)> = delta_mn.
for m=1:nz-1
  Gvec(:,m)=Gvec(:,m)/sqrt((Gvec(:,m).*drho)'*Gvec(:,m)/F);
end
Gvec = [zeros(1,nz-1); real(Gvec); zeros(1,nz-1)];     

% Calculate the F-modes from the G-modes
Fvec = diag(1./dz(1:nz))*(Gvec(1:nz,:)-Gvec(2:nz+1,:))*diag(-1./kd(1:nz-1)); 

% Normalize F-modes so that <Fm*Fn> = delta_mn and F(0) > 0
for m = 1:nz-1
  Fvec(:,m)=Fvec(:,m)/sqrt((Fvec(:,m).*dz)'*Fvec(:,m));
  if Fvec(1,m)<0; Fvec(:,m) = -Fvec(:,m); Gvec(:,m) = -Gvec(:,m); end
end
Fvec = [ones(nz,1) real(Fvec)]; 

% Include barotropic modes and eigenvalue
kd = [0; kd(:)];  

end