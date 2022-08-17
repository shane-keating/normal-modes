% FUNCITON [rho_out,max_err] = stable_strat(rho_in,tol)
%
% Stabilize a given density profile RHO_IN to return a strictly
% monotonically increasting profile RHO_OUT. Optional output MAX_ERR is the
% maximum error in the stablized profile, i.e. max(|rho_in - rho_out|).
% Optional input TOL sets the minimum density step (default = 1E-4 kg/m^3).

function [rho_out,max_err] = stable_strat(rho_in,tol)

  % check inputs and columnize data
  if nargin<2; tol = 1e-4; end
  rho_in = squeeze(rho_in(:)); rho_out = nan(size(rho_in)); 

  % find density steps that are smaller than a given tolerance
  drho = diff(rho_in); ii = drho<=tol; 

  % replace negative/zero density steps with TOL, then rebuild the
  % density profile by summing each density step from the top down. 
  if any(ii) && sum(~ii) > 1
    drho(ii) = tol; 
    rho_out = rho_in(1) + cumsum([0; drho]); 
    max_err = max(abs(rho_out-rho_in));
  else 
    rho_out = rho_in; 
    max_err = nan; 
  end
  
end

