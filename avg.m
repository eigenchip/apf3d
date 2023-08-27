function avg = avg(n,BC,period,kB)
% Helper function for function build_A.
% Returns avg, an averaging operator as an (n x n) 2d array. The operator
% can be applied along x, y or z depending on its usage in build_A.
% Inputs:
%   n:
%       number of locations. Options are nxEx, nxEy, nxEz, nyEx, nyEy, 
%       nyEz, nzEx, nzEy, nzEz, nxHx, nxHy, nxHz, nyHx, nyHy, nyHz, nzHx, 
%       nzHy and nzHz. 
%   BC:
%       string specifying the boundary conditions in the direction along
%       which avg operates. Options are 'periodic' and 'Bloch'. TODO: other
%       boundary conditions
%   period (optional):
%       required if BC is 'periodic' or 'Bloch'. Period of the system in 
%       the direction along which avg operates.
%   kB (optional):
%       required if BC is 'periodic' or 'Bloch'. Bloch wave number.
%
% Outputs:
%   avg:
%       (n x n) 2d array representing an averaging operator along one
%       dimension

if strcmpi(BC,'periodic') || strcmpi(BC,'Bloch')
    % Averaging operator
    avg = spdiags(1/2*[ones(n,1),ones(n,1),exp(1i*kB*period)*ones(n,1)],[1,0,1-n],n,n);

% TODO: elseif other boundary conditions
end


end