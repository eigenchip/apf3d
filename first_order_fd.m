function diff = first_order_fd(n,BC,kB,period)
% Helper function for function build_A.
% Returns diff, a first-order finite difference operator as an (n x n) 2d 
% array. The operator can act on x (diffx), y (diffy) or z (diffz) 
% depending on its usage in build_A. The finite difference is forward on 
% the E grid, and backward on the staggered H grid, according to Yee's 
% scheme (1966).
% 
% Inputs:
%   n:
%       number of x (nxE or nxH), y (nyE or nyH) or z (nzE or nzH)
%       locations
%   BC:
%       string specifying the boundary conditions in x, y or z. Options are
%       'Bloch' and 'periodic' (TODO: implement other boundary conditions).
%   kB (optional):
%       required if BC is 'Bloch' or 'periodic'. Specifies the Bloch number
%       in x (kxB), y (kyB) or z (kzB). If BC is 'periodic', kB = 0.
%   period:
%       period of the system in x (T), y (W) or z (P)
%   
% Outputs:  
%   diff:
%       (nx*ny*nz x nx*ny*nz) 2d array discretizing the first-order derivative


if strcmpi(BC,'Bloch') || strcmpi(BC,'periodic')
    % First-order finite forward difference operator in one dimension:
    % Boundary conditions are embedded in diff
    diff = spdiags([ones(n,1),-ones(n,1),exp(1i*kB*period)],[1,0,1-n],n,n);

% TODO: elseif other boundary conditions
end


end