function S = apf3d(solve_for,solve_for_direction,xBC,yBC,zBC,kBx,kBy,kBz,dx,D,lambda,eps_bg,eps_or_inv_eps,T,W,P,C_string)
%UNTITLED4 Summary of this function goes here
% Inputs:
%   C:
%       string option is 'transpose(B)'. TODO: implement other options,
%       e.g. numeric matrix.
%
%
%   D (optional):
% Outputs:
%   S:
%       

% Part 1: General initializations

% Start time of general initializations
t0 = clock;

% Set nx, ny, nz and components of epsilon or inverse epsilon
if strcmpi(solve_for,'E')
    if strcmpi(solve_for_direction,'x')
        % zx component of eps as a 3d array:
        eps_zx = eps_or_inv_eps{1};
        [nx,ny,nz] = size(eps_zx);
    elseif strcmpi(solve_for_direction,'y')
        eps_zy = eps_or_inv_eps{2};
        [nx,ny,nz] = size(eps_zy);
    else
        eps_zz = eps_or_inv_eps{3};
        [nx,ny,nz] = size(eps_zz);
    end
else
    if strcmpi(solve_for_direction,'x')
        [nx,ny,nz] = size(eps_or_inv_eps{2});                                      % TO VERIFY: UNLIKE MESTI, INV_EPS_XY,...,YZ ALL HAVE SAME SIZE %
    elseif strcmpi(solve_for_direction,'y')
        [nx,ny,nz] = size(eps_or_inv_eps{1});
    else
        [nx,ny,nz] = size(eps_or_inv_eps{1});
    end
end

% TODO: PMLs
% mesti.m l.678-781, 842-854
    
% Assign default value to the optional argument D
if ~exist('D','var')
    D = [];
end

% The following returns channels and B:
[channels,B] = build_B(T,W,nx,ny,xBC,yBC,n0,m0,dx,lambda,eps_bg,kBx,kBy,solve_for,source_type,source_shape,source,x0,y0);

% End time of general initializations
% Start time to build C
t1 = clock;
general_initializations_time = etime(t1,t0);

% Matrix C
if strcmpi(C_string,'transpose(B)')
    C = transpose(B);
% TODO: elseif C ~= transpose(B)
end

% End time to build C
t2 = clock;
build_C_time = etime(t2,t1);
fprintf('time: %7.3f secs\n', build_C_time);


% Part 2: Build A and solve S

% Build A
% The following also determines whether A is symmetric
[A,is_symmetric_A] = build_A(dx,nx,ny,nz,lambda,xBC,yBC,zBC,solve_for,solve_for_direction,kBx,kBy,kBz,T,W,P, eps_or_inv_eps);

% End time for building A
t3 = clock;
build_A_time = etime(t3,t2);
fprintf('time for building A: %7.3f secs\n', build_A_time);

% Solve S
[S,build_K_time,solving_time] = solve_S(A,is_symmetric_A,B,C_string);

% Start time for solving
t1 = clock;

% Add the constant coefficient from B and C
if strcmpi(xBC,'periodic') || strcmpi(xBC,'Bloch') || strcmpi(yBC,'periodic') || strcmpi(yBC,'Bloch')
    S = (channels.coeff)*S;
%TODO: elseif other boundary conditions
end

% Check size of D
if size(D,1) ~= size(C,1); error('size(D,1) and size(C,1) must be equal, but size(D,1) = %d and size(C,1) = %d.', size(D,1), size(C,1)); end
if size(D,2) ~= size(B,2); error('size(D,2) and size(B,2) must be equal, but size(D,2) = %d and size(B,2) = %d.', size(D,2), size(B,2)); end

% Subtract the contribution of E_in from the scattering matrix
S = S-D;

% End time for solving
t2 = clock;
solving_time = solving_time + etime(t2,t1);

% Total time
total_time = etime(t2,t0);
fprintf('Total elapsed time: %7.3f secs\n', total_time);


end