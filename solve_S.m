function [S,build_K_time,solving_time] = solve_S(A,is_symmetric_A,B,C_string)
% Helper function for function apf3d. Returns S as an (M_in x M_in) 2d 
% array. Prints information to the standard output or elsewhere depending
% on user's compiler.
% Follows https://github.com/complexphoton/MESTI.m/blob/main/src/
% mesti_matrix_solver.m
% Inputs:
%   A:
%       (nx*ny*nz x nx*ny*nz) sparse matrix representing the discretized 
%       wave operator.
%   B:
%       (nx*ny*nz x M_in) sparse matrix representing the input profiles
%       projected onto the orthonormal basis of channels of the system.
%   C_string:
%       string specifying whether C == transpose(B). Options are
%       'transpose(B)' and TODO: implement option 'not transpose(B)'.
%       (M_out x nx*ny*nz) sparse matrix for projecting the output profiles
%       onto the orthonormal basis of channels. C == transpose(B). TODO: 
%       implement the case C ~= transpose(B)
%
% Outputs:
%   S:
%       (M_out x M_in) scattering matrix of the system. When C ==
%       transpose(B), M_out == M_in. M_in is the number of input channels
%       used to decompose the input profiles (column vectors b). M_out is
%       the number of output channels used to decompose the output profiles
%       (E_scattered).

% Part 1: General initializations
% Start time for general initializations
t0 = clock;

% Check if the file zmumps.m exists in MATLAB's search path
if ~(exist('zmumps','file')); error('File zmumps not found in the current directory.'); end

% Check sizes of matrices A and B
if size(A,1) ~= size(A,2); error('size(A,1) and size(A,2) must be equal, but size(A,1) == %d and size(A,2) == %d.',size(A,1),size(A,2)); end
if size(A,2) ~= size(B,1); error('size(A,2) and size(B,1) must be equal, but size(A,2) == %d and size(B,1) == %d.',size(A,2),size(B,1)); end

% Determine whether matrix K = [A,B;C,0] is symmetric
if is_symmetric_A && strcmpi(C_string,'transpose(B)')
    is_symmetric_K = true;
else
    is_symmetric_K = false;
end
% TODO: Assess symmetry of A when PMLs are included

% Numbers of nonzero elements in A, B and C
fprintf('nnz(A) = %.3g; nnz(B) = %.3g\nnnz(C) = %.3g\n', nnz(A), nnz(B), nnz(C));

% End time for general initializations 
t2 = clock;
general_initializations_time = etime(t2,t0);


% Part 2: MUMPS
% Start time for building K
t1 = clock;

% Set matrices C and D
% Since size(B,2) == size(C,1), M_tot := M_in == M_out
M_tot = size(B,2);
D = sparse(M_tot,M_tot);
K = [[A;transpose(B)],[B;D]];

% Row vector of the locations of the Schur block in A. This vector is used 
% by MUMPS
indices_schur = nx*ny*nz + (1:M_tot);

% End time for building K
t2 = clock;
build_K_time = etime(t2,t1);
fprintf('time to build K and compute Schur complement of K: %7.3f secs\n', build_K_time); 

% Schur complement of K
[id, ordering_time] = schur_complement(nnzA,nnzB,nnzC,nnzS,is_symmetric_A,K,is_symmetric_K,indices_schur,general_initializations_time,build_K_time,factorization_time,solving_time);

% Start time for computing the ordering in MUMPS, using MATLAB's AMD
% See [http://s3.amazonaws.com/researchcompendia_prod/articles/
% b33f9bd20a0ad693c43394e8df1ea1f5-2013-12-23-02-56-16/p381-r_amestoy-
% algorithm.pdf].
t1 = clock;

% Schur complement == -C*inv(A)*B
S = -(id.SCHUR);

% End time for computing the ordering
id.JOB = -2;
[~] = zmumps(id);   
t2 = clock;

% With MUMPS, the time of the factorization is 0 second. So the total 
% factorization time is simply (t2 - t1).
factorization_time = etime(t2,t1);

% Since X = inv(A)*B is not solved, solving_time = 0
solving_time = 0;

% Needed for MUMPS:
str_ordering = 'AMD';

% Total time (general initializations + build K + AMD ordering +
% factorization + solving)
t2 = clock;
total_time = etime(t2,t0);


% TODO: include all INFOG (MUMPS's error messages) (l.883-914)


end