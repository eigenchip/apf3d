function [id,ordering_time] = schur_complement(nnzA,nnzB,nnzC,nnzS,K,is_symmetric_K,indices_schur,general_initializations_time,build_K_time,factorization_time,solving_time)
%
% Follows https://github.com/complexphoton/MESTI.m/blob/main/src/
% mesti_matrix_solver.m
%
%
%
%
%
%
%
%
%
%
%
%
%

% Size of K
N = size(K,1);

id = initmumps;

if is_symmetric_K
    % K symmetric and not positive definite
    id.SYM = 2;
else
    % K not symmetric
    id.SYM = 0;
end

% File zmumps.m must exist in MATLAB's search path
id = zmumps(id);

% Print to standard output
id.ICNTL(3) = 6;

% Do not store L and U factors
id.ICNTL(31) = 1;

% Start time for AMD ordering
% See [http://s3.amazonaws.com/researchcompendia_prod/articles/
% b33f9bd20a0ad693c43394e8df1ea1f5-2013-12-23-02-56-16/p381-r_amestoy-
% algorithm.pdf].
t1 = clock;

% Analysis/ordering...
id.JOB = 1;
% ...Done by AMD
id.ICNTL(7) = 0;
% Run zmumps()
id = zmumps(id,K);

% Check error message from MUMPS regarding the ordering
% MUMPS_error_message is a MUMPS helper function
if id.INFOG(1) < 0; error(MUMPS_error_message(id.INFOG)); end

% End time for ordering
t2 = clock;
ordering_time = etime(t2,t1);
fprintf('time for AMD ordering: %7.3f secs\n', ordering_time);

% Start time for factorization
t1 = clock;
id.JOB = 2;
id = zmumps(id,K);

% Check error message from MUMPS regarding factorization
if id.INFOG(1) < 0; error(MUMPS_error_message(id.INFOG)); end


% Directly from [https://github.com/complexphoton/MESTI.m/blob/main/src/
% mesti_matrix_solver.m#L865]:
% What the MEX function zmumpsmex() returns is not schur but its transpose, 
% probably because C is row major while MATLAB is column major. When A is 
% not symmetric, line 77 of the MATLAB interface zmumps.m attempts to undo 
% the transpose by returning schur', which would have worked if ' were to 
% be transpose. But ' is conjugate transpose. So we need to conjugate it
% When A is symmetric, MUMPS only returns the lower triangular part and the
% diagonal of the Schur complement (see MUMPS userguide) line 79 of the 
% MATLAB interface zmumps.m attempts to undo the transpose and to complete 
% the other half by returning triu(schur)+tril(schur',-1), which would have
% worked if ' were to be transpose. But ' is conjugate transpose. So the 
% lower triangular part (excluding the diagonal) has the wrong sign in its 
% imaginary part. Now we need to fix that.
if is_symmetric_K
    id.SCHUR = triu(id.SCHUR) + conj(tril(id.SCHUR,-1));
else
    id.SCHUR = conj(id.SCHUR);
end

% End time for factorization
t2 = clock;
factorization_time = etime(t2,t1);
fprintf('time for factorization: %7.3f secs\n', factorization_time);


% Lines 883-914 from [https://github.com/complexphoton/MESTI.m/blob/main/src/mesti_matrix_solver.m#L883C1-L914C4]:
% Interpret some of the error messages from MUMPS
function msg = MUMPS_error_message(INFOG)

fprintf('\n');
for nn = 1:length(INFOG)
    fprintf('INFOG(%d) = %d\n', nn, INFOG(nn));
end

% Interpret some of the error values; look at MUMPS user guide for complete listing
switch INFOG(1)
    case -1; err_msg = sprintf('An error occurred on processor %d', INFOG(2));
    case -2; err_msg = sprintf('NNZ (or NZ) = %d is out of range', INFOG(2));
    case -3; err_msg = 'MUMPS was called with an invalid value for JOB';
    case -4; err_msg = sprintf('Error in user-provided permutation array PERM_IN at position %d', INFOG(2));
    case -5; err_msg = sprintf('Not enough memory for real workspace allocation during analysis; INFOG(2) = %d', INFOG(2));
    case -6; err_msg = sprintf('Matrix is singular in structure; structural rank %d', INFOG(2));
    case -7; err_msg = sprintf('Not enough memory for integer workspace allocation during analysis; INFOG(2) = %d', INFOG(2));
    case -8; err_msg = 'Integer workarray IS too small for factorization; should increase ICNTL(14) and call again';
    case -9; err_msg = 'Real/complex workarray S too small for factorization; should increase ICNTL(14) and call again';
    case -10; err_msg = 'Matrix is numerically singular';
    case -13; err_msg = sprintf('Not enough memory for real/complex workspace allocation during factorization; INFOG(2) = %d; estimated memory usage = %d MB; actual memory allocation = %d MB', INFOG(2), INFOG(17), INFOG(19));
    case -14; err_msg = 'Integer workarray IS too small for solution; should increase ICNTL(14) and call again';
    case -15; err_msg = 'Integer workarray IS too small for iterative refinement and/or error analysis; should increase ICNTL(14) and call again';
    case -16; err_msg = sprintf('N = %d is out of range', INFOG(2));
    case -17; err_msg = 'Internal send buffer too small; should increase ICNTL(14) and call again';
    case -20; err_msg = 'Internal reception buffer too small; should increase ICNTL(14) and call again';
    case -44; err_msg = sprintf('The solve stage cannot be performed because the factors or part of the factors are not available; ICNTL(31) = %d', INFOG(2));
    otherwise; err_msg = 'See MUMPS user guide';
end
msg = sprintf('MUMPS failed with INFOG(1)=%d, INFOG(2)=%d: %s.', INFOG(1), INFOG(2), err_msg);

end

end