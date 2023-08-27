function [channels,B] = build_B(T,W,nx,ny,xBC,yBC,n0,m0,dx,lambda,eps_bg,kBx,kBy,solve_for,source_type,source_shape,source,x0,y0)
% Returns channels as an object with fields describing the propagating
% channels of the system and returns B as an (nx*ny*nz x M_in) 2d array 
% representing the input profiles projected on the orthonormal basis.
% Inputs:
%   x0:
%       optional
%   y0:
%       optional
%
%
%   source_type:
%       type of source specified by the user. Options are 'field values' and
%       'incidence angles' TODO: implement the option 'incident angles'
%   source_shape:
%       shape of the user-specified source. Options are 'point (x0,y0,0)',
%       'line (x,y0,0)', 'line (x0,y,0)', 'plane (x,y,0)'.
%   source:
%       if source_type is 'field', then source is an ...or... array
%       size(source,2) must be M_in
%   Remarks:
%       any backpropagation to the edges of the simulation domain (for sources) should be handled by the user;
%       if user wants to vary omega or lambda, they need to call mesti for every iteration of lambda or omega
%       user should call [channels,~] = build_B(...) to get helpful information about the orthonormal propagating channels, to build their own inputs
%       in order to vary the location of the source, user should call build_B anew, since build_B fixes the location once.
%   T:
%       period of the system in x
%   W:
%       period of the system in y
%   nx:
%       number of pixels in x
%   ny:
%       number of pixels in y
%   xBC:
%   yBC:
%   n0:
%   m0:
% Outputs:
%   channels:
%       object with fields specifying properties of the channels
%           channels.list_kxdx: 
%           channels.list_kydy:
%           channels.list_kzdz:
%           channels.list_kxdx_prop:
%               (Nprop x 1) column vector of propagating wave numbers kx
%           channels.list_kydy_prop:
%               (Nprop x 1) column vector of propagating wave numbers ky
%           channels.list_kzdz_prop:
%               (Nprop x 1) column vector of propagating wave numbers kz
%           channels.Nprop:
%               number of propagating channels in the orthonormal basis
%           channels.u:
%               function with handles kxdx and kydy. When kxdx and kydy are scalars, the
%               output is an (nx*ny x 1) column vector. When kxdx and kydy are
%               row vectors, the output is an (nx*ny x M_in) 2d array.
%           channels.indices_prop:
%           channels.indices_change_sign:
%           channels.flux_normalization_prop:
%
%
%   B: 2d array of size (nx*ny*nz) x M_in


% Part 1: Build the orthonormal basis of propagating channels

if strcmpi(xBC,'Bloch') && strcmpi(yBC,'Bloch') || strcmpi(xBC,'periodic') && strcmpi(yBC,'periodic')
    % Index of the transverse mode such that kxdx = kBx*T/nx or kydy = kBy*W/ny, is the median index. It is the same index for kx and ky.
    if mod(nx,2) == 1
        index_zero_kx = round((nx+1)/2);
    else
        index_zero_kx = round(nx/2);
    end
    
    % Wave numbers kx and ky
    channels.list_kxdx = kBx*T/nx + 2*pi/nx*((1:nx)-index_zero_kx);
    channels.list_kydx = kBy*W/ny + 2*pi/ny*((1:ny)-index_zero_kx);

    % Transverse modes satisfying periodic or Bloch periodic boundary conditions:
    % Create an (nx x ny) meshgrid with square spacing dx^2 for the domain of function u.
    [X,Y] = meshgrid(1:dx:T,1:dx:W);
    % reshape() permutes the y entries before x, so in order to permute x before y, we transpose u(x,y).
    % Since length(kxdx) == length(kydy) == M_in, reshape(...,[],M_in) returns an (nx*ny x M_in) 2d array.
    channels.u = @(kxdx,kydy) reshape(transpose(1/sqrt(nx)*exp(1i*kxdx*(X-n0)) * 1/sqrt(ny)*exp(1i*kydy*(Y-m0))),[],length(kxdx));

    % Total wave number squared with grid size normalization
    tot_wave_num_squared = (2*pi/lambda *dx)^2;

    % Wave numbers kz
    kzdz_squared = tot_wave_num_squared - (channels.list_kxdx).^2 - (channels.list_kydy).^2;
    channels.list_kzdz = sqrt(kzdz_squared);
    
    % Extracting the propagating channels:
    % Indices of the propagating channels in channels.list_kzdz
    channels.indices_prop = find(imag(kzdz)==0);
    % When tot_wave_num_squared is complex, singularities have to be handled so that kxdx and kydy 
    % still correspond to propagating waves. Let z^2:= tot_wave_num_squared/dx^2 = (omega/c)^2 *eps_bg
    % and c:= |kx^2 + ky^2|. Then the kz wavenumber function is of the form sqrt(z^2 - c) with singularities
    % at z = +c = +(kx^2 + ky^2) and z = -c = -(kx^2 + ky^2). For Bloch periodic boundary conditions
    % in x,y and z, circumventing these points corresponds to changing sign of tot_wave_num_squared.
    if ~isreal(tot_wave_num_squared)
        channels.indices_change_sign = find(real(kzdz_squared)<0 && imag(kzdz_squared)<0);
        channels.list_kzdz(channels.indices_change_sign) = - channels.list_kzdz(channels.indices_change_sign);
    end
    % Wave numbers of propagating channels
    channels.list_kxdx_prop = channels.list_kxdx(channels.indices_prop);
    channels.list_kydy_prop = channels.list_kydy(channels.indices_prop);
    channels.list_kzdz_prop = channels.list_kzdz(channels.indices_prop);
    % Number of propagating channels
    channels.Nprop = length(channels.indices_prop);
    % Flux normalization of propagating channels
    if strcmpi(solve_for,'E')
        % Column vector of flux normalization constants
        channels.flux_normalization_prop = sqrt(sin(channels.list_kzdz_prop));
        % Diagonal matrix of flux normalization constants
        channels.flux_normalization_prop_diag = spdiags(channels.flux_normalization_prop,0,size(channels.flux_normalization_prop),size(channels.flux_normalization_prop));
    else
        channels.flux_normalization_prop = sqrt(sin(channels.list_kzdz_prop))/eps_bg;
        channels.flux_normalization_prop_diag = spdiags(channels.flux_normalization_prop,0,size(channels.flux_normalization_prop),size(channels.flux_normalization_prop))/eps_bg;
    end

% TODO: elsif other boundary conditions
end


% Part 2: Build the matrix B of projected input profiles

% Reformat the user-specified source
if strcmpi(source_type, 'field values')
    % size(source,2) == M_in
    M_in = size(source,2);
    
    % Initialize an (nx*ny x M_in) sparse matrix
    source_reformat = sparse(nx*ny,M_in);
    
    if strcmpi(source_shape, 'point (x0,y0,0)')
        % size(source) == [1,M_in]
        % Row index of the nonzero pixel in the reformatted (nx*ny x M_in) source
        nonzero_pixel_row_index = nx*(y0-1)+x0;
        % Insert the nonzero row in the sparse matrix
        source_reformat(nonzero_pixel_row_index,:) = source;
    elseif strcmpi(source_shape, 'line (x,y0,0)')
        % size(source) == [nx,M_in]
        % Row index of the first nonzero pixel in the reformatted (nx*ny x M_in) source
        nonzero_pixels_first_row_index = size(source,1)*(y0-1)+1;
        % Insert the nx nonzero rows in the sparse matrix
        source_reformat(nonzero_pixels_first_row_index:(nonzero_pixels_first_row_index-1),:) = source;
    elseif strcmpi(source_shape, 'line (x0,y,0)')
        % size(source) == [ny,M_in]
        % Indices of the nonzero rows
        nonzero_row_indices = x0 + nx*((1:size(source,1))-1);
        % Create row vectors i,j and s of the nonzero locations for sparse(i,j,s)
        i = repmat(nonzero_row_indices,1,M_in);
        j = repelem(1:M_in,size(source,1));
        s = reshape(source,1,numel(source));
        % Assemble the sparse matrix
        source_reformat = sparse(i,j,s,nx*ny,M_in);    
    % Plane source (x,y,0):
    else
        % size(source) == [nx*ny,M_in]
        % source is already in the desired format
        source_reformat = source;
    end

% TODO: elseif other types of user-specified sources, ex: list of plane wave incidence angles
end


% (nx*ny x Nprop) matrix of propagating channels
u_prop = channels.u(channels.list_kxdx_prop,channels.list_kydy_prop);

% Project the user-specified source onto the matrix of propagating channels, using discrete complex inner product
% Resulting matrix has size (Nprop x M_in)
source_proj = u_prop'*source_reformat;

% Compute the nonzero block for matrix B
% Resulting matrix has size (nx*ny x M_in)
nonzero_block = channels.flux_normalization_prop_diag*u_prop*source_proj;

% Constant coefficient to be added at the end (improves C = transpose(B)):
if strcmpi(xBC,'Bloch') && strcmpi(yBC,'Bloch') || strcmpi(xBC,'periodic') && strcmpi(yBC,'periodic')
    channels.coeff = -2*i * -2*i;                                                                                   % NOT SURE %
% TODO: elseif other boundary conditions
end

% Pad with zeros to get matrix B
zero_block = zeros(nx*ny*(nz-1),M_in);
B = [nonzero_block;zero_block];


end