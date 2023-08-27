function [A,is_symmetric_A] = build_A(dx,nx,ny,nz,lambda,xBC,yBC,zBC,solve_for,solve_for_direction,kBx,kBy,kBz,T,W,P, eps_or_inv_eps)
% Builds the FDFD 3D wave operator from the wave equation in either E or H.
% Inputs:
%   dx:
%       grid size in x. dx == dy == dz.
%   nxE:                                                                    % TO BE ADDRESSED: I HAVE NOT IMPLEMENTED YEE'S SCHEME ANYWHERE YET%
%       number of x_n locations on the E grid
%   nyE:
%       number of y_m locations on the E grid
%   nzE:
%       number of z_p locations on the E grid
%   nxH:
%       number of x_{n-0.5} locations on the H grid
%   nyH:
%       number of y_{m-0.5} locations on the H grid
%   nzH:
%       number of z_{p-0.5} locations on the H grid 
%   lambda:
%       vacuum wavelength of the inputs
%   xBC:
%       boundary condition in x
%   yBC:
%       boundary condition in y
%   zBC:
%       boundary condition in z
%   solve_for:   
%       string specifying whether E or H is solved
%   kBx (optional):
%       optional, required for Bloch periodic boundary conditions in x
%   kBy (optional):
%       optional, required for Bloch periodic boundary conditions in y
%   kBz (optional):
%       optional, required for Bloch periodic boundary conditions in z
%   eps_or_inv_eps:
%       cell array which entries are 3d arrays of different sizes, each
%       representing one component of the permittivity or inverse 
%       permittivity tensor.
%       To ensure that the expression "eps_or_inv_eps(:)" (in assigning the
%       wave operators) permutes the x values first, then the y values and
%       the z values last, size(eps_or_inv_eps) must have the form 
%       [nx_{eps_or_inv_eps},ny_{eps_or_inv_eps},nz_{eps_or_inv_eps}].
%
%       if solve_for is 'E'
%           if solve_for_direction is 'x'
%               zx component is required
%           if solve_for_direction is 'y'
%               zy component is required
%           if solve_for_direction is 'z'
%               zz component is required
%
%           eps_zx = eps_or_inv_eps{1}, size(eps_zx) == [nxEx,nyEx,nzEx]
%           eps_zy = eps_or_inv_eps{2}, size(eps_zy) == [nxEy,nyEy,nzEy]
%           eps_zz = eps_or_inv_eps{3}, size(eps_zz) == [nxEz,nyEz,nzEz]
%
%       if solve_for is 'H'
%           if solve_for_direction is 'x'
%               xy, xz, yy and yz components are required
%           if solve_for_direction is 'y'
%               xx, xz, yx and yz components are required
%           if solve_for_direction is 'z'
%               xx, xy, yx and yy components are required
%               
%           inv_eps_xx = eps_or_inv_eps{1}, size(inv_eps_xx) == [nxEx,nyEx,nzEx]   TO BE CORRECTED.
%           inv_eps_xy = eps_or_inv_eps{2}, size(inv_eps_xy) == [nxEy,nyEx,nzEy]   .
%           inv_eps_xz = eps_or_inv_eps{3}, size(inv_eps_xz) == [nxEx,nyEx,nzEz]   .
%           inv_eps_yx = eps_or_inv_eps{4}, size(inv_eps_yx) == [nx,ny,nz]
%           eps_or_inv_eps{4} = eps_or_inv_eps{2}, by symmetry of the inverse permittivity tensor
%           inv_eps_yy = eps_or_inv_eps{5}, size(inv_eps_yy) == [nxEy,nyEz,nzEy]   .
%           inv_eps_yz = eps_or_inv_eps{6}, size(inv_eps_yz) == [nxEy,nyEz,nzEz]   .
%   
% Outputs:
%   A:
%       (nxE*nyE*nzE x nxE*nyE*nzE) (when solve_for == 'E') or (nxH*nyH*nzH x nxH*nyH*nzH) (when solve_for == 'H')
%       2d array representing the FDFD 3D wave operator for either E or H.


% Total wave number with grid size normalization
% dx == dy == dz
tot_wave_num = 2*pi/lambda *dx;

% Following Lin, HC., Wang, Z. & Hsu, C.W. Fast multi-source nanophotonic
% simulations using augmented partial factorization. Nat Comput Sci 2, 
% 815–822 (2022): We use outer Kronecker product to augment the diff
% operators to the remaining dimensions.

% 1D first-order finite difference operators
diffx = first_order_fd(nx,xBC,kBx,T);
diffy = first_order_fd(ny,yBC,kBy,W);
diffz = first_order_fd(nz,zBC,kBz,P);

if strcmpi(solve_for, 'E')
    % wave equations, vector form: 
    % [-(d/dx)^2 -(d/dy)^2 -(d/dz)^2 -(omega/c)^2*[eps_xx , eps_xy , eps_xz ; eps_yx , eps_yy , eps_yz ; 
    % eps_zx , eps_zy , eps_zz]]*dx*dy*dx * [Ex Ey Ez] = [bx by bz]

    % (d/dx)^2 == -diffx*diffx
    % Denote "(x)" the Kronecker outer product and "I" identity.
    % Augment (d/dx)^2 operator to 3D: Iz (x) Iy (x) (d/dx)^2, because
    % in the linear ordering of the nx*ny*nz pixels in matrices A and
    % B, z values permute the slowest and x values permute the fastest.
    diffx2 = kron(speye(nz) , kron(speye(ny),-diffx*diffx));
    % Augment (d/dy)^2 operator to 3D: Iz (x) (d/dy)^2 (x) Ix
    diffy2 = kron(speye(nz) , kron(-diffy*diffy,speye(nx)));
    % Augment (d/dz)^2 operator to 3D: (d/dz)^2 (x) Iy (x) Ix
    diffz2 = kron(kron(-diffz*diffz,speye(ny)) , speye(nx));

    if strcmpi(solve_for_direction,'x')
        eps = eps_or_inv_eps{1};
    elseif strcmpi(solve_for_direction,'y')
        eps = eps_or_inv_eps{2};
    else
        eps = eps_or_inv_eps{3};
    end
    
    % Discrete wave operator
    A = -diffx2 -diffy2 -diffz2 -spdiags((tot_wave_num^2)*eps(:),0,nx*ny*nz,nx*ny*nz);

else
    % wave equations, vector form:
    % [d/dy*(inv_eps_zy)*d/dz -d/dy*(inv_eps_zz)*d/dy -d/dz*(inv_eps_yy)*d/dz +d/dz*(inv_eps_yz)*d/dy ,
    % -d/dy*(inv_eps_zx)*d/dz +d/dy*(inv_eps_zz)*d/dx +d/dz*(inv_eps_yx)*d/dz -d/dz*(inv_eps_yz)*d/dx ,  
    % d/dy*(inv_eps_zx)*d/dy -d/dy*(inv_eps_zy)*d/dx -d/dz*(inv_eps_yx)*d/dy +d/dz*(inv_eps_yy)*d/dx ;  
    % d/dz*(inv_eps_xy)*d/dz -d/dz*(inv_eps_xz)*d/dy -d/dx*(inv_eps_zy)*d/dz +d/dx*(inv_eps_zz)*d/dy ,  
    % -d/dz*(inv_eps_xx)*d/dz +d/dz*(inv_eps_xz)*d/dx +d/dx*(inv_eps_zx)*d/dz -d/dx*(inv_eps_zz)*d/dx ,
    % d/dz*(inv_eps_xx)*d/dy -d/dz*(inv_eps_xy)*d/dx -d/dx*(inv_eps_zx)*d/dy +d/dx*(inv_eps_zy)*d/dx ;
    % d/dx*(inv_eps_yy)*d/dz -d/dx*(inv_eps_yz)*d/dy -d/dy*(inv_eps_xy)*d/dz +d/dy*(inv_eps_xz)*d/dy ,
    % -d/dx*(inv_eps_yx)*d/dz +d/dz*(inv_eps_yz)*d/dx +d/dy*(inv_eps_xx)*d/dz -d/dy*(inv_eps_xz)*d/dx ,
    % d/dx*(inv_eps_yx)*d/dy -d/dx*(inv_eps_yy)*d/dx -d/dy*(inv_eps_xx)*d/dy +d/dy*(inv_eps_xy)*d/dx
    % - (omega/c)^2*I]*dx*dy*dz * [Hx Hy Hz] = [bx by bz].
        
    % Set eps_or_inv_eps{4}
    eps_or_inv_eps{4} = eps_or_inv_eps{2};

    if strcmpi(solve_for_direction,'x')
        % Reformat xy, xz, yy and yz components of inverse permittivity
        inv_eps_xy = spdiags(eps_or_inv_eps{2}(:),0,nx*ny*nz,nx*ny*nz);
        inv_eps_xz = spdiags(eps_or_inv_eps{3}(:),0,nx*ny*nz,nx*ny*nz);
        inv_eps_yy = spdiags(eps_or_inv_eps{5}(:),0,nx*ny*nz,nx*ny*nz);
        inv_eps_yz = spdiags(eps_or_inv_eps{6}(:),0,nx*ny*nz,nx*ny*nz);
    elseif strcmpi(solve_for_direction,'y')
        % Reformat xx, xz, yx and yz components of inverse permittivity
        inv_eps_xx = spdiags(eps_or_inv_eps{1}(:),0,nx*ny*nz,nx*ny*nz);
        inv_eps_xz = spdiags(eps_or_inv_eps{3}(:),0,nx*ny*nz,nx*ny*nz);
        inv_eps_yx = spdiags(eps_or_inv_eps{4}(:),0,nx*ny*nz,nx*ny*nz);
        inv_eps_yz = spdiags(eps_or_inv_eps{6}(:),0,nx*ny*nz,nx*ny*nz);
    else
        % Reformat xx, xy, yx and yy components of inverse permittivity
        inv_eps_xx = spdiags(eps_or_inv_eps{1}(:),0,nx*ny*nz,nx*ny*nz);
        inv_eps_xy = spdiags(eps_or_inv_eps{2}(:),0,nx*ny*nz,nx*ny*nz);
        inv_eps_yx = spdiags(eps_or_inv_eps{4}(:),0,nx*ny*nz,nx*ny*nz);
        inv_eps_yy = spdiags(eps_or_inv_eps{5}(:),0,nx*ny*nz,nx*ny*nz);
    end

    % Averaging operators
    avgx = avg(nx,xBC,T,kBx);
    avgy = avg(ny,yBC,W,kBy);
    avgz = avg(nz,zBC,P,kBz);

    % 3D-augmented diff operators
    diffx = kron(speye(nz),kron(speye(ny),diffx));
    diffy = kron(speye(nz),kron(diffy,speye(nx)));
    diffz = kron(kron(diffz,speye(ny)),speye(nx));
        
    % Due to Yee's staggered grids, when diff operators of different 
    % directions are multiplied, they must be averaged in the other
    % directions. [Lin, HC., Wang, Z. & Hsu, C.W. Fast multi-source 
    % nanophotonic simulations using augmented partial factorization. Nat 
    % Comput Sci 2, 815–822 (2022)]
    % 3D-augmented and averaged diff operators
    diffx_avg = kron(avgz,kron(avgy,diffx));
    diffy_avg = kron(avgz,kron(diffy,avgx));
    diffz_avg = kron(kron(diffz,avgy),avgx);

    % Discrete wave operators
    if strcmpi(solve_for_direction, 'x') 
        A = diffx_avg*inv_eps_yy*-diffz_avg -diffx_avg*inv_eps_yz*-diffy_avg -diffy_avg*inv_eps_xy*-diffz_avg +diffy*inv_eps_xz*-diffy -(tot_wave_num^2)*speye(nx*ny*nz);
    elseif strcmpi(solve_for_direction,'y')
        A = -diffx_avg*inv_eps_yx*-diffz_avg +diffz_avg*inv_eps_yz*-diffx_avg +diffy_avg*inv_eps_xx*-diffz_avg -diffy_avg*inv_eps_xz*-diffx_avg -(tot_wave_num^2)*speye(nx*ny*nz);
    else
        % According to Lin, HC., Wang, Z. & Hsu, C.W. Fast multi-source 
        % nanophotonic simulations using augmented partial factorization. 
        % Nat Comput Sci 2, 815–822 (2022), it is better to add
        % -d/dx(inv_eps_yx)d/dy +d/dy(inv_eps_xy)d/dx separately to A, as 
        % most entries cancel:
        A = -diffx*inv_eps_yy*-diffx -diffy*inv_eps_xx*-diffy -(tot_wave_num^2)*speye(nx*ny*nz);
        A = A + (diffx_avg*inv_eps_yx*-diffy_avg + diffy_avg*inv_eps_xy*-diffx_avg);
    end

end

% Determine whether A is symmetric
is_symmetric_A = issymmetric(A);

% TODO: PMLs


end