%AMPLE 1.1: A Material Point Learning Environment
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   27/08/2020
% Description:
% Large deformation elasto-plastic (EP) material point method (MPM) code
% based on an updated Lagrangian (UL) descripition of motion with a 
% quadrilateral background mesh. 
%
%--------------------------------------------------------------------------
% See also:
% SETUPGRID             - analysis specific information
% ELEMMPINFO            - material point-element information
% DETEXTFORCE           - external forces
% DETFDOFS              - mesh unknown degrees of freedom
% LINSOLVE              - linear solver
% DETMPS                - material point stiffness and internal force
% UPDATEMPS             - update material points
% POSTPRO               - post processing function including vtk output
%--------------------------------------------------------------------------
clear;
addpath('constitutive','functions','plotting','setup');        
[lstps,g,mpData,mesh] = setupGrid_collapse;                                          % setup information
NRitMax = 10000; tol = 1e-3;                                                   % Newton Raphson parameters
[nodes,nD] = size(mesh.coord);                                              % number of nodes and dimensions
[nels,nen] = size(mesh.etpl);                                               % number of elements and nodes/element
nDoF = nodes*nD;                                                            % total number of degrees of freedom
nmp  = length(mpData);                                                      % number of material points
lstp = 0;                                                                   % zero loadstep counter (for plotting function)
uvw  = zeros(nDoF,1);                                                       % zeros displacements (for plotting function)
run postPro;                                                                % plotting initial state & mesh
tic
delete output/*
dt_scale = 0.25;
damping_scale = 0.01;
ficticous_mass = 1e3;
damping_factors = zeros(nmp,1);
est_dts = zeros(nmp,1);
est_ghost = zeros(nmp,1);
for i=1:nmp
  E = mpData(i).E;
  nu = mpData(i).nu;
  %rho = mpData(i).mpM/mpData(i).vp;
  rho = ficticous_mass / prod(mesh.h);
  p_mod = E/((1+nu)*(1-nu));
  est_dts(i) = min(mesh.h) * sqrt(rho/p_mod);
  est_ghost(i) = p_mod * 1e-3;
  damping_factors(i) = pi/2 * sqrt(mpData(i).E / (prod(mesh.h)*rho));
end
ghost_factor = min(est_ghost);
dt_est = min(est_dts);
damping = damping_scale*min(damping_factors);
for lstp=1:lstps                                                            % loadstep loop
  fprintf(1,'\n%s %4i %s %4i\n','loadstep ',lstp,' of ',lstps);             % text output to screen (loadstep)
  [mesh,mpData] = elemMPinfo(mesh,mpData);                                  % material point - element information
  fext = detExtForce(nodes,nD,g,mpData);                                    % external force calculation (total)
  fext = fext*lstp/lstps;                                                   % current external force value
  oobf = fext;                                                              % initial out-of-balance force
  fErr = 1;                                                                 % initial error
  frct = zeros(nDoF,1);                                                     % zero the reaction forces
  fdamp = zeros(nDoF,1);                                                    % zero the reaction forces
  uvw_prev_prev  = zeros(nDoF,1);                                                     % zero the displacements
  uvw_prev  = zeros(nDoF,1);                                                     % zero the displacements
  uvw  = zeros(nDoF,1);                                                     % zero the displacements
  vi  = zeros(nDoF,1);                                                      % zero the velocity
  ai  = zeros(nDoF,1);                                                      % zero the accelerations
  Mi  = ficticous_mass*eye(nDoF);                                           % Assemble ficticous mass matrix
  %Mi = detM(mpData,mesh);
  %Mi = diag(sum(Mi));
  R  = zeros(nDoF,1);                                                       % zero current residual
  R_n  = zeros(nDoF,1);                                                     % zero previous residaul
  fd   = detFDoFs(mesh);                                                    % free degrees of freedom
  dt = dt_scale  * dt_est;                                                  % Compute current dt
  NRit = 0;                                                                 % zero the iteration counter
  Kt   = 0;                                                                 % zero global stiffness matrix
  Ke_f = 0;
  Ke_n = 0;
  Ke = 0;
  Kj = ghostPenalty(mesh,mpData);                                         % Ghost penalty stabilisation 
  while (fErr > tol) && (NRit < NRitMax) || (NRit < 2)                      % global equilibrium loop
    %Take the linear solver out
    %We still need internal force, however Ks contribution removed from
    %detmps
    %uvw_prev_prev(fd) = uvw_prev(fd);
    uvw_prev(fd) = uvw(fd);
    fghost = -ghost_factor*Kj*uvw;
    ai(fd)        = Mi(fd,fd)\(fghost(fd));            % calculate acceleration vector
    vi(fd)        = vi(fd)+dt*ai(fd);                                       % update nodal velocity vector
    uvw(fd)  = uvw(fd)+dt*vi(fd);

    [fint,mpData] = detMPs(uvw,mpData);                                     % internal force
    fdamp(fd) = damping * Mi(fd,fd)*vi(fd);                                 %Compute the damping force
    
    ai(fd)        = Mi(fd,fd)\(fext(fd) - fint(fd) - fdamp(fd));            % calculate acceleration vector
    vi(fd)        = vi(fd)+dt*ai(fd);                                       % update nodal velocity vector
    uvw(fd)  = uvw(fd)+dt*vi(fd);                                           % Compute new trial displacments
    uvw(mesh.bc(:,1))=mesh.bc(:,2);
    oobf(fd) = (fext(fd)-fint(fd) + fghost(fd));                                         % out-of-balance force
    oobf(mesh.bc(:,1))=0;
    %Store last residual        
    R_n(fd) = R(fd);
    %Compute current residual
    R(fd) = fext(fd)-fint(fd)-fdamp(fd)+fghost(fd);
    R(mesh.bc(:,1))=0;
    % Compute KE to find critical damping
    Ke = sum(ficticous_mass*sqrt(vi(1:2:end).^2 + vi(2:2:end).^2));
    % Kinetic damping -> zero velocity when we pass through a peak KE
    if (Ke_n > Ke_f) && (Ke<Ke_n)
        fprintf(1,"Peak found\n");
        vi(fd) = vi(fd)*0;
        %interp_value = 
        uvw(fd) = uvw_prev(fd);
        %uvw_prev(fd) = uvw_prev(fd);
        Ke_n = 0;
        Ke_f = 0;
        Ke = 0;
        R(fd) = R(fd)*0;
        R_n(fd) = R_n(fd)*0;
    end
    % Update n-2,n-1 KE
    Ke_f = Ke_n;
    Ke_n = Ke;
    %Compute critical damping numerator and denominator
    num = vi'*(R_n-R);
    denom = abs(dt*Ke);
    if denom > 0
        damp_est = (2*sqrt(abs(num/denom)));
        %Update damping with estimated critical damping
        damping = damping_scale*damp_est;
    end
    oobf_nodal = sqrt(oobf(1:2:end).^2 + oobf(2:2:end).^2)/norm(fext(fd)+eps);
    fErr = norm(oobf(fd))/norm(fext(fd)+eps);                                                                   % normalised oobf error
    NRit = NRit+1;                                                          % increment the NR counter
    fprintf(1,'%s %2i %s %8.3e\n','  iteration ',NRit,' NR error ',fErr);   % text output to screen (NR error)
    if mod(NRit-1,50) == 0
        run postProIter
    end
  end
  mpData = updateMPs(uvw,mpData);                                           % update material points
  run postPro;                                                              % Plotting and post processing 
end
toc
