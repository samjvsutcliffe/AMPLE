function [M] = detM(mpData,mesh)

%Global mass matrix determination
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   23/01/2019
% Description:
% Function to determine the global mass matrix based on material point 
% information.
% 
%--------------------------------------------------------------------------
% [M] = DETM(mpData)
%--------------------------------------------------------------------------
% Input(s):
% mpData - material point structured array. The following fields are
%          required by the function:
%           - Svp   : basis functions (1,nn)
%           - nIN   : background mesh nodes associated with the MP (1,nn)
%           - mpM   : material point mass(1)
%           - SMe   : number stiffness matrix entries
% mesh   - mesh structured array.  The following fields are required by the
%          function:
%           - coord : nodal coordinates (nodes,nD)
%--------------------------------------------------------------------------
% Ouput(s);
% M      - consistent global mass matrix, sparse (nDoF,nDoF)
%--------------------------------------------------------------------------
% See also:
% 
%--------------------------------------------------------------------------

nmp   = length(mpData);                                                     % number of material points
npCnt = 0;                                                                  % counter for the number of entries in M
tnSMe = sum([mpData.nSMe]);                                                 % total number of mass matrix entries
mrow  = zeros(tnSMe,1); mcol=mrow; mval=mrow;                               % zero the stiffness information

[nodes,nD] = size(mesh.coord);                                              % number of nodes and dimensions

for mp=1:nmp                                                                % material point loop
    nIN = mpData(mp).nIN;                                                   % nodes associated with the material point 
    Svp = mpData(mp).Svp;                                                   % basis functions
    mpM = mpData(mp).mpM;                                                   % material point mass
    nn  = size(Svp,2);                                                      % no. dimensions & no. nodes
    N   = zeros(nD,nD*nn);                                                  % zero the shape function matrix
    if nD==1                                                                % 1D case
        N = Svp;                                                            % shape function matrix
    elseif nD==2                                                            % 2D case (plane strain & stress)
        N(1,1:nD:end)=Svp;                                                  % shape function matrix
        N(2,2:nD:end)=Svp;
    else                                                                    % 3D case
        N(1,1:nD:end)=Svp;                                                  % shape function matrix
        N(2,2:nD:end)=Svp;
        N(3,3:nD:end)=Svp;
    end
    Mp = N.'*N*mpM;                                                         % material point mass matrix
    
    ed  = repmat((nIN-1)*nD,nD,1)+repmat((1:nD).',1,nn);                    % degrees of freedom of nodes (matrix form)
    ed  = reshape(ed,1,nn*nD);                                              % degrees of freedom of nodes (vector form)
    
    npDoF = (size(ed,1)*size(ed,2))^2;                                      % no. entries in kp
    nnDoF = size(ed,1)*size(ed,2);                                          % no. DoF in kp                        
    mrow(npCnt+1:npCnt+npDoF) = repmat(ed.',nnDoF,1);                       % row position storage
    mcol(npCnt+1:npCnt+npDoF) = repmat(ed  ,nnDoF,1);                       % column position storage
    mval(npCnt+1:npCnt+npDoF) = Mp;                                         % mass value storage
    npCnt = npCnt+npDoF;                                                    % number of entries in M
end
nDoF = nodes*nD;                                                            % number of degrees of freedom
M = sparse(mrow,mcol,mval,nDoF,nDoF);                                       % form the global mass matrix

end