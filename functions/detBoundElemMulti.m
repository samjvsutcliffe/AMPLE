function [mesh] = detBoundElemMulti(mesh,mpData)

%Boundary element determination
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   11/06/2021
% Description:
% Function to determine the bounadry elements and the faces attached to the
% boundary elements. 
%
%--------------------------------------------------------------------------
% [mesh] = DETBOUNDELEMMULTI(mesh)
%--------------------------------------------------------------------------
% Input(s):
% mesh   - mesh structured array. Function requires: 
%           - coord : coordinates of the grid nodes (nodes,nD)
%           - etpl  : element topology (nels,nen) 
%           - ftpl  : face-element interactions
% mpData - material point structured array.  Function requires:
%           - bdy   : intger to identify different bodies
%           - lp    : domain lengths 
%           - mpC   : material point coordinates
%--------------------------------------------------------------------------
% Ouput(s);
% mesh   - mesh structured array. Function modifies:
%           - bE    : bounadry elements flag
%           - bF    : boundary face flag 
%--------------------------------------------------------------------------
% See also:
%
%--------------------------------------------------------------------------
    
[nels,~] = size(mesh.etpl);                                                 % number of elements
[~,nD]   = size(mesh.coord);                                                % number of dimensions
nmp = length(mpData);                                                       % number of material points 

%bdy   = [mpData.bdy];                                                       % material point body information
nBdy  = 1;%max(bdy);                                                           % number of bodies in the analysis 

mpCall = reshape([mpData.mpC],nD,nmp).';                                    % material point coordinates
lpAll  = reshape([mpData.lp],nD,nmp).';                                     % material point domain lengths

ftpl  = mesh.ftpl;                                                          % face topology
nface = length(ftpl);                                                       % number of faces

a  = false(nface,1);                                                        % zero boundary face flag
bEf = false(nels,1);                                                        % bounday element storage
for j = 1:nBdy
    mpC  = mpCall(:,:);                                                % current body material points
    lp   = lpAll(:,:);                                                 % domain lengths for current body
    eInA = zeros(nels,1);                                                   % zero elements taking part in the analysis
    nmp  = length(mpC);
    for mp = 1:nmp
        eIN  = elemForMP(mesh,mpC(mp,:),lp(mp,:));                          % elements connected to the material point
        eInA(eIN) = 1;                                                      % identify elements in the analysis
    end
    A  = eInA(ftpl);                                                        % elements attached to faces in/out
    B  = ne(A(:,1),A(:,2));                                                 % faces between in/out elements
    bF = (1:nface).';                                                       % list, 1 to number of faces
    bF = bF(B);                                                             % face numbers of faces between in and out elements
    bE = zeros(size(bF)); 
    for i = 1:length(bF)                                                    % loop over in/out faces
        fE = ftpl(bF(i),:);                                                 % elements attached to face
        if eInA(fE(1))==1
            bEf(fE(1)) = true;                                              % set boundary element = element 1
            bE(i) = fE(1);
        else
            bEf(fE(2)) = true;                                              % set boundary element = element 2
            bE(i) = fE(2); 
        end
    end
    bE = bE(bE>0);
    bE = unique(bE); 
    for i = 1:nface
        fE = ftpl(i,:);                                                     % elements attached to face
        b     = sum(eInA(fE))==2;                                           % both elements in analysis (yes, b1=1; no, b=0)
        a(i)  = max(a(i),(sum(bE==fE(1))+sum(bE==fE(2)))*b)>=1;             % is face attached to boundary element
    end
end
mesh.bE = bEf;                                                              % store bounadry elements flags
mesh.bF = a;                                                                % store bounadry face flag
end
