function [Kghost] = bulkghost(mesh)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
mesh = detBoundElem(mesh);                                                  % determine boundary element information 
bE = mesh.bE;
nbE = size(bE,1);
ngp = 4;
coord = mesh.coord;
[nodes,nD] = size(coord);
etpl = mesh.etpl;
nn = size(etpl,2);
neDoF=(nn*nD)^2;
krow=zeros(neDoF*nbE,1); kcol=krow; kval=krow;
[wgp,dNr,xsi,eta]=dershapefunc2D(ngp,nn);
for el = 1 : nbE
    nIN  = nodesForMP(etpl,bE(el)).';                                     % unique list of nodes associated with elements
    ed  = repmat((nIN-1)*nD,nD,1)+repmat((1:nD).',1,nn);                 % degrees of freedom of nodes (matrix form)
    ed  = reshape(ed,1,nn*nD);                                              % degrees of freedom of nodes (vector form)
    JT  = dNr*coord(etpl(bE(el),:),:);
    if nD == 1
        h = coord(etpl(bE(el),4),2)-coord(etpl(bE(el),1),1);
    elseif nD == 2
        hx = coord(etpl(bE(el),4),1)-coord(etpl(bE(el),1),1);
        hy = coord(etpl(bE(el),2),2)-coord(etpl(bE(el),1),2);
        h = hx*hy;
    else
        % hx =
        % hy =
        % hz = 
        % h = hx*hy*hz;
    end
    kghst = zeros(nn*nD);
    for gp = 1 : ngp 
        indx = nD*gp-(nD-1:-1:0);
        detJ = det(JT(indx,:));
        [N,N0] = shapefunc2D(xsi,eta,gp,nn);
        NN  = zeros(nD,nD*nn);                                                  % zero the shape function matrix
        NN0 = zeros(nD,nD*nn);
        if nD==1                                                                % 1D case
            NN  = N;                                                            % shape function matrix
            NN0 = N0;
        elseif nD==2                                                            % 2D case (plane strain & stress)
            NN(1,1:nD:end)  = N;                                                  % shape function matrix
            NN(2,2:nD:end)  = N;
            NN0(1,1:nD:end) = N0;                                                  % shape function matrix
            NN0(2,2:nD:end) = N0;
        else                                                                    % 3D case
            NN(1,1:nD:end)=N;                                                  % shape function matrix
            NN(2,2:nD:end)=N;
            NN(3,3:nD:end)=N;
            NN0(1,1:nD:end)=N0;                                                  % shape function matrix
            NN0(2,2:nD:end)=N0;
            NN0(3,3:nD:end)=N0;
        end
       kghst = kghst + detJ*wgp(gp)*(NN.'*(NN-NN0))/(h); 
    end
    krow((el-1)*neDoF+1:el*neDoF)=reshape(ed.'*ones(1,nn*nD),neDoF,1);
    kcol((el-1)*neDoF+1:el*neDoF)=reshape(ones(nn*nD,1)*ed  ,neDoF,1);
    kval((el-1)*neDoF+1:el*neDoF)=reshape(kghst,neDoF,1);
end
Kghost=sparse(krow,kcol,kval,nodes*nD,nodes*nD);
end

function [wp,dNr,xsi,eta]=dershapefunc2D(ngp,nen)
if ngp==4
  g2=1/sqrt(3);
  gp(:,1)=[-1 -1 1 1].'*g2;
  gp(:,2)=[-1 1 1 -1].'*g2;
  wp=ones(4,1);
else
  g2=sqrt(3/5);
  gp(:,1)=[-1 -1 -1  0  0  0  1  1  1].'*g2;
  gp(:,2)=[-1  1  0 -1  1  0 -1  1  0].'*g2;
  wp1=[5/9 5/9 5/9 8/9 8/9 8/9 5/9 5/9 5/9].';
  wp2=[5/9 5/9 8/9 5/9 5/9 8/9 5/9 5/9 8/9].';
  wp=wp1.*wp2;  
end
xsi=gp(:,1); eta=gp(:,2); r2=ngp*2;
if nen==8
  dNr(1:2:r2,1)=1/4*(2*xsi+eta).*(1-eta);
  dNr(1:2:r2,2)=-1/2*(1-eta.^2);
  dNr(1:2:r2,3)=1/4*(2*xsi-eta).*(1+eta); 
  dNr(1:2:r2,4)=-xsi.*(1+eta);
  dNr(1:2:r2,5)=1/4*(1+eta).*(2*xsi+eta);
  dNr(1:2:r2,6)=1/2*(1-eta.^2);
  dNr(1:2:r2,7)=1/4*(1-eta).*(2*xsi-eta); 
  dNr(1:2:r2,8)=-xsi.*(1-eta);
  dNr(2:2:r2+1,1)=1/4*(1-xsi).*(xsi+2*eta);
  dNr(2:2:r2+1,2)=-eta.*(1-xsi);
  dNr(2:2:r2+1,3)=1/4*(1-xsi).*(2*eta-xsi); 
  dNr(2:2:r2+1,4)=1/2*(1-xsi.^2);
  dNr(2:2:r2+1,5)=1/4*(1+xsi).*(2*eta+xsi);
  dNr(2:2:r2+1,6)=-eta.*(1+xsi);
  dNr(2:2:r2+1,7)=1/4*(1+xsi).*(2*eta-xsi);
  dNr(2:2:r2+1,8)=-1/2*(1-xsi.^2);
else
  dNr(1:2:r2,1)=-1/4*(1-eta);
  dNr(1:2:r2,2)=-1/4*(1+eta);
  dNr(1:2:r2,3)= 1/4*(1+eta);
  dNr(1:2:r2,4)= 1/4*(1-eta);
  dNr(2:2:r2+1,1)=-1/4*(1-xsi);
  dNr(2:2:r2+1,2)= 1/4*(1-xsi);
  dNr(2:2:r2+1,3)= 1/4*(1+xsi);
  dNr(2:2:r2+1,4)=-1/4*(1+xsi);
end
end

function [N,N0]=shapefunc2D(xsi,eta,gp,nen)
N=zeros(nen,1);
N0=ones(nen,1)/4; 
if nen==8
  N(1)=1/4*(1-xsi(gp))*(1-eta(gp))*(-xsi(gp)-eta(gp)-1);
  N(2)=1/2*(1-xsi(gp))*(1-eta(gp)^2);
  N(3)=1/4*(1-xsi(gp))*(1+eta(gp))*(-xsi(gp)+eta(gp)-1);
  N(4)=1/2*(1-xsi(gp)^2)*(1+eta(gp));
  N(5)=1/4*(1+xsi(gp))*(1+eta(gp))*(xsi(gp)+eta(gp)-1);
  N(6)=1/2*(1+xsi(gp))*(1-eta(gp)^2);
  N(7)=1/4*(1+xsi(gp))*(1-eta(gp))*(xsi(gp)-eta(gp)-1);
  N(8)=1/2*(1-xsi(gp)^2)*(1-eta(gp));
else
  N(1)=(1-xsi(gp))*(1-eta(gp))/4;
  N(2)=(1-xsi(gp))*(1+eta(gp))/4;
  N(3)=(1+xsi(gp))*(1+eta(gp))/4;
  N(4)=(1+xsi(gp))*(1-eta(gp))/4;
end
end
