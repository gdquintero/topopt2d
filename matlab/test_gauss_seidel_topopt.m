function [uc,ugs,uerr] = test_gauss_seidel_topopt

%Initial parameters:
load cantilever_beam_60x30;
str = cantilever_beam_60x30;
x = 0.4*ones(1800,1);
p = 3;
xmin = 0.001;
rmin = 2.5;

%Computing some constants:
n = str.nelem;
nndes = 2*str.nx*str.ny;
nonempty = 1:str.nelem;
opfil = 1;

%Computing the neighborhood of each finite element:
[N,numviz,weigh,wi] = get_neighborhood(str,rmin,[],nonempty,n,opfil);

%Applying the filter:
xfil = apply_filter(str,x,N,weigh,wi,0.0,opfil,[],nonempty,n,numviz,xmin);

%Computing the nodal equivalent forces vector of the structure:
f = vetfneq(str);

%Computing the element stiffness matrix:
kel = elem_stiff(str);
kelc = kel(:);

%Computing the global stiffness matrix, and applying
%the boundary conditions of the structure:
[row,col,inic] = get_rows_cols(str);
pos = get_pos(str,inic);
K = global_stiff(str,xfil,p,row,col,inic,kelc,pos);
linap = put_supports(str);
ap = speye(2*str.nnos);
ap1 = ap(linap,:);
ap2 = ap(:,linap);
K(linap,:) = ap1;
K(:,linap) = ap2;

%Solving the linear system Ku = f using Cholesky factorization:
R = chol(K);
w = (R')\f;
uc = R\w;

%Solving the linear system Ku = f using Gauss-Seidel method:
%u0 = uc + 0.1*rand(nndes,1);
u0 = uc;
ugs = gauss_seidel(K,f,u0,10000,1e-20);

%Computing the error between uc and ugs:
uerr = norm(ugs-uc,inf);