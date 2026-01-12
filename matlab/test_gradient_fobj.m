function [u,fobj,df0dx,vol,dfdx] = test_gradient_fobj

%Initial parameters:
load str3x3;
str = str3x3;
x = 0.4*ones(9,1);
p = 3;
xmin = 0.001;
rmin = 1.5;
frmax = 0.4;


%Computing some constants:
n = str.nelem;
nonempty = 1:str.nelem;
opfil = 1;
V = str.b*str.h*str.e;
Vmax = V*frmax;
velem = (V/n)*ones(n,1);

%Computing the neighborhood of each finite element:
[N,numviz,weigh,wi,gradxnew] = get_neighborhood(str,rmin,[],nonempty,n,opfil);

%Applying the filter:
xfil = apply_filter(str,x,N,weigh,wi,0.0,opfil,[],nonempty,n,numviz,xmin);

%Computing the nodal equivalent forces vector of the structure:
f = vetfneq(str);

%Computing the element stiffness matrix:
kel = elem_stiff(str);
kelc = kel(:);

%Computing the matrix that contains the indexes of the degrees of freedom associated to
%each element:
lin = linel(str);

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

%Solving the linear system Ku = f:
u = K\f;

%Computing the objective function:
fobj = f'*u;

%Computing the gradient of the objective function:
df0dx = grad_compliance_mean_density(str,xfil,p,u,kel,gradxnew,numviz,nonempty,n,N,lin);

%Computing the objective function:
vol = (velem'*xfil-Vmax)/Vmax;

%Computing the gradient of the volume constraint:
dfdx = grad_volume_mean_density(str,velem,N,weigh,wi,nonempty,n,Vmax,numviz);

