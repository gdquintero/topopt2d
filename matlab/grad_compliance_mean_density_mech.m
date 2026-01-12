function [df0dx,gdesl] = grad_compliance_mean_density_mech(str,xfil,p,u,kel,gradxnew,numviz,nonempty,kmax,N,lin,lambdain,lambdaout)

%Definindo algumas constantes:
nelem = str.nelem;

%Inicializando o vetor gradiente da funcao objetivo:
df0dx = zeros(nelem,1);
gdesl = zeros(nelem,1);
gradxdensobj = zeros(nelem,1);
gradxdensdesl = zeros(nelem,1);

%Calculando o gradiente da funcao objetivo:
for i = 1:kmax
    m1 = nonempty(i);
    lin2 = lin(:,m1);
    uel = u(lin2);
    gradxdensobj(m1) = lambdaout(lin2)'*(p*(xfil(m1)^(p-1))*kel*uel);
    gradxdensdesl(m1) = lambdain(lin2)'*(p*(xfil(m1)^(p-1))*kel*uel);
end

for i = 1:kmax
    m1 = nonempty(i);
    ind = numviz(m1);
    viz = N(1:ind,m1);
    df0dx(m1) = gradxnew(m1,viz)*gradxdensobj(viz);
    gdesl(m1) = gradxnew(m1,viz)*gradxdensdesl(viz);
end

%df0dx = -df0dx;
   
   






