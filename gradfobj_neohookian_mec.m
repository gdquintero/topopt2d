function [df0dx,gdesl] = gradfobj_neohookian_mec(estr,xfil,p,u,kel,lambdain,lambdaout,opfil,gradxnew,numviz,nvazio,N,coordno,kmax,large)

%Definindo algumas constantes:
ny = estr.ny;
nx = estr.nx;
nelem = (nx-1)*(ny-1);

%Inicializando o vetor gradiente da funcao objetivo:
df0dx = zeros(nelem,1);
gdesl = zeros(nelem,1);
if opfil == 1
   gradxdensobj = zeros(nelem,1);
   gradxdensdesl = zeros(nelem,1);
else
   gradxdensobj = [];
   gradxdensdesl = [];
end

%Calculando o gradiente da funcao objetivo:
if (large == 1)
    if (opfil == 0 || opfil == 2) %opfil = 0 (sem filtro), opfil = 2 (seno hiperbolico)
        for i = 1:nelem
            m1 = nvazio(i);
            lin = linel(m1,nx,ny);
            uel = u(lin);
            crd = [coordno(lin(2)/2,1) coordno(lin(2)/2,2);
                coordno(lin(4)/2,1) coordno(lin(4)/2,2);
                coordno(lin(6)/2,1) coordno(lin(6)/2,2);
                coordno(lin(8)/2,1) coordno(lin(8)/2,2)];
            fint = rhs_elem_neohookian(estr,uel,crd);
            df0dx(m1) = lambdaout(lin)'*(p*(xfil(m1)^(p-1))*fint);
            gdesl(m1) = lambdain(lin)'*(p*(xfil(m1)^(p-1))*fint);
        end
    elseif opfil == 1 %densidades ponderadas
        for i = 1:nelem
            m1 = nvazio(i);
            lin = linel(m1,nx,ny);
            uel = u(lin);
            crd = [coordno(lin(2)/2,1) coordno(lin(2)/2,2);
                coordno(lin(4)/2,1) coordno(lin(4)/2,2);
                coordno(lin(6)/2,1) coordno(lin(6)/2,2);
                coordno(lin(8)/2,1) coordno(lin(8)/2,2)];
            fint = rhs_elem_neohookian(estr,uel,crd);
            gradxdensobj(m1) = lambdaout(lin)'*(p*(xfil(m1)^(p-1))*fint);
            gradxdensdesl(m1) = lambdain(lin)'*(p*(xfil(m1)^(p-1))*fint);
        end
        for i = 1:nelem
            m1 = nvazio(i);
            ind = numviz(m1);
            viz = N(1:ind,m1);
            df0dx(m1) = gradxnew(m1,viz)*gradxdensobj(viz);
            gdesl(m1) = gradxnew(m1,viz)*gradxdensdesl(viz);
        end
    end
else
    if (opfil == 0 || opfil == 2) %opfil = 0 (sem filtro), opfil = 2 (seno hiperbolico)
        for i = 1:kmax
            m1 = nvazio(i);
            lin = linel(m1,nx,ny);
            uel = u(lin);
            df0dx(m1) = lambdaout(lin)'*(p*(xfil(m1)^(p-1))*kel*uel);
            gdesl(m1) = lambdain(lin)'*(p*(xfil(m1)^(p-1))*kel*uel);
        end
    elseif opfil == 1 %densidades ponderadas
        for i = 1:kmax
            m1 = nvazio(i);
            lin = linel(m1,nx,ny);
            uel = u(lin);    
            gradxdensobj(m1) = lambdaout(lin)'*(p*(xfil(m1)^(p-1))*kel*uel);
            gradxdensdesl(m1) = lambdain(lin)'*(p*(xfil(m1)^(p-1))*kel*uel);
        end
        for i = 1:kmax
            m1 = nvazio(i);
            ind = numviz(m1);
            viz = N(1:ind,m1);
            df0dx(m1) = gradxnew(m1,viz)*gradxdensobj(viz);
            gdesl(m1) = gradxnew(m1,viz)*gradxdensdesl(viz);
        end
    end    
end
   






