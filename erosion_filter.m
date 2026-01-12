function [xer, gradxer] = erosion_filter(str,x,beta,N,nonempty,kmax,numviz)
%************************************************************************************************************
%[xer, gradxer] = erosion_filter(str,x,beta,N,nonempty,kmax,numviz)
%
%This function applies the erosion filter to the original element
%densities.
%
%INPUT PARAMETERS:
%-----------------
%str: Structure with the problem data.
%x: Original element densities vector.
%beta: Parameter used in the erosion filter.
%N: Matrix where each column contains the neighbors of each finite element. 
%nonempty: Vector with the indexes of finite elements that are nonempty.
%kmax: Number of nonempty elements.
%numviz: Vector that contains the number of neighbors of each finite element.
%
%OUTPUT PARAMETERS:
%------------------
%xer = Filtered element densities vector.
%gradxer: Matrix where each line contains the coefficients associated to
%         the chain rule for the objective function and volume constraint derivatives,
%         when the erosion filter is used. 
%************************************************************************************************************

xer = zeros(str.nelem,1); 
gradxer = sparse(str.nelem,str.nelem);
for i = 1:kmax
    m1 = nonempty(i);
    viz = numviz(m1);
    ind = N(1:viz,m1);
    soma1 = sum(exp(beta*(1-x(ind))));
    soma2 = sum(ones(viz,1));
    xer(m1) = 1 - log(soma1/soma2)/beta;
    gradxer(m1,ind) = exp(beta*(1-x(ind)))/soma1;
end



