function [xdil,gradxdil] = dilation_filter(str,x,beta,N,nonempty,kmax,numviz)
%************************************************************************************************************
%[xdil,gradxdil] = dilation_filter(str,x,beta,N,nonempty,kmax,numviz)
%
%This function applies the dilation filter to the original element
%densities.
%
%INPUT PARAMETERS:
%-----------------
%str: Structure with the problem data.
%x: Original element densities vector.
%beta: Parameter used in dilation filter.
%N: Matrix where each column contains the neighbors of each finite element. 
%nonempty: Vector with the indexes of the finite elements that are nonempty.
%kmax: Number of nonempty elements.
%numviz: Vector that contains the number of neighbors of each finite element.
%
%OUTPUT PARAMETERS:
%------------------
%xdil = Filtered element densities vector.
%gradxdil: Matrix where each line contains the coefficients associated to
%          the chain rule for the objective function and volume constraint derivatives,
%          when the dilation filter is used. 
%************************************************************************************************************

xdil = zeros(str.nelem,1); 
gradxdil = sparse(str.nelem,str.nelem);
for i = 1:kmax
    m1 = nonempty(i);
    viz = numviz(m1);
    ind = N(1:viz,m1);
    soma1 = sum(exp(beta*x(ind)));
    soma2 = sum(ones(viz,1));
    xdil(m1) = log(soma1/soma2)/beta;
    gradxdil(m1,ind) = exp(beta*x(ind))/soma1;
end

