function gradxnew = graddens(str,N,weigh,wi,nonempty,kmax,numviz)
%************************************************************************************************************
%gradxnew = graddens(str,N,weigh,wi,nonempty,kmax,numviz)
%
%This function returns the matrix where each line contains the coefficients associated to
%the chain rule for the objective function and volume constraint gradient vectors, 
%when the mean density filter is used. 
%
%INPUT PARAMETERS:
%-----------------
%str: Structure with the problem data.
%N: Matrix where each column contains the neighbors of each finite element.
%weigh: Matrix where each row contains the weights of each finite element,
%       associated to the mean density filter or to the sinh filter.
%wi: Vector where each element is the sum of the finite element weights,
%    associated to the mean density filter or to the sinh filter.
%nonempty: Vector with the indexes of the elements of the domain that are nonempty.
%kmax: Number of nonempty elements.
%numviz: Vector that contains the number of neighbors of each finite element.
%
%OUTPUT PARAMETERS:
%------------------
%gradxnew: Matrix where each line contains the coefficients associated to
%          the chain rule for the objective function and volume constraint derivatives,
%          when the mean density filter is used. 
%************************************************************************************************************

gradxnew = sparse(str.nelem,str.nelem);
for k = 1:kmax
    m1 = nonempty(k);
    viz = numviz(m1);
    ind = N(1:viz,m1);
    gradxnew(m1,ind) = weigh(m1,1:viz)'./wi(ind);
end
       

