function lin = linel(str)
%**************************************************************************
%lin = linel(str)
%
%This function returns the indexes of the degrees of freedom associated to
%each element.
%
%INPUT PARAMETER:
%-----------------
%str = Structure with the problem data.
%
%OUTPUT PARAMETER:
%-----------------
%lin = matrix of order 8 x nelem where the j-th column contains
%      the indexes of the degrees of freedom associated to the j-th
%      element.
%**************************************************************************

n = str.nelem;
ny = str.ny;
k = ny-1;
lin = zeros(8,n);
for i = 1:n
    lin(2,i) = 2*((floor((i-1)/k))*ny+mod(i-1,k)+1); 
    lin(4,i) = lin(2,i) + 2*ny;
    lin(6,i) = lin(4,i) + 2;
    lin(8,i) = lin(2,i) + 2;
    lin(1,i) = lin(2,i) - 1;
    lin(3,i) = lin(4,i) - 1;
    lin(5,i) = lin(6,i) - 1;
    lin(7,i) = lin(8,i) - 1;
end