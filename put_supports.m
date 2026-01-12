function linap = put_supports(str)
%**************************************************************************
%linap = put_supports(str)
%
%This function returns the vector with the indexes associated to the
%nodal degrees of freedom where the nodal supports are applied.
%
%INPUT PARAMETER:
%----------------
%str: Structure with the problem data.
%
%OUTPUT PARAMETER:
%-----------------
%linap = Vector with the indexes associated to the
%        nodal degrees of freedom where the nodal supports are applied.
%**************************************************************************

linap = zeros(2*length(str.supp.node),1);
nlap = 0;
for i = 1:length(str.supp.node)
    if str.supp.ix(i) == 1
       nlap = nlap + 1;
       linap(nlap) = str.supp.node(i)*2-1;
    end
    if str.supp.iy(i) == 1
       nlap = nlap + 1;
       linap(nlap) = str.supp.node(i)*2;
    end
end
linap = linap(1:nlap,1);












