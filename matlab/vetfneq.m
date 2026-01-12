function fn = vetfneq(str)
%************************************************************************************************************
%fn = vetfneq(str)
%
%This function returns the vector of nodal equivalent forces of the
%structure.
%
%INPUT PARAMETER:
%----------------
%str: Structure with the problem data.
%
%OUTPUT PARAMETER:
%-----------------
%fn = Nodal equivalent forces vector of the structure.
%************************************************************************************************************

fn = zeros(2*str.nnos,1);
for j = 1:length(str.force.node)
    fn(str.force.node(j)*2-1) = fn(str.force.node(j)*2-1) + str.force.Fx(j);
    fn(str.force.node(j)*2) = fn(str.force.node(j)*2) + str.force.Fy(j);
    for i = 1:length(str.supp.node)
        if str.supp.ix(i) == 1
           fn(str.supp.node(i)*2-1) = 0;
        end
        if str.supp.iy(i) == 1
           fn(str.supp.node(i)*2) = 0;
        end
    end     
end 

% fn = zeros(2*str.nnos, length(str.force.node));
% for j = 1:length(str.force.node)
%     fn(str.force.node(j)*2-1,j) = fn(str.force.node(j)*2-1,j) + str.force.Fx(j);
%     fn(str.force.node(j)*2,j) = fn(str.force.node(j)*2,j) + str.force.Fy(j);
%     for i = 1:length(str.supp.node)
%         if str.supp.ix(i) == 1
%            fn(str.supp.node(i)*2-1,j) = 0;
%         end
%         if str.supp.iy(i) == 1
%            fn(str.supp.node(i)*2,j) = 0;
%         end
%     end     
% end 



