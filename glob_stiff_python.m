function K = glob_stiff_python(str,x,p,kel,lin,linap)

nndes = 2*str.nx*str.ny;
K = zeros(nndes,nndes);
n = str.nelem;

for k = 1:n   
    for i = 1:8
        for j = 1:8
            K(lin(i,k),lin(j,k)) = K(lin(i,k),lin(j,k)) + (x(k)^p)*kel(i,j);
        end
    end    
end

K = sparse(K);

nsupp = length(linap);

for i = 1:nsupp  
    for j = 1:nndes       
        if (linap(i) == j)       
           K(linap(i),j) = 1.0;      
        else       
           K(linap(i),j) = 0.0;
           K(j,linap(i)) = 0.0;      
        end   
    end  
end

% for i = 1:nsupp
%     for j = 1:nsupp
%         if (linap(i) == linap(j))
%            K(linap(i),linap(j)) = 1.0;
%         else
%            K(linap(i),linap(j)) = 0.0;
%         end
%     end
% end
        
     