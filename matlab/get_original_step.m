function s = get_original_step(n,splin)

s = zeros(n,1);

for i = 1:n
    s(i) = splin(i) + splin(i+n) + splin(i+2*n);
end

                       
