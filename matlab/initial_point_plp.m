function sp0 = initial_point_plp(n,s,breakpoints0,breakpoints1)

sp0 = zeros(3*n,1);

for i = 1:n
    
    if (s(i) <= breakpoints0(i))
       sp0(i) = s(i);
       sp0(i+n) = 0.0;
       sp0(i+2*n) = 0.0;
    elseif (s(i) > breakpoints0(i) && s(i) <= breakpoints1(i))
           sp0(i) = breakpoints0(i);
           sp0(i+n) = s(i)-breakpoints0(i);
           sp0(i+2*n) = 0.0;
    else
       sp0(i) = breakpoints0(i);
       sp0(i+n) = breakpoints1(i)-breakpoints0(i);
       sp0(i+2*n) = s(i)-breakpoints1(i);
    end
    
end


