function f0valplp = get_f0valplp(n,s,df0dxplin,bk,breakpoints0,breakpoints1,t0,t1,t2)

f0valplp = 0.0;
for i = 1:n
    if (s(i) < breakpoints0(i))
       f0valplp = f0valplp + df0dxplin(i)*s(i) - 0.5*bk(i)*((t0(i))^2);
    elseif (s(i) >= breakpoints0(i) && s(i) < breakpoints1(i))
           f0valplp = f0valplp + df0dxplin(i+n)*s(i) - 0.5*bk(i)*((t1(i))^2);
    else
       f0valplp = f0valplp + df0dxplin(i+2*n)*s(i) - 0.5*bk(i)*((t2(i))^2);
    end
end

                       
