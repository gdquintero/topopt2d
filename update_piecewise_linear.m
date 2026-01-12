function [breakpoints0,breakpoints1,t0,t1,t2,df0dxplin,dfdxplin] = update_piecewise_linear(n,gamma,x,xmin,df0dx,bk,dfdx)

t0 = zeros(n,1);
t1 = t0;
t2 = t0;
breakpoints0 = t0;
breakpoints1 = t0;
df0dxplin = zeros(3*n,1);
dfdxplin = df0dxplin';

for i = 1:n 
    
    %Unconstrained minimizer of the quadratic approximation of the
    %objective function:
    if (bk(i) > 0.0)
       v = -df0dx(i)/bk(i);
       v2 = 2*v;
    else
       v2 = 1e30;
    end
    
    %Bounds for the interpolation region:
    low1 = min(0.5*gamma, x(i)-xmin);
    upp1 = min(0.5*gamma, 1-x(i));
    low2 = min(gamma, x(i)-xmin);
    upp2 = min(gamma, 1-x(i));
    
    %Updating the bounds for the interpolation region:
    if (bk(i) > 0.0)
       if (v2 > upp1)
          if (v2 < upp2)
             upp = v2;
          else
             upp = upp2;
          end
          low = low1;
       elseif (v2 < -low1)
           if (v2 > -low2)
              low = -v2;
           else
              low = low2;
           end
           upp = upp1;
       else
          low = low1;
          upp = upp1;
       end
    else
       low = low1;
       upp = upp1;
    end
    
    %Getting the interpolation points of the separable quadratic term of
    %the separable quadratic approximation of the objective function:
    if ((abs(x(i)-xmin) <= 1e-6) && (abs(1-x(i)) > 1e-6))
       t0(i) = 0.0;
       t1(i) = 0.3*upp;
       t2(i) = upp;
    elseif ((abs(1-x(i)) <= 1e-6) && (abs(x(i)-xmin) > 1e-6))
           t0(i) = -low;
           t1(i) = -0.3*low;
           t2(i) = 0.0;
    elseif ((abs(x(i)-xmin) > 1e-6) && (abs(1-x(i)) > 1e-6))
           t0(i) = -low;
           t1(i) = 0.0;
           t2(i) = upp;
    else
           t0(i) = 0.0;
           t1(i) = 0.0;
           t2(i) = 0.0;
    end
    
    %Getting the breakpoints of the piecewise linear approximation of the
    %quadratic term of the separable quadratic approximation of the
    %objective function:
    breakpoints0(i) = 0.5*(t0(i)+t1(i));
    breakpoints1(i) = 0.5*(t1(i)+t2(i));
    
    %Getting the slopes of the piecewise linear approximation of the
    %quadratic approximation of the objective function:
    df0dxplin(i) = df0dx(i) + bk(i)*t0(i);
    df0dxplin(i+n) = df0dx(i) + bk(i)*t1(i);
    df0dxplin(i+2*n) = df0dx(i) + bk(i)*t2(i);
    
    %Getting the vector used in the linear approximatino of the volume
    %constraint used in the piecewise linear subproblem:
    dfdxplin(i) = dfdx(i);
    dfdxplin(i+n) = dfdx(i);
    dfdxplin(i+2*n) = dfdx(i);
         
end   



    

