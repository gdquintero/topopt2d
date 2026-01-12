function [xstar,iter,itint,f0val,f0val2,beta,time] = linseq_mech(str1,str2,str3,x,xmin,p,delta,frmax,tol,tolcount,maxit,...
                                                                 prnt,opfil,rmin,beta,hole,empty)
%***********************************************************************************************************************
%[xstar,iter,itint,f0val,f0val2,beta,time] = linseq_mech(str1,str2,str3,x,xmin,p,delta,frmax,tol,tolcount,maxit,...
%                                                        prnt,opfil,rmin,beta,hole,empty)
%
%This function solves the topology optimization for compliant mechanisms, 
%following the Nishiwaki et al. formulation* and using a globally convergent version of 
%the Sequential Linear Programming (SLP) algorithm**, proposed by Gomes and Senne.
%
%*Nishiwaki, S.; Frecker, M. I.; Seungjae, M.; Kikuchi, N.; Topology
%optimization of compliant mechanisms using the homogenization method.
%International Journal for Numerical Methods in Engineering, 42, 1998, p.
%535-559.
%
%**Gomes, F. A. M.; Senne, T. A.; An SLP algorithm and its application to
%topology optimization. Computational and Applied Mathematics, Vol. 30, N. 1, pp. 53-89, 2011.
%
%
%INPUT PARAMETERS:
%-----------------
%str1: Structure with the problem data, considering only the applied traction  
%      and the respective boundary conditions in the first load case (Nishiwaki et al. formulation*).
%str2: Structure with the problem data, considering only the dummy load  
%      and the respective boundary conditions in the first load case (Nishiwaki et al. formulation*).
%str3: Structure with the problem data, considering only the reaction force
%      and the respective boundary conditions in the second load case (Nishiwaki et al. formulation*).
%x: Initial point for the SLP algorithm (initial element densities).
%xmin: Small positive value that the element densities can assume, to avoid numerical instabilities.
%p: Penalty parameter of the SIMP model for the element densities,
%   or penalty parameter of the sinh filter for the volume constraint. 
%delta: Initial trust region radius for the SLP algorithm.
%frmax: Volume fraction of the rectangular domain that the optimal
%       structure can contain. 
%tol: Minimum value for the step norm obtained by the SLP algorithm (stopping criterion).
%tolcount: Number of times that the step norm obtained by the SLP algorithm
%          needs to reach a value less or equal than the value of the "tol"
%          parameter before stopping the algorithm (stopping criterion).
%maxit: Maximum number of iterations that the algorithm can perform (stopping criterion).
%prnt: Number of times that the outer iteration data and the topology will be
%      printed on the screen. 
%opfil: Parameter that indicates what filter is used.
%       opfil = 0 (no filter)
%       opfil = 1 (mean density filter)
%       opfil = 2 (dilation filter)
%       opfil = 3 (erosion filter)
%       opfil = 4 (sinh filter)
%rmin: Filter radius.
%beta: Initial value of the parameter used in the dilation and erosion
%      filters. If another filter is used, set beta = 0. 
%hole: Parameter that indicates if the rectangular domain contains some hole.
%      If hole = 1, the domain has a rectangular hole.
%      If hole = 0, the domain doesn't have any hole. 
%empty: Vector that contains the indexes of the elements that belongs to
%       the rectangular hole. If the domain doesn't have any hole, set
%       empty = [].
%
%OUTPUT PARAMETERS:
%------------------
%xstar: Approximated optimal solution of the topology optimization problem.
%iter: Total number of outer iterations performed by the SLP algorithm to find
%      the approximated optimal solution of the topology optimization
%      problem.
%itint: Total number of inner iterations performed by the SLP algorithm to find
%       the approximated optimal solution of the topology optimization
%       problem.
%f0val: Approximated optimal value of the objective function.
%f0val2: Approximated optimal value of the objective function, obtained setting p = 1 and without using any filter
%        at the approximated optimal solution.
%beta: Parameter used in the dilation and erosion filters. 
%time: Total time taken by the SLP algorithm to find the approximated
%      optimal solution of the topology optimization problem. 
%***********************************************************************************************************************

%Starting to count the time:
tic; 

%Computing some constants:
str = str1; 
n = str.nelem;
sinhp = sinh(p);
spinv = 1/sinhp;
V = str.b*str.h*str.e;
velem = (V/str.nelem)*ones(str.nelem,1);
nonempty = 1:str.nelem;
if hole == 1
   nonempty(empty) = [];
   kmax = length(nonempty);
   velem(empty) = 0;
   Vdom = sum(velem);
   x(empty) = xmin;
else
   kmax = n;
   Vdom = V;
end
Vmax = frmax*Vdom;
csnhp = spinv/Vmax;
lpoptions = optimset('Display','off');
lpoptions2 = optimset('Display','off','LargeScale','off','Simplex','on');
snorm = 1e20;
iter = 0;
itint = 0;
thold1 = 1;
thold2 = 1;
thmax = 1;   
c = [zeros(n,1); 1; 1];
count = 0;
delta0 = delta;

%Computing the neighborhood of each finite element:
if opfil ~= 0
   [N,numviz,weigh,wi,gradxnew] = get_neighborhood(str,rmin,empty,nonempty,kmax,opfil);
else
   N = [];
   numviz = [];
   weigh = [];
   wi = [];
   gradxnew = [];
end

%Applying the filter:
[xfil,gradxdil,gradxer] = apply_filter(str,x,N,weigh,wi,beta,opfil,empty,nonempty,kmax,numviz,xmin);

%Computing the nodal equivalent forces vectors of the compliant mechanism:
fa = vetfneq(str1);
fb = vetfneq(str2);
fc = vetfneq(str3);

%Computing the element stiffness matrix:
kel = elem_stiff(str);
kelc = kel(:);

%Choosing the appropriate value of the penalty parameter to compute
%the objective function value. If the sinh filter is used, this penalty
%parameter is always equal to 1. 
if opfil == 4
   pp = 1;
else
   pp = p;
end

%Computing the global stiffness matrix, and applying
%the boundary conditions of the compliant mechanism:
[row,col,inic] = get_rows_cols(str);
Ka = global_stiff(str,xfil,pp,row,col,inic,kelc);
Kc = Ka;
linap1 = put_supports(str1);
linap2 = put_supports(str3);
ap = speye(2*str.nnos);
ap11 = ap(linap1,:);
ap12 = ap(:,linap1);
ap21 = ap(linap2,:);
ap22 = ap(:,linap2);
Ka(linap1,:) = ap11;
Ka(:,linap1) = ap12;
Kc(linap2,:) = ap21;
Kc(:,linap2) = ap22;
perm_a = symamd(Ka);
perm_c = symamd(Kc);
fa2 = fa(perm_a);
fb2 = fb(perm_a);
fc2 = fc(perm_c);

%Computing the initial value of the objective function and of the volume
%constraint:
[f0val,b,xnew,xvol,ua,ub,uc,ubKua,ucKuc] = obj_func_mech(str,xfil,p,fa,fc,...
fa2,fb2,fc2,Ka,Kc,perm_a,perm_c,velem,Vmax,opfil,N,weigh,wi,empty,nonempty,kmax,numviz,sinhp,xmin);
f0valold = f0val;

%Computing the gradient vectors of the objective function and of the volume
%constraint:
if opfil == 0
   A = velem'/Vmax;
   df0dx = grad_func_mech(n,str,xfil,xvol,p,ua,ub,uc,kel,velem,Vmax,gradxnew,gradxdil,gradxer,...
                                                 numviz,nonempty,kmax,N,opfil,ubKua,ucKuc,csnhp);
elseif opfil == 1
       A = grad_volume_mean_density(str,velem,N,weigh,wi,nonempty,kmax,Vmax,numviz)';
       df0dx = grad_func_mech(n,str,xfil,xvol,p,ua,ub,uc,kel,velem,Vmax,gradxnew,gradxdil,gradxer,...
                                                     numviz,nonempty,kmax,N,opfil,ubKua,ucKuc,csnhp);
else
   [df0dx,A] = grad_func_mech(n,str,xfil,xvol,p,ua,ub,uc,kel,velem,Vmax,gradxnew,gradxdil,gradxer,...
                                                     numviz,nonempty,kmax,N,opfil,ubKua,ucKuc,csnhp);
end

%Computing the lower and upper bounds for the LP problem:
sL = max(-delta, xmin-x);
sL(empty) = 0.0;
sU = min(delta, 1-x);
sU(empty) = 0.0;

%Showing the initial topology of the structure:
if prnt > 0
   figure;
   get_picture(xnew,str.nely,str.nelx);
   iprn = 0;
end

%Solving the topology optimization problem while the stopping criterion is not satisfied:
while ( (count < tolcount) && (iter <= maxit) )
  
      %Solving the LP problem that approximates the original one around the
      %current point:
      [s,fval,flag2] = linprog(df0dx,[],[],A,-b,sL,sU,zeros(n,1),lpoptions); 
      if flag2 ~= 1
         A2 = [A -1 1];
         sL2 = [max(-0.8*delta, xmin-x); 0; 0];
         sU2 = [min(0.8*delta, 1-x); inf; inf];
         [s,~,flag] = linprog(c,[],[],A2,-b,sL2,sU2,zeros(n+2,1),lpoptions2);
         s = s(1:n);
         fval = df0dx'*s;
      else
         flag = flag2;
      end
      snorm = norm(s,inf);
      xold = x; 
      x = x + s;
      deltaold = delta;
      grfold = df0dx;
      Aold = A;
      bold = b;
     
      %Applying the filter:
      [xfil,gradxdil,gradxer] = apply_filter(str,x,N,weigh,wi,beta,opfil,empty,nonempty,kmax,numviz,xmin);
      
      %Computing the global stiffness matrix, and applying
      %the boundary conditions of the structure:
      Ka = global_stiff(str,xfil,pp,row,col,inic,kelc);
      Kc = Ka;
      Ka(linap1,:) = ap11;
      Ka(:,linap1) = ap12;
      Kc(linap2,:) = ap21;
      Kc(:,linap2) = ap22;
      
      %Computing the objective function and the volume
      %constraint at the current point:
      [f0val,b,xnew,xvol,ua,ub,uc,ubKua,ucKuc] = obj_func_mech(str,xfil,p,fa,fc,...
      fa2,fb2,fc2,Ka,Kc,perm_a,perm_c,velem,Vmax,opfil,N,weigh,wi,empty,nonempty,kmax,numviz,sinhp,xmin);
      vol = sum(xnew);
      
      %Computing the predicted reduction of the objective function (predopt) and the 
      %predicted reduction of the infeasibility (predfsb):
      predopt = -fval;
      predfsb = abs(bold) - abs(Aold*s + bold);
      
      %Updating the penalty parameter associated to the merit function:
      thkmin = min(thold1,thold2);
      thklarge = (1 + (1e6/((iter+1)^(1.1))))*thkmin;
      if predopt > 0.5*predfsb
         thksup = 1;
      else
         thksup = (0.5*predfsb)/(predfsb-predopt);
      end
      theta = min(thksup,thklarge);
      theta = min(theta,thmax);
      thold2 = thold1;
      thold1 = theta;
      
      %Computing the predicted reduction of the merit function:
      pred = theta*predopt + (1-theta)*predfsb;
      
      %Computing the actual reduction of the merit function:
      aredopt = f0valold - f0val;
      aredfsb = abs(bold) - abs(b);
      ared = theta*aredopt + (1-theta)*aredfsb;
      
      %Testing if the new point will be accepted or rejected:
      if ared < 0.1*pred
            
         delta = max(0.25*snorm,0.1*deltaold);
         thmax = theta;
         x = xold;
         f0val = f0valold;
         df0dx = grfold;
         A = Aold;
         b = bold;
         
         %Updating the lower and the upper bounds for the next LP problem:
         sL = max(-delta, xmin-x);
         sL(empty) = 0.0;
         sU = min(delta, 1-x);
         sU(empty) = 0.0;
         
         %Counting the number of inner iterations:
         itint = itint + 1;
 
      else
          
         if (ared >= 0.5*pred) 
            delta = min(2.5*deltaold,1.0-xmin);
         else
            delta = delta0;
         end
         thmax = 1.0;
         
         %Updating the lower and the upper bounds for the next LP problem:
         sL = max(-delta, xmin-x);
         sL(empty) = 0.0;
         sU = min(delta, 1-x);
         sU(empty) = 0.0;
        
         %Updating the objective function value:
         f0valold = f0val;
         
         %Counting the number of inner iterations:
         itint = itint + 1;
         
         %Counting the number of outer iterations:
         iter = iter + 1;
            
         %Updating the gradient vectors of the objective function and of
         %the volume constraint:
         if opfil >= 2
            [df0dx,A] = grad_func_mech(n,str,xfil,xvol,p,ua,ub,uc,kel,velem,Vmax,gradxnew,gradxdil,gradxer,...
                                                              numviz,nonempty,kmax,N,opfil,ubKua,ucKuc,csnhp);
         else
            df0dx = grad_func_mech(n,str,xfil,xvol,p,ua,ub,uc,kel,velem,Vmax,gradxnew,gradxdil,gradxer,...
                                                          numviz,nonempty,kmax,N,opfil,ubKua,ucKuc,csnhp);
         end
               
      end
               
      %Showing the data about the inner iteration and the topology of the
      %compliant mechanism:
      if (prnt>0)
         if (mod(iter,prnt)==0)
            get_picture(xnew,str.nely,str.nelx);
            pause(1e-6);
         end
         if (iprn==0),
            disp([]);
            disp(' iter     ||s||          delta           f0val           ared            pred           fracvol         volume        flag');
            disp('----- ------------- --------------- --------------- --------------- --------------- --------------- ---------------- -----');
         end
         iprn = mod(+1,20);
         fprintf('%5u %+13.5E %+15.7E %+15.7E %+15.7E %15.7E %15.7E %5u',iter,snorm,deltaold,f0val,ared,pred,vol/Vdom,vol,flag);
         fprintf('\n');
      end
  
      %Counting the number of times that the norm step associated to the
      %SLP algorithm reaches the minimum value given by the user:
      if (snorm <= tol), 
        count = count+1;
      else
        count = 0;
      end
  
end

%Approximated optimal solution for the topology optimization problem:
xstar = xnew;

%Approximated value of the objective function without any filter:
if (opfil == 0 || opfil == 1) 
   if p == 3      
      opfil = 0; 
      Ka = global_stiff(str,xstar,1,row,col,inic,kelc);
      Kc = Ka;
      Ka(linap1,:) = ap11;
      Ka(:,linap1) = ap12;
      Kc(linap2,:) = ap21;
      Kc(:,linap2) = ap22;  
      f0val2 = obj_func_mech(str,xstar,1,fa,fc,...
               fa2,fb2,fc2,Ka,Kc,perm_a,perm_c,velem,Vmax,opfil,N,weigh,wi,empty,nonempty,kmax,numviz,sinhp,xmin);        
   else      
      f0val2 = f0val;      
   end
end

if (opfil == 2 || opfil == 3) 
   if (p == 3 && beta == 1.6)              
      opfil = 0;
      Ka = global_stiff(str,xstar,1,row,col,inic,kelc);
      Kc = Ka;
      Ka(linap1,:) = ap11;
      Ka(:,linap1) = ap12;
      Kc(linap2,:) = ap21;
      Kc(:,linap2) = ap22;  
      f0val2 = obj_func_mech(str,xstar,1,fa,fc,...
               fa2,fb2,fc2,Ka,Kc,perm_a,perm_c,velem,Vmax,opfil,N,weigh,wi,empty,nonempty,kmax,numviz,sinhp,xmin);         
   else     
      f0val2 = f0val;     
   end
end

if opfil == 4 
   if p == 6       
      opfil = 0; 
      Ka = global_stiff(str,xstar,1,row,col,inic,kelc);
      Kc = Ka;
      Ka(linap1,:) = ap11;
      Ka(:,linap1) = ap12;
      Kc(linap2,:) = ap21;
      Kc(:,linap2) = ap22;  
      f0val2 = obj_func_mech(str,xstar,1,fa,fc,...
               fa2,fb2,fc2,Ka,Kc,perm_a,perm_c,velem,Vmax,opfil,N,weigh,wi,empty,nonempty,kmax,numviz,sinhp,xmin);         
   else    
      f0val2 = f0val;
   end
end

%Showing the optimal topology of the structure and some data:
if (prnt>=0)
   if (prnt==0)
      disp([]);
      disp(' iter     ||s||          delta           f0val           ared            pred           fracvol         volume        flag');
      disp('----- ------------- --------------- --------------- --------------- --------------- --------------- ---------------- -----');
      fprintf('%5u %+13.5E %+15.7E %+15.7E %+15.7E %15.7E %15.7E %5u',iter,snorm,deltaold,f0val,ared,pred,vol/Vdom,vol,flag);
      fprintf('\n');
      figure;
   end      
   get_picture(xstar,str.nely,str.nelx);
end

%Total time taken by the SLP algorithm to find the approximated optimal
%solution of the topology optimization problem:
time = toc;

