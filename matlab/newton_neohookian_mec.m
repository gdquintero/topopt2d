function [f0val,ustar,lambdain,lambdaout] = newton_neohookian_mec(estr,x,p,u,f,ein,eout,dout,kout,coordno,lin,col,inic,linap,ap1,ap2,tol,maxit)
%[f0val,ustar,lambdain,lambdaout] = newton_neohookian_mec(estr,x,p,u,f,ein,eout,dout,kout,coordno,lin,col,inic,linap,ap1,ap2,tol,maxit)

%Inicializando algumas variaveis:
iter = 0; 
fractol = 1.0;

%Matriz de rigidez global associada aos grandes deslocamentos:
K = rigglobal_grandesl_neohookian(estr,x,p,lin,col,inic,u,coordno);

%Aplicando a rigidez da mola:
K(dout,dout) = K(dout,dout) + kout;
      
%Aplicando apoios:
K(linap,:) = ap1;
K(:,linap) = ap2;
      
%Fatoracao LDLt de K:
[L,D,P] = ldl(K);

%Impressao de alguns dados:
disp('iter        tols         tolrhs');
disp('-----  ------------- -------------');

%Montagem do lado direito do sistema linear:
rhs = rhs_newton_neohookian(estr,x,p,u,f,coordno);
rhs(dout) = rhs(dout) + kout*u(dout);
rhs(linap) = 0.0;
tolrhs = norm(rhs,2);
    
while ( (iter <= maxit) && (tolrhs > tol) && (fractol <= 10.0) ) 
          
      %Obtendo a direcao de Newton:
      s = calculau_ldl(-rhs,L,D,P);
         
      %Atualizando a solucao:
      u = u + s; 
          
      %Atualizando o residuo:
      rhs = rhs_newton_neohookian(estr,x,p,u,f,coordno);
      rhs(dout) = rhs(dout) + kout*u(dout);
      rhs(linap) = 0.0;
         
      %Matriz de rigidez global associada aos grandes deslocamentos:
      K = rigglobal_grandesl_neohookian(estr,x,p,lin,col,inic,u,coordno);
      
      %Aplicando a rigidez da mola:
      K(dout,dout) = K(dout,dout) + kout;

      %Aplicando apoios:
      K(linap,:) = ap1;
      K(:,linap) = ap2;
         
      %Fatoracao LDLt de K:
      [L,D,P] = ldl(K);
                    
      %Numero de iteracoes:
      iter = iter + 1;
         
      %Norma do passo:
      tols = norm(s,2);
        
      %Norma do residuo:
      tolrhs = norm(rhs,2);
               
      %Impressao de alguns dados:
      fprintf('%5u %13.5E %+13.5E\n',iter,tols,tolrhs);
               
end

%Solucao obtida pelo metodo de Newton ou pelo metodo de comprimento de arco:
ustar = u;

%Calculando o valor da funcao objetivo (compliance):
f0val = f'*ustar;

%Obtendo o vetor de multiplicadores de Lagrange para o 
%calculo do gradiente da funcao objetivo pelo metodo adjunto:
lambdaout = calculau_ldl(-eout,L,D,P);
lambdain = calculau_ldl(-ein,L,D,P);












