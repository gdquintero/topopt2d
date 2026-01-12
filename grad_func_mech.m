function [df0dx,gdesl] = grad_func_mech(str,xfil,p,u,kel,gradxnew,numviz,nonempty,kmax,N,lin,lambdain,lambdaout)

[df0dx,gdesl] = grad_compliance_mean_density_mech(str,xfil,p,u,kel,gradxnew,numviz,nonempty,kmax,N,lin,lambdain,lambdaout);

