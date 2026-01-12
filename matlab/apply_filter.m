function [xfil,gradxdil,gradxer] = apply_filter(str,x,N,weigh,wi,beta,opfil,empty,nonempty,kmax,numviz,xmin)
%function [xfil,xdens] = apply_filter(str,x,N,weigh,wi,beta,opfil,empty,nonempty,kmax,numviz,xmin)
%************************************************************************************************************
%[xfil,gradxdil,gradxer] = apply_filter(str,x,N,weigh,wi,beta,opfil,empty,nonempty,kmax,numviz,xmin)
%
%This function applies the filter to the densities.
%
%INPUT PARAMETERS:
%-----------------
%str: Structure with the problem data.
%x: Original element densities vector.
%N: Matrix where each column contains the neighbors of each finite element.
%weigh: Matrix where each row contains the weights of each finite element,
%       associated to the mean density filter or to the sinh filter.
%       If the dilation or the erosion filter are used, "weigh" is an empty
%       matrix.
%wi: Vector where each element is the sum of the finite element weights,
%    associated to the mean density filter or to the sinh filter.
%    If the dilation or erosion filter are used, "wi" is an empty
%    vector.
%beta: Parameter used only in dilation and erosion filters. 
%opfil: Parameter that indicates what filter is used.
%       opfil = 1 (mean density filter)
%       opfil = 2 (dilation filter)
%       opfil = 3 (erosion filter)
%       opfil = 4 (sinh filter)
%empty: Vector with the indexes of the finite elements that are always empty.
%nonempty: Vector with the indexes of the finite elements that are nonempty.
%kmax: Number of nonempty elements.
%numviz: Vector that contains the number of neighbors of each finite element.
%xmin: Small positive value that the element densities can assume, to avoid numerical instabilities. 
%
%OUTPUT PARAMETERS:
%------------------
%xfil: Filtered element densities vector. 
%gradxdil: Matrix where each line contains the coefficients associated to
%          the chain rule for the objective function derivatives,
%          when the dilation filter is used. 
%          If another filter is used, "gradxdil" is an empty matrix.
%gradxer: Matrix where each line contains the coefficients associated to
%         the chain rule for the objective function derivatives,
%         when the erosion filter is used. 
%         If another filter is used, "gradxer" is an empty matrix.
%************************************************************************************************************

if opfil == 1
   xdens = mean_density_filter(str,x,weigh,wi,nonempty,kmax,numviz,N);
   xfil = xdens;
   xfil(empty) = xmin;
   xdens(empty) = xmin;
elseif opfil == 2
   xdens = mean_density_filter(str,x,weigh,wi,nonempty,kmax,numviz,N);
   xfil = 1-exp(-beta*xdens)+xdens*exp(-beta);
   xfil(empty) = xmin;
   xdens(empty) = xmin;
end
gradxdil = [];
gradxer = [];

% if (opfil == 0 || opfil == 4) %No filter (opfil == 0), Sinh filter (opfil == 4)
%    xfil = x;
%    xfil(empty) = xmin;
%    gradxdil = [];
%    gradxer = [];
% elseif opfil == 1 %Mean density filter (opfil == 1)
%    xvol = mean_density_filter(str,x,weigh,wi,nonempty,kmax,numviz,N);
%    xfil = xvol;
%    xfil(empty) = xmin;
%    gradxdil = [];
%    gradxer = [];
% elseif opfil == 2 %Dilation filter (opfil == 2)
%    [xdil,gradxdil] = dilation_filter(str,x,beta,N,nonempty,kmax,numviz);
%    xfil = xdil;
%    xfil(empty) = xmin;
%    gradxer = [];
% elseif opfil == 3 %Erosion filter (opfil == 3)
%    [xer,gradxer] = erosion_filter(str,x,beta,N,nonempty,kmax,numviz);
%    xfil = xer;
%    xfil(empty) = xmin;
%    gradxdil = [];
% end
 





 