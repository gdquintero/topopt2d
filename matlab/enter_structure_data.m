function str = enter_structure_data
%**************************************************************************
%str = enter_structure_data
%
%This function reads the problem data (structure or compliant mechanism).
%
%OUTPUT PARAMETER:
%-----------------
%str: Structure with the problem data.
%**************************************************************************

%**************************************************************************
%Reading the structure parameters
%**************************************************************************

str.b = input('Basis: ');
  if str.b <= 0
     disp('Error: this value has to be positive.')
     return;
  end
  
str.h = input('Height: ');
  if str.h <= 0
     disp('Error: this value has to be positive.')
     return;
  end 


str.nelx = input('Number of elements in horizontal: '); 
  if str.nelx <= 0
     disp('Error: this value has to be positive.')
     return;
  end
str.nely = input('Number of elements in vertical: ');
  if str.nely <= 0
     disp('Error: this value has to be positive.')
     return;
  end

str.e = input('Thickness: ');
    if str.e <= 0
       disp('Error: this value has to be positive.');
       return;
    end
    
str.nu = input('Poisson coefficient: ');
    if ( (str.nu <= 0) || (str.nu > 1) )
       disp('Error: this value varies between 0 and 1. ')
       return; 
    end
    
str.E = input('Young modulus: ');
    if str.E <= 0
       disp('Error: this value has to be positive. ');
       return;
    end
disp(' ');

be = str.b/str.nelx; %element basis
he = str.h/str.nely; %element height
str.be2 = be/2; %element basis medium point
str.he2 = he/2; %element height medium point

str.nx = str.nelx + 1; %number of nodes on the horizontal
str.ny = str.nely + 1; %number of nodes on the vertical

str.nnos = str.nx*str.ny; %total number of nodes 
str.nelem = str.nelx*str.nely; %total number of elements

%**************************************************************************
%Reading the applied loads
%**************************************************************************

i = 1;
j = 1; 
while 1    
    disp(['Load ',num2str(i)]);
    disp(' ');
    val = input('x-component of load: ');
    if (isempty(val))      
       break;
    else
      loads.Qx(i) = val;
    end
    val = [];
    while (isempty(val))
       val = input('y-component of load: ');
    end
    loads.Qy(i) = val;  
    op = input('Choose the load distribution (1 = concentrated, 2 = horizontal, 3 = vertical): ');
    if (op == 1) %Concentrated load
       loads.x0(i) = input('x-coordinate of load application: ');
       if ( (loads.x0(i) < 0) || (loads.x0(i) > str.b) )
          disp('Error: the load doesn´t reach the structure.');
          continue;
       end
       loads.y0(i) = input('y-coordinate of load application: ');
       if ( (loads.y0(i)<0) || (loads.y0(i)>str.h) )
          disp('Error: the load doesn´t reach the structure.');
          continue;
       end
       loads.xf(i) = loads.x0(i);
       loads.yf(i) = loads.y0(i);
       nox = floor((loads.x0(i) + str.be2)/be) + 1 ;
       noy = floor((loads.y0(i) + str.he2)/he) + 1 ;
       force.Fx(j) = loads.Qx(i);
       force.Fy(j) = loads.Qy(i);
       force.node(j) = (nox-1)*str.ny + noy;
       j = j+1;
    elseif (op == 2) %Horizontal distributed load
       loads.x0(i) = input('Initial point of horizontal distributed load: ');
       if ( (loads.x0(i) < 0) || (loads.x0(i) > str.b) )
          disp('Error: the load doesn´t reach the structure.');
          continue;
       end
       loads.xf(i) = input('Final point of horizontal distributed load: ');
       if ( (loads.xf(i) <= loads.x0(i)) || (loads.xf(i) > str.b) )
          disp('Error: the load doesn´t reach the structure.');
          continue;
       end          
       loads.y0(i) = input('Vertical distance between the horizontal load and the basis of the structure: ');
       if ( (loads.y0(i) < 0) || (loads.y0(i) > str.h))
          disp('Error: the load doesn´t reach the structure.');
          continue;
       end
       loads.yf(i) = loads.y0(i);
       
       %Numbering the nodes related to the initial and to the final points
       %of horizontal load distribution:
       nox1 = floor((loads.x0(i) + str.be2)/be) + 1 ;
       nox2 = floor((loads.xf(i) + str.be2)/be) + 1 ;
       noy = floor((loads.y0(i) + str.he2)/he) + 1 ;  
       if nox1 == nox2 %Case where the load distribution is between two consecutive elements
          force.Fx(j) = (loads.xf(i) - loads.x0(i))*loads.Qx(i);
          force.Fy(j) = (loads.xf(i) - loads.x0(i))*loads.Qy(i);
          force.node(j) = (nox1-1)*str.ny + noy;
          j = j+1;
       else
          %Initial element of horizontal load distribution:
          force.Fx(j) = loads.Qx(i)*(nox1*be - str.be2 - loads.x0(i));
          force.Fy(j) = loads.Qy(i)*(nox1*be - str.be2 - loads.x0(i));
          force.node(j) = (nox1-1)*str.ny + noy;
          j = j+1;
          %Intermediate elements of horizontal load distribution:
          for k = (nox1+1):(nox2-1)
              force.Fx(j) = loads.Qx(i)*be;
              force.Fy(j) = loads.Qy(i)*be;
              force.node(j) = force.node(j-1) + str.ny;
              j = j+1;
          end
          %Final element of horizontal load distribution:
          force.Fx(j) = loads.Qx(i)*(loads.xf(i) - (nox2-1)*be + str.be2);
          force.Fy(j) = loads.Qy(i)*(loads.xf(i) - (nox2-1)*be + str.be2);
          force.node(j) = (nox2-1)*str.ny + noy;
          j = j+1;
       end   
    elseif (op == 3) %Vertical distributed load
       loads.y0(i) = input('Initial point of vertical distributed load: ');
       if ( (loads.y0(i) < 0) || (loads.y0(i) > str.h) )
          disp('Error: the load doesn´t reach the structure.');
          continue;
       end
       loads.yf(i) = input('Final point of vertical distributed load: ');
       if ( (loads.yf(i) <= loads.y0(i) ) || (loads.yf(i) > str.h) )
          disp('Error: the load doesn´t reach the structure.');
          continue;
       end
       loads.x0(i) = input('Horizontal distance between the horizontal load and the basis of the structure: ');
       if ( (loads.x0(i) < 0) || (loads.x0(i) > str.b) )
          disp('Error: the load doesn´t reach the structure.');
          continue;
       end
       loads.xf(i) = loads.x0(i);
       
       %Numbering the nodes related to the initial and to the final points
       %of vertical load distribution:
       noy1 = floor((loads.y0(i) + str.he2)/he) + 1 ;
       noy2 = floor((loads.yf(i) + str.he2)/he) + 1 ;
       nox = floor((loads.x0(i) + str.be2)/be) + 1 ;  
       if noy1 == noy2 %Case where the load distribution is between two consecutive elements
          force.Fx(j) = (loads.yf(i) - loads.y0(i))*loads.Qx(i);
          force.Fy(j) = (loads.yf(i) - loads.y0(i))*loads.Qy(i);
          force.node(j) = (nox-1)*str.ny + noy1;
          j = j+1;
       else
          %Initial element of vertical load distribution:
          force.Fx(j) = loads.Qx(i)*(noy1*he - str.he2 - loads.y0(i));
          force.Fy(j) = loads.Qy(i)*(noy1*he - str.he2 - loads.y0(i));
          force.node(j) = (nox-1)*str.ny + noy1;
          j = j+1;
          %Intermediate elements of vertical load distribution:
          for k = (noy1+1):(noy2-1)
              force.Fx(j) = loads.Qx(i)*he;
              force.Fy(j) = loads.Qy(i)*he;
              force.node(j) = force.node(j-1) + 1;
              j = j+1;
          end
          %Final element of vertical load distribution:
          force.Fx(j) = loads.Qx(i)*(loads.yf(i) - (noy2-1)*he + str.he2);
          force.Fy(j) = loads.Qy(i)*(loads.yf(i) - (noy2-1)*he + str.he2);
          force.node(j) = (nox-1)*str.ny + noy2;
          j = j+1;
       end
    end
    disp(' ');
    i = i+1;
end

str.loads = loads;
str.force = force;

%**************************************************************************
%Reading the applied supports (boundary conditions)
%**************************************************************************

i = 1;
while 1
    disp(['Support ',num2str(i)]);
    disp(' ');
    val = input('x-coordinate of support: ');
        if (isempty(val))
          break;
        else
          supp.x(i) = val;
        end
        nox = floor((supp.x(i)+str.be2)/be)+1;
        if ( (nox > str.nx) || (nox < 1) )
           disp('Error: the support doesn´t reach the structure.');
           continue;
        end
    supp.y(i) = input('y-coordinate of support: ');
        noy = floor((supp.y(i)+str.he2)/he)+1;
        if ( (noy > str.ny) || (noy < 1) )
           disp('Error: the support doesn´t reach the structure.');
           continue;
        end
    supp.node(i) = (nox-1)*str.ny+noy;   
    val = input('Does the support prevent the movement in horizontal (1 = yes, 0 = no)? ');
    if (val == 0)
        supp.ix(i) = 0;
    else
        supp.ix(i) = 1;
    end
    val = input('Does the support prevent the movement in vertical (1 = yes, 0 = no)? ');
    if (val == 0)
        supp.iy(i) = 0;
    else
        supp.iy(i) = 1;
    end
    ifix = i;
    val = input('Is the support repeated in horizontal (1 = yes, 0 = no)? ');
    if (val == 1)
        val = input('Support length (negative -> the support is on the left of initial point; positive -> the support is on the right of initial point): ');
        if (val ~= 0)
            nox2 = floor((supp.x(i)+val+str.be2)/be)+1;
            nox2 = max(min(nox2,str.nx),1);
            steps = sign(val);
            for j=steps:steps:nox2-nox
               i = i+1;
               supp.x(i) = supp.x(ifix)+j*be;
               supp.y(i) = supp.y(ifix);
               supp.node(i) = supp.node(ifix)+j*str.ny;
               supp.ix(i) = supp.ix(ifix);
               supp.iy(i) = supp.iy(ifix);
            end
        end
    end
    val = input('Is the support repeated in vertical (1 = yes, 0 = no)? ');
    if (val == 1)
        val = input('Support length (negative -> the support is below the initial point; positive -> the support is above the initial point): ');
        if (val ~= 0)
            noy2 = floor((supp.y(i)+val+str.he2)/he)+1;
            noy2 = max(min(noy2,str.ny),1);
            steps = sign(val);
            for j=steps:steps:noy2-noy
               i = i+1;
               supp.x(i) = supp.x(ifix);
               supp.y(i) = supp.y(ifix)+j*he;
               supp.node(i) = supp.node(ifix)+j;
               supp.ix(i) = supp.ix(ifix);
               supp.iy(i) = supp.iy(ifix);
            end
        end
    end
    disp(' ');
    i = i+1;
end

str.supp = supp; 

disp(' ');
