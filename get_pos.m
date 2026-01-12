function pos = get_pos(str,inic)
%**************************************************************************
%pos = get_pos(str,inic)
%
%This function returns the positions of the nonzero elements of the 
%global stiffness matrix.
%
%INPUT PARAMETERS:
%-----------------
%str = Structure with the problem data.
%inic = Vector where each element represents the number minus 1 of the  
%       first nonzero element of each column of the  global stiffness 
%       matrix.
%
%OUTPUT PARAMETER:
%-----------------
%pos = matrix of order 64 x nelem that contains the positions of the
%      nonzero elements of the global stiffness matrix. 
%**************************************************************************

%**************************************************************************
%Defining some constants:
%**************************************************************************

ny = str.ny;
nx = str.nx;
nely = str.ny-1;
nelem = (nx-1)*(ny-1);
nndes = nx*ny;
ind1 = [1;2;5;6;7;8;3;4]; 
ind2 = [1;2;7;8;9;10;3;4];
ind3 = [3;4;9;10;11;12;5;6];
ind4 = [5;6;9;10;11;12;7;8];
ind5 = [7;8;13;14;15;16;9;10];
ind6 = [9;10;15;16;17;18;11;12];
pos = zeros(64,nelem);

%**************************************************************************
%Obtaining the global stiffness matrix columns associated to the
%first column of elements of a rectangular domain:
%**************************************************************************

pos(:,1) = [inic(1)+ind1;inic(2)+ind1;...
            inic(2*ny+1)+ind1;inic(2*ny+2)+ind1;...
            inic(2*ny+3)+ind2;inic(2*ny+4)+ind2;...
            inic(3)+ind2;inic(4)+ind2];
   

for i = 2:nely-1
    pos(:,i) = [inic(2*i-1)+ind3;inic(2*i)+ind3;...
                inic(2*(i+ny)-1)+ind3;inic(2*(i+ny))+ind3;...
                inic(2*(i+ny)+1)+ind2;inic(2*(i+ny)+2)+ind2;...
                inic(2*i+1)+ind2;inic(2*i+2)+ind2];
end

pos(:,nely) = [inic(2*ny-3)+ind3;inic(2*ny-2)+ind3;...
               inic(4*ny-3)+ind3;inic(4*ny-2)+ind3;...
               inic(4*ny-1)+ind1;inic(4*ny)+ind1;...
               inic(2*ny-1)+ind1;inic(2*ny)+ind1];

%**************************************************************************
%Obtaining the global stiffness matrix columns associated to the
%intermediate columns of elements of a rectangular domain:
%**************************************************************************

for j = 2:(nx-2)
    
    pos(:,(j-1)*nely+1) = [inic(2*(j-1)*ny+1)+ind4;inic(2*(j-1)*ny+2)+ind4;...
                           inic(2*j*ny+1)+ind1;inic(2*j*ny+2)+ind1;...
                           inic(2*j*ny+3)+ind2;inic(2*j*ny+4)+ind2;...
                           inic(2*(j-1)*ny+3)+ind5;inic(2*(j-1)*ny+4)+ind5];
      
  
    for i = 2:nely-1
        pos(:,(j-1)*nely+i) = [inic(2*((j-1)*ny+i)-1)+ind6;inic(2*((j-1)*ny+i))+ind6;...
                               inic(2*(j*ny+i)-1)+ind3;inic(2*(j*ny+i))+ind3;...
                               inic(2*(j*ny+i)+1)+ind2;inic(2*(j*ny+i)+2)+ind2;...
                               inic(2*((j-1)*ny+i)+1)+ind5;inic(2*((j-1)*ny+i)+2)+ind5];
    end
  
    pos(:,j*nely) = [inic(2*((j-1)*ny+nely)-1)+ind6;inic(2*((j-1)*ny+nely))+ind6;...
                     inic(2*(j*ny+nely)-1)+ind3;inic(2*(j*ny+nely))+ind3;...
                     inic(2*(j*ny+nely)+1)+ind1;inic(2*(j*ny+nely)+2)+ind1;...
                     inic(2*((j-1)*ny+nely)+1)+ind4;inic(2*((j-1)*ny+nely)+2)+ind4];

end

%**************************************************************************
%Obtaining the global stiffness matrix columns associated to the
%last column of elements of a rectangular domain:
%**************************************************************************

pos(:,nelem-nely+1) = [inic(2*(nndes-2*ny)+1)+ind4;inic(2*(nndes-2*ny)+2)+ind4;...
                       inic(2*(nndes-ny)+1)+ind1;inic(2*(nndes-ny)+2)+ind1;...
                       inic(2*(nndes-ny)+3)+ind2;inic(2*(nndes-ny)+4)+ind2;...
                       inic(2*(nndes-2*ny)+3)+ind5;inic(2*(nndes-2*ny)+4)+ind5];

for i = 2:nely-1
    pos(:,nelem-nely+i) = [inic(2*(nndes-2*ny+i)-1)+ind6;inic(2*(nndes-2*ny+i))+ind6;...
                           inic(2*(nndes-ny+i)-1)+ind3;inic(2*(nndes-ny+i))+ind3;...
                           inic(2*(nndes-ny+i)+1)+ind2;inic(2*(nndes-ny+i)+2)+ind2;...
                           inic(2*(nndes-2*ny+i)+1)+ind5;inic(2*(nndes-2*ny+i)+2)+ind5];
end

pos(:,nelem) = [inic(2*(nndes-ny)-3)+ind6;inic(2*(nndes-ny)-2)+ind6;...
                inic(2*nndes-3)+ind3;inic(2*nndes-2)+ind3;...
                inic(2*nndes-1)+ind1;inic(2*nndes)+ind1;...
                inic(2*(nndes-ny)-1)+ind4;inic(2*(nndes-ny))+ind4];
   
