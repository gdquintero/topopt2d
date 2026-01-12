function [row,col,inic] = get_rows_cols(str)
%**************************************************************************
%[row,col,inic] = get_rows_cols(str)
%
%This function returns the vectors used to compute the global 
%stiffness matrix.
%
%INPUT PARAMETER:
%----------------
%str: Structure with the problem data.
%
%OUTPUT PARAMETERS:
%------------------
%row = Vector that contains the indexes of the rows associated to 
%      nonzero elements of the global stiffness matrix.
%col = Vector that contains the indexes of the columns associated to 
%      nonzero elements of the global stiffness matrix.
%inic = Vector where each element represents the number minus 1 of the 
%       first nonzero element of each column of the global stiffness 
%       matrix.
%**************************************************************************

%**************************************************************************
%Defining some constants:
%**************************************************************************

ny = str.ny;
nx = str.nx;
nndes = nx*ny;
ncol = 2*str.nx*str.ny;
nelmat = 64 + 48*(ny-2) + 48*(nx-2) + 36*(nx-2)*(ny-2);
inic = zeros(ncol+1,1);
row = zeros(nelmat,1);
col = zeros(nelmat,1);
ant = [-2*ny-3;-2*ny-2;-2*ny-1;-2*ny;-2*ny+1;-2*ny+2]; 
atu = [-3;-2;-1;0;1;2]; 
pro = [2*ny-3;2*ny-2;2*ny-1;2*ny;2*ny+1;2*ny+2];

%**************************************************************************
%Generating the data associated to the first column of the elements of a
%rectangular domain:
%**************************************************************************

nel = 8;
inic(2) = inic(1)+nel;
inic(3) = inic(2)+nel;
col(inic(1)+1:inic(1)+nel) = 1;
col(inic(2)+1:inic(2)+nel) = 2;
row(inic(1)+1:inic(1)+nel) = 2+[atu(3:6);pro(3:6)];
row(inic(2)+1:inic(2)+nel) = 2+[atu(3:6);pro(3:6)];

nel = 12;
for i = 2:ny-1
    inic(2*i) = inic(2*i-1)+nel;
    inic(2*i+1) = inic(2*i)+nel;
    col(inic(2*i-1)+1:inic(2*i-1)+nel) = 2*i-1;
    col(inic(2*i)+1:inic(2*i)+nel) = 2*i;
    row(inic(2*i-1)+1:inic(2*i-1)+nel) = 2*i+[atu;pro];
    row(inic(2*i)+1:inic(2*i)+nel) = 2*i+[atu;pro];
end

nel = 8;
inic(2*ny) = inic(2*ny-1)+nel;
inic(2*ny+1) = inic(2*ny)+nel;
col(inic(2*ny-1)+1:inic(2*ny-1)+nel) = 2*ny-1;
col(inic(2*ny)+1:inic(2*ny)+nel) = 2*ny;
row(inic(2*ny-1)+1:inic(2*ny-1)+nel) = 2*ny+[atu(1:4);pro(1:4)];
row(inic(2*ny)+1:inic(2*ny)+nel) = 2*ny+[atu(1:4);pro(1:4)];

%**************************************************************************
%Generating the data associated to the intermediate columns of the elements 
%of a rectangular domain:
%**************************************************************************

for j = 2:(nx-1)

    nel = 12;
    inic(2*(j-1)*ny+2) = inic(2*(j-1)*ny+1)+nel;
    inic(2*(j-1)*ny+3) = inic(2*(j-1)*ny+2)+nel;
    col(inic(2*(j-1)*ny+1)+1:inic(2*(j-1)*ny+1)+nel) = 2*(j-1)*ny+1;
    col(inic(2*(j-1)*ny+2)+1:inic(2*(j-1)*ny+2)+nel) = 2*(j-1)*ny+2;
    row(inic(2*(j-1)*ny+1)+1:inic(2*(j-1)*ny+1)+nel) = 2*(j-1)*ny+2+[ant(3:6);atu(3:6);pro(3:6)];
    row(inic(2*(j-1)*ny+2)+1:inic(2*(j-1)*ny+2)+nel) = 2*(j-1)*ny+2+[ant(3:6);atu(3:6);pro(3:6)];

    nel = 18;
    for i = 2:ny-1
        inic(2*(j-1)*ny+2*i) = inic(2*(j-1)*ny+2*i-1)+nel;
        inic(2*(j-1)*ny+2*i+1) = inic(2*(j-1)*ny+2*i)+nel;
        col(inic(2*(j-1)*ny+2*i-1)+1:inic(2*(j-1)*ny+2*i-1)+nel) = 2*(j-1)*ny+2*i-1;
        col(inic(2*(j-1)*ny+2*i)+1:inic(2*(j-1)*ny+2*i)+nel) = 2*(j-1)*ny+2*i;
        row(inic(2*(j-1)*ny+2*i-1)+1:inic(2*(j-1)*ny+2*i-1)+nel) = 2*(j-1)*ny+2*i+[ant;atu;pro];
        row(inic(2*(j-1)*ny+2*i)+1:inic(2*(j-1)*ny+2*i)+nel) = 2*(j-1)*ny+2*i+[ant;atu;pro];
    end

    nel = 12;
    inic(2*j*ny) = inic(2*j*ny-1)+nel;
    inic(2*j*ny+1) = inic(2*j*ny)+nel;
    col(inic(2*j*ny-1)+1:inic(2*j*ny-1)+nel) = 2*j*ny-1;
    col(inic(2*j*ny)+1:inic(2*j*ny)+nel) = 2*j*ny;
    row(inic(2*j*ny-1)+1:inic(2*j*ny-1)+nel) = 2*j*ny+[ant(1:4);atu(1:4);pro(1:4)];
    row(inic(2*j*ny)+1:inic(2*j*ny)+nel) = 2*j*ny+[ant(1:4);atu(1:4);pro(1:4)];

end

%**************************************************************************
%Generating the data associated to the last column of the elements of a 
%rectangular domain:
%**************************************************************************

nel = 8;
inic(2*nndes-2*ny+2) = inic(2*nndes-2*ny+1)+nel;
inic(2*nndes-2*ny+3) = inic(2*nndes-2*ny+2)+nel;
col(inic(2*nndes-2*ny+1)+1:inic(2*nndes-2*ny+1)+nel) = 2*nndes-2*ny+1;
col(inic(2*nndes-2*ny+2)+1:inic(2*nndes-2*ny+2)+nel) = 2*nndes-2*ny+2;
row(inic(2*nndes-2*ny+1)+1:inic(2*nndes-2*ny+1)+nel) = 2*nndes-2*ny+2+[ant(3:6);atu(3:6)];
row(inic(2*nndes-2*ny+2)+1:inic(2*nndes-2*ny+2)+nel) = 2*nndes-2*ny+2+[ant(3:6);atu(3:6)];

nel = 12;
for i = 2:ny-1
    inic(2*nndes-2*ny+2*i) = inic(2*nndes-2*ny+2*i-1)+nel;
    inic(2*nndes-2*ny+2*i+1) = inic(2*nndes-2*ny+2*i)+nel;
    col(inic(2*nndes-2*ny+2*i-1)+1:inic(2*nndes-2*ny+2*i-1)+nel) = 2*nndes-2*ny+2*i-1;
    col(inic(2*nndes-2*ny+2*i)+1:inic(2*nndes-2*ny+2*i)+nel) = 2*nndes-2*ny+2*i;
    row(inic(2*nndes-2*ny+2*i-1)+1:inic(2*nndes-2*ny+2*i-1)+nel) = 2*nndes-2*ny+2*i+[ant;atu];
    row(inic(2*nndes-2*ny+2*i)+1:inic(2*nndes-2*ny+2*i)+nel) = 2*nndes-2*ny+2*i+[ant;atu];
end

nel = 8;
inic(2*nndes) = inic(2*nndes-1)+nel;
inic(2*nndes+1) = inic(2*nndes)+nel;
col(inic(2*nndes-1)+1:inic(2*nndes-1)+nel) = 2*nndes-1;
col(inic(2*nndes)+1:inic(2*nndes)+nel) = 2*nndes;
row(inic(2*nndes-1)+1:inic(2*nndes-1)+nel) = 2*nndes+[ant(1:4);atu(1:4)];
row(inic(2*nndes)+1:inic(2*nndes)+nel) = 2*nndes+[ant(1:4);atu(1:4)];

