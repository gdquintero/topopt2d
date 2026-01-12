function get_picture(x,nely,nelx)
%**************************************************************************
%get_picture(x,nely,nelx)
%
%This function shows the structure or compliant mechanism topologies.
%
%INPUT PARAMETERS:
%-----------------
%x: Element densities vector.
%nely: Number of elements in vertical in a rectangular domain.
%nelx: Number of elements in horizontal in a rectangular domain.
%**************************************************************************

R = -reshape(x,nely,nelx);
colormap('gray');
imagesc(R);
axis xy;
