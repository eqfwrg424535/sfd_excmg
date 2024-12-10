%% script for creating the Cascadia model using the open data
% the interpolation method used here is griddatan('linear') 
% the semi-uniform grid used is 256*256*256
% the original data is of size 80*78*34
clc;
clear;

% read in the data
[x,y,z,rho,nzAir,type,origin,rotation] = read_WS3d_model('cascad_prior.ws');
%sum(x),sum(y),sum(z)
x = x-origin(1);
y = y-origin(2);
% build the original mesh

aax =cumsum(x)-sum(x)/2+x/2; % cell centers
bby =cumsum(y)-sum(y)/2+y/2;
ccz =cumsum(z)-z/2;
coord = zeros(size(x,1)*size(y,1)*size(z,1),3);
for k = 1:size(z,1)
    for j = 1:size(y,1)
       for i= 1:size(x,1)
           coord((k-1)*size(y,1)*size(x,1)+(j-1)*size(x,1)+i,1) = aax(i);
           coord((k-1)*size(y,1)*size(x,1)+(j-1)*size(x,1)+i,2) = bby(j);
           coord((k-1)*size(y,1)*size(x,1)+(j-1)*size(x,1)+i,3) = ccz(k);
       end
    end
 end
% define the new mesh
a = zeros(256,1);
a(65:192) = 6e3;
for i = 1:64
  a(192+i) = 6e3*1.02^i;
end
a(1:64) =a(256:-1:193);
b = a;
c = zeros(256,1);
for i = 1:256
    c(i) = 50*1.022^i;
end
aax2 =cumsum(a)-sum(a)/2+a/2; % cell centers
bby2 =cumsum(b)-sum(b)/2+b/2;
ccz2 =cumsum(c)-c/2;
[X Y Z] = meshgrid(aax2,bby2,ccz2);
xq = [X(:) Y(:) Z(:)];
% interpolation
Z1=griddatan(coord,rho(:),xq,'linear');

level = 7; 
eps = 1e-10;
solver = 1;
output_type = 2;
nair = 128;
topo = 0;
ext = [1.02 1.022];
abu = [65 192 65 192 0 0];
conduc_type =1;
% output
fid = fopen('cascad.fwd','w');
fprintf(fid,'#the Cascadia model, generated with Matlab\n');
fprintf(fid,'%d %e %d %d\n', level, eps,solver,output_type);
fprintf(fid,'%d %d %d %d %d\n', size(a,1),size(b,1),size(c,1),nair,topo);
for i = 1:size(a,1)
    if i~= size(a,1)
       fprintf(fid,'%.2e ',a(i));
    else
        fprintf(fid,'%.2e\n',a(i));
    end
end
for i = 1:size(b,1)
    if i~= size(b,1)
       fprintf(fid,'%.2e ',b(i));
    else
        fprintf(fid,'%.2e\n',b(i));
    end
end
for i = 1:size(c,1)
    if i~= size(c,1)
       fprintf(fid,'%.2e ',c(i));
    else
        fprintf(fid,'%.2e\n',c(i));
    end
end
fprintf(fid,'%f %f\n', ext);
fprintf(fid,'%d %d %d %d %d %d\n', abu);
fprintf(fid,'%d\n',conduc_type);
for i =1:size(Z1,1)
    fprintf(fid,'%.2e\n',Z1(i));
end

fclose(fid);