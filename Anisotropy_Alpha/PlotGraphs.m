
fileID = fopen('file1.txt','r');

formatSpec = '%f %f %f';
sizeA = [3 Inf];

A = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);

fileID = fopen('file2.txt','r');

formatSpec = '%f %f %f';
sizeA = [3 Inf];

B = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);


fileID = fopen('file3.txt','r');

formatSpec = '%f %f %f';
sizeA = [3 Inf];

C = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);

fileID = fopen('file4.txt','r');

formatSpec = '%f %f %f';
sizeA = [3 Inf];

D = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);




%D = zeros(3,size(B,2));
%E = zeros(3,size(C,2));
k = 1;
for i=1:1:size(A,2)
    for j=1:1:size(B,2)
    if(A(1,i) == B(1,j) && A(2,i) == B(2,j))
        E(1,k) = B(1,j);
        E(2,k) = B(2,j);
        E(3,k) = B(3,j) + A(3,i);
        k = k +1;
        break;
    end
    end
end

k = 1;
for i=1:1:size(E,2)
    for j=1:1:size(C,2)
    if(E(1,i) == C(1,j) && E(2,i) == C(2,j))
        F(1,k) = C(1,j);
        F(2,k) = C(2,j);
        F(3,k) = C(3,j) + E(3,i);
        k = k +1;
        break;
    end
    end
end

k = 1;
for i=1:1:size(F,2)
    for j=1:1:size(D,2)
    if(F(1,i) == D(1,j) && F(2,i) == D(2,j))
        G(1,k) = D(1,j);
        G(2,k) = D(2,j);
        G(3,k) = D(3,j) + F(3,i);
        k = k +1;
        break;
    end
    end
end


tri = delaunay(G(1,:),G(2,:)); %x,y,z column vectors
trisurf(tri,G(1,:),G(2,:),G(3,:));
axis ([-1 1 0 10 3e-25 100]);
x = G(1,:);
y = G(2,:);
z = G(3,:);
indexmin = find(min(z) == z);
xmin = x(indexmin);
ymin = y(indexmin);
zmin = z(indexmin);

   line(x(indexmin),y(indexmin),z(indexmin),...
         'marker','o',...
         'markersize',10,...
         'markeredgecolor',[0,0,0],...
         'markerfacecolor',[1,1,1],...
         'linestyle','none');

indexmax = find(max(z) == z);
xmax = x(indexmax);
ymax = y(indexmax);
zmax = z(indexmax);

   line(x(indexmax),y(indexmax),z(indexmax),...
         'marker','o',...
         'markersize',10,...
         'markeredgecolor',[0,0,0],...
         'markerfacecolor',[1,1,1],...
         'linestyle','none');


strmin = ['Minimum = ',num2str(xmin),',',num2str(ymin),',',num2str(zmin)];
text(xmin,ymin,zmin,strmin,'HorizontalAlignment','left');

strmax = ['Maximum = ',num2str(xmax),',',num2str(ymax),',',num2str(zmax)];
text(xmax,ymax,zmax,strmax,'HorizontalAlignment','right');