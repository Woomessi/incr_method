function [triangles, max_x, min_x, max_y, min_y] = reCons(triangles, size_tri)
%全局最值
max_x = max([triangles(:,1); triangles(:,4);triangles(:,7)])+1e-5;
min_x = min([triangles(:,1); triangles(:,4);triangles(:,7)])-1e-5;

max_y = max([triangles(:,2); triangles(:,5);triangles(:,8)])+1e-5;
min_y = min([triangles(:,2); triangles(:,5);triangles(:,8)])-1e-5;
%重构triangles的数据结构
triangles = [triangles,max(triangles(:,[1 4 7]),[],2),min(triangles(:,[1 4 7]),[],2),max(triangles(:,[ 2 5 8]),[],2),min(triangles(:,[2 5 8]),[],2),max(triangles(:,[3 6 9]),[],2),min(triangles(:,[3 6 9]),[],2)] ;
triangles_1 = zeros(size_tri,22);
for i = 1:size_tri
    x_centroid = (triangles(i,1)+triangles(i,4)+triangles(i,7))./3;
    y_centroid = (triangles(i,2)+triangles(i,5)+triangles(i,8))./3;
    z_centroid = (triangles(i,3)+triangles(i,6)+triangles(i,9))./3;
    %[A] 已知三角形三顶点求面积
    A = 0.5.*norm(cross([triangles(i,4)-triangles(i,1),triangles(i,5)-triangles(i,2),triangles(i,6)-triangles(i,3)],[triangles(i,7)-triangles(i,1),triangles(i,8)-triangles(i,2),triangles(i,9)-triangles(i,3)]));
    triangles_1(i,:) = [triangles(i,:),x_centroid, y_centroid,z_centroid,A];
end
triangles = triangles_1;
end