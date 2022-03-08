% the global unit is mm
clc
clear
%{
The generation of surface sections starts from origin of the
triangulated surface which is shifted to one corner of
the surface with X-axis oriented along the boundary of the
surface
%}
triangles = read_binary_stl_file('tb.STL');
size_tri = size(triangles,1);
d = 20;
w = 50;
D = w - d;
x_incr = w/10;
% 重构
[triangles, max_x, min_x, max_y, min_y] = reCons(triangles, size_tri);
% 确定初始三角形
[tri_ori, cor_ori] = triOri(size_tri, triangles, min_x);
x_ori = triangles(tri_ori,cor_ori);
y_ori = triangles(tri_ori,cor_ori+1);
% 生成surface_section
points_path = zeros(50,6);
i = 1;
% d是为了保证曲率变化大的区域也能被搜索到
while y_ori+D <= max_y+d
while x_ori+D <= max_x
[point_surface_section,normal_surface_section,surface_section] = surfaceSection(size_tri, triangles, x_ori, y_ori, D);
% 据曲面上的点沿法线方向偏移一个offset，法线反向
[point_surface_section,normal_surface_section] = offsetting_1(100,point_surface_section,normal_surface_section);
points_path(i,:) = [point_surface_section,normal_surface_section];

x_ori = x_ori+x_incr;
i = i+1;
end
y_ori = y_ori+d;
x_ori = triangles(tri_ori,cor_ori);
end
points_path(all(points_path==0,2),:)=[];

plot_stl(triangles)
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
hold on
quiver3(points_path(:,1),points_path(:,2),points_path(:,3),points_path(:,4),points_path(:,5),points_path(:,6),'r')
hold on
scatter3(points_path(:,1),points_path(:,2),points_path(:,3),10,[0.9290 0.6940 0.1250],'filled')
hold on