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
d = 43.46;
w = 100;
D = w - d;
x_incr = w/10;
%全局最值
max_x = max([triangles(:,1); triangles(:,4);triangles(:,7)])+1e-5;
min_x = min([triangles(:,1); triangles(:,4);triangles(:,7)])-1e-5;

max_y = max([triangles(:,2); triangles(:,5);triangles(:,8)])+1e-5;
min_y = min([triangles(:,2); triangles(:,5);triangles(:,8)])-1e-5;

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
%确定初始三角形
for i = 1:size_tri
    if abs(triangles(i,14) - min_x) <= 1e-4
        tri_ori = i;
    end
end
%确定初始顶点
if abs(triangles(tri_ori,1) - min_x) <= 1e-4
    cor_ori = 1;
else
    if abs(triangles(tri_ori,4) - min_x) <= 1e-4
        cor_ori = 4;
    else
        if abs(triangles(tri_ori,7) - min_x) <= 1e-4
            cor_ori = 7;
        end
    end
end
%生成surface_section
index_surface_section = zeros(size_tri,1).*NaN;
n = 1;

for i = 1:size_tri
if triangles(i,14)>=triangles(tri_ori,cor_ori) && triangles(i,16)>=triangles(tri_ori,cor_ori+1) && triangles(i,13)<=triangles(tri_ori,cor_ori)+D && triangles(i,15)<=triangles(tri_ori,cor_ori+1)+D
    index_surface_section(n) = i;
    n = n+1;
end
end

index_surface_section(isnan(index_surface_section))=[];
surface_section = triangles(index_surface_section,:);
size_surface_section = size(surface_section,1);
% hold on
% plot_stl(triangles)

%{
the triangles are sorted based on the relative distances of their centroids 
compared to a reference triangle which has lowest y coordinate value within
the sorted triangles
%}
[M,I] = min(surface_section(:,16));
dis = zeros(size_surface_section,1);
for i = 1:size_surface_section
dis(i) = norm(surface_section(i,19:21)-surface_section(I,19:21));
end
%{
Starting from the reference triangle,the areas of triangles are then 
successively added till their cumulative area does not exceed pi.*(0.5.*D).^2
%}
surface_section = [surface_section,dis];
surface_section = sortrows(surface_section,23);
i = 1;
sum_A = surface_section(i,22);
pi*(0.5.*D).^2
while sum_A <= pi*(0.5.*D).^2
    i = i+1;
    sum_A = sum_A + surface_section(i,22);
end
size_surface_section_1 = i-1;
surface_section = surface_section(1:size_surface_section_1,:);
%面积加权平均形心
point_surface_section = sum(surface_section(:,19:21).*surface_section(:,22))./sum(surface_section(:,22));
%取三角形上离此点之最近点
dis_1 = zeros(size_surface_section_1,1);
for i = 1:size_surface_section_1
dis_1(i) = norm(surface_section(i,19:21)-point_surface_section);
end
[M_1,I_1] = min(dis_1);
point_I_1 = [surface_section(I_1,1:3);surface_section(I_1,4:6);surface_section(I_1,7:9);surface_section(I_1,19:21)];
dis_2 = zeros(4,1);
for i = 1:4
    dis_2(i) = norm(point_I_1(i,:)-point_surface_section);
end
[M_2,I_2] = min(dis_2);
point_surface_section = point_I_1(I_2,:);

normal_surface_section = sum(surface_section(:,10:12).*surface_section(:,22))./sum(surface_section(:,22));

% plot_stl(surface_section)
% xlabel('x')
% ylabel('y')
% zlabel('z')
% axis equal
% hold on
% quiver3(point_surface_section(1),point_surface_section(2),point_surface_section(3),normal_surface_section(1),normal_surface_section(2),normal_surface_section(3),'r')
% hold on
% scatter3(point_surface_section(1),point_surface_section(2),point_surface_section(3),10,[0.9290 0.6940 0.1250],'filled')
% hold on

%next surface section
index_surface_section_1 = zeros(size_tri,1).*NaN;
n = 1;

for i = 1:size_tri
if triangles(i,14)>=triangles(tri_ori,cor_ori)+x_incr && triangles(i,16)>=triangles(tri_ori,cor_ori+1) && triangles(i,13)<=triangles(tri_ori,cor_ori)+D+x_incr && triangles(i,15)<=triangles(tri_ori,cor_ori+1)+D
    index_surface_section_1(n) = i;
    n = n+1;
end
end