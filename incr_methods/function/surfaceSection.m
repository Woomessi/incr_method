function [point_surface_section,normal_surface_section,surface_section] = surfaceSection(size_tri, triangles, x_ori, y_ori, D)
index_surface_section = zeros(size_tri,1).*NaN;
n = 1;

for i = 1:size_tri
if triangles(i,14)>=x_ori && triangles(i,16)+2>=y_ori && triangles(i,13)<=x_ori+D && triangles(i,15)<=y_ori+D+2
    index_surface_section(n) = i;
    n = n+1;
end
end

index_surface_section(isnan(index_surface_section))=[];
surface_section = [triangles(index_surface_section,:)];
size_surface_section = size(surface_section,1);
% hold on
% plot_stl(triangles)

%{
the triangles are sorted based on the relative distances of their centroids 
compared to a reference triangle which has lowest y coordinate value within
the sorted triangles
%}
[~,I] = min(surface_section(:,16));
dis = zeros(size_surface_section,1);
for i = 1:size_surface_section
dis(i) = norm(surface_section(i,19:21)-surface_section(I,19:21));
end
%{
Starting from the reference triangle,the areas of triangles are then 
successively added till their cumulative area does not exceed pi.*(0.5.*D).^2
%}
surface_section = [surface_section,dis,index_surface_section];
surface_section = sortrows(surface_section,23);
i = 1;
sum_A = surface_section(i,22);
while sum_A <= pi*(0.5.*D).^2 && i<size_surface_section
    i = i+1;
    sum_A = sum_A + surface_section(i,22);
end
if i == size_surface_section
    size_surface_section_1 = i;
else
size_surface_section_1 = i-1;
end
surface_section = surface_section(1:size_surface_section_1,:);
%面积加权平均形心
point_surface_section = sum(surface_section(:,19:21).*surface_section(:,22))./sum(surface_section(:,22));
%取三角形上离此点之最近点
dis_1 = zeros(size_surface_section_1,1);
for i = 1:size_surface_section_1
dis_1(i) = norm(surface_section(i,19:21)-point_surface_section);
end
[~,I_1] = min(dis_1);
point_I_1 = [surface_section(I_1,1:3);surface_section(I_1,4:6);surface_section(I_1,7:9);surface_section(I_1,19:21)];
dis_2 = zeros(4,1);
for i = 1:4
    dis_2(i) = norm(point_I_1(i,:)-point_surface_section);
end
[~,I_2] = min(dis_2);
point_surface_section = point_I_1(I_2,:);

normal_surface_section = sum(surface_section(:,10:12).*surface_section(:,22))./sum(surface_section(:,22));
end