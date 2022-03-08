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
h = 100;
v = 150;
beta = 2.2;
q_max = 0.025;

% d = 10;
% w = 25;
% D = w - d;
% x_incr = w/10;
% h = 10;
% v = 150;
% beta = 2.2;
% q_max = 0.005;

M = [];
thickness = zeros(size_tri,1);
% 重构
[triangles, max_x, min_x, max_y, min_y] = reCons(triangles, size_tri);
% 确定初始三角形
[tri_ori, cor_ori] = triOri(size_tri, triangles, min_x);
x_ori = triangles(tri_ori,cor_ori);
y_ori = triangles(tri_ori,cor_ori+1);
% 生成surface_section
% all_surface_section = cell(10,50);
% all_points_path = cell(10,50);
points_path = zeros(50,6);
i = 1;
j = 1;
k = 1;

if y_ori < 0.5*(max_y + min_y)
    % d是为了保证曲率变化大的区域也能被搜索到
    while y_ori+D <= max_y+d
        while x_ori+D <= max_x
            [point_surface_section,normal_surface_section,surface_section] = surfaceSection(size_tri, triangles, x_ori, y_ori, D);
            cell_surface_section = num2cell(surface_section);
            all_surface_section (i,j) = {cell_surface_section};
            % 据曲面上的点沿法线方向偏移一个offset，法线反向
            [point_surface_section,normal_surface_section] = offsetting_1(h,point_surface_section,normal_surface_section);
            points_path(k,:) = [point_surface_section,normal_surface_section];
            cell_points_path = num2cell(points_path(k,:));
            all_points_path(i,j) = {cell_points_path};

            x_ori = x_ori+x_incr;
            j = j+1;
            k = k+1;
        end
        size_points_path = size(all_points_path(i,:),2);
        m = zeros(1,size_points_path-1);
        for a = 1:size_points_path-1
            m(a) = norm(cell2mat(all_points_path{i,a+1}(1:3))-cell2mat(all_points_path{i,a}(1:3)));
        end

        M = [M;m];
        i = i+1;
        j = 1;
        y_ori = y_ori+d;
        x_ori = triangles(tri_ori,cor_ori);
    end
else
    while y_ori-D >= min_y-d
        while x_ori+D <= max_x+10
            [point_surface_section,normal_surface_section,surface_section] = surfaceSection_1(size_tri, triangles, x_ori, y_ori, D);
            cell_surface_section = num2cell(surface_section);
            all_surface_section (i,j) = {cell_surface_section};
            % 据曲面上的点沿法线方向偏移一个offset，法线反向
            [point_surface_section,normal_surface_section] = offsetting_1(h,point_surface_section,normal_surface_section);
            points_path(k,:) = [point_surface_section,normal_surface_section];
            cell_points_path = num2cell(points_path(k,:));
            all_points_path(i,j) = {cell_points_path};

            x_ori = x_ori+x_incr;
            j = j+1;
            k = k+1;
        end
        size_points_path = size(all_points_path(i,:),2);
        m = zeros(1,size_points_path-1);
        for a = 1:size_points_path-1
            m(a) = norm(cell2mat(all_points_path{i,a+1}(1:3))-cell2mat(all_points_path{i,a}(1:3)));
        end
        M = [M;m];
        i = i+1;
        j = 1;
        y_ori = y_ori-d;
        x_ori = triangles(tri_ori,cor_ori);
    end
end
points_path(all(points_path==0,2),:)=[];

thickness = computeThickness(all_surface_section, size_points_path, all_points_path, h, beta, q_max, w, M, v, thickness);
triangles = [triangles,thickness];

patchForm(size_tri, triangles);

quiver3(points_path(:,1),points_path(:,2),points_path(:,3),points_path(:,4),points_path(:,5),points_path(:,6),'r')
hold on
scatter3(points_path(:,1),points_path(:,2),points_path(:,3),10,[0.9290 0.6940 0.1250],'filled')
hold on