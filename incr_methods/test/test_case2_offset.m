% the global unit is mm
clc
clear
%{
The generation of surface sections starts from origin of the
triangulated surface which is shifted to one corner of
the surface with X-axis oriented along the boundary of the
surface
%}
triangles = read_binary_stl_file('tb_1.STL');
size_tri = size(triangles,1);
d = 20;
w = 50;
D = w - d;
x_incr = w/10;
h = 5;
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

while y_ori-D >= min_y-d
    while x_ori+D <= max_x+10
        index_surface_section = zeros(size_tri,1).*NaN;
        n = 1;

        for i_2 = 1:size_tri
            if triangles(i_2,14)>=x_ori && triangles(i_2,15)<=y_ori && triangles(i_2,13)<=x_ori+D && triangles(i_2,16)>=y_ori-D
                index_surface_section(n) = i_2;
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
        [~,I] = max(surface_section(:,15));
        dis = zeros(size_surface_section,1);
        for i_2 = 1:size_surface_section
            dis(i_2) = norm(surface_section(i_2,19:21)-surface_section(I,19:21));
        end
        %{
Starting from the reference triangle,the areas of triangles are then 
successively added till their cumulative area does not exceed pi.*(0.5.*D).^2
        %}
        surface_section = [surface_section,dis,index_surface_section];
        surface_section = sortrows(surface_section,23);
        i_2 = 1;
        sum_A = surface_section(i_2,22);
        while sum_A <= pi*(0.5.*D).^2 && i_2<size_surface_section
            i_2 = i_2+1;
            sum_A = sum_A + surface_section(i_2,22);
        end
        if i_2 == size_surface_section
            size_surface_section_1 = i_2;
        else
            size_surface_section_1 = i_2-1;
        end
        surface_section = surface_section(1:size_surface_section_1,:);
        %面积加权平均形心
        point_surface_section = sum(surface_section(:,19:21).*surface_section(:,22))./sum(surface_section(:,22));
        %取三角形上离此点之最近点
        dis_1 = zeros(size_surface_section_1,1);
        for i_2 = 1:size_surface_section_1
            dis_1(i_2) = norm(surface_section(i_2,19:21)-point_surface_section);
        end
        [~,I_1] = min(dis_1);
        point_I_1 = [surface_section(I_1,1:3);surface_section(I_1,4:6);surface_section(I_1,7:9);surface_section(I_1,19:21)];
        dis_2 = zeros(4,1);
        for i_2 = 1:4
            dis_2(i_2) = norm(point_I_1(i_2,:)-point_surface_section);
        end
        [~,I_2] = min(dis_2);
        point_surface_section = point_I_1(I_2,:);

        normal_surface_section = sum(surface_section(:,10:12).*surface_section(:,22))./sum(surface_section(:,22));

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
points_path(all(points_path==0,2),:)=[];

thickness = computeThickness(all_surface_section, size_points_path, all_points_path, h, beta, q_max, w, M, v, thickness);
triangles = [triangles,abs(thickness)];

patchForm(size_tri, triangles);

quiver3(points_path(:,1),points_path(:,2),points_path(:,3),points_path(:,4),points_path(:,5),points_path(:,6),'r')
hold on
scatter3(points_path(:,1),points_path(:,2),points_path(:,3),10,[0.9290 0.6940 0.1250],'filled')
hold on