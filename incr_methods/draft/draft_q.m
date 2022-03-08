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
% i:当前pass
i = 1;
% j:当前surfacesection
j = 1;
k = 1;
% d是为了保证曲率变化大的区域也能被搜索到
while y_ori+D <= max_y+d
    while x_ori+D <= max_x
        [point_surface_section,normal_surface_section,surface_section] = surfaceSection(size_tri, triangles, x_ori, y_ori, D);
        cell_surface_section = num2cell(surface_section);
        all_surface_section (i,j) = {cell_surface_section};
        % 据曲面上的点沿法线方向偏移一个offset，法线反向
        [point_surface_section,normal_surface_section] = offsetting_1(100,point_surface_section,normal_surface_section);
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
   
    % 遍历当前pass的所有surface_section
    for b = 1:size(all_surface_section(i,:),2)
        % 遍历当前surface_section的所有三角形
        for c = 1:size(all_surface_section{i,b},1)
            all1_q_1 = 0;
            % 遍历所有的路径点
            for d = 2:size_points_path
                % t:当前三角形的形心、法向量及索引信息
                t = cell2mat(all_surface_section{i,b}(c,[19:21,10:12,24]));
                % p:当前路径点的信息
                p = cell2mat(all_points_path{i,b});
                % l：路径点指向当前三角形形心的向量
                l = t(1:3)-p(1:3);
                cos_gamma = sum(l.*-t(4:6))./(norm(l).*norm(t(4:6)));
                cos_theta = abs(sum(l.*p(4:6))./(norm(l).*norm(p(4:6))));
                ref = p(1:3) + h.*p(4:6)./norm(p(4:6));
                s0 = [(p(6).*l(1).*(ref(3)-p(3)) + p(5).*l(1).*(ref(2)-p(2)) + p(1).*(p(5).*l(2)+p(6).*l(3)) + p(4).*ref(1).*l(1))./(p(4).*l(1)+p(5).*l(2)+p(6).*l(3)),(p(6).*l(2).*(ref(3)-p(3)) + p(4).*l(2).*(ref(1)-p(1)) + p(2).*(p(4).*l(1)+p(6).*l(3)) + p(5).*ref(2).*l(2))./(p(4).*l(1)+p(5).*l(2)+p(6).*l(3)),(p(5).*l(3).*(ref(2)-p(2)) + p(4).*l(3).*(ref(1)-p(1)) + p(3).*(p(4).*l(1)+p(5).*l(2)) + p(6).*ref(3).*l(3))./(p(4).*l(1)+p(5).*l(2)+p(6).*l(3))];
                double r
                r = norm(s0-ref);
                q_0 = Beta(r,beta,q_max,w);
                q_1 = m(d-1).*CDB(q_0,h,l,cos_gamma,cos_theta);
                all1_q_1 = all1_q_1 + q_1;
                thickness(t(7)) = all1_q_1./v;
            end
        end
    end
    i = i+1;
    j = 1;
    y_ori = y_ori+d;
    x_ori = triangles(tri_ori,cor_ori);
end
points_path(all(points_path==0,2),:)=[];

% a = cell2mat(all_points_path{1,1}(1:3))

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