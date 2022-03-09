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
M = [];
thickness = zeros(size_tri,1);
% 重构
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
% 确定初始三角形
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
x_ori = triangles(tri_ori,cor_ori);
y_ori = triangles(tri_ori,cor_ori+1);
% 生成surface_section
% all_surface_section = cell(10,50);
% all_points_path = cell(10,50);
points_path = zeros(5000,6);
i = 1;
j = 1;
k = 1;
% d是为了保证曲率变化大的区域也能被搜索到
while y_ori+D <= max_y+d
    while x_ori+D <= max_x
        index_surface_section = zeros(size_tri,1).*NaN;
        n = 1;
%         修改了surface_section的搜索规则
        for i_0 = 1:size_tri
            if triangles(i_0,14)>=x_ori && triangles(i_0,16)+2>=y_ori && triangles(i_0,13)<=x_ori+D && triangles(i_0,15)<=y_ori+D+2
                index_surface_section(n) = i_0;
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
        for i_0 = 1:size_surface_section
            dis(i_0) = norm(surface_section(i_0,19:21)-surface_section(I,19:21));
        end
        %{
Starting from the reference triangle,the areas of triangles are then 
successively added till their cumulative area does not exceed pi.*(0.5.*D).^2
        %}
        surface_section = [surface_section,dis,index_surface_section];
        surface_section = sortrows(surface_section,23);
        i_0 = 1;
        sum_A = surface_section(i_0,22);
        while sum_A <= pi*(0.5.*D).^2 && i_0<size_surface_section
            i_0 = i_0+1;
            sum_A = sum_A + surface_section(i_0,22);
        end
        if i_0 == size_surface_section
            size_surface_section_1 = i_0;
        else
            size_surface_section_1 = i_0-1;
        end
        surface_section = surface_section(1:size_surface_section_1,:);
        %面积加权平均形心
        point_surface_section = sum(surface_section(:,19:21).*surface_section(:,22))./sum(surface_section(:,22));
        %取三角形上离此点之最近点
        dis_1 = zeros(size_surface_section_1,1);
        for i_0 = 1:size_surface_section_1
            dis_1(i_0) = norm(surface_section(i_0,19:21)-point_surface_section);
        end
        [~,I_1] = min(dis_1);
        point_I_1 = [surface_section(I_1,1:3);surface_section(I_1,4:6);surface_section(I_1,7:9);surface_section(I_1,19:21)];
        dis_2 = zeros(4,1);
        for i_0 = 1:4
            dis_2(i_0) = norm(point_I_1(i_0,:)-point_surface_section);
        end
        [~,I_2] = min(dis_2);
        point_surface_section = point_I_1(I_2,:);

        normal_surface_section = sum(surface_section(:,10:12).*surface_section(:,22))./sum(surface_section(:,22));

        cell_surface_section = num2cell(surface_section);
        all_surface_section (i,j) = {cell_surface_section};
        % 据曲面上的点沿法线方向偏移一个offset，法线反向
        normal_surface_section = normal_surface_section ./ norm(normal_surface_section);
        normal_surface_section = h.* normal_surface_section;
        point_surface_section = [point_surface_section,1]';

        operator = [eye(3) normal_surface_section';[0 0 0 1]];
        point_surface_section = operator * point_surface_section;
        point_surface_section = point_surface_section(1:3,:)';
        normal_surface_section = -1 .* normal_surface_section;
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
points_path(all(points_path==0,2),:)=[];
% i = 1;
for i = 1:size(all_surface_section,1)
%     sum_q_1 = 0;
    % 遍历当前pass的所有surface_section
    % b = 1;
    for b = 1:size(all_surface_section(i,:),2)
%         sum_q_1 = 0;
        % 遍历当前surface_section的所有三角形
        for c = 1:size(all_surface_section{i,b},1)
            sum_q_1 = 0;
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
                q_1 = M(i,d-1).*CDB(q_0,h,l,cos_gamma,cos_theta);
                sum_q_1 = sum_q_1 + q_1;
                thickness(t(7)) = sum_q_1./v;
            end
        end
    end
end
triangles = [triangles,thickness];

F = [];
V = [];
F_0 = [1,2,3];
for i = 1:size_tri
    v = [triangles(i,[1:3]); triangles(i,[4:6]); triangles(i,[7:9])];
    V = [V;v];
    F = [F;F_0];
    F_0 = F_0+3;
end
col = triangles(:,23);
patch('Faces',F,'Vertices',V,'FaceVertexCData',col,'FaceColor','flat');
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
hold on
quiver3(points_path(:,1),points_path(:,2),points_path(:,3),points_path(:,4),points_path(:,5),points_path(:,6),'r')
hold on
scatter3(points_path(:,1),points_path(:,2),points_path(:,3),10,[0.9290 0.6940 0.1250],'filled')
hold on
% quiver3(triangles(:,19),triangles(:,20),triangles(:,21),triangles(:,10),triangles(:,11),triangles(:,12),'r')
% hold on