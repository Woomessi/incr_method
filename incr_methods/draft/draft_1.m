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
tri_size = size(triangles,1);
d = 43.46;
w = 100;
D = w - d;

% plot_stl(triangles);

%【1e-5】摄动法，避免切平面仅与两端三角面片的一个顶点相交，这会导致程序错误
max_x = max([triangles(:,1); triangles(:,4);triangles(:,7)])+1e-5;
min_x = min([triangles(:,1); triangles(:,4);triangles(:,7)])-1e-5;

max_y = max([triangles(:,2); triangles(:,5);triangles(:,8)])+1e-5;
min_y = min([triangles(:,2); triangles(:,5);triangles(:,8)])-1e-5;

triangles = [triangles,max(triangles(:,[1 4 7]),[],2),min(triangles(:,[1 4 7]),[],2),max(triangles(:,[ 2 5 8]),[],2),min(triangles(:,[2 5 8]),[],2),max(triangles(:,[3 6 9]),[],2),min(triangles(:,[3 6 9]),[],2)] ;

for i = 1:tri_size
    if abs(triangles(i,14) - min_x) <= 1e-4
        tri_ori = i;
    end
end

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

surface_section_index = zeros(tri_size,1).*NaN;
n = 1;

for i = 1:tri_size
if triangles(i,14)>=triangles(tri_ori,cor_ori) && triangles(i,16)>=triangles(tri_ori,cor_ori+1) && triangles(i,13)<=triangles(tri_ori,cor_ori)+D && triangles(i,15)<=triangles(tri_ori,cor_ori+1)+D
    surface_section_index(n) = i;
    n = n+1;
end
end


surface_section_index(isnan(surface_section_index))=[];
surface_section = triangles(surface_section_index,:);
plot_stl(surface_section)
% hold on
% plot_stl(triangles)




