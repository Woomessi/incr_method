function [tri_ori, cor_ori] = triOri(size_tri, triangles, min_x)
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
end