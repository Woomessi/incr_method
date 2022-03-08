function patchForm(size_tri, triangles)
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
end