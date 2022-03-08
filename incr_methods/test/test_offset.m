clc
clear
[point_surface_section,normal_surface_section] = offsetting_1(300,[1,2,3],[2,3,4]);
a = norm([1,2,3]-point_surface_section);
b = norm(normal_surface_section);
quiver3(point_surface_section(1),point_surface_section(2),point_surface_section(3),normal_surface_section(1),normal_surface_section(2),normal_surface_section(3),'r')
axis equal
hold on
scatter3(point_surface_section(1),point_surface_section(2),point_surface_section(3),10,[0.9290 0.6940 0.1250],'filled')
hold on
scatter3(1,2,3,10,[0.9290 0.6940 0.1250],'filled')
hold on
quiver3(1,2,3,2,3,4,'b')
hold on