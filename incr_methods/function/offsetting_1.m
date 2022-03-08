function [point_surface_section,normal_surface_section] = offsetting_1(offset,point_surface_section,normal_surface_section)
normal_surface_section = normal_surface_section ./ norm(normal_surface_section);
normal_surface_section = offset.* normal_surface_section;
point_surface_section = [point_surface_section,1]';

operator = [eye(3) normal_surface_section';[0 0 0 1]];
point_surface_section = operator * point_surface_section;
point_surface_section = point_surface_section(1:3,:)';
normal_surface_section = -1 .* normal_surface_section;
end