function thickness = computeThickness_1(all_surface_section, size_points_path, all_points_path, h, beta, q_max, w, M, v, thickness)
for i5 = 1:size(all_surface_section,1)
    % 遍历当前pass的所有surface_section
    for b5 = 1:size(all_surface_section(i5,:),2)
        % 遍历当前surface_section的所有三角形
        for c5 = 1:size(all_surface_section{i5,b5},1)
            all1_q_1 = 0;
            % 遍历所有的路径点
            for d5 = 2:size_points_path
                % t:当前三角形的形心、法向量及索引信息
                t = cell2mat(all_surface_section{i5,b5}(c5,[19:21,10:12,24]));
                % p:当前路径点的信息
                p = cell2mat(all_points_path{i5,b5});
                % l：路径点指向当前三角形形心的向量
                l = t(1:3)-p(1:3);
                cos_gamma = sum(l.*-t(4:6))./(norm(l).*norm(t(4:6)));
                cos_theta = abs(sum(l.*p(4:6))./(norm(l).*norm(p(4:6))));
                ref = p(1:3) + h.*p(4:6)./norm(p(4:6));
                s0 = [(p(6).*l(1).*(ref(3)-p(3)) + p(5).*l(1).*(ref(2)-p(2)) + p(1).*(p(5).*l(2)+p(6).*l(3)) + p(4).*ref(1).*l(1))./(p(4).*l(1)+p(5).*l(2)+p(6).*l(3)),(p(6).*l(2).*(ref(3)-p(3)) + p(4).*l(2).*(ref(1)-p(1)) + p(2).*(p(4).*l(1)+p(6).*l(3)) + p(5).*ref(2).*l(2))./(p(4).*l(1)+p(5).*l(2)+p(6).*l(3)),(p(5).*l(3).*(ref(2)-p(2)) + p(4).*l(3).*(ref(1)-p(1)) + p(3).*(p(4).*l(1)+p(5).*l(2)) + p(6).*ref(3).*l(3))./(p(4).*l(1)+p(5).*l(2)+p(6).*l(3))];
                double r
                r = norm(s0-ref);
                q_0 = Beta(r,beta,q_max,w);
                mat_M = cell2mat(M{1,i5});
                %%%%%%%%%%%%%%%%%%%%%%
                if d5-1 > size(mat_M,2)
                    break
                else
                %%%
                q_1 = mat_M(d5-1).*CDB(q_0,h,l,cos_gamma,cos_theta)./v;
                all1_q_1 = all1_q_1 + q_1;
                thickness(t(7)) = all1_q_1;
                end
            end
        end
    end
end
end