function [ n1f ] = selectNormalsVarol(p1, edges, n1s, n2s)
% Greedy algorimthm to select Normals based on neighboring Normals

npts = size(n1s(2).n,2); % number of points

nN = 2*ones(1,npts); % two normals for every point

for viewnum = 2: length(n1s)
    % labels to indicate where the algorithms have run
    runlabels = (nN == 1);        
    n1 = n1s(viewnum).n;
    n2 = n2s(viewnum).n;    
    % Start from the points around runlabels == 1
    % Find edgesum for runlabels = Max-1.

    ind = find(runlabels==0);
%     frs = 1.2; % initial fraction of required normals
    % Until all normals are selected loop:
    while sum(runlabels)<npts    
%         presumlabels = sum(runlabels);
        for i = ind(1):ind(end)
            % Find the neighbors of point i
            [~,pidx1,pidx2]=InvEdgeList(edges,i);
            vertices = [edges(pidx1,2);edges(pidx2,1)]; 
            vertices = round(vertices);
            % if enough points with a single normal:
%             if sum(runlabels(vertices)) > numel(vertices)/frs
                % Get the normal by selecting the best matching one
                % Stack dot products for the 1st normal:
                n1i = n1(:,i);
                n2i = n2(:,i);
                ngN1 = n1(:,vertices);
                ngN2 = n2(:,vertices);
%                 dotsum1 = sum(abs(n1i'*ngN1)+abs(n1i'*ngN2));
%                 dotsum2 = sum(abs(n2i'*ngN1)+abs(n2i'*ngN2));       

                dotsum1 = sum((n1i'*ngN1)+(n1i'*ngN2));
                dotsum2 = sum((n2i'*ngN1)+(n2i'*ngN2));   
                if dotsum1>dotsum2
                    n2(:,i) = [0;0;0]';
                else
                    n1(:,i) = n2(:,i);
                    n2(:,i) = [0;0;0]';
                end
                runlabels(i) = 1;
%             end        
        end
%         if presumlabels == sum(runlabels)
%             frs=frs+0.5;
%         else
%             frs= 1.2;
%         end
    end
    n1f(viewnum).n = n1;
end

end

