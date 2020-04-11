function [rj, tj ] = rayInt_qs2( ir, it, rc, tc, V)
%ray interesection test with triangles designed to meet the requirements
% of the 'queryset' (qs) function for D. Engiwirda's aabb BVH 

%   [RJ, TJ] = INCIRCLE(IR,IT,RC,TC) compute the pairwise in-
%   tersections between the rays RC(IR,:) and the triangles
%   TC(IT,:). V are vertices of the triangles. 
%   [RJ,TJ] are pairs of intersections, such that the ray
%   RJ(K) intersects with the  triangle TJ(K).

    li = cell(length(it),1);
    lj = cell(length(it),1);

    ry = rc(ir,:) ;
 
    for ii = +1 : length(ir)
        [~, pos, faceInds] = intersectLineMesh3d(ry(ii,:), V, tc(it,:));
        
        lj{ii} = faceInds(pos>0);
        li{ii} = ii * ones(length(lj{ii}),1);
    end

    rj = ir(vertcat(li{:})) ;
    tj = it(vertcat(lj{:})) ;




end

