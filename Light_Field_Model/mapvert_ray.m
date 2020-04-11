function [tm,im] = mapvert_ray(tr,rays)
%MAPVERT find the tree-to-vertex mappings.
%   [TM,IM] = MAPVERT(TR,PI) returns the tree-to-vertex and 
%   vertex-to-tree mappings for a given aabb-tree TR and a 
%   collection of query vertices PI.
%
%   The tree-to-item mapping TM is a structure representing
%   the intersection of the rays with the tree TR. TM.II 
%   is an M-by-1 array of tree indices and TM.LL is an 
%   M-by-1 cell array of item lists. Specifically, items in 
%   the list TM.LL{JJ} intersect with the node TM.II(JJ).
%
%   The item-to-tree mapping IM is a structure representing
%   the inverse mapping. IM.II is an N-by-1 array of item
%   indices and IM.LL is an N-by-1 cell array of node lists.
%   Specifically, nodes in the list IM.LL{JJ} intersect with
%   the item IM.II(JJ).
%
%   See also QUERYSET, MAPRECT, MAKETREE

%   Darren Engwirda : 2014 --
%   Email           : engwirda@mit.edu
%   Last updated    : 06/04/2017



%----------------------- call SCANTREE to do the actual work
    if (nargout == +1)
       [tm   ] = scantree(tr,rays,@partvert);     
    else
       [tm,im] = scantree(tr,rays,@partvert);
    end
    
end

function [j1,j2] = partvert(rays,b1,b2)
%PARTVERT partition points between boxes B1,B2 for SCANTREE.
%determines if rays intersect boxes B1 and/or B2;


    j1 = true(size(rays,1),1);
    j2 = true(size(rays,1),1);
    
    nrays = size(rays, 1);
    
    for i = 1:nrays
        t_b1 = rayBoxIntersection_mod(rays(i,1:3), rays(i,4:6), b1(1:3), b1(4:6));
        t_b2 = rayBoxIntersection_mod(rays(i,1:3), rays(i,4:6), b2(1:3), b2(4:6));
        if t_b1 == -1  %-1 indicates no interesction
            j1(i) = false;
        end
        if t_b2 == -1
            j2(i) = false;
        end
    end
    

end



