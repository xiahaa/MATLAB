function [r,p,y] = rot2euler(rot)

    p = -asin(rot(1,3));
    r = atan2(rot(2,3),rot(3,3));
    y = atan2(rot(1,2),rot(1,1));

end