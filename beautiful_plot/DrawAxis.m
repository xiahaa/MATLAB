%
% Copyright ? 2012, The Massachusetts Institute of Technology. All rights reserved. 
%
% THE LICENSOR EXPRESSLY DISCLAIMS ANY AND ALL WARRANTIES CONCERNING THIS 
% SOFIWARE AND DOCUMENTATION, INCLUDING ANY WARRANTIES OF MERCHANTABILITY, 
% FITNESS FOR ANY PARTICULAR PURPOSE, NON- INFRINGEMENT AND WARRANTIES OF 
% PERFORMANCE, AND ANY WARRANTY THAT MIGHT OTHERWISE ARISE FROM COURSE OF 
% DEALING OR USAGE OF TRADE. NO WARRANTY IS EITHER EXPRESS OR IMPLIED WITH 
% RESPECT TO THE USE OF THE SOFIWARE OR DOCUMENTATION. Under no circumstances 
% shall the Licensor be liable for incidental, special, indirect, direct or 
% consequential damages, or loss of profits, interruption of business, or 
% related expenses which may arise from use of Software or Documentation, 
% including but not limited to those resulting from defects in Software 
% and/or Documentation, or loss or inaccuracy of data of any kind. 
%
% This software is licensed under the "LIMITED RESEARCH LICENSE (SOURCE
% CODE)" as described in the included LICENSE.txt
%
% Please cite the paper below if you are using this software in your work:
% Brookshire, J.; Teller, S. Extrinsic Calibration from Per-Sensor Egomotion. 
%   Robotics: Science and Systems, 2012.
%
function DrawAxis(T, scale, cx, cy, cz)    

    origin = T * [0; 0; 0; 1];
    xaxis = T * [scale; 0; 0; 1];
    yaxis = T * [0; scale; 0; 1];
    zaxis = T * [0; 0; scale; 1];

    a = [origin xaxis];
    line(a(1,:), a(2,:), a(3,:), 'Color', cx, 'LineWidth', 3);
    a = [origin yaxis];
    line(a(1,:), a(2,:), a(3,:), 'Color', cy, 'LineWidth', 3);
    a = [origin zaxis];
    line(a(1,:), a(2,:), a(3,:), 'Color', cz, 'LineWidth', 3);
