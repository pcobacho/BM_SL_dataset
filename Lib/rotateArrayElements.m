function Cr=rotateArrayElements(C,theta,axis)

R_x = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)]; % rotation with respect to the x-axis
R_y = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)]; % rotation with respect to the y-axis
R_z = [cosd(theta) -sind(theta) 0; sind(theta) cosd(theta) 0; 0 0 1]; % rotation with respect to the z-axis

if strcmpi(axis, 'x')
    Cr = C' * R_x';
elseif strcmpi(axis, 'y')
    Cr = C' * R_y';
elseif strcmpi(axis, 'z')
    Cr = C' * R_z';
else
    disp('Invalid "axis" value');
end
