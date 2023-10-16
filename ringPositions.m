clear,clc,close all;

% Convert hexadecimal color code to RGB (range 0 to 1)
colors = [sscanf('FF8080', '%2x').'; sscanf('D2E0FB', '%2x').';...
    sscanf('CDFAD5', '%2x').';sscanf('FFCF96', '%2x').';...
    sscanf('F9F3CC', '%2x').';sscanf('BEADFA', '%2x').';...
    sscanf('EDB7ED', '%2x').'] / 255;

% Scenario Parameters
interSiteDist = 150;  % Inter site distance
cell_radio = 50;  % radio of each cell/sector

% Base station coordinates
gNBpos = [0,0];

% cell angles (degrees)
sector1_ang = 0;
sector2_ang = 120;
sector3_ang = 240;

% gNB positions coordinates
num_gNB = 7;  % number of gNBs
gNB_ang = linspace(0, 360, num_gNB);  % gNB position angles (in degrees)
gNB_ang = gNB_ang(1:end-1);
puntos_x = [0, gNBpos(1) + interSiteDist * cosd(gNB_ang)];
puntos_y = [0, gNBpos(2) + interSiteDist * sind(gNB_ang)];

figure;
hold on;

plot(0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
text(10, 10, 'Estación Base Original', 'FontSize', 12);

for i = 1:num_gNB
    % Coordinates of sector centers for current base station
    current_base_station_x = puntos_x(i);
    current_base_station_y = puntos_y(i);
    
    current_sector1_x = current_base_station_x + cell_radio * cosd(sector1_ang);
    current_sector1_y = current_base_station_y + cell_radio * sind(sector1_ang);

    current_sector2_x = current_base_station_x + cell_radio * cosd(sector2_ang);
    current_sector2_y = current_base_station_y + cell_radio * sind(sector2_ang);

    current_sector3_x = current_base_station_x + cell_radio * cosd(sector3_ang);
    current_sector3_y = current_base_station_y + cell_radio * sind(sector3_ang);

    % Drawing the base station and its sectors with custom color fills
    plot(current_base_station_x, current_base_station_y, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    draw_hexagon(current_sector1_x, current_sector1_y, cell_radio, colors(i,:));
    draw_hexagon(current_sector2_x, current_sector2_y, cell_radio, colors(i,:));
    draw_hexagon(current_sector3_x, current_sector3_y, cell_radio, colors(i,:));
end

% Print the points
plot(puntos_x, puntos_y, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'w');

% Graph config
axis equal;
title('Hexagonal Cells Scenario');
xlabel('X (m)');
ylabel('Y (m)');
grid on;

hold off;

% Function to draw a hexagon with center at (x, y) and radius r
function draw_hexagon(x, y, r, color)
    angles = linspace(0, 2*pi, 7);
    x_points = r * cos(angles) + x;
    y_points = r * sin(angles) + y;
    fill(x_points, y_points, color);
end
