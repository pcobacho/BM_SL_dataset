clear,clc,close all;

% Convertir el código de color hexadecimal a RGB (rango 0 a 1)
color_personalizado = sscanf('D2E0FB', '%2x').' / 255;

% Parámetros del escenario
radio_sector = 50;  % Radio del sector hexagonal en metros

% Coordenadas de la estación base
base_station_x = 0;
base_station_y = 0;

% Ángulos de los sectores en grados
angulo_sector1 = 0;
angulo_sector2 = 120;
angulo_sector3 = 240;

% Coordenadas de los centros de los sectores
sector1_x = base_station_x + radio_sector * cosd(angulo_sector1);
sector1_y = base_station_y + radio_sector * sind(angulo_sector1);

sector2_x = base_station_x + radio_sector * cosd(angulo_sector2);
sector2_y = base_station_y + radio_sector * sind(angulo_sector2);

sector3_x = base_station_x + radio_sector * cosd(angulo_sector3);
sector3_y = base_station_y + radio_sector * sind(angulo_sector3);

% Crear una figura
figure;
hold on;

% Dibujar los sectores hexagonales con relleno de color personalizado
draw_hexagon(sector1_x, sector1_y, radio_sector, color_personalizado);
draw_hexagon(sector2_x, sector2_y, radio_sector, color_personalizado);
draw_hexagon(sector3_x, sector3_y, radio_sector, color_personalizado);

% Etiquetar los sectores
text(sector1_x - 20, sector1_y + 20, 'Sector 1', 'FontSize', 12);
text(sector2_x + 20, sector2_y - 20, 'Sector 2', 'FontSize', 12);
text(sector3_x - 20, sector3_y - 20, 'Sector 3', 'FontSize', 12);

% Dibujar la estación base
plot(base_station_x, base_station_y, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); hold on;
text(base_station_x + 10, base_station_y + 10, 'gNB', 'FontSize', 12);

% Configuración de la gráfica
axis equal;
title('Escenario de Comunicaciones Móviles 2D');
xlabel('Distancia (metros)');
ylabel('Distancia (metros)');
grid on;

hold off;

% Función para dibujar un hexágono con centro en (x, y) y radio r
function draw_hexagon(x, y, r, color)
    angles = linspace(0, 2*pi, 7);
    x_points = r * cos(angles) + x;
    y_points = r * sin(angles) + y;
    fill(x_points, y_points, color);  % Rellenar con color personalizado
end
