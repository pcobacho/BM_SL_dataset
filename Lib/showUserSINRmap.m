function showUserSINRmap(sinrdB,userPos,num_users)

sinrdB = abs(sinrdB);

% mapa de colores
cMap_log = interp1( [0; 1/3; 3/4; 1], flip([0.5373    0.8784    0.4667; 1.0000    0.9176    0.4980; 0.9961    0.7294    0.3098; 0.9490    0.3843    0.4196]), linspace(0,1,256) ); %progresion color verde-amarillo-naranja-rojo logaritmica

% Crea una figura y establece los ejes
figure;
hold on;

% Dibuja los círculos coloreados según los valores de los usuarios
rango = max(sinrdB)-min(sinrdB);
for i = 1:num_users
    color_index = round(((sinrdB(i) - min(sinrdB)) / abs(rango)) * 254) + 1;  % Mapea el valor al índice del colormap
    color = colormap(cMap_log);
    scatter(userPos(1,i), userPos(2,i), 50, color(color_index, :), 'filled');
end

% Establece el colorbar
colorbar;
title('SINR Map');
xlabel('X Coordinate (m)'); ylabel('Y Coordinate (m)');
grid on

hold off;