function C = arrayElementCoordinates(numRows,numCols,lambda)

% Distancia entre puntos
distancia_x = 0.5*lambda; % Distancia en la dimensión x
distancia_y = 0.5*lambda; % Distancia en la dimensión y

% Crear coordenadas en el eje x e y centradas en (0,0)
y = ((1:numCols) - (numCols + 1) / 2) * distancia_x;  % Coordenadas en el eje x
z = ((1:numRows) - (numRows + 1) / 2) * distancia_y;  % Coordenadas en el eje y

% Crear malla de puntos
[Y, Z] = meshgrid(y, z);

% Reorganizar las coordenadas en una matriz de nxm filas y 2 columnas
C = [zeros(1,length(Y(:)));Y(:)'; Z(:)'];

% Graficar los puntos
% figure();
% plot(C(2, :), C(3, :), 'o');
% xlabel('Coordenada Y'); ylabel('Coordenada Z');
% ttl = sprintf(['Panel de ' num2str(numRows) 'x' num2str(numCols) ' puntos equidistantes centrado en (0,0)']);
% title(ttl); axis square; grid on;
