function scatPos = get_scatters_positions(numScat,scen_center,interSiteDist)

r=5*interSiteDist/3;

theta = 2 * pi * rand(1,numScat);
rho = r * sqrt(rand(1,numScat));

% generate random points
scatPos(1,:) = scen_center(1)+rho.*cos(theta); 
scatPos(2,:) = scen_center(2)+rho.*sin(theta); 
scatPos(3,:) = zeros(1,numScat);%randi(25,[1 numScat]);