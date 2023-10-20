function scatPos = get_scatters_positions(numScat,scen_center,interSiteDist)

r=5*interSiteDist/3;

theta = 2 * pi * rand(1,numScat);
rho = r * sqrt(rand(1,numScat));

scatPos = [scen_center(1)+rho.*cos(theta); scen_center(2)+rho.*sin(theta)]; % generate random points 