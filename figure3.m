c = [0;0];
G = [-2 0 1; 0 2 1];
E = [1 0 3; 0 1 1];
disp(2.*E)
pZ = polyZonotope(c, G, [], 2.*E);

figure; hold on;
plot(pZ, [1, 2], 'b', 'Splits', 10);
