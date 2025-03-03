% observables and system matrix A
load('Roessler');
samplingTime = 0.05;

% reachability parameter
params.R0 = zonotope(interval([-0.05;-8.45;-0.05],[0.05;-8.35;0.05]));
params.tFinal = 6;

% specification 
i = 3;
c = [-1 0 0]; d = 7.1; %d = 6.375 - 0.025 * i;
    
hs = halfspace(c,d);
spec = specification(hs,'safeSet');

% visualization
figure; hold on; box on;
t = 0:samplingTime:params.tFinal;
pZ = polyZonotope(g(taylm(params.R0)));
R = expm(A*t(80))*pZ;

R = polyZonotope(R.c, 3*R.G, R.GI, (R.E), R.id);
plot(R,[1,2],'k', LineWidth=2);
xlim([-7.4,-5.9]);
ylim([2.7,3.2]);
plot(halfspace(-spec.set.c,-spec.set.d),[1,2],'FaceColor','r', ...
                                       'FaceAlpha',0.5,'EdgeColor','none');
plot(zonotope(R),[1,2],'--k', LineWidth=2);
set(gcf,'Units','centimeters','OuterPosition', [0, 0, 12, 10]);
xlabel('x_1'); ylabel('x_2');

figure; hold on; box on;
temp = splitLongestGen(R);
xlim([-7.4,-5.9]);
ylim([2.7,3.2]);
plot(halfspace(-spec.set.c,-spec.set.d),[1,2],'FaceColor','r', ...
                                       'FaceAlpha',0.5,'EdgeColor','none');
plot(temp{1},[1,2],'c', LineWidth=2);
plot(zonotope(temp{1}),[1,2],'--c', LineWidth=2);
plot(temp{2},[1,2],'b', LineWidth=2);
plot(zonotope(temp{2}),[1,2],'--b', LineWidth=2);
set(gcf,'Units','centimeters','OuterPosition', [0, 0, 12, 10]);
xlabel('x_1'); ylabel('x_2');