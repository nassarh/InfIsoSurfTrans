%%
% figure
tiledlayout(1, 6, 'TileSpacing', 'none', 'Padding', 'none');

%% surface definition
% curvilinear coordinates
pts = 10; % ~pts per unit cell
uc = 4;  % ~number of unit cells
x = linspace(-1/2,1/2,1+2*uc*pts)-100*eps;
y = linspace(-uc/2,uc/2,1+2*uc*pts)-100*eps;

dx = x(2)-x(1);
dy = y(2)-y(1);

[x,y] = meshgrid(x,y);

% path & profile (graphs)
% df = sign(sin(x*2*pi))+1/2*sign(sin(2*x*2*pi)); % non-centrosymmetric
df = sign(sin(x*2*pi));
dg = sign(sin(y*2*pi));

% inclination
theta = 0;
s = sin(theta);
c = cos(theta);

% path & profile (adapted base)
alpha1 = df*0+sign(cos(x*2*pi));
alpha2 = df;

beta1 = c+s*dg;
beta2 = -s+c*dg;

% surface (canonical base)
f = antid(df,2)*dx;
g = antid(dg,1)*dy;

r1 = antid(alpha1,2)*dx;
r2 = y - s*f;
r3 = g + c*f;

% plot surface
nexttile
myplot(r1,r2,r3,pts);

% deflection nominal amplitude
d = 0.4;

%% stretching
% adapted base
u = antid(alpha2.^2./alpha1,2)*dx;
v = -antid(beta2.^2./beta1,1)*dy;
w = -antid(alpha2,2)*dx+antid(beta2,1)*dy;

% canonical base
t1 = u;
t2 = c*v - s*w;
t3 = s*v + c*w;

nexttile
myplot(r1+d*t1,r2+d*t2,r3+d*t3,pts)

%% twisting
% adapted base
u = -antid(alpha2,2)*dx.*antid(beta1,1)*dy+antid(beta1.*antid(beta2,1)*dy-beta2.*antid(beta1,1)*dy,1)*dy;
v = -antid(alpha1,2)*dx.*antid(beta2,1)*dy+antid(alpha1.*antid(alpha2,2)*dx-alpha2.*antid(alpha1,2)*dx,2)*dx;
w = antid(alpha1,2)*dx.*antid(beta1,1)*dy;

% canonical base
t1 = u;
t2 = c*v - s*w;
t3 = s*v + c*w;

nexttile
myplot(r1+d*t1,r2+d*t2,r3+d*t3,pts)

%% bending
% adapted base
% solution s
us = -antid(beta2,1)*dy.*antid(alpha2.^2./alpha1,2)*dx-antid(alpha2.*antid(alpha2.^2./alpha1,2)*dx-alpha2.^2./alpha1.*antid(alpha2,2)*dx,2)*dx;
vs = -antid(alpha2,2)*dx.*antid(beta2.^2./beta1,1)*dy-antid(beta2.*antid(beta2.^2./beta1,1)*dy-beta2.^2./beta1.*antid(beta2,1)*dy,1)*dy;
ws = antid(alpha2,2)*dx.*antid(beta2,1)*dy+antid(beta1.*antid(beta2.^2./beta1,1)*dy-beta2.*antid(beta2,1)*dy,1)*dy+antid(alpha1.*antid(alpha2.^2./alpha1,2)*dx-alpha2.*antid(alpha2,2)*dx,2)*dx;

% solution p
up = -antid(beta1,1)*dy.*antid(alpha2.^2./alpha1,2)*dx;
vp = -antid(alpha2,2)*dx.*antid(beta2,1)*dy-antid(beta2.*antid(beta2,1)*dy-beta2.^2./beta1.*antid(beta1,1)*dy,1)*dy+antid(alpha1.*antid(alpha2.^2./alpha1,2)*dx-alpha2.*antid(alpha2,2)*dx,2)*dx;
wp = antid(alpha2,2)*dx.*antid(beta1,1)*dy+antid(beta1.*antid(beta2,1)*dy-beta2.*antid(beta1,1)*dy,1)*dy;

% solution p bar
ub = -antid(beta2,1)*dy.*antid(alpha2,2)*dx-antid(alpha2.*antid(alpha2,2)*dx-alpha2.^2./alpha1.*antid(alpha1,2)*dx,2)*dx+antid(beta1.*antid(beta2.^2./beta1,1)*dy-beta2.*antid(beta2,1)*dy,1)*dy;
vb = -antid(alpha1,2)*dx.*antid(beta2.^2./beta1,1)*dy;
wb = antid(beta2,1)*dy.*antid(alpha1,2)*dx+antid(alpha1.*antid(alpha2,2)*dx-alpha2.*antid(alpha1,2)*dx,2)*dx;

% solution twisting
ut = -antid(alpha2,2)*dx.*antid(beta1,1)*dy+antid(beta1.*antid(beta2,1)*dy-beta2.*antid(beta1,1)*dy,1)*dy;
vt = -antid(alpha1,2)*dx.*antid(beta2,1)*dy+antid(alpha1.*antid(alpha2,2)*dx-alpha2.*antid(alpha1,2)*dx,2)*dx;
wt = antid(alpha1,2)*dx.*antid(beta1,1)*dy;

% solution (out of plane) bending
u = c*us+s*up;
v = c*vs+s*vp;
w = c*ws+s*wp;

% canonical base
t1 = u;
t2 = c*v - s*w;
t3 = s*v + c*w;

nexttile
myplot(r1+d*t1,r2+d*t2,r3+d*t3,pts)

% the other solution (in plane) bending
u = -(c^2-s^2)*ub-2*c*s*ut;
v = -(c^2-s^2)*vb-2*c*s*vt;
w = -(c^2-s^2)*wb-2*c*s*wt;

% canonical base
t1 = u;
t2 = c*v - s*w;
t3 = s*v + c*w;

nexttile
myplot(r1+d*t1,r2+d*t2,r3+d*t3,pts)

% solution (in plane) bending
u = -s*us+c*up;
v = -s*vs+c*vp;
w = -s*ws+c*wp;

% canonical base
t1 = u;
t2 = c*v - s*w;
t3 = s*v + c*w;

nexttile
myplot(r1+d*t1,r2+d*t2,r3+d*t3,pts)

% check isometry
% compute (linearized) change in metric
clc
% stretch in y
i = 1;
max(abs(diff(r1,1,i).*diff(t1,1,i)+diff(r2,1,i).*diff(t2,1,i)+diff(r3,1,i).*diff(t3,1,i)),[],"all")
% stretch in x
i = 2;
max(abs(diff(r1,1,i).*diff(t1,1,i)+diff(r2,1,i).*diff(t2,1,i)+diff(r3,1,i).*diff(t3,1,i)),[],"all")
% shear
max(abs(diff(r1(:,1:end-1),1,1).*diff(t1(1:end-1,:),1,2)+diff(r2(:,1:end-1),1,1).*diff(t2(1:end-1,:),1,2)+diff(r3(:,1:end-1),1,1).*diff(t3(1:end-1,:),1,2)+...
    diff(t1(:,1:end-1),1,1).*diff(r1(1:end-1,:),1,2)+diff(t2(:,1:end-1),1,1).*diff(r2(1:end-1,:),1,2)+diff(t3(:,1:end-1),1,1).*diff(r3(1:end-1,:),1,2)),[],"all")


%% save figure
exportgraphics(gcf,'my_figure.png', 'Resolution', 900)

%% auxiliary function
function ad = antid(alpha,axis)
    ad = cumsum(alpha, axis);
    ad = ad - ad((size(alpha,1)*size(alpha,2)+1)/2);
end

function mesh = myplot(r1,r2,r3,pts)
    % plot parameters
    FaceAlpha=0.8;
    LineWidth=1;
    FaceColor=[0.8,0.6,1];
    EdgeAlpha=0.1;
    
    % plot surface and secondary mesh
    mesh = surf(r1,r2,r3,FaceAlpha=FaceAlpha,LineWidth=LineWidth,FaceColor=FaceColor,EdgeAlpha=EdgeAlpha);
    
    % plot primary mesh
    hold on
    plot3(r1(:,1:pts:end),r2(:,1:pts:end),r3(:,1:pts:end),'k','LineWidth',1)
    plot3(r1(1:pts:end,:)',r2(1:pts:end,:)',r3(1:pts:end,:)','k','LineWidth',1)
    hold off
    
    % set axes
    axis equal
    set(gca, 'XTickLabel', {});
    set(gca, 'YTickLabel', {});
    set(gca, 'ZTickLabel', {});
    
    % axis off

    a = 2;
    xlim([-a,a])
    ylim([-a,a])
    zlim([-a,a])
end
