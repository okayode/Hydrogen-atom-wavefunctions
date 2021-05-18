clear;
clc;

fprintf('\n#################################################################################\n')
fprintf('                         Hydrogen atom wavefunctions                                 ')
fprintf('\n#################################################################################\n')
fprintf('                         Kayode Olumoyin                                           \n')
fprintf('                         COMS 7100                                                 \n')
fprintf('                         Project 4                                                 \n')
fprintf('\n#################################################################################\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long e % using the highest matlab precision

set(0,'DefaultFigureWindowStyle','docked');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

fprintf('\n#################################################################################\n')
fprintf('Part 1: The Radial Function R(r) \n')
fprintf('\n#################################################################################\n')

%% Part 1: The Radial Function R(r)
k1 = 1;
figure
nmax = 4;
syms x r
for n = 1:nmax
    for l = 0:n-1
    alp = simplify(assoLaguerrePolynomial(n-l-1,2*l+1,x));
    rf = simplify(radialFunc(n,l,r));
    nI = round(int((r*radialFunc(n,l,r))^2,r,[0 Inf]));
    rradialFunc2 = (r*radialFunc(n,l,r))^2;
    fprintf('Part 1 \n')
    fprintf('Radial function: \n')
    fprintf('n = %i, l = %i\n', n,l)
    fprintf('n-l-1 = %i, 2l+1 = %i\n',n-l-1,2*l+1)
    fprintf('associated Laguerre polynomial = %s\n', char(alp))
    fprintf('function Rnl(r) = %s\n', char(rf))
    fprintf('Normalization Integral = %s\n', char(nI))
    fprintf('#############################################################################\n')
    rr = linspace(0,11*n,1000);
    yy = subs(radialFunc(n,l,r),{r},{rr});
    yyy = subs(rradialFunc2,{r},{rr});
    
    subplot(2,5,k1)
    plot(rr,yy, rr,yyy)
    pbaspect([1 1 1])
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    title(['n = ',num2str(n), ', l = ', num2str(l)])
    k1 = k1 + 1;
    
    sgtitle('The Radial Function R(r)')    
    lgd = legend({'Rnl(r)','r^2Rnl(r)^2'},'Orientation','horizontal','Location','southoutside');
    lgd.FontSize = 5.0; 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n#################################################################################\n')
fprintf('Part 2: The angular function Y(theta, phi) \n')
fprintf('\n#################################################################################\n')

%% Part 2: % The angular function Y(theta, phi) (1)
k2 = 1;
figure
syms x y z p t
lmax = nmax-3;
for l = 0:lmax
    for m = -l:l 
        fprintf('Part 2 \n')
        fprintf('Angular function: \n')
        fprintf('l = %i, m = %i\n', l,m)
        cartY = simplify(cartRealsphHarm(l,m,x,y,z));
        fprintf('cartesian function = %s\n', char(cartY))
        int_1 = int(realSpHarm100(l,m,t,p).^2.*sin(t),t,[0 pi]);
        normInt = round(int(int_1,p,[0 2*pi]));
        fprintf('Normalization Integral in spherical frame = %s\n', char(normInt))
        fprintf('###########################################################################\n')
      
        xR = linspace(-1.5,1.5,100);
        yR = linspace(-1.5,1.5,100);
        zR = linspace(-1.5,1.5,100);
        [xx,yy,zz] = meshgrid(xR,yR,zR);
        %ym = cartRealsphHarm(l,m,xx,yy,zz);
        [theta, phi, nr] = mycart2sph(xx, yy, zz);
        ym = realSpHarm100(l,m,theta,phi).*exp(-nr);
        
        subplot(2,3,k2);
        p1 = patch(isosurface(xx,yy,zz,ym,-0.15));
        isonormals(xx,yy,zz,ym,p1);
        set(p1,'FaceColor','blue','EdgeColor','none');
        
        p2 = patch(isosurface(xx,yy,zz,ym,0.15));
        isonormals(xx,yy,zz,ym,p2);
        set(p2,'FaceColor','red','EdgeColor','none'); 
        
        title(['l = ',num2str(l), ', m = ', num2str(m)])
        k2 = k2 + 1;
        q = 1.2;       
        xlim([-q q])
        ylim([-q q])
        zlim([-q q])        
        xlabel('x axis')
        ylabel('y axis')
        zlabel('z axis')
        daspect([1 1 1]);
        axis square
        gca;
        view(3)
        camlight 
        lighting gouraud
    end
end

%% Part 2: % The angular function Y(theta, phi) (2)
k22 = 1;
syms x y z p t
figure
l = 2;
for m = -l:l
    fprintf('Part 2 \n')
    fprintf('Angular function: \n')
    fprintf('l = %i, m = %i\n', l,m)
    cartY = cartRealsphHarm(l,m,x,y,z);
    fprintf('cartesian function = %s\n', char(cartY))
    int_1 = int(realSpHarm100(l,m,t,p).^2.*sin(t),t,[0 pi]);
    normInt = round(int(int_1,p,[0 2*pi]));
    fprintf('Normalization Integral in spherical frame = %s\n', char(normInt))
    fprintf('#############################################################################\n')

    xR = linspace(-1.5,1.5,100);
    yR = linspace(-1.5,1.5,100);
    zR = linspace(-1.5,1.5,100);
    [xx,yy,zz] = meshgrid(xR,yR,zR);
    %ym = cartRealsphHarm(l,m,xx,yy,zz);
    [theta, phi, nr] = mycart2sph(xx, yy, zz);
    ym = realSpHarm100(l,m,theta,phi).*exp(-nr);

    subplot(2,3,k22);
    p1 = patch(isosurface(xx,yy,zz,ym,-0.15));
    isonormals(xx,yy,zz,ym,p1);
    set(p1,'FaceColor','blue','EdgeColor','none');

    p2 = patch(isosurface(xx,yy,zz,ym,0.15));
    isonormals(xx,yy,zz,ym,p2);
    set(p2,'FaceColor','red','EdgeColor','none'); 

    title(['l = ',num2str(l), ', m = ', num2str(m)])

    k22 = k22 + 1; 
    q = 1.5;
    xlim([-q q])
    ylim([-q q])
    zlim([-q q]) 
    xlabel('x axis')
    ylabel('y axis')
    zlabel('z axis')
    daspect([1 1 1]);
    axis square
    gca;
    view(3)
    camlight 
    lighting gouraud
end

%% Part 2: % The angular function Y(theta, phi) (3)
k23 = 1;
syms x y z p t
figure
l = 3;
for m = -l+1:l-1
    fprintf('Part 2 \n')
    fprintf('Angular function: \n')
    fprintf('l = %i, m = %i\n', l,m)
    cartY = cartRealsphHarm(l,m,x,y,z);
    fprintf('cartesian function = %s\n', char(cartY))
    int_1 = int(realSpHarm100(l,m,t,p).^2.*sin(t),t,[0 pi]);
    normInt = round(int(int_1,p,[0 2*pi]));
    fprintf('Normalization Integral in spherical frame = %s\n', char(normInt))
    fprintf('###########################################################################\n')

    xR = linspace(-1.5,1.5,100);
    yR = linspace(-1.5,1.5,100);
    zR = linspace(-1.5,1.5,100);
    [xx,yy,zz] = meshgrid(xR,yR,zR);
    %ym = cartRealsphHarm(l,m,xx,yy,zz);
    [theta, phi, nr] = mycart2sph(xx, yy, zz);
    ym = realSpHarm100(l,m,theta,phi).*exp(-nr);

    subplot(2,3,k23);
    p1 = patch(isosurface(xx,yy,zz,ym,-0.15));
    isonormals(xx,yy,zz,ym,p1);
    set(p1,'FaceColor','blue','EdgeColor','none');

    p2 = patch(isosurface(xx,yy,zz,ym,0.15));
    isonormals(xx,yy,zz,ym,p2);
    set(p2,'FaceColor','red','EdgeColor','none'); 

    title(['l = ',num2str(l), ', m = ', num2str(m)])

    k23 = k23 + 1;
    q = 2.0;
    xlim([-q q])
    ylim([-q q])
    zlim([-q q])        
    xlabel('x axis')
    ylabel('y axis')
    zlabel('z axis')
    daspect([1 1 1]);
    axis square
    gca;
    view(3)
    camlight 
    lighting gouraud
end

%% Part 2: % The angular function Y(theta, phi) (4)
k24 = 1;
syms x y z p t
figure
l = 3;
m = [-3 3];
for ii = 1:2
    fprintf('Part 2 \n')
    fprintf('Angular function: \n')
    fprintf('l = %i, m = %i\n', l,m(ii))
    cartY = cartRealsphHarm(l,m(ii),x,y,z);
    fprintf('cartesian function = %s\n', char(cartY))
    int_1 = int(realSpHarm100(l,m(ii),t,p).^2.*sin(t),t,[0 pi]);
    normInt = round(int(int_1,p,[0 2*pi]));
    fprintf('Normalization Integral in spherical frame = %s\n', char(normInt))
    fprintf('###########################################################################\n')

    xR = linspace(-1.5,1.5,100);
    yR = linspace(-1.5,1.5,100);
    zR = linspace(-1.5,1.5,100);
    [xx,yy,zz] = meshgrid(xR,yR,zR);
    %ym = cartRealsphHarm(l,m(ii),xx,yy,zz);
    [theta, phi, nr] = mycart2sph(xx, yy, zz);
    ym = realSpHarm100(l,m(ii),theta,phi).*exp(-nr);

    subplot(2,3,k24);
    p1 = patch(isosurface(xx,yy,zz,ym,-0.15));
    isonormals(xx,yy,zz,ym,p1);
    set(p1,'FaceColor','blue','EdgeColor','none');

    p2 = patch(isosurface(xx,yy,zz,ym,0.15));
    isonormals(xx,yy,zz,ym,p2);
    set(p2,'FaceColor','red','EdgeColor','none'); 

    title(['l = ',num2str(l), ', m = ', num2str(m(ii))])
    k24 = k24 + 1;
    q = 1.4;
    xlim([-q q])
    ylim([-q q])
    zlim([-q q])        
    xlabel('x axis')
    ylabel('y axis')
    zlabel('z axis')
    daspect([1 1 1]);
    axis square
    gca;
    view(3)
    camlight 
    lighting gouraud
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Part 3: % The complete hydrgen wavefunction

fprintf('\n#################################################################################\n')
fprintf('Part 3: The complete hydrgen wavefunction \n')
fprintf('\n#################################################################################\n')

%% 
syms x y z p t r
k3 = 1;
%figure
for n = 1:nmax
    for l = 0:n-1
        for m = -l:l
            fprintf('Part 3 \n')
            fprintf('n = %i, l = %i, m = %i\n',n, l,m)
            rf = radialFunc(n,l,r);
            fprintf('radial function = %s\n', char(rf))
            ay = cartRealsphHarm(l,m,x,y,z);
            fprintf('Angular function = %s\n', char(ay))
            psii = cartCompWaveFunc(n,l,m,x,y,z);
            fprintf('Total function = %s\n', char(psii))
            int_1 = int(conj(compWaveFunc(n,l,r,m,t,p))...
                .*compWaveFunc(n,l,r,m,t,p)...
                .*sin(t).*r.^2,p,[0 2*pi]);
            int_2 = int(int_1,t,[0 pi]);
            normInt = round(int(int_2,r,[0 Inf]));
            fprintf('Normalization check = %s\n', char(normInt))
            fprintf('#######################################################################\n')            
            q = 3*n^2+4;
            xR = linspace(-q,q,200);
            yR = linspace(-q,q,200);
            zR = linspace(-q,q,200);
            [xxx,yyy,zzz] = meshgrid(xR,yR,zR);
            % psi = cartCompWaveFunc(n,l,m,xxx,yyy,zzz);
            [theta, phi, nr] = mycart2sph(xxx, yyy, zzz);
            psi = compWaveFunc(n,l,nr,m,theta,phi);
            
            figure
            subplot(1,3,1);
            contourslice(xxx,yyy,zzz,psi,[],[],0,20);
            xlim([-q q])
            ylim([-q q])
            xlabel('x axis')
            ylabel('y axis')
            %view ([0 0 90])
            view(0,90)
            daspect([1 1 1]);
            
            subplot(1,3,2);
            contourslice(xxx,yyy,zzz,psi,[],0,[],20);
            xlim([-q q])
            zlim([-q q])
            xlabel('x axis')
            zlabel('z axis')
            %view ([0 90 0])
            view(0,0)
            daspect([1 1 1]);
            
            subplot(1,3,3);
            contourslice(xxx,yyy,zzz,psi,0,[],[],20);
            ylim([-q q])
            zlim([-q q])
            ylabel('y axis')
            zlabel('z axis')
            %view ([90 0 0])
            view(90,0)
            daspect([1 1 1]);
                   
            sgtitle(['n = ',num2str(n), ', l = ',num2str(l), ', m = ', num2str(m)])           
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% functions starts here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% The Radial Function R(r)
function R = radialFunc(n,l,r)
a = factorial(n-l-1);
b = 2*n*factorial(n+l)^3;
N = sqrt(a/b); 
L =  assoLaguerrePolynomial(n-l-1,2*l+1,2*r/n);
R = (2/n)^(l+(3/2))*N*(r.^l).*exp(-r/n).*L;
end

%% associated Laguerre polynomial
function L = assoLaguerrePolynomial(p,q,r)
L = 0;
for i = 0:p
    L = L + factorial(p+q)*nchoosek(p+q,p-i)/factorial(i)*(-r).^i;
end
end

%% The angular function Y(theta, phi)
function Y = angularFunc(l,m,theta,phi)
a = (2*l+1)*factorial(l-m);
b = 4*pi*factorial(l+m);
c = sqrt(a/b);

% Plm = legendre(l,cos(theta));
% if l ~= 0
%     Plm = reshape(Plm(m+1,:,:),size(phi));
% end
PP = assoLegendrePolynomial(l,m,cos(theta));
Y = (-1)^(m+abs(m))*c.*PP.*exp(1i*m*phi);
% Y = (-1)^(m+abs(m))*c.*Plm.*exp(1i*m*phi);
end

%% associated Legendre polynomial 
% function L = assoLegendrePolynomial(l,m,x)
% L = zeros(size(x));
% nC = floor(1/2*l - 1/2*abs(m));
% for n = 0:nC
%     PC = (-1)^n*nchoosek(l-2*n,abs(m))*nchoosek(l,n)*nchoosek(2*l-2*n,l);
%     L = L + PC*x.^(l - 2*n - abs(m));
% end
% L = (-1)^m*(1-x.^2).^(abs(m)/2).*(factorial(abs(m))/2^l*L);
% end

function L = assoLegendrePolynomial(l,m,x)
L = 0;
    for r = 0 : floor(1/2*l - 1/2*abs(m))
        L = L + (-1)^r*nchoosek(l - 2*r, abs(m))*nchoosek(l, r)*...
            nchoosek(2*l - 2*r, l)*x.^(l - 2*r - abs(m));
    end
    L = (1 - x.^2).^(abs(m)/2).*(factorial(abs(m))/2^l*L);
end

%% The Real Spherical harmonic function (eqn 96 in proj4 notes)
function y = realSpHarm96(l,m,theta,phi)
if m < 0 
    Y = angularFunc(l,abs(m),theta,phi);
    y = sqrt(2)*imag(Y);
elseif m == 0
    Y = angularFunc(l,0,theta,phi);
    y = Y;
elseif m > 0
    Y = angularFunc(l,abs(m),theta,phi);
    y = sqrt(2)*real(Y);
end
end

%% The Real Spherical harmonic function (eqn 100 in proj4 notes)
function y = realSpHarm100(l,m,theta,phi)

PP = assoLegendrePolynomial(l,abs(m),cos(theta));
if m == 0
    a = (2*l+1)*factorial(l-abs(m));
    b = 2*(1+1)*pi*factorial(l+abs(m));
    Nlm = sqrt(a/b);
else
    a = (2*l+1)*factorial(l-abs(m));
    b = 2*(1)*pi*factorial(l+abs(m));
    Nlm = sqrt(a/b);
end

kk = (-1)^abs(m)*Nlm*PP;
if m > 0 
    y = kk.*cos(abs(m).*phi);
elseif m == 0
    y = kk;
elseif m < 0
    y = kk.*sin(abs(m).*phi);
end
end

%% cartRealSphHarm

function cartY = cartRealsphHarm(l,m,x,y,z)
r = sqrt(x.^2+y.^2+z.^2);
theta = acos(z./r);
phi = atan2(y,x); 
% cartY = realSpHarm100(l,m,theta,phi);
cartY = realSpHarm100(l,m,theta,phi);
end

function [theta, phi, r] = mycart2sph(x, y, z)
r = sqrt(x.^2+y.^2+z.^2);
theta = acos(z./r); 
%phi = atan(y./x);
phi  = atan2(y,x);     
end

%% The complete hydrogen wavefunction

function psi = compWaveFunc(n,l,r,m,theta,phi)
% y = realSpHarm100(l,m,theta,phi);
y = realSpHarm100(l,m,theta,phi);
R = radialFunc(n,l,r);
% psi = conj(R).*R.*y.*y;
psi = R.*y;
end

function psi = cartCompWaveFunc(n,l,m,x,y,z)
r = sqrt(x.^2+y.^2+z.^2);
theta = acos(z./r);
phi   = atan2(y,x); 
% y = realSpHarm100(l,m,theta,phi);
y = realSpHarm100(l,m,theta,phi);
R = radialFunc(n,l,r);
% psi = conj(R).*R.*y.*y;
psi = R.*y;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% functions ends here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




