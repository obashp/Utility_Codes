%Generate Co-ordinates for rotated cartesian grid.
dx = 1.0;
dy = dx;

tCrd = [[0,0]; [1,0]; [1,1]; [0,1]];
Crd = [dx*tCrd(:,1) dy*(tCrd(:,2))];% - tCrd(:,1))];

%Generate Shape functions and Derivatives for Gauss-Quadrature domain
nQuad = 4;
isqrt3 = 1.0/sqrt(3);
gW = [1, 1, 1, 1]';
gP = [-isqrt3, -isqrt3; isqrt3 -isqrt3; isqrt3 isqrt3; -isqrt3 isqrt3];

Psi1=@(xi,eta) 0.25*(1.0-xi).*(1.0-eta);
Psi2=@(xi,eta) 0.25*(1.0+xi).*(1.0-eta);
Psi3=@(xi,eta) 0.25*(1.0+xi).*(1.0+eta);
Psi4=@(xi,eta) 0.25*(1.0-xi).*(1.0+eta);

dPsi11=@(xi,eta) -0.25*(1.0-eta);
dPsi21=@(xi,eta) 0.25*(1.0-eta);
dPsi31=@(xi,eta) 0.25*(1.0+eta);
dPsi41=@(xi,eta) -0.25*(1.0+eta);

dPsi12=@(xi,eta) -0.25*(1.0-xi);
dPsi22=@(xi,eta) -0.25*(1.0+xi);
dPsi32=@(xi,eta) 0.25*(1.0+xi);
dPsi42=@(xi,eta) 0.25*(1.0-xi);

JPsi=@(xi,eta) [[dPsi11(xi,eta), dPsi21(xi,eta), dPsi31(xi,eta), dPsi41(xi,eta)]; ...
                [dPsi12(xi,eta), dPsi22(xi,eta), dPsi32(xi,eta), dPsi42(xi,eta)]];

%Visualize Co-ordinates and Gauss Points in Physical domain
figure()
CrdGP =   Crd(1,:).*Psi1(gP(:,1),gP(:,2)) ...
        + Crd(2,:).*Psi2(gP(:,1),gP(:,2)) ...
        + Crd(3,:).*Psi3(gP(:,1),gP(:,2)) ...
        + Crd(4,:).*Psi4(gP(:,1),gP(:,2));
vertices = [Crd(1,:);Crd(2,:);Crd(3,:);Crd(4,:);Crd(1,:)];
plot(vertices(:,1), vertices(:,2));
hold on;
plot(CrdGP(:,1),CrdGP(:,2),'x');
hold off;

%Compute Jacobians at each of the quadrature points
Jac = zeros(nQuad,1);
for i = 1:nQuad
   Jac(i) = det(JPsi(gP(i,1),gP(i,2))*Crd);
end

%Function to Evaluate integral of
F=@(x,y) (x.^2).*y;

%Evaluate integral using Gauss Quadrature
FGP = F(CrdGP(:,1), CrdGP(:,2));
S = sum(gW.*FGP.*Jac);









