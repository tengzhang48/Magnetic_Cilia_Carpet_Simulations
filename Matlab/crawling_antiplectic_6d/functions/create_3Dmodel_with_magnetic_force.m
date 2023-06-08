function [lattice_model] = create_3Dmodel_with_magnetic_force(obj)

% lattice model
lattice_model=model();
mu = obj.mu;
lambda = obj.lambda;
rho = obj.rho;
coord = obj.coord;
connect = obj.connect;
id = obj.id;


nnode = size(obj.coord,1); 
Ne= size(obj.connect,1);                                                        % total element number

m = zeros(nnode,1);                                                              % store the mass of each node
ks = zeros(Ne,28);                                                              % total spring constant
kbond = zeros(28*Ne,1);
ms = zeros(Ne,8);                                                               % total mass at element node
v0 = zeros(Ne,1);                                                               % initial volume/area
bond_left = zeros(Ne,28);                                                       % store left bond nodes
bond_right = zeros(Ne,28);                                                      % store right bond nodes
n_left =[ones(1,7),2*ones(1,6),3*ones(1,5),4*ones(1,4),5*ones(1,3),6,6,7];
n_right =[2:8,3:8,4:8,5:8,6:8,7:8,8];

% Magnetic force
km = zeros(Ne,24);                                                           % store the stiffness of the magnetic force
fm = zeros(nnode,3); 

%%%
Nint = 8;
[xi,w] = Gaussian_point_weight(Nint,8,3);

% Loop over the element
for ii = 1:Ne
    
    node = connect(ii,1:8);
    xa = coord(node,1);
    ya = coord(node,2);
    za = coord(node,3);
    kse = zeros(8,8);    % store coefficent of springs
    me = zeros(8,8);     % store coefficent of springs
    kme = zeros(8,3);  
    
    % Loop over the Gaussian Points
    for k = 1:Nint
        z1 = xi(1,k);
        z2 = xi(2,k);
        z3 = xi(3,k);
        [dNdx,yita] = shape_function_derivative(z1,z2,z3,xa,ya,za);
        kx = dNdx*transpose(dNdx);
        kse = kse +kx*yita*w(k);
        Nme = shape_function(z1,z2,z3);                                        % shape function
        me = me + connect(ii,10)*Nme*Nme'*yita*w(k);
        kme = kme + dNdx*yita*w(k);
    end
    
    [ve,~,~,~] = pressure_volume(connect(ii,9),connect(ii,11),1,xa,ya,za);
    ms(ii,:) = sum(me);
    
    
    % lumped mass matrix
    for jj = 1:28
    ks(ii,jj) = -kse(n_left(jj),n_right(jj));    % bond stiffness
    % factor 2 is for lammps convention
    kbond(Ne*(jj-1)+ii) = connect(ii,9)*ks(ii,jj)/2;   
    end
    bond_left(ii,:) = node(n_left);    % 
    bond_right(ii,:) = node(n_right);   % node 2
    
    %
    v0(ii) = ve;     % initial area
    km(ii,1:8) = kme(:,1);                                                  % dNdx_1
    km(ii,9:16) = kme(:,2);
    km(ii,17:end) = kme(:,3);
    
end

bond = [bond_left(:),bond_right(:)];                                        % bond connection
                                                        % bond stiffness
% clear the old data
clear bond_left bond_right;


% loop over element to add the mass to each node
% loop over element to calculate the magnetic force
beta =  90/180*pi;
beta0 =  7/8*pi;
beta1 =  -30/180*pi;
gamma = 0/180*pi;
B_a = [cos(gamma),0,sin(gamma)]';  
B_a = B_a*40.0; 

for ii = 1:Ne
    
    node = connect(ii,1:8); % node   % calculate the center of the element
    xyz_c = sum(coord(node,:))/8;
    
    % Residual mangetic is a function of position
    if xyz_c(3) <=4
    x_l=fix((xyz_c(1)-6)/4); 
    theta = -30*x_l/180*pi;
    alpha = beta+theta;
    B_r = [cos(alpha),0,sin(alpha)]';
    B_r = B_r*40;
    % m(node) = m(node) + ms(ii,:)';
    kme = [km(ii,1:8);km(ii,9:16);km(ii,17:end)];
    BBe = B_a*B_r'; % 1/mu0*B_a*B_r has the unit of N/m^2
    fme = (1.0*BBe/1.000E+03)*kme;      % normlized by stress unit
    fm(node,1) =  fm(node,1) + fme(1,:)';
    fm(node,2) =  fm(node,2) + fme(2,:)';
    fm(node,3) =  fm(node,3) + fme(3,:)';
    end
    
    % Residual mangetic is a function of position
%     x_l=fix((xyz_c(1)-10)/15); % for 6X1 array
%     beta=beta0-pi/6*x_l; % for 6X1 array
%     theta = (xyz_c(1)-10-15*x_l)/10*1.75*pi; % for 6X1 array
%     alpha = beta - theta;
%     B_r = [cos(alpha),0,sin(alpha)]';
%     B_r = B_r*40;
    %m(node) = m(node) + ms(ii,:)';
%     kme = [km(ii,1:8);km(ii,9:16);km(ii,17:end)];
%     BBe = B_a*B_r'; % 1/mu0*B_a*B_r has the unit of N/m^2
     % f=-du/dx (there is a double negative signs here)
%     fme = (1.0*BBe/0.1)*kme;    % normlized by stress unit
%     fm(node,1) =  fm(node,1) + fme(1,:)';
%     fm(node,2) =  fm(node,2) + fme(2,:)';
%     fm(node,3) =  fm(node,3) + fme(3,:)';
end

% combine the same bond by adding their stiffness together
bond_1 = bond;
% switch the bond index to let the first one be the smaller one
l = bond(:,1) > bond(:,2);
bond(l,1) = bond_1(l,2);
bond(l,2) = bond_1(l,1);
% find the unique bond, ic store the bond number
[bond_lammps,~,ic] = unique(bond,'rows');
kbond_lammps = zeros(length(bond_lammps),1);
% loop over the original bonds and add the stiffness of bonds with the same
% number together
for ii = 1:length(bond)
    bb = ic(ii);
    kbond_lammps(bb) = kbond_lammps(bb) + kbond(ii);
end


% clear the old data
clear bond kbond;

% Identify the distinguished bond types
tol = 1e-3;
[k_bond,~,ty_bond] = uniquetol(kbond_lammps,tol);

% Identify the distinguished dbond types
% tol = 1e-9;
% [improper_v,~,ty_improper] = uniquetol(v0,tol);
improper_v = v0;
% v_dbond is the initial volume of the cube

% loop over element to add the mass to each node
for ii = 1:Ne
    node = connect(ii,1:8);                                                   % node
    m(node) = m(node) + ms(ii,:)';
end


% Identify the distinguished atom types
tol = 1e-1;
[lattice_model.M_atom,~,ty_atom] = uniquetol(m,tol);


x_atom = coord;
clear coord;
angle = [];
ty_angle = [];

lattice_model.N_i = [length(x_atom),length(bond_lammps),0,0,size(connect,1)];
lattice_model.T_i = [length(lattice_model.M_atom),length(k_bond),0,0,length(improper_v)];

lattice_model.data_atom = [id,ty_atom,x_atom];
clear ty_atom x_atom;
lattice_model.data_bond = [ty_bond,bond_lammps];
clear ty_bond bond_lammps;
lattice_model.data_angle = [ty_angle,angle];
clear ty_angle angle;
lattice_model.data_dbond = [];
lattice_model.data_improper = connect(:,1:4);
clear ty_improper;


% Coeffs
bond_t = [1:lattice_model.T_i(2)]'; % Bond type serial
lattice_model.bond_coeffs = [k_bond,0*bond_t]; % the initial length is zero

lattice_model.dbond_coeffs = []; % Dihedral type serial

improper_mu = connect(:,9);
improper_lambda = connect(:,11);
improper_2 = connect(:,5:8);
lattice_model.improper_coeffs = [improper_mu,improper_lambda,improper_v,improper_2];
lattice_model.magnetic_forces = fm;

clear connect;


end

function N = shape_function(z1,z2,z3)
% quadratic shape function in reduced space

N1 = 1/8*(1-z1)*(1-z2)*(1-z3);
N2 = 1/8*(1+z1)*(1-z2)*(1-z3);
N3 = 1/8*(1+z1)*(1+z2)*(1-z3);
N4 = 1/8*(1-z1)*(1+z2)*(1-z3);
N5 = 1/8*(1-z1)*(1-z2)*(1+z3);
N6 = 1/8*(1+z1)*(1-z2)*(1+z3);
N7 = 1/8*(1+z1)*(1+z2)*(1+z3);
N8 = 1/8*(1-z1)*(1+z2)*(1+z3);
N = [N1; N2; N3; N4;N5;N6;N7;N8];

end

function [dNdx,yita] = shape_function_derivative(z1,z2,z3,xa,ya,za)

% derivative of shape function in reduced space

dNdz = [[ -((z2 - 1)*(z3 - 1))/8, -(z1/8 - 1/8)*(z3 - 1), -(z1/8 - 1/8)*(z2 - 1)]
    [  ((z2 - 1)*(z3 - 1))/8,  (z1/8 + 1/8)*(z3 - 1),  (z1/8 + 1/8)*(z2 - 1)]
    [ -((z2 + 1)*(z3 - 1))/8, -(z1/8 + 1/8)*(z3 - 1), -(z1/8 + 1/8)*(z2 + 1)]
    [  ((z2 + 1)*(z3 - 1))/8,  (z1/8 - 1/8)*(z3 - 1),  (z1/8 - 1/8)*(z2 + 1)]
    [  ((z2 - 1)*(z3 + 1))/8,  (z1/8 - 1/8)*(z3 + 1),  (z1/8 - 1/8)*(z2 - 1)]
    [ -((z2 - 1)*(z3 + 1))/8, -(z1/8 + 1/8)*(z3 + 1), -(z1/8 + 1/8)*(z2 - 1)]
    [  ((z2 + 1)*(z3 + 1))/8,  (z1/8 + 1/8)*(z3 + 1),  (z1/8 + 1/8)*(z2 + 1)]
    [ -((z2 + 1)*(z3 + 1))/8, -(z1/8 - 1/8)*(z3 + 1), -(z1/8 - 1/8)*(z2 + 1)]];

%

dxdz = [xa'*dNdz;ya'*dNdz;za'*dNdz];
yita = det(dxdz);                       % Jacobian of the mapping

% derivative of shape function in physical space
dNdx = dNdz/dxdz;

end

function [rse,fse,kse] = bond_length_derivative(kb,xye)

% rse is the bond length
% fse is the force/drivative of the spring length squre respect to the nodal position
% r^2 ---> r*dr/dx

ke = [eye(3),-eye(3);-eye(3),eye(3)];
rse = sqrt(xye'*ke*xye);                                                % spring length
fse = -kb*ke*xye;                                                       % forces (corresponds to energy with a 1/2 factor)
kse = ke*kb;                                                            % stiffness matrix

end


function [ve,pe,fpe,kpe] = pressure_volume(mu,lambda,v0,xa,ya,za)

% ve is the volume
% dse is the drivative of the area respect to the nodal position
% s0 is the initial area
% the functional energy form is -mu*lnJ + lambda/2*(lnJ)^2
% pe is the enrgy associated with the area change
% fpe is the force applied on each nodes due to area change
% kpe is the stiffness matrix

kab = [ 0     0     0     0     0     0     0     0
    0     0    -1    -1     1     1     0     0
    0     1     0    -1     0     0     0     0
    0     1     1     0    -1     0     0    -1
    0    -1     0     1     0    -1     0     1
    0    -1     0     0     1     0     0     0
    0     0     0     0     0     0     0     0
    0     0     0     1    -1     0     0     0];

nc = [  1 4 3 2 5 8 7 6;
    2 1 4 3 6 5 8 7;
    3 2 1 4 7 6 5 8;
    4 3 2 1 8 7 6 5;
    5 8 7 6 1 4 3 2;
    6 5 8 7 2 1 4 3;
    7 6 5 8 3 2 1 4;
    8 7 6 5 4 3 2 1];

ddve = zeros(24,24);

for i = 1:4
    ci = nc(:,i);
    ddve(9:16,17:24) = ddve(9:16,17:24) + kab(ci,ci)*xa(i);     % ep_231
    ddve(17:24,9:16) = ddve(17:24,9:16) - kab(ci,ci)*xa(i);     % ep_321
    ddve(17:24,1:8) = ddve(17:24,1:8) + kab(ci,ci)*ya(i);       % ep_312
    ddve(1:8,17:24) = ddve(1:8,17:24) - kab(ci,ci)*ya(i);       % ep_132
    ddve(1:8,9:16) = ddve(1:8,9:16) + kab(ci,ci)*za(i);         % ep_123
    ddve(9:16,1:8) = ddve(9:16,1:8) - kab(ci,ci)*za(i);         % ep_213
    
    cci = nc(:,i+4);
    ddve(9:16,17:24) = ddve(9:16,17:24) - kab(cci,cci)*xa(i+4);     % ep_231
    ddve(17:24,9:16) = ddve(17:24,9:16) + kab(cci,cci)*xa(i+4);     % ep_321
    ddve(17:24,1:8) = ddve(17:24,1:8) - kab(cci,cci)*ya(i+4);       % ep_312
    ddve(1:8,17:24) = ddve(1:8,17:24) + kab(cci,cci)*ya(i+4);       % ep_132
    ddve(1:8,9:16) = ddve(1:8,9:16) - kab(cci,cci)*za(i+4);         % ep_123
    ddve(9:16,1:8) = ddve(9:16,1:8) + kab(cci,cci)*za(i+4);         % ep_213
    
end

xyz = [xa;ya;za];

dve = ddve*xyz;

ve = dve'*xyz;                                                      % volume
dve = dve/24;                                                       % dV/dx_i^a
ddve = ddve/12;                                                     % dV/dx_i^ax_j^b
ve = ve/72;

J = ve/v0;

if J <0
    aa = 0;
end
pe = v0*(-mu*log(J)+lambda/2*(log(J)).^2);
pp = v0/J*(-mu+lambda*log(J));                                          % pressure due to area change (dpe/dJ)
dpp = v0/J^2*(mu + lambda*(1-log(J)));                                  % derivative of pressure respect to J

% calculate force
fpe = -pp*dve/v0;                                                       % fpe = -dpe/dJ*dJ/dx_i^a and dJ/dx_i^a = dve/v0
% calculate stiffness
kpe = pp*ddve/v0 + dpp*dve*dve'/v0^2;                                   % the order is node first and than x_1 to x_3

end

function [xi,w] = Gaussian_point_weight(Nint,shape,dim)

w = [1.,1.,1.,1.,1.,1.,1.,1.];                                  % Gaussian weights
xi = zeros(3,8);                                                % Gaussian points
x1D = [-1/sqrt(3),1/sqrt(3)];

for kk = 1:2
    for jj = 1:2
        for ii = 1:2
            nn = 4*(kk-1) + 2*(jj-1) + ii;
            xi(1,nn) = x1D(ii);
            xi(2,nn) = x1D(jj);
            xi(3,nn) = x1D(kk);
        end
    end
end

end