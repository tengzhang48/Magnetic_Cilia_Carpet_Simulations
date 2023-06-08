function [lattice_model] = create_2Dmodel_with_magnetic_force(obj)

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
ks = zeros(Ne,6);                                                              % total spring constant
kbond = zeros(6*Ne,1);
ms = zeros(Ne,4);                                                               % total mass at element node
s0 = zeros(Ne,1);                                                               % initial volume/area
bond_left = zeros(Ne,6);                                                       % store left bond nodes
bond_right = zeros(Ne,6);                                                      % store right bond nodes

% Magnetic force
km = zeros(Ne,8);                                                           % store the stiffness of the magnetic force
fm = zeros(nnode,2); 

%%%
Nint = 4;
[xi,w] = Gaussian_point_weight(Nint,4,2);

% Loop over the element
for ii = 1:Ne
    
    node = connect(ii,1:4);
    xa = coord(node,1);
    ya = coord(node,2);
    kse = zeros(4,4);    % store coefficent of springs
    me = zeros(4,4);     % store coefficent of springs
    se = 0;  
    kme = zeros(4,2);
    
    % Loop over the Gaussian Points
    for k = 1:Nint
        z1 = xi(1,k);
        z2 = xi(2,k);
        [dNdx,yita] = shape_function_derivative(z1,z2,xa,ya);
        kx = dNdx*transpose(dNdx);
        kse = kse +kx*yita*w(k);
        se = se + yita*w(k);
        Nme = shape_function(z1,z2);                                        % shape function
        me = me + connect(ii,6)*Nme*Nme'*yita*w(k);
        kme = kme + dNdx*yita*w(k);
        
    end
    
    ms(ii,:) = sum(me);  % lumped mass matrix
    
    % factor 2 is for lammps convention
    ks(ii,:) = -[kse(1,2),kse(1,3),kse(1,4),kse(2,3),kse(2,4),kse(3,4)]*connect(ii,5)/2;   % bond stiffness
    
    bond_left(ii,:) = [node(1),node(1),node(1),node(2),node(2),node(3)];    % node 1
    bond_right(ii,:) = [node(2),node(3),node(4),node(3),node(4),node(4)];   % node 2
    
    %
    s0(ii) = se;     % initial area
    km(ii,1:4) = kme(:,1);                                                  % dNdx_1
    km(ii,5:8) = kme(:,2);
    
end

bond = [bond_left(:),bond_right(:)]; 
kbond = ks(:); % bond connection
                                                        % bond stiffness
% clear the old data
clear bond_left bond_right;

% loop over element to add the mass to each node
% loop over element to calculate the magnetic force
beta0 =  7/8*pi;
gamma = 0/180*pi;
B_a = [cos(gamma),sin(gamma)]';  
B_a = B_a*40.0; 

for ii = 1:Ne
    
    node = connect(ii,1:4); % node   % calculate the center of the element
    xyz_c = sum(coord(node,:))/4;
    
    x_l=fix((xyz_c(1)-10)/15); % for 6X1 array
    beta=beta0-pi/6*x_l; % for 6X1 array
    theta = (xyz_c(1)-10-15*x_l)/10*1.85*pi; % for 6X1 array
    % Residual mangetic is a function of position
    alpha = beta - theta;
    B_r = [cos(alpha),sin(alpha)]';
    B_r = B_r*40;
    %m(node) = m(node) + ms(ii,:)';
    kme = [km(ii,1:4);km(ii,5:8)];
    BBe = B_a*B_r'; % 1/mu0*B_a*B_r has the unit of N/m^2
    % f=-du/dx (there is a double negative signs here)
    fme = (1.0*BBe/0.126)*kme;    % normlized by stress unit
    fm(node,1) =  fm(node,1) + fme(1,:)';
    fm(node,2) =  fm(node,2) + fme(2,:)';
   
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
%tol = 1e-9;
%[improper_v,~,ty_improper] = uniquetol(s0,tol);
improper_v = s0;
% v_dbond is the initial volume of the cube

% loop over element to add the mass to each node
for ii = 1:Ne
    node = connect(ii,1:4);                                                   % node
    m(node) = m(node) + ms(ii,:)';
end


% Identify the distinguished atom types
tol = 1e-3;
[lattice_model.M_atom,~,ty_atom] = uniquetol(m,tol);


x_atom = coord;
ty_molecule = connect(:,7); 
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
%lattice_model.data_improper = [ty_improper, connect(:,1:4)];
lattice_model.data_improper = connect(:,[1 3 4 2]);
clear ty_improper;


% Coeffs
bond_t = [1:lattice_model.T_i(2)]'; % Bond type serial
lattice_model.bond_coeffs = [k_bond,0*bond_t]; % the initial length is zero

lattice_model.dbond_coeffs = []; % Dihedral type serial

improper_mu = connect(:,5);
improper_lambda = connect(:,7);

lattice_model.improper_coeffs = [improper_mu,improper_lambda,improper_v];
lattice_model.magnetic_forces = fm;

clear connect;


end

function N = shape_function(z1,z2)

% quadratic shape function in reduced space

N1 = 1/4*(1-z1)*(1-z2);
N2 = 1/4*(1+z1)*(1-z2);
N3 = 1/4*(1+z1)*(1+z2);
N4 = 1/4*(1-z1)*(1+z2);
N = [N1; N2; N3; N4];

end

function [dNdx,yita] = shape_function_derivative(z1,z2,xa,ya)

% derivative of shape function in reduced space

dNdz = [ z2/4 - 1/4,   z1/4 - 1/4;
    1/4 - z2/4, - z1/4 - 1/4;
    z2/4 + 1/4,   z1/4 + 1/4;
    -z2/4 - 1/4,   1/4 - z1/4];
%
dxdz = [xa'*dNdz;ya'*dNdz];
yita = det(dxdz);                                                       % Jacobian of the mapping

% derivative of shape function in physical space
dNdx = dNdz/dxdz;

end

function [se,dse] = area_value_derivative(xa,ya,ddse)

% se is the area
% dse is the drivative of the area respect to the nodal position

xye = [xa,ya];
xye = xye';
xye = xye(:);                                                           % re-organize the coordinate
se = 1/2*xye'*ddse*xye;                                                 % area
dse = ddse*xye;

end

function [pe,fpe,kpe] = pressure_area(mu,lambda,s0,xa,ya,ddse)


% s0 is the initial area
% the functional energy form is -mu*lnJ + lambda/2*(lnJ)^2
% pe is the enrgy associated with the area change
% fpe is the force applied on each nodes due to area change
% kpe is the stiffness matrix

[se,dse] = area_value_derivative(xa,ya,ddse);
J = se/s0;
pe = s0*(-mu*log(J)+lambda/2*(log(J)).^2);
pp = s0/J*(-mu+lambda*log(J));                                          % pressure due to area change (dpe/dJ)
dpp = s0/J^2*(mu + lambda*(1-log(J)));                                  % derivative of pressure respect to J

% calculate force
fpe = -pp*dse/s0;                                                       % fpe = dpe/dJ*dJ/dx_i^a and dJ/dx_i^a = dse/s0
% calculate stiffness
kpe = pp*ddse/s0 + dpp*dse*dse'/s0^2;

end

function [rse,fse,kse] = bond_length_derivative(kb,xye)

% rse is the bond length
% fse is the force/drivative of the spring length squre respect to the nodal position
% r^2 ---> r*dr/dx

ke = [eye(2),-eye(2);-eye(2),eye(2)];
rse = sqrt(xye'*ke*xye);                                                % spring length
fse = -kb*ke*xye;                                                       % forces (corresponds to energy with a 1/2 factor)
kse = ke*kb;                                                            % stiffness matrix

end

function [xi,w] = Gaussian_point_weight(Nint,shape,dim)

w = [1.,1.,1.,1.];                                                      % Gaussian weights
xi = zeros(2,4);                                                        % Gaussian points
xi(1,1) = -1/sqrt(3);
xi(2,1) = xi(1,1);
xi(1,2) = -xi(1,1);
xi(2,2) = xi(1,1);
xi(1,3) = xi(1,1);
xi(2,3) = -xi(1,1);
xi(1,4) = -xi(1,1);
xi(2,4) = -xi(1,1);

end