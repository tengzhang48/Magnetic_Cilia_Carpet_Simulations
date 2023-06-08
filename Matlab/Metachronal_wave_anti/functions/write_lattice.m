% write the information into a data file for lammps simulation

function write_lattice(fn,Mclass, b_s)

N_atom = Mclass.N_i(1);                % number of atom
N_bond = Mclass.N_i(2);                % number of bond
N_angle = Mclass.N_i(3);               % number of angle
N_dihedral = Mclass. N_i(4);            % number of dihedral
N_improper = Mclass.N_i(5);            % number of improper

T_atom = Mclass.T_i(1);                % type of atom
T_bond = Mclass.T_i(2);                % type of bond
T_angle = Mclass.T_i(3);               % type of angle
T_dihedral = Mclass.T_i(4);            % type of dihedral
T_improper = Mclass.T_i(5);            % type of improper

xlo = b_s(1,1); xhi = b_s(1,2);
ylo = b_s(2,1); yhi = b_s(2,2);
zlo = b_s(3,1); zhi = b_s(3,2);

% open the file with write permission
dim = size(Mclass.data_atom,2)-2;

if dim == 3
filename = [fn,'-cube-lattice.lam']; % create file name
else
filename = [fn,'-quad-lattice.lam'];
end
fid = fopen(filename, 'w');

% LAMMPS Description

fprintf(fid,'LAMMPS Description\n\n');

% atom, bond, angles, dihedrals, impropers

fprintf(fid,'%d\t', N_atom);
fprintf(fid,'atoms\n');
fprintf(fid,'%d\t', N_bond);
fprintf(fid,'bonds\n');
fprintf(fid,'%d\t', N_angle);
fprintf(fid,'angles\n');
fprintf(fid,'%d\t', N_dihedral);
fprintf(fid,'dihedrals\n');
fprintf(fid,'%d\t', N_improper);
fprintf(fid,'impropers\n\n');

fprintf(fid,'%d\t', T_atom);
fprintf(fid,'atom types\n');
fprintf(fid,'%d\t', T_bond);
fprintf(fid,'bond types\n');
fprintf(fid,'%d\t', T_angle);
fprintf(fid,'angle types\n');
fprintf(fid,'%d\t', T_dihedral);
fprintf(fid,'dihedral types\n');
fprintf(fid,'%d\t', T_improper);
fprintf(fid,'improper types\n\n');

% size of box

fprintf(fid,'%10.8f\t%10.8f\t', xlo,xhi);
fprintf(fid,'xlo xhi\n');
fprintf(fid,'%10.8f\t%10.8f\t', ylo,yhi);
fprintf(fid,'ylo yhi\n');
fprintf(fid,'%10.8f\t%10.8f\t', zlo,zhi);
fprintf(fid,'zlo zhi\n\n');

% Masses
fprintf(fid,'Masses\n\n');
t_atom = [1:T_atom]';
fprintf(fid,'%d %10.8f\n',[t_atom, Mclass.M_atom]');
fprintf(fid,'\n');

% Bond coeffs
fprintf(fid,'Bond Coeffs\n\n');
t_bond = [1:T_bond]';
fprintf(fid,'%d %12.9f %12.9f\n',[t_bond, Mclass.bond_coeffs]');
fprintf(fid,'\n');

% improper angles coeffs
fprintf(fid,'Improper Coeffs\n\n');
i_improper = [1:N_improper]';
if dim ==3
    fprintf(fid,'%d %12.9f %12.9f %12.9f %d %d %d %d\n',[i_improper,Mclass.improper_coeffs]');
else
    fprintf(fid,'%d %12.9f %12.9f %12.9f\n',[i_improper,Mclass.improper_coeffs]');
end
fprintf(fid,'\n');

% Atoms
% molecular style

fprintf(fid,'Atoms\n\n');
i_atom = [1:N_atom]';
if dim == 3
    data_atom = [i_atom,Mclass.data_atom];
else
    z_coord = ones(N_atom,1);
    data_atom = [i_atom,Mclass.data_atom,0.5*z_coord];
end
fprintf(fid,'%d %d %d %15.10e %15.10f %15.10f\n',data_atom');


% bonds
if N_bond > 0
    
    fprintf(fid,'\nBonds\n\n');
    i_bond = [1:N_bond]';
    data_bond = [i_bond,Mclass.data_bond];
    fprintf(fid,'%d %d %d %d\n',data_bond');

end


% angles
if N_angle > 0
    fprintf(fid,'\nAngles\n\n');
    i_angle = [1:N_angle]';
    data_angle = [i_angle,Mclass.data_angle];
    fprintf(fid,'%d %d %d %d %d\n',data_angle');
  
end
% dihedral angles

if N_dihedral > 0
    fprintf(fid,'\nDihedrals\n\n');
    i_dihedral = [1:N_dihedral]';
    data_dihedral = [i_dihedral,Mclass.data_dbond];
    fprintf(fid,'%d %d %d %d %d %d\n',data_dihedral');
end

% N_improper

if N_improper > 0
    fprintf(fid,'\nImpropers\n\n');
    i_improper = [1:N_improper]';

    data_improper = [i_improper,i_improper,Mclass.data_improper];

    fprintf(fid,'%d %d %d %d %d %d\n',data_improper');
end


fclose(fid);