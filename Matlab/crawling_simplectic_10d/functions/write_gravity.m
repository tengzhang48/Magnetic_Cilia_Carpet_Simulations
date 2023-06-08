% write the information into a data file for lammps simulation

function write_gravity(Mclass, gravity)

N_atom = Mclass.N_i(1);                % number of atom
T_atom = Mclass.T_i(1);                % type of atom

mg = zeros(N_atom,2);

for ii=1:N_atom
    mg(ii,1) = ii;
    atom_type = Mclass.data_atom(ii,2);
    mg(ii,2) = Mclass.M_atom(atom_type)*gravity;
end



fid = fopen('gravity.txt', 'w');
fprintf(fid,'%d\t', N_atom);
fprintf(fid,'\n');
fprintf(fid,'%d %15.10f\n',mg');

fclose(fid);


