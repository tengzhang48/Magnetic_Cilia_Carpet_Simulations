% write the information into a data file for lammps simulation

function write_magnetic_forces(Mclass)

N_atom = Mclass.N_i(1);                % number of atom
fm = Mclass.magnetic_forces;
[numRows,numCols] = size(fm);
fid = fopen('robot_m06_br0ba0_forcex.txt', 'w');
fprintf(fid,'%d\n',length(fm));
i_atom = [1:N_atom]';
data_mforcex = [i_atom,fm(:,1)];
fprintf(fid,'%d %15.10f\n',data_mforcex');
fclose(fid);
fid = fopen('robot_m06_br0ba0_forcey.txt', 'w');
fprintf(fid,'%d\n',length(fm));
data_mforcey = [i_atom,fm(:,2)];
fprintf(fid,'%d %15.10f\n',data_mforcey');
fclose(fid);
if numCols == 3
fid = fopen('robot_m06_br0ba0_forcez.txt', 'w');
fprintf(fid,'%d\n',length(fm));
data_mforcez = [i_atom,fm(:,3)];
fprintf(fid,'%d %15.10f\n',data_mforcez');
fclose(fid);
end






