function [sum_model] = assemble_model_with_magnetic_force(S_array)
%ASSEMBLE_MODEL Summary of this function goes here
%   Detailed explanation goes here
sum_solid=solid();
numofsolid=size(S_array,2);
id_offset = 0;

for i=1:numofsolid
    S_array(i).id = ones(size(S_array(i).coord,1),1)*i;
    sum_solid.id = [sum_solid.id;S_array(i).id];
    sum_solid.coord = [sum_solid.coord;S_array(i).coord];
    S_array(i).offset(id_offset);
    sum_solid.connect = [sum_solid.connect;S_array(i).connect];
    id_offset = id_offset + size(S_array(i).coord,1);
end

dim = size(sum_solid.coord,2);

if dim == 3
%sum_model=create_3Dmodel(sum_solid);
sum_model=create_3Dmodel_with_magnetic_force(sum_solid); 
elseif dim == 2
%sum_model=create_2Dmodel(sum_solid); 
sum_model=create_2Dmodel_with_magnetic_force(sum_solid); 
else 
error('Dimension in the solid class is wrong')
end

