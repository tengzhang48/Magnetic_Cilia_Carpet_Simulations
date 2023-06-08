classdef model < handle  %Inheritate from a handle class in order to change properties using methods
    
    properties
        N_i = zeros(1,5);
        T_i = zeros(1,5);
        M_atom
        data_atom
        data_bond
        bond_coeffs
        data_angle
        data_dbond
        dbond_coeffs
        data_improper
        improper_coeffs
        magnetic_forces
    end
    
    methods
        function obj = model()
        
        end
             
    end
end
