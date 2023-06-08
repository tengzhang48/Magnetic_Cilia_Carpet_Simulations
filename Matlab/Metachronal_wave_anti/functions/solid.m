classdef solid < handle  %Inheritate from a handle class in order to change properties using methods
    
    properties
        mu
        rho
        lambda
        coord
        connect
        id
    end
    
    methods
        %Construct an instance of this class
        function obj = solid(varargin)
            
            if nargin == 0
            
            elseif nargin == 13
                obj.mu = varargin{1};
                obj.rho = varargin{2};
                obj.lambda = varargin{3};
                start_id = varargin{10};
                [obj.coord, obj.connect]= generateBox3d(varargin{4},varargin{5},varargin{6},varargin{7},varargin{8},varargin{9});
                numofelements = size(obj.connect,1);
                murow = ones(numofelements,1);
                rhorow = ones(numofelements,1);
                lambdarow = ones(numofelements,1);
                for i=1:numofelements
                    if i<=start_id
                    murow(i)= varargin{1};
                    rhorow(i) = varargin{2};
                    lambdarow(i) = varargin{3};
                    else
                    murow(i)= varargin{11};
                    rhorow(i) = varargin{12};
                    lambdarow(i) = varargin{13};    
                    end
                end
                obj.connect = [obj.connect,murow,rhorow,lambdarow];  
            elseif nargin == 11
                obj.mu = varargin{1};
                obj.rho = varargin{2};
                obj.lambda = varargin{3};
                start_id = varargin{8};
                [obj.coord, obj.connect]= generateBox2d(varargin{4},varargin{5},varargin{6},varargin{7});
                numofelements = size(obj.connect,1);
                murow = ones(numofelements,1);
                rhorow = ones(numofelements,1);
                lambdarow = ones(numofelements,1);
                for i=1:numofelements
                    if i<=start_id
                    murow(i)= varargin{1};
                    rhorow(i) = varargin{2};
                    lambdarow(i) = varargin{3};
                    else
                    murow(i)= varargin{9};
                    rhorow(i) = varargin{10};
                    lambdarow(i) = varargin{11};    
                    end
                end
                obj.connect = [obj.connect,murow,rhorow,lambdarow];   
            elseif nargin == 10
                obj.mu = varargin{1};
                obj.rho = varargin{2};
                obj.lambda = varargin{3};
                start_id = varargin{6};
                [obj.coord, obj.connect]= read_mesh(varargin{4},varargin{5});
                numofelements = size(obj.connect,1);
                murow = ones(numofelements,1);
                rhorow = ones(numofelements,1);
                lambdarow = ones(numofelements,1);
                 for i=1:numofelements
                    if i<=start_id
                    murow(i)= varargin{1};
                    rhorow(i) = varargin{2};
                    lambdarow(i) = varargin{3};
                    else
                    murow(i)= varargin{7};
                    rhorow(i) = varargin{8};
                    lambdarow(i) = varargin{9};    
                    end
                end
                obj.connect = [obj.connect,murow,rhorow,lambdarow];    
            elseif nargin == 9
                obj.mu = varargin{1};
                obj.rho = varargin{2};
                obj.lambda = varargin{3};
                [obj.coord, obj.connect]= generateBox3d(varargin{4},varargin{5},varargin{6},varargin{7},varargin{8},varargin{9});
                numofelements = size(obj.connect,1);
                newrow = ones(numofelements,1);
                obj.connect = [obj.connect,varargin{1}*newrow,varargin{2}*newrow,varargin{3}*newrow];
            elseif nargin == 7
                obj.mu = varargin{1};
                obj.rho = varargin{2};
                obj.lambda = varargin{3};
                [obj.coord, obj.connect]= generateBox2d(varargin{4},varargin{5},varargin{6},varargin{7});
                numofelements = size(obj.connect,1);
                newrow = ones(numofelements,1);
                obj.connect = [obj.connect,varargin{1}*newrow,varargin{2}*newrow,varargin{3}*newrow];
            elseif nargin == 5
                obj.mu = varargin{1};
                obj.rho = varargin{2};
                obj.lambda = varargin{3};
                [obj.coord, obj.connect]= read_mesh(varargin{4},varargin{5});
                numofelements = size(obj.connect,1);
                newrow = ones(numofelements,1);
                obj.connect = [obj.connect,varargin{1}*newrow,varargin{2}*newrow,varargin{3}*newrow];
            else
                error('Solid constructor do not support %d input parameters',nargin);
            end
        end
        
        % To move all the meshes along a vector (dx,dy,dz)
        function move(obj,dx,dy,dz)
            [numRows, numCols] = size(obj.coord);
            for rr = 1 : numRows
                obj.coord(rr,1) = obj.coord(rr,1)+dx;
                obj.coord(rr,2) = obj.coord(rr,2)+dy;
                if numCols == 3
                obj.coord(rr,3) = obj.coord(rr,3)+dz;
                end
            end
        end
        
        function scale(obj,rx,ry,rz)
            [numRows, numCols] = size(obj.coord);
            for rr = 1 : numRows
                obj.coord(rr,1) = obj.coord(rr,1)*rx;
                obj.coord(rr,2) = obj.coord(rr,2)*ry;
                if numCols == 3
                obj.coord(rr,3) = obj.coord(rr,3)*rz;
                end
            end
        end
        
        % To give the element ID a offset
        function offset(obj,id_offset)
                numCols = size(obj.connect,2);
                if numCols == 11
                    obj.connect(:,1:8) = obj.connect(:,1:8)+id_offset;
                elseif numCols == 7
                    obj.connect(:,1:4) = obj.connect(:,1:4)+id_offset;
                end
        end
             
    end
end

function [coord, connect] = generateBox3d(Lx,Ly,Lz,N1,N2,N3)

% The mesh is extruded along the positive z diection

nnode = (N1+1)*(N2+1)*(N3+1);
Ne = N1*N2*N3;
coord = zeros(nnode,3);
connect = zeros(Ne,8);


dx3 = Lz/N3;
x3 = [0:N3]*dx3;

[coord2, connect2] = generateBox2d(Lx,Ly,N1,N2);
nn2 = size(coord2,1);
ne2 = size(connect2,1);
x30 = ones(nn2,1);
for ii = 1:N3
    local_n = (1:nn2) + (ii-1)*nn2;
    coord(local_n,1:2) = coord2;
    coord(local_n,3) = x30*x3(ii);
    local_e = (1:ne2) + (ii-1)*ne2;
    connect(local_e,1:4) = connect2 + (ii-1)*nn2;
    connect(local_e,5:8) = connect2 + ii*nn2;
end
ii = N3 + 1;
local_n = (1:nn2) + (ii-1)*nn2;
coord(local_n,1:2) = coord2;
coord(local_n,3) = x30*x3(ii);

end

% create a 2 d mesh
function [coord, connect] = generateBox2d(Lx,Ly,N1,N2)

% N---
% 2
% 1  2 N

dx1 = Lx/N1;
dx2 = Ly/N2;
Nx1 = [1:N1];
Nx2 = [1:N2];
x1 = [0:N1]*dx1;
x2 = [0:N2]*dx2;
[xy1,xy2] = meshgrid(x1,x2);
coord = [xy1(:),xy2(:)];
[Nxy1,Nxy2] = meshgrid(Nx1,Nx2);
connect = zeros(N1*N2,4);
Nxy = [Nxy1(:),Nxy2(:)];
connect(:,1) = (Nxy(:,1)-1).*(N2+1) + Nxy(:,2);
connect(:,2) = (Nxy(:,1)).*(N2+1) + Nxy(:,2);
connect(:,3) = connect(:,2) + 1;
connect(:,4) = connect(:,1) + 1;

end

% read coord and connectivity from .csv files for coordinates and
% connectivity table
function [coord, connect] = read_mesh(coordinate_file, connectivity_file)

if ischar(coordinate_file)
    coord = csvread(coordinate_file,1,1);
else
    coord = coordinate_file;
end

if ischar(coordinate_file)
    connect = csvread(connectivity_file,1,1);
else
    connect = connectivity_file;
end


end



