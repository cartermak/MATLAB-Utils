% -------------------------------------------------------------------------
% Copyright (c) 2020 Carter Mak
% https://github.com/cartermak/MATLAB-Utils
% -------------------------------------------------------------------------
% Description:
%	Object type to store and process a 2D bar finite element model.
% -------------------------------------------------------------------------
classdef FEM
	properties (SetAccess = private,GetAccess = public)
		N;           % Number of 2D nodes
		bars;        % Vector of BarElement objects
		nodes;       % Array of node positions
		forces;      % Array of external nodal forces
		constraints; % Logical array of whether node is constrained
		K_global;    % Global stiffness matrix
		K_reduced;   % Reduced stiffness matrix
		u;           % Vector of nodal displacements
	end
	
	methods (Access = public)
        
		function obj = FEM(N)
            % -------------------------------------------------------------
            % Inputs:
            %     N - number of nodes
            % -------------------------------------------------------------
			obj.N = N;
			obj.bars = BarElement.empty(0,1);
			obj.nodes = zeros(N,2);
			obj.forces = zeros(N,2);
			obj.constraints = false(N,2);
			obj.K_global = zeros(2*N,2*N);
			obj.u = zeros(2*obj.N,1);
        end
        
		function obj = addNode(obj,id,x,y)
            % -------------------------------------------------------------
            % Description:
            %     Define the position of a node in the model
            % Inputs:
            %     id - Node index
            %     x  - horizontal position coordinate
            %     y  - vertical position coordinate
            % -------------------------------------------------------------
			obj.nodes(id,:) = [x,y];
		end
		
		function obj = addBar(obj,bar_id,E,A,id_i,id_j)
            % -------------------------------------------------------------
            % Description:
            %     Add a bar connecting two nodes to the model
            % Inputs:
            %     bar_id - Index of bar to define
            %     E      - Young's Modulus
            %     A      - Bar cross-sectional area
            %     id_i   - Nodal index of first connection
            %     id_j   - Nodal index of second connection
            % -------------------------------------------------------------
			
			% Add bar to the vector of bar elements
			obj.bars(bar_id) = BarElement(...
				E,A,id_i,id_j,...
				obj.nodes(id_i,:),...
				obj.nodes(id_j,:));
			
			% Add element stiffness to global stiffness matrix
			obj.K_global = obj.K_global + obj.getPaddedGlobal(bar_id);
		end
		
		function obj = constrain(obj,id,dim)
            % -------------------------------------------------------------
            % Description:
            %     Define a kinematic zero-displacement constraint on a node
            %     in a given dimenison
            % Inputs:
            %     id  - Index of node to constrain
            %     dim - Index of dimension to constrain (1 for x, 2 for y)
            % -------------------------------------------------------------
			obj.constraints(id,dim) = true;
		end
		
		function obj = addForce(obj,id,dim,f)
            % -------------------------------------------------------------
            % Description:
            %     Add an external load to a given node in a given direction
            % Inputs:
            %     id  - Index of node where force is applied
            %     dim - Index of direction of force (1 for x, 2 for y)
            % -------------------------------------------------------------
			obj.forces(id,dim) = f;
		end
		
		function obj = solve(obj)
            % -------------------------------------------------------------
            % Description:
            %     Solve the finite element problem and store the results
            % -------------------------------------------------------------
			
			% Unwrap matrices and extract vectors of constraints and forces
			cons = obj.constraints';
			cons = cons(:);
			unknowns = ~cons;
			f = obj.forces';
			f = f(:);
			
			% Apply constraints to form reduced system
			obj.K_reduced = obj.K_global(unknowns,unknowns);
			f_reduced = f(unknowns);
			
			% Solve reduced system
			u_reduced = obj.K_reduced\f_reduced;
			
			% Populate full output vector
			obj.u(unknowns) = u_reduced;

		end
		
		function printResults(obj)
            % -------------------------------------------------------------
            % Description:
            %     Print calculated results of the finite element problem
            % -------------------------------------------------------------
            
			for id = 1:length(obj.bars)
				fprintf('Bar %u Element Stiffness Matrix:\n',id)
				disp(obj.bars(id).K_element_global)
			end

			fprintf('Global Master Stiffness Matrix:\n')
			disp(obj.K_global)

			fprintf('Global Reduced Stiffness Matrix:\n')
			disp(obj.K_reduced)

			fprintf('Vector of Nodal Displacements:\n')
			disp(obj.u)
		end
	end
	
	methods (Access = private)
		function K_e_pad = getPaddedGlobal(obj,bar_id)
            % -------------------------------------------------------------
            % Description:
            %     Get the padded global stiffness matrix for the bar with
            %     index `bar_id`
            % -------------------------------------------------------------
			
			% Initialize output matrix
			K_e_pad = zeros(2*obj.N,2*obj.N);
			
			% Get element stiffness matrix in global coordinates, unpadded
			K_e = obj.bars(bar_id).K_element_global;
			
			% Get node indices
			id_i = obj.bars(bar_id).id_i;
			id_j = obj.bars(bar_id).id_j;
			
			% Vector defining indices where we have interest
			idx = [2*id_i-1,2*id_i,2*id_j-1,2*id_j];
			
			% Fill padded matrix
			K_e_pad(idx,idx) = K_e;
			
		end
	end
end