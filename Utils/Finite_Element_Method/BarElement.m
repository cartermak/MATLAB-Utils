% -------------------------------------------------------------------------
% Copyright (c) 2020 Carter Mak
% -------------------------------------------------------------------------
% Description:
%	Object type to store and process data for a 2D bar element in a finite
%	element model.
% -------------------------------------------------------------------------

classdef BarElement
	
	properties(SetAccess = private, GetAccess = public)
		E;                % Young's Modulus
		A;                % Cross-sectional Area
		L;                % Length
		id_i;             % First connection index
		id_j;             % Second connection index
		pos_id_i;         % Position vector of node i
		pos_id_j;         % Position vector of node j
		phi;              % Angle of bar w.r.t. global coordinate frame
		K_element_global; % Element stiffness matrix in global coordinates
	end
	
	methods
		function obj = BarElement(E,A,id_i,id_j,pos_id_i,pos_id_j)
            % -------------------------------------------------------------
            % Description:
            %     Constructor to instantiate and define parameters of a bar
            %     element
            % Inputs:
            %     E        - Young's Modulus
            %     A        - Cross-sectional Area
            %     id_i     - First connection index
            %     id_j     - Second connection index
            %     pos_id_i - Position vector of node i
            %     pos_id_j - Position vector of node j
            % -------------------------------------------------------------
            
            % Set bar properties
			obj.E = E;
			obj.A = A;
			obj.id_i = id_i;
			obj.id_j = id_j;
			obj.pos_id_i = pos_id_i;
			obj.pos_id_j = pos_id_j;
            
            % Calculate length using coorinates of endpoints
			obj.L = sqrt(...
				(pos_id_j(2)-pos_id_i(2))^2 ...
			 +(pos_id_j(1)-pos_id_i(1))^2);
			obj.phi = atan2(...
				pos_id_j(2) - pos_id_i(2),...
				pos_id_j(1) - pos_id_i(1));
			
            % Calculate element stiffness matrix in global coordinates
			c = cos(obj.phi);
			s = sin(obj.phi);
			
			obj.K_element_global = (E*A/obj.L)*([...
				c^2,s*c,-c^2,-s*c;...
				s*c,s^2,-s*c,-s^2;...
				-c^2,-s*c,c^2,s*c;...
				-s*c,-s^2,s*c,s^2;...
				]);
			
		end
	end
end