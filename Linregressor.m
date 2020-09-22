classdef Linregressor
% classdef Linregressor
% -------------------------------------------------------------------------
% Carter Mak
% -------------------------------------------------------------------------
% Description:
% Linear regression class for calculating and managing a linear regression
% of input data, along with associated errors.
% -------------------------------------------------------------------------

	properties (SetAccess = private)
		N     % Number of samples in regression
		x     % Raw x-values
		y     % Raw y-values
		m     % Regression slope
		b     % Regression y-intercept
		sig_m % Uncertainty in m
		sig_b % Uncertainty in b
		SSE_y % Absolute uncertainty in y-values (sum of squared error)
		Q     % Q matrix of uncertainties
	end
	
	methods (Access = public) % Publically accessible methods
		function obj = Linregressor(x,y) % Constructor
			% -------------------------------------------------------------
			% Inputs:
			%	x - vector of independent variable values
			%	y - vector of dependent variable values
			%
			% Outputs:
			%	obj - Linregressor object
			% -------------------------------------------------------------
			
			% Add raw values to class object as column vectors
			obj.x = x(:);
			obj.y = y(:);
			
			% Check length of input values
			if length(obj.x) ~= length(obj.y)
				error('Input vectors x and y must have same length')
			end
			
			% Calculate regression
			obj = calcReg(obj);
		end
		
		function [yNew,sig_yNew] = regOutput(obj,xi) % Get estimated y-value(s) from regression
			% -------------------------------------------------------------
			% Inputs:
			%	obj - Linregressor object
			%	xi  - Independent var. index
			%
			% Outputs:
			%	yNew     - Regression output at xi
			%	sig_yNew - Uncertainty in yNew
			% -------------------------------------------------------------
			Ni = length(xi); % Number of input x-values
			
			% Preallocate
			yNew = zeros(Ni,1);
			sig_yNew = zeros(Ni,1);
			
			% Loop over input x-values, get y-value with uncertainty
			for i = 1:Ni
				yNew(i) = getOutput(obj,xi(i));
				sig_yNew(i) = getError(obj,xi(i));
			end
		end
		
		function xNew = regInput(obj,yi) % Get estimated x-value(s) from regression
			% -------------------------------------------------------------
			% Inputs:
			%	obj - Linregressor object
			%	yi  - Dependent var. index
			%
			% Outputs:
			%	xNew     - Inverse regression output at yi
			% -------------------------------------------------------------
			
			% Number of input y-values
			Ni = length(yi);
			
			% Preallocate
			xNew = zeros(Ni,1);
			
			% Loop over input y-values, get estmiated x-values
			for i = 1:Ni
				xNew(i) = getInput(obj,yi(i));
			end
		end
		
	end
	
	methods (Access = private) % Private class methods
		
		function obj = calcReg(obj) % Calculates linear regression parameters
			% -------------------------------------------------------------
			% Inputs:
			%	obj - Linregressor object
			%
			% Outputs:
			%	obj - Linregressor object (with m,b,sig_m,sig_b,Q,SSE_y)
			% -------------------------------------------------------------
			
			% Get number of input value-pairs
			obj.N = length(obj.x);
			
			% A-matrix (ones-padded for non-zero y-intercept)
			A = [obj.x,ones(obj.N,1)];
			
			% Regression coefficients
			t = (A'*A)\(A'*obj.y);
			obj.m = t(1);
			obj.b = t(2);
			
			% Get absolute error in y-values
			obj.SSE_y = obj.sseError;
			
			% Weighting matrix
			W = (1/(obj.SSE_y.^2))*speye(obj.N);
			
			% Q matrix (output coefficient uncertainties)
			obj.Q = inv(A'*W*A);
			obj.sig_m = sqrt(obj.Q(1));
			obj.sig_b = sqrt(obj.Q(4));
			
		end
		
		function sig_yNew = getError(obj,xi) % Get dominant error at given x-value
			% -------------------------------------------------------------
			% Inputs:
			%	obj - Linregressor object
			%	xi  - Independent var. index
			%
			% Outputs:
			%	sig_yNew - Uncertainty in yNew
			% -------------------------------------------------------------
			
			% Get both absolute and extrapolated errors at given point
			sig_yNew = zeros(1,2);
			sig_yNew(1) = obj.SSE_y;
			sig_yNew(2) = extrapError(obj,xi);
			
			% Take the higher error value
			sig_yNew = max(sig_yNew);
		end
		
		function sig_yNew = sseError(obj) % Get SSE error
			% -------------------------------------------------------------
			% Inputs:
			%	obj - Linregressor object
			%
			% Outputs:
			%	sig_yNew - SSE error in obj regression y-values
			% -------------------------------------------------------------
			
			% Error function (Eq. 8.15 from Taylor textbook, pg. 187)
			sig_yNew = sqrt((1/(obj.N-2))*sum((obj.y-obj.b-obj.m*obj.x).^2));
		end
		
		function sig_yNew = extrapError(obj,xi) % Get extrapolated error at a given x-value
			% -------------------------------------------------------------
			% Inputs:
			%	obj - Linregressor object
			%	xi  - Independent var. index
			%
			% Outputs:
			%	sig_yNew - Extrapolation error in regression at xi
			% -------------------------------------------------------------
			
			% Extrapolated error formula from class
			sig_yNew = sqrt([xi,1]*obj.Q*[xi;1]);
		end
		
		function yNew = getOutput(obj,xi) % Get y-value from regression
			% -------------------------------------------------------------
			% Inputs:
			%	obj - Linregressor object
			%	xi  - Independent var. index
			%
			% Outputs:
			%	yNew - Regression value at xi
			% -------------------------------------------------------------
			
			% y = m*x + b
			yNew = obj.b + obj.m*xi;
		end
		
		function xNew = getInput(obj,yi) % Get x-value from regression
			% -------------------------------------------------------------
			% Inputs:
			%	obj - Linregressor object
			%	yi  - Dependent var. index
			%
			% Outputs:
			%	xNew - Inverse regression output at yi
			% -------------------------------------------------------------
			
			% x = (y - b)/m
			xNew = (yi - obj.b)/obj.m;
		end
	end
end
