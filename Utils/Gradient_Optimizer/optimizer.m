function x = optimizer(g_fun,x0)
% function x = optimizer(g_fun,x0)
% -------------------------------------------------------------------------
% Copyright (c) 2020 Carter Mak
% https://github.com/cartermak/MATLAB-Utils
% -------------------------------------------------------------------------
% MIT License
% 
% Copyright (c) 2020 Carter Mak
% 
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the
% following conditions:
% 
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.
% -------------------------------------------------------------------------
% Description: 
%	Gradient/Curvature-based numerical optimization. Stop condition is
%	local gradient of less than 0.1% of the gradient at initial condition
%	in all function arguments.
% -------------------------------------------------------------------------
% Inputs:
%	g_fun - function to be optimized
%	x0    - Initial condition vector
%   
% Ouputs: 
%   x     - State vector at local extremum

% Assumptions: 1% of inital value will be a sufficent step size
% -------------------------------------------------------------------------

%% Initialization
x0 = x0(:); % Ensure x0 is a column vector
x = x0;

steps = 0.01*x0; % Use 1% of initial value as step size
steps(~steps) = 0.01; % If IC is zero, use 0.01 as step size

maxIts = 1e3; % Set max number of iterations
counter = 0;

currGrad = gradient(g_fun,x,steps);
	currCurv = curvature(g_fun,x,steps);
stopCondition = abs(0.001*currGrad); % Stop condition

%% Loop

while any(abs(currGrad)>abs(stopCondition))
	% Extrapolate extremum based on local curvatures
	d = currGrad./(currCurv.*sqrt(1 + currGrad.^(-2)));
	
	% Shift output
	x = x - d;
	
	% Get local gradient and curvature
	currGrad = gradient(g_fun,x,steps);
	currCurv = curvature(g_fun,x,steps);
	
	% Check max iteration condition
	counter = counter+1;
	if counter>maxIts
		warning('Max number of iterations exceeded.')
		break;
	end

end

end

function gradient = gradient(g_fun,x,steps)
% -------------------------------------------------------------------------
% Description: Returns discrete-valued gradient at x, using a radius on
% each variable as defined by steps.
% -------------------------------------------------------------------------


	dim = length(x);
	gradient = zeros(dim,1);
	
	for i = 1:dim
		% Temp vars
		delta = steps(i);
		xCurr = x;
		% Interrogate g_fun at nearby points
		xCurr(i) = xCurr(i) + 0.5*delta;
		val1 = g_fun(xCurr);
		xCurr(i) = xCurr(i) - delta;
		val2 = g_fun(xCurr);
		% Calculate slope, centered about x
		gradient(i) = (val1-val2)/delta;
	end
end

function curvature = curvature(g_fun,x,steps)
% -------------------------------------------------------------------------
% Description: Returns discrete-valued curvature at x, using a radius on
% each variable as defined by steps.
% -------------------------------------------------------------------------
	dim = length(x);
	curvature = zeros(dim,1);
	
	for i = 1:dim
		% Temp vars
		delta = steps(i);
		xCurr = x;
		% Interrogate g_fun above, below, and at x (in each dimension)
		val2 = g_fun(xCurr);
		xCurr(i) = xCurr(i) + 0.5*delta;
		val3 = g_fun(xCurr);
		xCurr(i) = xCurr(i) - delta;
		val1 = g_fun(xCurr);
		% Get slopes
		slope1 = 2*(val2-val1)/delta;
		slope2 = 2*(val3-val2)/delta;
		% Curvature
		curvature(i) = 2*(slope2-slope1)/delta;
	end
end