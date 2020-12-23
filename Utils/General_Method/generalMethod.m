function [val,err] = generalMethod(f,inVals,inErrs)
% function [val,err] = generalMethod(f,inVals,inErrs)
% -------------------------------------------------------------------------
% Copyright (c) 2020 Carter Mak
% https://github.com/cartermak/MATLAB-General-Method
% -------------------------------------------------------------------------
% Description:
% Pass an anonymous function handle and cell arrays of input values with
% errors. Returns output value with error estimated using the general
% method based on partial derivatives.
% -------------------------------------------------------------------------
% Inputs:
%	f:      function handle
%	inVals: cell array of input values to f, in order
%	inErrs: Cell array of uncertainties associated with each input value,
%	in order. If no uncertainty, use 0.
% 
% Outputs:
%	val: Vector of output values
%	err: Vector of uncertainties in output values
% -------------------------------------------------------------------------

%% Initialization
% Convert f to a symbolic function
g = sym(f);

% Get list of input variables
symVars = symvar(g);

% Number of input variables
nVars = length(symVars);

% Check for correct number of input variables and errors
if (nVars ~= length(inVals) || nVars ~= length(inErrs))
	error('Incorrect number of input values and/or errors.')
end

% Number of samples is length of longest input variable list
nSamples = 0;
for i = 1:length(inVals)
	if length(inVals{i})>nSamples
		nSamples = length(inVals{i});
	end
end

% Check that all variables and errors have the same length. Inputs should
% either have length Nsamples or 1 (implies same value for all samples)
for i = 1:length(inVals)
	
	% Check all variables. Extend if necessary.
	if length(inVals{i})==1
		inVals{i} = repmat(inVals{i},nSamples);
	end
	
	% Check all error values. Extend if necessary.
	if length(inErrs{i})==1
		inErrs{i} = repmat(inErrs{i},nSamples);
	end
	
	% If the length of any inputs or error values isn't gucci, throw error
	if nSamples~=length(inVals{i}) || nSamples~=length(inErrs{i})
		error('Inconsistent number of samples')
	end
	
end

% Preallocate outputs
val = zeros(nSamples,1);
err = val;

%% Calculus

% Preallocate (ish. They're cell arrays. It doesn't really matter.)
derivs = cell(nVars,1);
derivInputs = derivs;

% Calculate all partials
for i = 1:nVars
	currVar = symVars(i);
	currDiff = diff(g,currVar);
	derivs{i} = matlabFunction(currDiff);
	[~,derivInputs{i}] = intersect(symVars,symvar(currDiff),'stable');
end

%% Calculations

% Preallocate array of intermediate error terms
intermedErrs = zeros(nVars,1);
sampleVals = zeros(nVars,1);

% Loop over samples
for i = 1:nSamples
	
	% Store input variables to a temporary vector, sampleVals
	for j = 1:nVars
		sampleVals(j) = inVals{j}(i);
	end
	
	% Calculate output value for sample i
	tmpCell = num2cell(sampleVals);
	val(i) = f(tmpCell{:});
	
	% Loop over input variables for error calculation
	for j = 1:nVars
		currInputs = sampleVals(derivInputs{j});
		tmpCell = num2cell(currInputs);
		currFunc = derivs{j};
		intermedErrs(j) = currFunc(tmpCell{:})*inErrs{j}(i);
	end
	err(i) = sqrt(sum(intermedErrs.^2,'omitnan'));
end
