clearvars;close all;clc;

% Define parameters
% These can also vary for different members
E = 500;
A = 0.5;

% Instantiate finite element model with 4 nodes
fem = FEM(4);

% Define node positions
fem = fem.addNode(1,0,0);
fem = fem.addNode(2,5,0);
fem = fem.addNode(3,10,0);
fem = fem.addNode(4,5+5/sqrt(2),5/sqrt(2));

% Define bars connecting nodes
fem = fem.addBar(1,E,A,1,2);
fem = fem.addBar(2,E,A,2,3);
fem = fem.addBar(3,E,A,2,4);

% Define kinematic constraints
fem = fem.constrain(1,2);
fem = fem.constrain(2,2);
fem = fem.constrain(3,1);
fem = fem.constrain(3,2);
fem = fem.constrain(4,1);
fem = fem.constrain(4,2);

% Add external load(s)
fem = fem.addForce(1,1,-200);

% Solve finite element problem and print results
fem = fem.solve();
fem.printResults();