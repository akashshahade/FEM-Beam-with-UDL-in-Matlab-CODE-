%-------------------------------------------------------------------%
% Software by -  Akash S. Shahade
% Title - Program for FE Analysis of Beam.
%-------------------------------------------------------------------%
clc;clear;
fprintf(' \nFEM ANALYSIS OF BEAM. \n\n');

%------------------------------------------------------%
% STEP 01 - PRE-PROCESSING
%------------------------------------------------------%

P = -12*10^3;                        % 12kN UDL acting in -ve Y direction
E=2*10^11;                           % Young's Modulus
moi = 4*10^(-6);                     % Moment of Inertia
L=1;                                 % Length
ei=E*moi;

ele_nod=[1 2;2 3];                   % Elements are connected with these nodes
nod_coor=[0 0;1 0;2 0];              % Node coordinates
num_ele=size(ele_nod,1);             % Number of elements
ele_dof=[1 2 3 4;3 4 5 6];           % D.O.F  associated with Nodes
num_nod=3;                           % Number of Nodes
dof = 2;                             % D.O.F per node

displacement = zeros(dof*num_nod,1);    % Zero Matrix for Displacement
force = zeros(dof*num_nod,1);           % Zero Matrix for Force
stiffness = zeros(dof*num_nod);         % Zero Matrix for Stiffness

%------------------------------------------------------%
% Stiffness matrix calculation & ASSEMBLY
%------------------------------------------------------%

for e=1:num_ele                               % For 1 to Number of elements
    
    k = ((ei)/(L^3))*[12 6*L -12 6*L;...      % Stiffness matrix Calculation
        6*L 4*(L^2) -6*L 2*(L^2);...
        -12 -6*L 12 -6*L;...
        6*L 2*(L^2) -6*L 4*(L^2)];

   
% extract the rows of ele_dof (for each element e)
ele_dof_vec=ele_dof(e,:);

    for i=1:4
        for j=1:4
                                              % Assembly of Global Matrix
  stiffness(ele_dof_vec(1,i),ele_dof_vec(1,j))=...
  stiffness(ele_dof_vec(1,i),ele_dof_vec(1,j))+k(i,j);

        end
    end
end


f2 = [(P*L)/2; (P*(L^2))/12; (P*L)/2; -(P*(L^2))/12];
force = [0;0;f2];

fprintf('Global Stiffness Matrix: \n');
disp(stiffness);
fprintf('\n Global Load Vector: \n');
disp(force);
fprintf('\n------------------------------------------------\n');

%------------------------------------------------------%
% Boundary Conditions
%------------------------------------------------------%

fixed_dof = [1 2 3 5];               % Fixed D.O.F.
k=stiffness;
k(fixed_dof,:)=[];                   % Eliminating Rows
k(:,fixed_dof)=[];                   % Eliminating Columns

f=force;
f(fixed_dof,:)=[];                   % Eliminating Rows

%------------------------------------------------------%
% STEP 02 - SOLVE
%------------------------------------------------------%

q = k\f ;

%------------------------------------------------------%
% STEP 03 - POST-PROCESSING
%------------------------------------------------------%

displacement=[0;0;0;q(1);0;q(2)];               % Displacement Vector

reaction = stiffness*displacement - force;      % Calculate Reaction Forces

Node = [1;2;3;4;5;6];
qxy = {'w1';'Theta1';'w2';'Theta2';'w3';'Theta3'};
Q = displacement;
T=table(Node,qxy,Q);
disp(T);

fprintf('\n------------------------------------------------\n');

F = {'R1x';'R1y';'R2x';'R2y';'R3x';'R3y'};
Reaction_N = reaction;
J = table(F,Reaction_N);
disp(J);
fprintf('END OF PROGRAM.\n');
