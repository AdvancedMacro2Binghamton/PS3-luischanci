%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Binghamton University                         %
% PhD in Economics                              %
% ECON634 Advanced Macroeconomics               %
% Fall 2017                                     %
% Luis Chancí (lchanci1@binghamton.edu)         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

% PARAMETERS
beta  = 0.9932; sigma = 1.5; b     = 0.5;    
y_s   = [1, b]; PI    = [.97 .03; .5 .5]; 

% ASSET VECTOR
num_a = 1000;  a = linspace(-2, 5, num_a); 

% INITIAL GUESS FOR q
q_min = 0.98; q_max = 1.1;

aggsav = 1;
while (abs(aggsav) >= 0.01)
    q_guess = (q_min + q_max) / 2;
    c  = bsxfun(@plus, bsxfun(@minus,a',q_guess*a),permute(y_s,[1 3 2]));
    u  = (c.^(1-sigma))./(1-sigma); u(c<0)=-Inf;
    v0 = zeros(2, num_a);
    e1 = 1;
    while e1 >1e-06
          v = u + beta * repmat(permute((PI*v0),[3 2 1]),[num_a 1 1]);
          [vfn,idx] = max(v,[],2);        
          e1 = max(max(abs(permute(vfn,[3 1 2]) - v0)));
          v0 = permute(vfn,[3 1 2]);
    end
    pol_indx = permute(idx,[3 1 2]);
    g  = a(pol_indx);
    Mu = ones(size(g))/numel(g);
    e2 = 1; 
    while e2>=1e-06;
        [emp_ind, a_ind, mass] = find(Mu);
        MuNew  = zeros(size(Mu));
        for ii = 1:length(emp_ind)
            apr_ind = pol_indx(emp_ind(ii), a_ind(ii));   
            MuNew (:, apr_ind) = MuNew (:, apr_ind) + ...
                  (PI(emp_ind(ii), :) * mass(ii))';
        end
        e2 = max(max(abs(MuNew  - Mu)));
        Mu = MuNew;
    end
    aggsav = sum(sum(Mu.*g));
    if aggsav>0
       q_min = q_guess; else; q_max = q_guess; 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lorenz Curve and Gini Coeff.
f  = reshape(Mu',[numel(Mu) 1]);
y  = [reshape(bsxfun(@plus,repmat(a, [2,1]),y_s')',[2*num_a 1]),...% Wealth
	  reshape(repmat(y_s',[1 num_a])',[2*num_a 1])];                % Income
L    = []; G = [];
for i=1:2; 
    s = cumsum(sortrows([f,f.*y(:,i),y(:,i)],3));
    L = [L,bsxfun(@rdivide,s,s(end,:))*100];
    G = [G,1 - sum((s(1:end-1,2)+s(2:end,2)).*diff(s(:,1)))];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Welfare
c_bar  = 1*0.9434+0.5*.0566;                       % Using invariant Pr.
W_FB   = ((c_bar^(1-sigma))/(1-sigma))*inv(1-beta);% Welfare First Best
lambda = ((v0.^(-1)).*W_FB).^(1/(1-sigma))-1;      % Consumption equivalent
Frac   = sum(sum((lambda>0).*Mu));                 % Fraction W_FB > v(s,a)
WG     = sum(sum(lambda.*Mu));
%%%%%%%%%%%%%%    Now we can Plot  %%%%%%%%%%%%%%