clc
clear
%defining the objective function
function f = springObjective(x)
x1 = x(1);
x2 = x(2);
x3 = x(3);
f = x1^2*x2*(x3+2);
end
%Creating constraints function
function [c,ceq] =springconstraints(x)
x1 = x(1);
x2 = x(2);
x3 = x(3);
g1 = 1 - (x2^3*x3)/(71785*x1^4);
g2 = (4*x2^2 - x1*x2)/(12566*(x2*x1^3 -x1^4))+ 1/(5108*x1^2)- 1;
g3 = 1 - (140.45*x1)/(x2^2*x3);
g4 = (x1 + x2)/1.5 - 1;
c = [g1;g2;g3;g4];
ceq = [];
end
%defining limits
Ib = [0.05 0.25 2];
ub = [2 1.3 15];
%employing FMINCON Optimizer
%first defining initial guess
x0 = [0.5 0.5 5];
tic
opts = optimoptions('fmincon','Display','off','Algorithm','sqp');
[x,fval] = fmincon(@springObjective,x0,[],[],[],[],Ib,ub,[],opts);
disp('Best values of x:')
disp(x)
disp('Miimum objective value:')
disp(fval)
time_fmincon = toc;

%employing ga Optimizer
%first defining number of variables 
nvars = 3;
tic 
[x,fval] = ga(@springObjective,nvars,[],[],[],[],Ib,ub,@springconstraints);
disp('Best values of x:')
disp(x)
disp('Miimum objective value:')
disp(fval)
time_ga = toc;

%employing particleswarm
%using penalty technique since particleswarm doesnt accept nonlinear
%constraints directly
function f = springPenalty(x) 
f = springObjective(x);
[c,~] = springconstraints(x);
penalty = sum(max(c,0).^2);
f = f + 1e6*penalty;
end
nvars = 3;
tic 
[x,fval] = particleswarm(@springObjective,nvars,Ib,ub);
disp('Best values of x:')
disp(x)
disp('Miimum objective value:')
disp(fval)
time_particleswarm = toc;

%Employing patternsearch optimizer
x0 = [0.5 0.5 5];
tic 
[x,fval] = patternsearch(@springObjective,x0,[],[],[],[],Ib,ub,[],opts);
disp('Best values of x:')
disp(x)
disp('Miimum objective value:')
disp(fval)
time_patternsearch = toc;

%employing the surrogateopt optimizer
tic 
[x,fval] = surrogateopt(@springObjective,Ib,ub);
disp('Best values of x:')
disp(x)
disp('Miimum objective value:')
disp(fval)
time_surrogateopt = toc;

%employing simulannealbnd optimizer
x0 = [0.5 0.5 5];
tic 
[x,fval] = simulannealbnd(@springPenalty,x0,Ib,ub);
disp('Best values of x:')
disp(x)
disp('Miimum objective value:')
disp(fval)
time_simulannealbnd= toc;

%graph of the times for all optimizers
Optimizers = ['FMINCON','GA','PARTICLE SWARM','PATTERN SEARCH','SURROGATEOPT','SIMULANNAEALBND'];
TIMES = [time_fmincon time_ga time_particleswarm time_patternsearch time_surrogateopt time_simulannealbnd];
bar(TIMES)
xticks(1:numel(Optimizers))
xticklabels(Optimizers)
xtickangle(45)
set(gca,'XtickLabel',Optimizers)
ylabel('TIMES(seconds)')
title('OPTIMIZATION SOLVER TIME COMPARISON')
grid on
save
