% The algorithm is the Genetic Algorithm (GA) that uses Bipolar Mating Tendency (BMT) as a selection method. So, the algorithm is called "GA-BMT".
% BMT takes ten parameters as inputs and returns three parameters as an outputs:

% --------------------------- Inputs --------------------------------------
% popsize                       : Size of the population
% maxIter                       : Number of generations
% t_size                        : Tournament Size
% bipolarity value              : Probability value that decides who is going to mate;
%                                 the best or the worse individual
% pc                            : Crossover Rate
% pm                            : Mutation Probability
% elitism_rate                  : What percentage of the best members of the population we need to keep
% GC(Gene_constraints)          : This parameter has 2*n parameters, where n is the #gene, which are;
%                                       ul : upper limit of the range 
%                                       ll : lower limit of the range
%                                       GC must be like "[ ul1 ul2
%                                                         ll1  ll2  ]" if every individual has 2 genes.
% f_name                        : The name of objective function which must
%                                 be entered like 'function name'.
% -------------------------------------------------------------------------

% ------------------------ Outputs ----------------------------------------
% mincost                       : The found minimum fitness value
% value                         : The individual having mincost value
% bests                         : All bests throughout iterations
% -------------------------------------------------------------------------

% As an example, you can call the algorithm like;
% [mincost,value,bests]=GA_BMT(100,100,4,0.25,0.7,0.05,0.05,[5 5;-5 -5],'dejong',1)
