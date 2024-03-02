%  --------------- BIPOLAR MATING TENDENCY --------------------------------

% How to cite the article?

% Gençal, M.C., Oral, M. "Bipolar Mating Tendency: Harmony Between the Best and the Worst Individuals",
% Arabian Journal for Science and Engineering, vol. 47, no. 2, Springer Science and Business Media LLC,
% Sept. 2021, pp. 1849–71, doi:10.1007/s13369-021-06105-5.

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

% The following algorithm is the Genetic Algorithm (GA) that uses Bipolar
% Mating Tendency (BMT) as a selection method. So, the algorithm is called
% "GA-BMT".

function [mincost,value,bests]=GA_BMT(popsize, maxIter,t_size,bipolarity_value,pc,pm,elitism_rate,GC,f_name,seed)

rng(seed);
% As an example, you can call the algorithm like;
% [mincost,value,bests]=GA_BMT(100,100,4,0.25,0.7,0.05,0.05,[5 5;-5 -5],'dejong',1)

n=size(GC,2);

%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%  
pop=ones(popsize,n);

% By using uls and lls in GC, we can create the initial population
for i=1:n
    pop(:,i)=unifrnd(GC(2*i),GC(2*i-1),popsize,1);
end

% Calculate fitness values of each individual
cost=feval(f_name,pop(:,1),pop(:,2));
% Depends on n, you can write this equality as follows;
% cost=feval(f_name,pop(:,1),pop(:,2),...,pop(:,n) )

% Create the matrix "bests" which is one of our outputs
bests=zeros(maxIter,1);

% For plotting
% iter=1:maxIter;

count=maxIter;
while maxIter>0
    
    %%%%%%%% Selection %%%%%%%%
    
    [first_mate,second_mate]=bmt(pop,t_size,cost,bipolarity_value);
    
    %%%%%%%% CrossOver %%%%%%%%
    
    % rpc decides what individual is going to mate
    rpc=rand(popsize,1);
    
    % Individuals not attend the crossover
    ind_co=find(rpc>pc);
    
    % Compute offspring
    newpop=first_mate./2+second_mate./2+(1/2+0.1).*abs(second_mate - first_mate).*(2*rand-1);
    newpop(ind_co,:)=first_mate(ind_co,:);

    %%%%%%%% Mutation %%%%%%%%
    
    % Mutate the individual, if applicable, by creating new point for it.
    mut_pop=ones(popsize,n);
    for i=1:n
        mut_pop(:,i)=unifrnd(GC(2*i),GC(2*i-1),popsize,1);
    end
    
    % rpm decides what individual is going to mutate
    rpm=rand(popsize,1);
    mut_value=rpm<=pm;
    
    % mgpr,number of mutated gene probability, value to decide how many genes are going to mutate
    mgpr=0.95;
    mpm=rand(popsize,1);
    indices=mpm<=mgpr;
    
    % The matrix temp is used in the case of mutating all genes
    temp=ones(popsize,1);
    temp(indices)=0;
    
    % mg is a matrix that helping us to choose which gene we're going to
    % mutate in the case of choosing 1 gene.
    mg=rand(popsize,n);
    for j=1:n
        if mod(j,2)==0
            mg(:,j)=round(mg(:,j));
        else
            mg(:,j)=~round(mg(:,j));
        end
    end
    
    % After the mutation procedure, mut_matrix assists us to create the new
    % population.
    mut_matrix=zeros(popsize,n);
    for i=1:popsize
        if mut_value(i)==0
            mut_matrix(i,:)=zeros(1,n);
        else
            if temp(i)==1
                mut_matrix(i,:)=ones(1,n);
            else
                mut_matrix(i,:)=mg(i,:);
            end
        end
    end
    newpop=newpop.*~mut_matrix+mut_pop.*mut_matrix;
    
    maxIter=maxIter-1;
    
    %%%%%%%% Elitism %%%%%%%%
    
    [~,in1]=sort(cost);
    n_elit=round(size(in1,1)*elitism_rate);
    if n_elit<1
        n_elit=1;
    end
    part1=pop(in1(1:n_elit),:);
    cost_new=feval(f_name,newpop(:,1),newpop(:,2));
    [~,in2]=sort(cost_new);
    part2=newpop(in2(1:popsize-n_elit),:);
    pop=[part1;part2];
    
    % Finding the best for each iteration
    index=count-maxIter;
    cost=feval(f_name,pop(:,1),pop(:,2));
    [val,~]=min(cost);
    bests(index)=val;
end
last_cost=feval(f_name,pop(:,1),pop(:,2));
[mincost,in]=min(last_cost);
value=pop(in,:);

% For plotting
% figure, semilogy(iter,bests,'r');
% xlabel('Iterations');
% ylabel('Best Value');


%%%%%%%%% Bipolar Mating Tendency %%%%%%%%%%

function [first_mate,second_mate]=bmt(pop,tournamentsize,cost,equality_value)
attend1=zeros(tournamentsize,1);
attend2=zeros(tournamentsize,1);
popsize=size(pop,1);
first_mate=ones(popsize,2);
second_mate=ones(popsize,2);
rand_value=rand(popsize,1);
for j=1:popsize
    for k=1 : tournamentsize
        attend1(k)=ceil((popsize-1)*rand+1);
        attend2(k)=ceil((popsize-1)*rand+1);
    end
    winner1=find(cost==min(cost(attend1)));
    winner2_min=find(cost==min(cost(attend2)));
    winner2_max=find(cost==max(cost(attend2)));
    % If we have, more than 1 winner
    first_mate(j,:)=pop(winner1(1),:);
    if rand_value(j)<=equality_value
        second_mate(j,:)=pop(winner2_min(1),:);
    else
        second_mate(j,:)=pop(winner2_max(1),:);
    end
end
