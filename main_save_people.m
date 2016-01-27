% There is a normal network, which gets perturbed by a calamity.
% Rescuers explores the new perturbed network trying to rescue people.
% Rescuers start from hospital (node 1) and must find and bring back 
% injured people before it is too late.
% Different strategies of exploration and rescuing are tested.
% Goal: find out if there is a superior strategy or a mix thereof.


clc
close all
clear



%% Create Networks


from_nodes = [6 1 2 2 3 4 4 5 5 6 1 7 2 8 9 10 11 11];
to_nodes =   [2 6 3 5 4 1 6 3 7 3 5 11 9 11 2 8 10 2];

% Use discrete weights, so we can consider them a single simulation steps.
Wt = [4 10 1 12 14 8 8 9 7 7 9 3 4 5 6 7 8 2]; %weights of time
Nt = sparse(from_nodes, to_nodes, Wt); %sparse matrix that is the network

%view(biograph(Nt,[],'ShowWeights','on')) %display the network with weights

% Perturb the network.
Ws = Wt * (1+rand()); 
Ns = sparse(from_nodes, to_nodes, Ws); %sparse matrix that is the network

%view(biograph(Ns,[],'ShowWeights','on')) %display the network with weights

% graphallshortestpaths(Ns) %evaluate total time spent from every node of the network to each other node

% Create biograph object (used to compute shortest paths).
Obj_s = biograph(Ns); 
Obj_t = biograph(Nt); 



% Nr. Rescuers
N = 10;

% Nr. Rounds = the time to live (ttl). Injured people need to be found,
% and brought to hospital in less than this time, otherwise they die.
% This is a simplication. All injured people have the same time to live.
T = 200;

% Nr. of rescued people in the simulation.
RESCUED = 0;
TOTAL_INJURED = 0; % will be init later.

% Multiplier for how many injured per node maximum.
INJ_MULTIPLIER = 20;

% Under normal conditions, the time to reach a node from another one.
normal_sp = graphallshortestpaths(Nt);

% Under a natural disaster the network was perturbed, and rescuers do
% not know for sure how much it takes to move between nodes.
% The values of this matrix have to be learnt by the agents (rescuers).
perturbed_sp = graphallshortestpaths(Ns);

% Each rescuer build an own image of the network after perturbation.
% It starts as a copy of the normal network, equal for all rescuers.
% It might get updated during the exploration or not, depending on
% strategy.
networks = repmat(normal_sp, N);

% Number of nodes in the network.
nNodes = length(Nt);

% Node indexes.
nodes = 1:nNodes;

z_max=10;

R_matrix=nan(z_max,1);
D_matrix=nan(z_max,1);
strategy_1_matrix=nan(z_max,1);
strategy_2_matrix=nan(z_max,1);
strategy_3_matrix=nan(z_max,1);




for z=1:z_max
% Place injured people randomly on the network.
injured = (randn(nNodes,1) < 0.5) .* randi(INJ_MULTIPLIER, nNodes, 1);

% Hypotesis: there are always injured in the network. 
% if (injured == zeros(size(injured)))
% id=randi(nNodes,1,1);
% injured(id)=randi(INJ_MULTIPLIER,nNodes,1);
%   if id==1
% injured(2)=7; % set chosen value 
%   end
% end

num_injured=injured;

% Place Hospital in the network (node 1);
injured(1) = 0;

TOTAL_INJURED = sum(injured);

% Current position in the network for the rescuers.
% (At the beginning all in position 1 - hospital).
positions = ones(N,1);

% Reformat network data for use in the main simulation loop.
% If a rescuer is in node X and wants to reach node Y, which node should
% he / she go to next, according to the shortest path algorithm?
next_move = zeros(nNodes);
for x = 1 : nNodes
    for y = 1 : nNodes
        if (x ~= y) 
            [dist_t, path_t, pred_t] = shortestpath(Obj_t, x,y);
            next_move(x,y) = path_t(2);
        end
    end
end

next_move_perturbed = zeros(nNodes);
for x = 1 : nNodes
    for y = 1 : nNodes
        if (x ~= y) 
            [dist_t, path_t, pred_t] = shortestpath(Obj_s, x,y);
            next_move_perturbed(x,y) = path_t(2);
        end
    end
end

next_move_BFS = zeros(nNodes);
for x = 1 : nNodes
    for y = 1 : nNodes
        if (x ~= y) 
            [dist_t, path_t, pred_t] = shortestpath(Obj_s, x,y,'Method','BFS');
            next_move_BFS(x,y) = path_t(2);
        end
    end
end



% Where a rescuer is heading to. 0 means no target.
target_nodes = zeros(N, 1);

% Whether a rescuer is carrying an injured.
carrying_injured = zeros(N, 1);

% Current waiting time in node.
% To move from one node to another, rescuers need to wait the value
% of the wait in the connection between two nodes.
waiting_in_node = zeros(N,1);

% Available strategies.
avail_strategies = [ 1 2 3];

% Nr. available strategies each player each round.
nr_strategies = length(avail_strategies);

% Currently chosen strategy by each player.
strategies = ones(N,1);
strategy_1=0;
strategy_2=0;
strategy_3=0;
strategy_salvation=zeros(N,nr_strategies);
salvation_1=zeros(N,1);
salvation_2=zeros(N,1);

% Initial propensities for each strategy of each player.
propensities = ones(N, nr_strategies);

% Probabilities to choose a strategy. (all equally probable at t=0).
probabilities = ones(N, nr_strategies)*0.25;

% Records payoffs and routes for each player at each iteration.
payoff_matrix = zeros(N,T);
route_choice = zeros(N,T);

% Number of simulation
n_simulation=1000;
C=zeros(n_simulation,2);
c1=0;
c2=0;
D=zeros(n_simulation,1);
R=zeros(n_simulation,1);


DEAD=0;
RESCUED=0;

%% Model the strategy choices of each player 
% For each round, each player decides which strategy to choose
% depending on the previous payoff and propensity. 
% Probability is given by x(i,j)/sum(x(i,j))

    

%% Routh choice and interaction


for t = 1 : T
        
        
        for player = 1 : N
            
            % Reset variables for each player.
            injured_found = 0;
            injured_rescued = 0;
            
            %if (num_injured ~= zeros(size(num_injured))) % if there are still injured to rescue
            
                %if TOTAL_INJURED>RESCUED                    % NO DEAD < 0
                if (sum(injured)>0)
           
            if (target_nodes(player) == 0)        
            % Randomly draw a strategy according to the prob distribution.
            curStrategy = randsample(avail_strategies, ...
                1, true, probabilities(player,:));
            
            strategies(player) = curStrategy;
            end
            
            %Find a target node

            % STRATEGY 1
            if curStrategy==1
                
               if (target_nodes(player) == 0)
                % Notice: +1 is to exclude the hospital as a target.
                target_nodes(player) = selectTarget(nodes, ...
                                                    [positions(player)]);
                                                
                [next_node, time_to_node] = enterNewNode(player, ...
                                                         positions, ...
                                                         next_move, ...
                                                         target_nodes, ...             
                                                         perturbed_sp);
                                                     
                                                     
                % Set waiting time in current node.
                waiting_in_node(player) = time_to_node;
                                                     
               end
            end
            
            % STRATEGY 2
            if curStrategy==2
                
             if (target_nodes(player) == 0)
                % Notice: +1 is to exclude the hospital as a target.
             
              
                target_nodes(player) = selectTarget2(find(injured), ...
                                                    [positions(player)],Obj_s);
                                                
                [next_node, time_to_node] = enterNewNode(player, ...
                                                         positions, ...
                                                         next_move, ...
                                                         target_nodes, ...             
                                                         perturbed_sp);
                                                     
                                                     
                % Set waiting time in current node.
                waiting_in_node(player) = time_to_node;
                                                     
              end
            end
            
             % STRATEGY 3
            if curStrategy==3
                
               if (target_nodes(player) == 0)
                % Notice: +1 is to exclude the hospital as a target.
                target_nodes(player) = selectTarget3(nodes, ...
                                                    [positions(player)],Obj_s);
                                                
                [next_node, time_to_node] = enterNewNode3(player, ...
                                                         positions, ...
                                                         next_move_BFS, ...
                                                         target_nodes, ...             
                                                         perturbed_sp);
                                                     
                                                     
                % Set waiting time in current node.
                waiting_in_node(player) = time_to_node;
                                                     
               end
            end
            
 
            % Entering a new node.
            if (waiting_in_node(player) <= 0)
                
            
            % Evaluate node.
            % There are still injureds in the node evaluated  
            
            [next_node, time_to_node] = enterNewNode(player, ...
                                                         positions, ...
                                                         next_move, ...
                                                         target_nodes, ...             
                                                         perturbed_sp);

                % Move player in new position.
                positions(player) = next_node;
            
                
                % Player is in hospital
                if (next_node == 1)
                    
                    % Transit or rescued successful?
                    if (carrying_injured(player) == 1)
                        injured_rescued = 1;
                        carrying_injured(player) = 0;
                        RESCUED = RESCUED + 1;

                    % Count strategy when injured are saved
                    if curStrategy==1
                        strategy_1=strategy_1+1;
                    elseif curStrategy==2
                        strategy_2=strategy_2+1;
                    else
                        strategy_3=strategy_3+1;
                    end
                    
                         % Display injured list on the command window
                         if injured(next_node) == 0                            
                            injured;
                         end
                        
                        
                        % Reset target node and time_to_node.
                        % They will be initialized at next iteration.
                        target_nodes(player) = 0;
                        time_to_node = 0;
                    end
                 end
                
                % Evaluate node.
                % Limits to 1 person rescued at the same time.
                if (injured(next_node) > 0 && ...
                        carrying_injured(player) == 0)
                    
                    % If an injured is found, remove it,
                    % change target_node to hospital,
                    % and update other variables.
                    injured(next_node) = injured(next_node) - 1;
                    carrying_injured(player) = 1;
                    target_nodes(player) = 1;
                    injured_found = true;
                    [next_node, time_to_node] = enterNewNode(player, ...
                                                         positions, ...
                                                         next_move, ...
                                                         target_nodes, ...             
                                                         perturbed_sp);
                else
                    if (next_node == target_nodes(player))
                        % We reached target node, but did not find an
                        % injured person. Reset!
                        target_nodes(player) = 0;
                        time_to_node = 0;
                    end
                end
                
                % Set waiting time in current node.
                waiting_in_node(player) = time_to_node;
            else
                % Reduce waiting time.
                waiting_in_node(player) = waiting_in_node(player) - 1;    
            end
            
            end  
            

                
            % PAYOFF PART 
            % Compute payoff each player.
            % TODO: Need to decide how to do it exactly. 
            % Might be strategy dependent or not.
            payoff = 0;
            if (injured_found)
                payoff = payoff + 50;
            elseif (injured_rescued)
                payoff = payoff + 100;
            else
                payoff = payoff - 2;
            end
            
            % Updates propensities.
            if payoff >= 0
                propensities(player,curStrategy) = ...
                propensities(player,curStrategy) + payoff;  
            else
                not_choosen_strategies = ...
                setdiff(avail_strategies, curStrategy);
                propensities(player,not_choosen_strategies) = ...
                propensities(player,not_choosen_strategies) - payoff;
            end
            
           
        end
             
        end
   
        
%end

% Display values: RESCUED, DEAD, strategy_1, strategy_2
RESCUED;
DEAD = TOTAL_INJURED - RESCUED;
strategy_1;
strategy_2;
strategy_3;
strategy_salvation=[salvation_1, salvation_2];



R_matrix(z)=RESCUED;
D_matrix(z)=DEAD;
strategy_1_matrix(z)=strategy_1;
strategy_2_matrix(z)=strategy_2;
strategy_3_matrix(z)=strategy_3;

end
% Plot: to display how many dead and rescued people, and to display how many times a strategy
% (1 or 2) was chosen in the code

STRATEGY_matrix=[strategy_1_matrix, strategy_2_matrix,strategy_3_matrix];
s1=std(strategy_1_matrix);
s2=std(strategy_2_matrix);
s3=std(strategy_3_matrix);

mean(strategy_1_matrix);
mean(strategy_2_matrix);
mean(strategy_3_matrix);

x=[1:z];
m1=ones(size(x))*mean(strategy_1_matrix);
m2=ones(size(x))*mean(strategy_2_matrix);
m3=ones(size(x))*mean(strategy_3_matrix);

figure
x=[1:z];
plot(x,strategy_1_matrix,'r',x,strategy_2_matrix+30,'b',x,strategy_3_matrix+60,'g')
legend('strategy1', 'strategy2', 'strategy3')
ylabel('Number of strategy that rescued people')
xlabel('Times')
title('Simulation and result for strategies')

figure
x=[1:z];
plot(x,strategy_1_matrix,'r',x,strategy_2_matrix+30,'b',x,strategy_3_matrix+60,'g',x,m3+60,'k',x,m1,'k',x,m2+30,'k')
legend('strategy1','strategy2', 'strategy3')
ylabel('Number of strategy that rescued people')
xlabel('Times')
title('Simulation and result for strategies')

figure
plot(x,R_matrix,'r',x,D_matrix,'b')
ylabel('Number of people')
xlabel('Times')
title('Simulation and result for rescued and dead')
legend('rescued people', 'dead people')

figure
hist(strategy_1_matrix)
h= findobj(gca,'Type','patch');
h.FaceColor=[0 0.5 0.5];
h.EdgeColor = 'w';


figure
hist(STRATEGY_matrix)
legend('strategy1','strategy2', 'strategy3')
title('Numbers of rescued people in each strategy')
ylabel('Number of times')
xlabel('Number of rescued people')

figure
subplot(2,1,1)
plot(1,RESCUED,'*g',2,DEAD,'*r')
legend('rescued','dead')
subplot(2,1,2)
plot(1,strategy_1,'*g',2,strategy_2,'*r', 3, strategy_3, '*b')
legend('strategy 1', 'strategy 2','strategy 3')



