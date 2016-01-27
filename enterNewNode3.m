function [next, time_to_node] = enterNewNode3(p, positions, next_move_BFS, ...
                                             target_nodes, perturbed_sp)
    
    % Which node should I go next according to shortest paths.
    % TODO: this could be decided according to the strategy
    % currently played by player.
    next = next_move_BFS(positions(p), target_nodes(p));

    % How much time do I need to next node?
    time_to_node = perturbed_sp(positions(p), next);
    
    
end