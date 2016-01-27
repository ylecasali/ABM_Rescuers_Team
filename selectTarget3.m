function target = selectTarget3( nodes, excluded_nodes,Obj_s )

available_nodes = setdiff(nodes, excluded_nodes);


selection = available_nodes(randi(length(available_nodes))); %random selection
%of the nodes from the available nodes

[dist_t_3, path_t_3, pred_t_3] = shortestpath(Obj_s,excluded_nodes, available_nodes,'Method','BFS');
idx=find(dist_t_3==min(dist_t_3));
%idx_2=randi(idx)
idx_2=randi([1,length(idx)],1,1);
idx_3=idx(idx_2);
target = available_nodes(idx_3);

end

