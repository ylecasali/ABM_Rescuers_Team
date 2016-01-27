function target = selectTarget2( nodes, excluded_nodes,Obj_s)

   available_nodes = setdiff(nodes, excluded_nodes);
   %target = available_nodes(randi(length(available_nodes)));
   [dist_t_2, path_t_2, pred_t_2] = shortestpath(Obj_s,excluded_nodes, available_nodes);
   idx=find(dist_t_2==min(dist_t_2));
   target = available_nodes(idx);
   
end
