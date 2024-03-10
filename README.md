# communities
A toolkit to support analysis of communities in networks, based on R+igraph.
Contents:

1- benchmark networks: 
  make_ring_of_cliques()
  make_LFR()
  
2- enrich_communities: 
 add_names_from_vids(), 
 add_size(), 
 sort_by_size(), 

3- Visualize communities: 
 expanded_circle_layout(), 
 make_comms_network()
 
4- Analyse communities: 
 mixing_parameter_partition(), 
 mixing_parameter_comms(), 
 plot_mixing_par_quality_check)(), 
 check_stability_repeated(), 
 check_input_ordering_bias(), 
 assess_uncertainty_by_pairwise_comparison)()

5- Compare communities: 
 comms_evolution_over_time(C1, C2) - returns: continue/split/merge/created/estinguished.. 


