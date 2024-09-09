## Check's the LAST frame of a simulation, is there a network connecting all the particles?
## Extract data to CSV and run primary network analysis on the results 
## of a DPD bimodal colloid simulation
## NOTE: requires matching Fortran module
## NOTE: this code assumes 1 solvent type (typeid=0) and 
##		2 colloid types (typeid=1, typeid=2)
##
## - Extracts the list of all bonded colloid pairs (i.e network edges) from
##     GSD to CSV files for the LAST frame
## - Calculates primary network analysis metrics for the LAST frames
## 	- number of connected components 
##	- average degree (AKA average coordination number in gels)
##      - largest connected component (LCC)
##      - average clustering coefficient
##      - average square clustering coefficient 
##
## (Rob Campbell)

######### MODULE LIBRARY
import numpy as np
import gsd.hoomd
import math
import module
import pandas as pd
import networkx as nx
import os
from statistics import mean


######### INPUT PARAMETERS
filename = '../Gelation-angularrepulsion-BD.gsd'

edge_output = 'frame-edges/edgelist' # + <#>.csv in f90

colloid_type_id = 0

# create "frame-edges" subfolder if it doesn't exit
if os.path.exists('frame-edges') == False:
  os.mkdir('frame-edges')

######### DEFINE SIM CHECKS

#######
## edgelistCSV
# create a CSV file of all the bonded particle pairs (i.e. edges of the network) for the LAST frame
# i,j position for each bond/edge
def edgelistCSV_last_calc(filename, edge_output):
  traj = gsd.hoomd.open(filename, 'r')
  nframes = 1
  lastframe = len(traj)-1
  colloids = np.where(traj[-1].particles.typeid == [colloid_type_id])[0]
  ncolloids = len(colloids)
  radii = 0.5*traj[-1].particles.diameter[colloids]
  rcut = 0.1
  lbox = traj[-1].configuration.box[:3]
  inv_lbox = 1/lbox

  # get the xyz position of all colloids in the last frame    
  allpos = np.zeros((1,ncolloids,3))
  allpos[0,:,:] = traj[-1].particles.position[colloids] 
 
  module.edgelist_calc(nframes,ncolloids,radii,allpos,lbox,rcut,edge_output)

  os.rename(edge_output+'0.csv', edge_output+str(lastframe)+'.csv')

  print("Edgelist calculation complete for last frame of the simulation")
#######


#######
## networkx analysis
# use the networkx package to calculate for the LAST frame:
#   - the number of connected components
#   - average degree (AKA average coordination number)
#   - largest connected component (LCC)
#   - average clustering coefficient
#   - average square clustering coefficient
def primary_networkx_last_calc(filename, edge_output):
  # get the number of particles
  traj = gsd.hoomd.open(filename, 'r')
  nframes = 1
  lastframe = len(traj)-1
  colloids = np.where(traj[-1].particles.typeid == [colloid_type_id])[0]
  ncolloids = len(colloids)

  # create empty dataframe for all network data  
  network_frames_df = pd.DataFrame()
  
  '''
  # get length of sim from CSV data
  dir_path = r'frame-edges'
  # NOTE: this counts ALL CSV files in this directory
  nframes = 0
  # Iterate directory
  for filename in os.listdir(dir_path):
    # Check if the file has a CSV extension
    if filename.endswith('.csv'):
      nframes += 1
  '''

  # import data from the last frame
  filepath = edge_output+str(lastframe)+'.csv'
  # import CSV data
  lastframe_df = pd.read_csv(filepath)
  # rename colums as needed
  lastframe_df = lastframe_df.rename(columns={"i": "source", "j": "target"})
  lastframe_df.insert(loc=0, column='frame', value=lastframe)
 
  # network analysis

  df = lastframe_df[["source", "target"]]

  # create the network from edge list
  g = nx.from_pandas_edgelist(df)

  # if a node is not in the network, add it
  for particle in range(ncolloids):
    if (  not(   g.has_node(particle)   )  ):
      g.add_node(particle)

  # calculate the total number of edges      
  nedges = nx.number_of_edges(g)

  # calculate the average degree of the network, i.e. avg contact number in gels
  avg_degree = 2 * nedges / ncolloids

  # number of connected components
  n_cc = nx.number_connected_components(g)

  # return the indices of the nodes of the largest connected components
  lcc_nodes = max(nx.connected_components(g), key=len)

  # calculate the size of the largest connected component
  lcc_size = len(lcc_nodes)

  # calculate the avg clustering coefficient of the network
  avg_clustering_coeff = nx.average_clustering(g)

  # calculate the avg square clustering coefficient of the network
  square_clustering_per_node = nx.square_clustering(g)
  avg_sq_clustering_coeff = mean(square_clustering_per_node[jj] for jj in square_clustering_per_node)

  # compile outputs
  data={'frame'                  :[lastframe],
        'n_components'           :[n_cc], 
        'lcc_size'               :[lcc_size],
        'avg_degree'             :[avg_degree],
        'avg_clustering_coeff'   :[avg_clustering_coeff],
        'avg_sq_clustering_coeff':[avg_sq_clustering_coeff]}

  res_df = pd.DataFrame(data)


  # write all data to one CSV file
  res_df.to_csv('gelcheck-lastframe.csv',index = False)
  print("Primary network analysis complete for last frame")
#######


######### RUN CHECKS ON A SIMULATION

if __name__ == '__main__':
  edgelistCSV_last_calc(filename, edge_output)	
  primary_networkx_last_calc(filename, edge_output)	
