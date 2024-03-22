## Extract data to CSV and run primary network analysis on the results 
## of a BD colloid simulation
## NOTE: requires matching Fortran module
## NOTE: this code assumes 1 colloid type (typeid=1)
##
"""
## This code performs the following analyses:
## - Extracts colloid position data to from GSD to CSV files for all frames
## - Extracts the list of all bonded colloid pairs (i.e network edges) from
##     GSD to CSV files for all frames
## - Calculates primary network analysis metrics for all frames
## 	- number of connected components 
##	- average degree (AKA average coordination number in gels)
##      - largest connected component (LCC)
##      - average clustering coefficient
##      - average square clustering coefficient 
"""
## (Rob Campbell)


########################
""" MODULE LIBRARY """
########################
import numpy as np
import gsd.hoomd
import math
import fortranmod as module 
import pandas as pd
import networkx as nx
import os
import re
from statistics import mean


##########################
""" INPUT PARAMETERS """
##########################
# path to sim data
filename = '../Gelation-BD.gsd'

# filepath to folder where data files will be created
data_outpath = 'data'

###########################
""" DEFINE SIM CHECKS """
###########################

# create "data" subfolder if it doesn't exit
if os.path.exists(data_outpath) == False:
  os.mkdir(data_outpath)

#######
## posCSV
"""
# create a CSV of particle (i.e. node) position information for each frame
# (tag, x, y, z, typeID, radius)
"""
def posCSV_calc(filename):
  # open the simulation GSD file as "traj" (trajectory)
  traj = gsd.hoomd.open(filename, 'r')

  # set path and filename
  dir_path = data_outpath+'/frame-pos'
  pos_output = dir_path+'/positions_frame' # + <#>.csv in python loop

  # check for existing CSV data 
  if os.path.exists(dir_path) == False:
    os.mkdir(dir_path)
  if os.path.exists(dir_path) == True:
    # NOTE: this counts ALL CSV files in this directory
    nframes = 0
    # set the pattern for files ending in <number>.csv
    pattern = re.compile(r'positions_frame\d+\.csv$')
    # Iterate directory
    for filename in os.listdir(dir_path):
      # Check if the file has a CSV extension
      if pattern.match(filename):
        nframes += 1

    if nframes == len(traj):
      print('position data CSV files already seem to exist for all frames. Not creating new CSV files.')
      return

  nframes = len(traj)
  colloids = np.where(traj[-1].particles.typeid == [0])[0]
  ncolloids = len(colloids)
  typeid = traj[-1].particles.typeid[colloids]
  radii = 0.5*traj[-1].particles.diameter[colloids]
  rcut = 0.1
  lbox = traj[-1].configuration.box[:3]

  for f in range(nframes):
    rpos = traj[f].particles.position[colloids]

    # make a data frame to export to CSV
    df_pos = pd.DataFrame()
    df_pos['tag'] = colloids
    df_pos['x'] = rpos[:,0]
    df_pos['y'] = rpos[:,1]
    df_pos['z'] = rpos[:,2] 
    df_pos['typeID'] = typeid
    df_pos['radius'] = radii

    df_pos.to_csv(pos_output+str(f)+'.csv', index=False)
  print("Position data saved to CSV for "+str(nframes)+" frames")
#######


#######
## edgelistCSV
"""
# create a CSV file of all the bonded particle pairs (i.e. edges of the network) for each frame
# i,j position for each bond/edge
"""

def edgelistCSV_calc(filename):
  # open the simulation GSD file as "traj" (trajectory)
  traj = gsd.hoomd.open(filename, 'r')
  
  # set path and filename
  edge_dir_path = data_outpath+'/frame-edges'
  edge_output = edge_dir_path+'/edgelist' # + <#>.csv in python loop
  
  # check for existing CSV data 
  if os.path.exists(edge_dir_path) == False:
    os.mkdir(edge_dir_path)
  if os.path.exists(edge_dir_path) == True:
    # NOTE: this counts ALL CSV files in this directory
    nframes = 0
    # set the pattern for files ending in <number>.csv
    pattern = re.compile(r'edgelist\d+\.csv$')
    # Iterate directory
    for filename in os.listdir(edge_dir_path):
      # Check if the file has a CSV extension
      if pattern.match(filename):
        nframes += 1
    
    if nframes == len(traj):
      print('edgelist data CSV files already seem to exist for all frames. Not creating new CSV files.')
      return

  nframes = len(traj)
  colloids = np.where(traj[-1].particles.typeid == [0])[0]
  ncolloids = len(colloids)
  radii = 0.5*traj[-1].particles.diameter[colloids]
  rcut = 0.1
  lbox = traj[-1].configuration.box[:3]
  inv_lbox = 1/lbox

  # create an array of xyz positon of all colloids in all frames    
  allpos = np.zeros((nframes,ncolloids,3))
  for i in range(0,nframes):
    allpos[i,:,:] = traj[i].particles.position[colloids] 
 
  module.edgelist_calc(nframes,ncolloids,radii,allpos,lbox,rcut,edge_output)
  print("Edgelist calculation complete for "+str(nframes)+" frames")
#######


#######
## networkx analysis
"""
# use the networkx package to calculate for all frames:
#   - the number of connected components
#   - average degree (AKA average coordination number)
#   - largest connected component (LCC)
#   - average clustering coefficient
#   - average square clustering coefficient
"""

def primary_networkx_calc(filename):

  edge_output = data_outpath+'/frame-edges/edgelist' # + <#>.csv in f90

  # get the number of colloids
  traj = gsd.hoomd.open(filename, 'r')
  nframes = len(traj)
  colloids = np.where(traj[-1].particles.typeid == [0])[0]
  ncolloids = len(colloids)

  # create empty dataframe for all network data  
  network_frames_df = pd.DataFrame()
  
  # import all data into one dataframe
  frame_dfs = []
  for frame in range(nframes):
    # loop through all frames
    filepath = edge_output+str(frame)+'.csv'
    # import CSV data
    df = pd.read_csv(filepath)
    # rename colums as needed
    df = df.rename(columns={"i": "source", "j": "target"})
    df.insert(loc=0, column='frame', value=frame)
    frame_dfs.append(df)

  alledge_df = pd.concat(frame_dfs, ignore_index=True) 

  # network analysis
  for frame in range(nframes):
    df = alledge_df[alledge_df['frame'] == frame][["source", "target"]]

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
    data={'frame'                  :[frame],
          'n_components'           :[n_cc], 
          'lcc_size'               :[lcc_size],
          'ncolloids'              :[ncolloids],
          'avg_degree'             :[avg_degree],
          'avg_clustering_coeff'   :[avg_clustering_coeff],
          'avg_sq_clustering_coeff':[avg_sq_clustering_coeff]}

    res_df = pd.DataFrame(data)
    network_frames_df = pd.concat([network_frames_df, res_df])

  # write all data to one CSV file
  network_frames_df.to_csv(data_outpath+'/networkx-allframes.csv',index = False)
  print("Primary network analysis complete")
#######


####################################
""" RUN CHECKS ON A SIMULATION """
####################################

if __name__ == '__main__':
  posCSV_calc(filename)	
  edgelistCSV_calc(filename)	
  primary_networkx_calc(filename)	
