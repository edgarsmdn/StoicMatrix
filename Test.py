"""
Reaction network stoichiometry matrix generation Test-File
-------------------------------------------------
Author: Edgar Ivan Sanchez Medina
Email: sanchez@mpi-magdeburg.mpg.de
"""
from ReactionNetwork import ReactionNetwork

main_reactions = [
    'A+B>C+D',
    'F>A',
    'D+E>F'
    ]


network_id = '001'

network_1 = ReactionNetwork(main_reactions) # Initialize reacion network
network_1.write_subnetworks(network_id)     # Write subnetworks for each main reaction 
network_1.write_SMatrix(network_id)         # Write Stoichiometric Matrix

