# StoicMatrix
Generation of stoichiometric matrix from a reaction network

## Reaction network: list of strings

A reaction network class must be initialized with a list of main reactions represented as strings: `A+B>C+D`.

## Write sub-networks

Once a reaction network is constructed, the function `.write_subnetworks(network_id)` can be used to write all the sub-networks for each main reaction.
Within these sub-networks, an enzime is denoted by `E#`, where `#` is replaced by the number of the corresponding main reaction. This function will create a set of .txt files containing each of the sub-networks. The `network_id` is just a string identifier for the current network that can be set arbitrary by the user.

## Write stoichiometric matrix

To obtain the stoichiometric matrix, one must use the function `.write_SMatrix(network_id)` which will create a .csv file containing the stoichiometric matrix of the reaction network. The `network_id` is just a string identifier for the current network that can be set arbitrary by the user.

## Example

An example of how to use it is within the `Test.py` file.

## Limitations

Only reactions with one or two reactants and products (without considering the enzime) are tested.

## License

This repository contains a the following [license](https://github.com/edgarsmdn/StoicMatrix/blob/main/LICENSE)
