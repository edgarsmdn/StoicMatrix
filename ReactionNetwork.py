"""
Reaction network stoichiometry matrix generation
-------------------------------------------------
Author: Edgar Ivan Sanchez Medina
Email: sanchez@mpi-magdeburg.mpg.de
"""
import os
import numpy as np
import pandas as pd


class ReactionNetwork():
    
    def __init__(self, main_reactions):
        
        self.reactions = main_reactions
        self.n_subnetworks = len(main_reactions)
        
    def write_SMatrix(self, network_id):
        name_folder = 'Network_' + network_id
        network = []
        for i in range(self.n_subnetworks):
            subnetwork = self.get_subnetwork(self.reactions[i], i+1)
            network += subnetwork
        SMatrix = self.get_SMatrix(network)
        SMatrix.to_csv(name_folder+'/SMatrix.csv')
        
    def write_subnetworks(self, network_id):
        name_folder = 'Network_' + network_id
        if not os.path.exists(name_folder):
            os.makedirs(name_folder)
        for i in range(self.n_subnetworks):
            subnetwork = self.get_subnetwork(self.reactions[i], i+1)
            
            with open(name_folder+'/Subnetwork_'+str(i+1)+'.txt', 'w') as f:
                f.write('Subnetwork_'+str(i+1))
                f.write('\nMain reaction: '+ str(self.reactions[i])+'\n')
                f.write('-'*40)
                f.write('\n')
                for item in subnetwork:
                    f.write("%s\n" % item)
    def get_SMatrix(self, network):
        species = []
        for reaction in network:
            reactants, products = reaction.rsplit('>')
            reactants_indv = reactants.rsplit('+')
            products_indv = products.rsplit('+')
            for r in reactants_indv:
                if r not in species:
                    species.append(r)
            for p in products_indv:
                if p not in species:
                    species.append(p)
                    
        # initialized matrix
        SMatrix_dict = {}
        for indx, s in enumerate(species):
            SMatrix_dict[s] = indx
        SMatrix = np.zeros((len(species), len(network)))
            
        for j, reaction in enumerate(network):
            reactants, products = reaction.rsplit('>')
            reactants_indv = reactants.rsplit('+')
            products_indv  = products.rsplit('+')
            for r in reactants_indv:
                SMatrix[SMatrix_dict[r], j] = -1
            for p in products_indv:
                SMatrix[SMatrix_dict[p], j] = 1
         
        SMatrix_df = pd.DataFrame(SMatrix, index=species, columns=network)
        
        return SMatrix_df
            
    def get_subnetwork(self, reaction, subnetwork_id):
        reactants, products = reaction.rsplit('>')
        self.reactants_indv = reactants.rsplit('+')
        self.products_indv  = products.rsplit('+')
        
        self.n_reactants = len(self.reactants_indv)
        self.n_products = len(self.products_indv)
        
        self.enzime = 'E'+str(subnetwork_id)
        
        # 1. Biding reactions
        bidings = self.get_biding_reactions() 
        # 2. Enzime complex + other reactants reactions
        if self.n_reactants > 1:
            Complex_reactions, unique_reactants = self.get_complex_reactions(bidings) 
        else:
            Complex_reactions, unique_reactants = [], self.reactants_indv
        # 3. ReactantsComplex to ProductComplex reactions
        prod_complex, ProductComplex = self.get_reactcomplex2prodcomplex(unique_reactants)
        # 4. Unbiding reactions
        unbidings = self.get_unbiding_reactions(ProductComplex)
        
        subnetwork = bidings + Complex_reactions + prod_complex + unbidings
        return subnetwork
    
    def get_biding_reactions(self):
        bidings = []
        for r in self.reactants_indv:
            biding = r+'+'+ self.enzime +'>'+r+self.enzime
            bidings.append(biding)
        return bidings
    
    def get_complex_reactions(self, bidings):
        enzimecomplex_OtherReactants = []
        enzime = self.enzime
        for b in bidings:
            EnzimeComplex = b.rsplit('>')[1]
            reactant_in_complex = EnzimeComplex.rsplit(enzime)[0]
            other_reactants = [r for r in self.reactants_indv if reactant_in_complex != r]
            for i in range(self.n_reactants-1):
                other_reactant = other_reactants[i]
                enzime_complex_other = other_reactant +'+'+ EnzimeComplex +'>'+ other_reactant+EnzimeComplex
                enzimecomplex_OtherReactants.append(enzime_complex_other)         
        # Replace name of complexes being equal
        unique_reactants = []
        for r in enzimecomplex_OtherReactants:
            complete_complex = r.rsplit('>')[1]
            reacts_bided = complete_complex.rsplit(enzime)[0]
            reacts_bided_sorted = ''.join(sorted(reacts_bided))
            if reacts_bided_sorted not in unique_reactants:
                unique_reactants.append(reacts_bided_sorted)
        enzimecomplex_OtherReactants_unique = []
        for r in enzimecomplex_OtherReactants:
            complete_complex = r.rsplit('>')[1]
            r = r.replace(complete_complex, unique_reactants[0]+enzime)
            enzimecomplex_OtherReactants_unique.append(r)
        return enzimecomplex_OtherReactants_unique, unique_reactants
    
    def get_reactcomplex2prodcomplex(self, unique_reactants):
        prod_complex = []
        products_sorted = ''.join(sorted(self.products_indv))
        ReactantComplex = unique_reactants[0]+self.enzime
        ProductComplex  =  ReactantComplex.replace(unique_reactants[0], products_sorted)
        prod_complex.append(ReactantComplex+'>'+ProductComplex)
        
        return prod_complex, ProductComplex
    
    def get_unbiding_reactions(self, ProductComplex):
        unbidings = []
        prods_unbid = [x for x in ProductComplex.rsplit('E')[0]]
        for i in range(self.n_products):
            new_complex = ProductComplex.replace(prods_unbid[i],'')
            unbidings.append(ProductComplex+'>'+prods_unbid[i] +'+'+new_complex)
            if self.n_products > 1:
                unbidings.append(new_complex+'>'+prods_unbid[i-1]+'+'+self.enzime)
        return unbidings
    

  
    
        
        
        

