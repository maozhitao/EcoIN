# coding: utf-8

#!pip uninstall networkx
#!pip install networkx==1.9.1


# import package

'''
Created on 2019/01/08

@author: Mr.Mao
'''
import networkx as nx
import re
import matplotlib.pyplot as plt
import copy
import os
import shutil
import glob

# get_met_reaction_list,delete currency metabolites
def get_met_reaction_list(met_file,reaction_file):    
    metlist=[]
    reactionlist=[]
    for eachdata in open(met_file):   
        #metdata=eachdata.split('\r\n')
        metdata=eachdata.split('\n')     
        metlist.append(metdata[0])    
    for eachdata in open(reaction_file):   
        #metdata=eachdata.split('\r\n')
        metdata=eachdata.split('\n') 
        reactionlist.append(metdata[0])
         
    return [metlist,reactionlist]

def initial_network(G,total_file):
    #regulator
    network_list=[]
    for eachdata in open(total_file):   
        network_list.append(eachdata)           
    for eachline in network_list:
        #tm=eachline.split('\r\n') 
        tm=eachline.split('\n')
        data=tm[0].split('\t')        
        if re.search('cpi',data[3]):
            G.add_node(data[0], nodetype='met')
            G.add_node(data[1], nodetype='gene')
            if data[4]=='0':
                G.add_edge(data[0], data[1], edgetype='cpi',data_att=data[2])    
            elif data[4]=='1':
                G.add_edge(data[0], data[1], edgetype='cpi',data_att=data[2])
                G.add_edge(data[1], data[0], edgetype='cpi',data_att=data[2]) 
        elif re.search('reaction_gene',data[3]):
            G.add_node(data[0], nodetype='gene')
            G.add_node(data[1], nodetype='reaction') 
            if data[4]=='0':
                G.add_edge(data[0], data[1], edgetype='reaction_gene',data_att=data[2])    
            elif data[4]=='1':
                G.add_edge(data[0], data[1], edgetype='reaction_gene',data_att=data[2])
                G.add_edge(data[1], data[0], edgetype='reaction_gene',data_att=data[2])     
        elif re.search('ppi',data[3]):
            G.add_node(data[0], nodetype='gene')
            G.add_node(data[1], nodetype='gene')   
            if data[4]=='0':
                G.add_edge(data[0], data[1], edgetype='ppi',data_att=data[2])    
            elif data[4]=='1':
                G.add_edge(data[0], data[1], edgetype='ppi',data_att=data[2])
                G.add_edge(data[1], data[0], edgetype='ppi',data_att=data[2])     
        elif re.search('tf-gene',data[3]):
            G.add_node(data[0], nodetype='gene')
            G.add_node(data[1], nodetype='gene')   
            if data[4]=='0':
                G.add_edge(data[0], data[1], edgetype='tf_gene',data_att=data[2])    
            elif data[4]=='1':
                G.add_edge(data[0], data[1], edgetype='tf_gene',data_att=data[2])
                G.add_edge(data[1], data[0], edgetype='tf_gene',data_att=data[2]) 
        elif re.search('leader_peptide',data[3]):
            if data[2]=='Binding;':
                G.add_node(data[0], nodetype='met')
                G.add_node(data[1], nodetype='leader')   
                if data[4]=='0':
                    G.add_edge(data[0], data[1], edgetype='leader_peptide',data_att=data[2])    
                elif data[4]=='1':
                    G.add_edge(data[0], data[1], edgetype='leader_peptide',data_att=data[2])
                    G.add_edge(data[1], data[0], edgetype='leader_peptide',data_att=data[2])                 
            elif data[2]=='Attenuation;':
                G.add_node(data[0], nodetype='leader')
                G.add_node(data[1], nodetype='gene')   
                if data[4]=='0':
                    G.add_edge(data[0], data[1], edgetype='leader_peptide',data_att=data[2])    
                elif data[4]=='1':
                    G.add_edge(data[0], data[1], edgetype='leader_peptide',data_att=data[2])
                    G.add_edge(data[1], data[0], edgetype='leader_peptide',data_att=data[2])                 
                         
        elif re.search('rna_regulation',data[3]):
            G.add_node(data[0], nodetype='gene')
            G.add_node(data[1], nodetype='gene')   
            if data[4]=='0':
                G.add_edge(data[0], data[1], edgetype='rna_regulation',data_att=data[2])    
            elif data[4]=='1':
                G.add_edge(data[0], data[1], edgetype='rna_regulation',data_att=data[2])
                G.add_edge(data[1], data[0], edgetype='rna_regulation',data_att=data[2])                 
        elif re.search('sigma_factor_regulation',data[3]):
            G.add_node(data[0], nodetype='gene')
            G.add_node(data[1], nodetype='gene')   
            if data[4]=='0':
                G.add_edge(data[0], data[1], edgetype='sigma_factor_regulation',data_att=data[2])    
            elif data[4]=='1':
                G.add_edge(data[0], data[1], edgetype='sigma_factor_regulation',data_att=data[2])
                G.add_edge(data[1], data[0], edgetype='sigma_factor_regulation',data_att=data[2])                 
        elif re.search('signal_transduction',data[3]):  
            if re.search('[c]',data[0]):
                #print data[0]
                G.add_node(data[0], nodetype='met')
                G.add_node(data[1], nodetype='gene')   
                if data[4]=='0':
                    G.add_edge(data[0], data[1], edgetype='Signal_transduction',data_att=data[2])    
                elif data[4]=='1':
                    G.add_edge(data[0], data[1], edgetype='Signal_transduction',data_att=data[2])
                    G.add_edge(data[1], data[0], edgetype='Signal_transduction',data_att=data[2])                 
            else:
                G.add_node(data[0], nodetype='gene')
                G.add_node(data[1], nodetype='gene')   
                if data[4]=='0':
                    G.add_edge(data[0], data[1], edgetype='Signal_transduction',data_att=data[2])    
                elif data[4]=='1':
                    G.add_edge(data[0], data[1], edgetype='Signal_transduction',data_att=data[2])
                    G.add_edge(data[1], data[0], edgetype='Signal_transduction',data_att=data[2])                    
        elif re.search('Protein_modification',data[3]):
            if re.search('[c]',data[0]):
                #print data[0]
                G.add_node(data[0], nodetype='met')
                G.add_node(data[1], nodetype='gene')   
                if data[4]=='0':
                    G.add_edge(data[0], data[1], edgetype='Protein_modification',data_att=data[2])    
                elif data[4]=='1':
                    G.add_edge(data[0], data[1], edgetype='Protein_modification',data_att=data[2])
                    G.add_edge(data[1], data[0], edgetype='Protein_modification',data_att=data[2])                 
            else:
                G.add_node(data[0], nodetype='gene')
                G.add_node(data[1], nodetype='gene')   
                if data[4]=='0':
                    G.add_edge(data[0], data[1], edgetype='Protein_modification',data_att=data[2])    
                elif data[4]=='1':
                    G.add_edge(data[0], data[1], edgetype='Protein_modification',data_att=data[2])
                    G.add_edge(data[1], data[0], edgetype='Protein_modification',data_att=data[2])                                 
        elif re.search('substrate',data[2]):
            G.add_node(data[0], nodetype='met')
            G.add_node(data[1], nodetype='reaction')
            if data[4]=='0':
                G.add_edge(data[0], data[1], edgetype='reaction_met',data_att=data[2])    
            elif data[4]=='1':
                G.add_edge(data[0], data[1], edgetype='reaction_met',data_att=data[2])
                G.add_edge(data[1], data[0], edgetype='reaction_met',data_att=data[2]) 
        elif re.search('product',data[2]):
            G.add_node(data[1], nodetype='met')
            G.add_node(data[0], nodetype='reaction')
            if data[4]=='0':
                G.add_edge(data[0], data[1], edgetype='reaction_met',data_att=data[2])    
            elif data[4]=='1':
                G.add_edge(data[0], data[1], edgetype='reaction_met',data_att=data[2])
                G.add_edge(data[1], data[0], edgetype='reaction_met',data_att=data[2]) 
                                 
        else:
            pass
    return G

# initial regulator network
def initial_regulator_network(G,total_file):
    #regulator
    network_list=[]
    for eachdata in open(total_file):   
        network_list.append(eachdata)           
    for eachline in network_list:
        #tm=eachline.split('\r\n') 
        tm=eachline.split('\n')
        data=tm[0].split('\t')        
        if re.search('cpi',data[3]):
            G.add_node(data[0], nodetype='met')
            G.add_node(data[1], nodetype='gene')
            if data[4]=='0':
                G.add_edge(data[0], data[1], edgetype='cpi',data_att=data[2])    
            elif data[4]=='1':
                G.add_edge(data[0], data[1], edgetype='cpi',data_att=data[2])
                G.add_edge(data[1], data[0], edgetype='cpi',data_att=data[2]) 
        elif re.search('reaction_gene',data[3]):
            G.add_node(data[0], nodetype='gene')
            G.add_node(data[1], nodetype='reaction') 
            if data[4]=='0':
                G.add_edge(data[0], data[1], edgetype='reaction_gene',data_att=data[2])    
            elif data[4]=='1':
                G.add_edge(data[0], data[1], edgetype='reaction_gene',data_att=data[2])
                G.add_edge(data[1], data[0], edgetype='reaction_gene',data_att=data[2])     
        elif re.search('ppi',data[3]):
            G.add_node(data[0], nodetype='gene')
            G.add_node(data[1], nodetype='gene')   
            if data[4]=='0':
                G.add_edge(data[0], data[1], edgetype='ppi',data_att=data[2])    
            elif data[4]=='1':
                G.add_edge(data[0], data[1], edgetype='ppi',data_att=data[2])
                G.add_edge(data[1], data[0], edgetype='ppi',data_att=data[2])     
        elif re.search('tf-gene',data[3]):
            G.add_node(data[0], nodetype='gene')
            G.add_node(data[1], nodetype='gene')   
            if data[4]=='0':
                G.add_edge(data[0], data[1], edgetype='tf_gene',data_att=data[2])    
            elif data[4]=='1':
                G.add_edge(data[0], data[1], edgetype='tf_gene',data_att=data[2])
                G.add_edge(data[1], data[0], edgetype='tf_gene',data_att=data[2]) 
        elif re.search('leader_peptide',data[3]):
            if data[2]=='Binding;':
                G.add_node(data[0], nodetype='met')
                G.add_node(data[1], nodetype='leader')   
                if data[4]=='0':
                    G.add_edge(data[0], data[1], edgetype='leader_peptide',data_att=data[2])    
                elif data[4]=='1':
                    G.add_edge(data[0], data[1], edgetype='leader_peptide',data_att=data[2])
                    #G.add_edge(data[1], data[0], edgetype='leader_peptide',data_att=data[2])                 
            elif data[2]=='Attenuation;':
                G.add_node(data[0], nodetype='leader')
                G.add_node(data[1], nodetype='gene')   
                if data[4]=='0':
                    G.add_edge(data[0], data[1], edgetype='leader_peptide',data_att=data[2])    
                elif data[4]=='1':
                    G.add_edge(data[0], data[1], edgetype='leader_peptide',data_att=data[2])
                    G.add_edge(data[1], data[0], edgetype='leader_peptide',data_att=data[2])                 
                         
        elif re.search('rna_regulation',data[3]):
            G.add_node(data[0], nodetype='gene')
            G.add_node(data[1], nodetype='gene')   
            if data[4]=='0':
                G.add_edge(data[0], data[1], edgetype='rna_regulation',data_att=data[2])    
            elif data[4]=='1':
                G.add_edge(data[0], data[1], edgetype='rna_regulation',data_att=data[2])
                G.add_edge(data[1], data[0], edgetype='rna_regulation',data_att=data[2])                 
        elif re.search('sigma_factor_regulation',data[3]):
            G.add_node(data[0], nodetype='gene')
            G.add_node(data[1], nodetype='gene')   
            if data[4]=='0':
                G.add_edge(data[0], data[1], edgetype='sigma_factor_regulation',data_att=data[2])    
            elif data[4]=='1':
                G.add_edge(data[0], data[1], edgetype='sigma_factor_regulation',data_att=data[2])
                G.add_edge(data[1], data[0], edgetype='sigma_factor_regulation',data_att=data[2])                 
        elif re.search('signal_transduction',data[3]):  
            if re.search('[c]',data[0]):
                #print data[0]
                G.add_node(data[0], nodetype='met')
                G.add_node(data[1], nodetype='gene')   
                if data[4]=='0':
                    G.add_edge(data[0], data[1], edgetype='Signal_transduction',data_att=data[2])    
                elif data[4]=='1':
                    G.add_edge(data[0], data[1], edgetype='Signal_transduction',data_att=data[2])
                    G.add_edge(data[1], data[0], edgetype='Signal_transduction',data_att=data[2])                 
            else:
                G.add_node(data[0], nodetype='gene')
                G.add_node(data[1], nodetype='gene')   
                if data[4]=='0':
                    G.add_edge(data[0], data[1], edgetype='Signal_transduction',data_att=data[2])    
                elif data[4]=='1':
                    G.add_edge(data[0], data[1], edgetype='Signal_transduction',data_att=data[2])
                    G.add_edge(data[1], data[0], edgetype='Signal_transduction',data_att=data[2])                    
        elif re.search('Protein_modification',data[3]):
            if re.search('[c]',data[0]):
                #print data[0]
                G.add_node(data[0], nodetype='met')
                G.add_node(data[1], nodetype='gene')   
                if data[4]=='0':
                    G.add_edge(data[0], data[1], edgetype='Protein_modification',data_att=data[2])    
                elif data[4]=='1':
                    G.add_edge(data[0], data[1], edgetype='Protein_modification',data_att=data[2])
                    G.add_edge(data[1], data[0], edgetype='Protein_modification',data_att=data[2])                 
            else:
                G.add_node(data[0], nodetype='gene')
                G.add_node(data[1], nodetype='gene')   
                if data[4]=='0':
                    G.add_edge(data[0], data[1], edgetype='Protein_modification',data_att=data[2])    
                elif data[4]=='1':
                    G.add_edge(data[0], data[1], edgetype='Protein_modification',data_att=data[2])
                    G.add_edge(data[1], data[0], edgetype='Protein_modification',data_att=data[2])                               
                                 
        else:
            pass
    return G

# initial regulonDB_tf network
def regulonDB_tfnetwork(G,total_file):
    #regulator
    network_list=[]
    for eachdata in open(total_file):   
        network_list.append(eachdata)           
    for eachline in network_list: 
        tm=eachline.split('\n')
        data=tm[0].split('\t')        
        if re.search('tf-gene',data[3]):
            G.add_node(data[0], nodetype='gene')
            G.add_node(data[1], nodetype='gene')   
            if data[4]=='0':
                G.add_edge(data[0], data[1], edgetype='tf_gene',data_att=data[2])    
            elif data[4]=='1':
                G.add_edge(data[0], data[1], edgetype='tf_gene',data_att=data[2])
                G.add_edge(data[1], data[0], edgetype='tf_gene',data_att=data[2])         
        else:
            pass
    return G
    
# initial metabolic network including reaction_enz
def met_network(G,total_file):
    #regulator
    network_list=[]
    for eachdata in open(total_file):   
        network_list.append(eachdata)           
    for eachline in network_list:  
        tm=eachline.split('\n')
        data=tm[0].split('\t')        
        if re.search('substrate',data[2]):
            G.add_node(data[0], nodetype='met')
            G.add_node(data[1], nodetype='reaction')
            if data[4]=='0':
                G.add_edge(data[0], data[1], edgetype='reaction_met',data_att=data[2])    
            elif data[4]=='1':
                G.add_edge(data[0], data[1], edgetype='reaction_met',data_att=data[2])
                G.add_edge(data[1], data[0], edgetype='reaction_met',data_att=data[2]) 
        elif re.search('reaction_gene',data[3]):
            G.add_node(data[0], nodetype='gene')
            G.add_node(data[1], nodetype='reaction') 
            if data[4]=='0':
                G.add_edge(data[0], data[1], edgetype='reaction_gene',data_att=data[2])    
            elif data[4]=='1':
                G.add_edge(data[0], data[1], edgetype='reaction_gene',data_att=data[2])
                G.add_edge(data[1], data[0], edgetype='reaction_gene',data_att=data[2])                   
        elif re.search('product',data[2]):
            G.add_node(data[1], nodetype='met')
            G.add_node(data[0], nodetype='reaction')
            if data[4]=='0':
                G.add_edge(data[0], data[1], edgetype='reaction_met',data_att=data[2])    
            elif data[4]=='1':
                G.add_edge(data[0], data[1], edgetype='reaction_met',data_att=data[2])
                G.add_edge(data[1], data[0], edgetype='reaction_met',data_att=data[2])                        
        else:
            pass 
    return G

#get met to reaction n length path from regulate network
def search_regulation_pair_length(G,unimet,unimet2,n):
    outdata=[] 
    outdata2=[]   
    for eachmet in unimet:
        for eachmet2 in unimet2:
            if (eachmet2 in G.nodes()) and (eachmet in G.nodes()):            
                try:      
                    path=nx.all_simple_paths(G, eachmet, eachmet2,n)
                    for eachpath in path:
                        if eachpath:
                            outstr1 = ''
                            outstr2 = ''
                            outstr3 = ''
                            i = 0
                            j = 0
                            for eachdata in eachpath:
                                outstr1 = outstr1 + '\t' + eachdata
                                outstr2 = outstr2 + '\t' + G.node[eachdata]['nodetype'] 
                            while i < len(eachpath) - 1:                 
                                #outstr3 = outstr3 + '\t' + G.edge[eachpath[i]][eachpath[i + 1]]['edgetype'] + '_' + G.edge[eachpath[i]][eachpath[i + 1]]['data_att']
                                outstr3 = outstr3 + '\t' + G.edges[eachpath[i],eachpath[i + 1]]['edgetype'] + '_' + G.edges[eachpath[i],eachpath[i + 1]]['data_att']
                                #if G.edge[eachpath[i]][eachpath[i + 1]]['data_att'] == 'Binding;':
                                if G.edges[eachpath[i],eachpath[i + 1]]['data_att'] == 'Binding;':
                                    j = j + 1
                                i = i + 1
                            if j < 3:
                                outstr = outstr1 + outstr2 + outstr3 + '\r\n'
                                outdata.append(outstr) 
                        if len(eachpath) == 3:                   
                            G4 = G.copy()
                            G4.remove_edge(eachpath[0], eachpath[1])
                            path = nx.all_shortest_paths(G4, eachmet, eachmet2)
                            for eachpath in path:
                                if eachpath and len(eachpath)>3: 
                                    outstr1 = ''
                                    outstr2 = ''
                                    outstr3 = ''
                                    i = 0
                                    j = 0
                                    for eachdata in eachpath:
                                        outstr1 = outstr1 + '\t' + eachdata
                                        outstr2 = outstr2 + '\t' + G4.node[eachdata]['nodetype'] 
                                    while i < len(eachpath) - 1:                 
                                        #outstr3 = outstr3 + '\t' + G4.edge[eachpath[i]][eachpath[i + 1]]['edgetype'] + '_' + G4.edge[eachpath[i]][eachpath[i + 1]]['data_att']
                                        outstr3 = outstr3 + '\t' + G4.edges[eachpath[i],eachpath[i + 1]]['edgetype'] + '_' + G4.edges[eachpath[i],eachpath[i + 1]]['data_att']
                                        #if G4.edge[eachpath[i]][eachpath[i + 1]]['data_att'] == 'Binding;':
                                        if G4.edges[eachpath[i],eachpath[i + 1]]['data_att'] == 'Binding;':
                                            j = j + 1
                                        i = i + 1
                                    if j < 3:
                                        outstr = outstr1 + outstr2 + outstr3 + '\r\n'
                                        outdata.append(outstr)  
                      
                except nx.NetworkXNoPath:
                    pass           
    for eachdi in list(set(outdata)):
        if re.search('ppi_Binding;\tppi_Binding', eachdi):
            pass
        else:
            outdata2.append(eachdi)
            #print(eachdi)
    return outdata2


#reaction，gene and metbolite annotation
def reaction_anontation(reaction,total_reaction_file):   
    outdata=''
    for eachdata in open(total_reaction_file):   
        #tm=eachdata.split('\r\n')
        tm=eachdata.split('\n')
        data=tm[0].split('\t')
        if data[0]==reaction:
            outdata= data[5]
    return outdata
               
def gene_anontation(gene,total_gene_file): 
    outdata=''  
    for eachdata in open(total_gene_file):   
        data=eachdata.split('\t')
        if data[0]==gene:
            outdata=data[2]+' , '+data[3]
    return outdata    

def met_anontation(met,total_met_file):  
    outdata='' 
    for eachdata in open(total_met_file):   
        data=eachdata.split('\t')
        if data[0]==met:
            outdata=data[1]
    return outdata   

def ano_nocircle(regulator_list,outfilename,total_reaction_file,total_gene_file,total_met_file):
    outFile2=open(outfilename,'w')
    for eachreg in regulator_list:
        htline=eachreg.split('\r\n')[0]  
        if re.search('reaction_gene_Catalysis',htline):
            data=htline.split('\tmet\t')
            firstdata=data[0].split('\t')
            data2=htline.split('\treaction\t')
            seconddata=data2[1].split('\t')
            i=0
            for eachdata in firstdata:  
                if eachdata:
                    if i<len(seconddata):
                        if seconddata[i]=='ppi_Binding':
                            direction_str=' <==> '
                        else:
                            direction_str=' --> '        
                        usedata=eachdata+'['+reaction_anontation(eachdata,total_reaction_file)+gene_anontation(eachdata,total_gene_file)+met_anontation(eachdata,total_met_file)+']('+seconddata[i]+')'+direction_str
                    else:
                        usedata=eachdata+'['+reaction_anontation(eachdata,total_reaction_file)+gene_anontation(eachdata,total_gene_file)+met_anontation(eachdata,total_met_file)+']'
                    outFile2.write(usedata) 
                    i=i+1  
            outFile2.write('\n') 
    outFile2.close()    

#feedback annotation and category
def regulation_classification(classification_file,infile,an_name):
    outfile=r'%s/%s_Regulation_of_enzyme_activity.txt'%(classification_file,an_name)
    outFile=open(outfile,'w')     
    outfile2=r'%s/%s_Transcriptional_regulation.txt'%(classification_file,an_name)
    outFile2=open(outfile2,'w')
    outfile3=r'%s/%s_Complex_regulation.txt'%(classification_file,an_name)
    outFile3=open(outfile3,'w')   
    
    for eachdata in open(infile):
        data=eachdata.split('-->')
        if len(data)==3:
            outFile.write(eachdata)
        elif len(data)==4:
            outFile2.write(eachdata)
        else:
            outFile3.write(eachdata)              
                 
    outFile.close()
    outFile2.close()
    outFile3.close()
    
#feedback annotation and category
def regulation_classification_new(classification_file,infile,an_name):
    outfile=r'%s/%s_Regulation_of_enzyme_activity.txt'%(classification_file,an_name)
    outFile=open(outfile,'w')     
    outfile2=r'%s/%s_Regulation_of_the_amount_of_enzyme.txt'%(classification_file,an_name)
    outFile2=open(outfile2,'w')  
    
    for eachdata in open(infile):
        data=eachdata.split('-->')
        if re.search('cpi_',data[-3]):
            outFile.write(eachdata)
        elif re.search('Signal_transduction_',data[-3]):        
            outFile.write(eachdata)
        elif re.search('Protein_modification_',data[-3]):        
            outFile.write(eachdata)
        else:
            outFile2.write(eachdata)
             
                 
    outFile.close()
    outFile2.close()

# create file to save results
def create_file(store_path):
    if os.path.exists(store_path):
        print("path exists")
        shutil.rmtree(store_path)
        os.mkdir(store_path)
    else:      
        os.mkdir(store_path)
        print(store_path) 
        
def run_feed_feedforwad(eachname,regulator_n,pathway_file,regulation_file,classification_file): 
    #initial_file  
    met_file_name=eachname+'_met'
    reaction_file_name=eachname+'_reaction'
    met_file=r'./superpathway/%s.txt'%(met_file_name)
    reaction_file=r'./superpathway/%s.txt'%(reaction_file_name)
    DG_core=nx.DiGraph()
    DG_core=initial_regulator_network(DG_core,regulation_file) 
    total_reaction_file=r'./reaction_information.txt'
    total_gene_file=r'./ecogene_information.txt'
    total_met_file=r'./metabolite_information.txt'
    
    data_ini=get_met_reaction_list(met_file,reaction_file)
    metlist=data_ini[0] 
    reactionlist=data_ini[1]
    eachname=met_file_name+'_2_'+reaction_file_name
    outfile2=r'%s%s_feed_feedforwad_ano.txt'%(pathway_file,eachname) 
 
    regulator_list=search_regulation_pair_length(DG_core,metlist,reactionlist,regulator_n)

    ano_nocircle(regulator_list,outfile2,total_reaction_file,total_gene_file,total_met_file)
    
    regulation_classification_new(classification_file,outfile2,eachname)
    print(eachname+' OK!')

def get_met_list(reaction_metfile,reaction_file):
    met_list=[]
    reactionlist=[]    
    for eachdata in open(reaction_file):   
        #metdata=eachdata.split('\r\n')
        reactiondata=eachdata.split('\n') 
        reactionlist.append(reactiondata[0])
        for eachrecord in open(reaction_metfile):
            reaction_met_record=eachrecord.split('\t')
            if re.search('product',reaction_met_record[2]):
                if reactiondata[0]==reaction_met_record[0]:
                    met_list.append(reaction_met_record[1])            
            elif re.search('substrate',reaction_met_record[2]):
                if reactiondata[0]==reaction_met_record[1]:
                    met_list.append(reaction_met_record[0])        
    return [met_list,reactionlist]  

def run_feed_feedforwad_new(eachname,regulator_n,pathway_file,regulation_file,classification_file,reaction_metfile,infile): 
    #initial_file    
    reaction_file_name=eachname+'_reaction'
    print(reaction_file_name)
    reaction_file=r'./ini_pathway190617/%s/%s.txt'%(infile,reaction_file_name)
    DG_core=nx.DiGraph()
    DG_core=initial_regulator_network(DG_core,regulation_file) 
    total_reaction_file=r'./reaction_information.txt'
    total_gene_file=r'./ecogene_information.txt'
    total_met_file=r'./metabolite_information.txt'
    
    data_ini=get_met_list(reaction_metfile,reaction_file)
    metlist=data_ini[0] 
    reactionlist=data_ini[1]
    print(metlist)
    print(reactionlist)
    eachname=reaction_file_name+'_%s'%infile
    outfile2=r'%s%s_feed_feedforwad_ano.txt'%(pathway_file,eachname) 
 
    regulator_list=search_regulation_pair_length(DG_core,metlist,reactionlist,regulator_n)

    ano_nocircle(regulator_list,outfile2,total_reaction_file,total_gene_file,total_met_file)
    
    regulation_classification_new(classification_file,outfile2,eachname)
    print(eachname+' OK!')
    