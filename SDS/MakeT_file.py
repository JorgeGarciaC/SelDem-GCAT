import sys

###############################################                                                                                                                       
# Author: Jorge García                          #                                                                                                                     
# Input: A vcf file                               #                                                                                                                   
# Output: A formated matrix (t_file) with       #                                                                                                                      
# at least one of each genotypes for each       #                                                                                                                      
# variant site                                  #                                                                                                                      
###############################################                                                                                                                        

# Input file                                                                                                                                                           
refalt_vcf = open(sys.argv[1], "r")

# Output file                                                                                                                                                          
t_file = open(sys.argv[2], "w")

# we need to go line by line:                                                                                                                                          
for line in refalt_vcf:
# if the line starts by # --> we do not print it                                                                                                                       
    if line.startswith("#"):
        pass
    else:
    # We change in order to have AA = 0, AD = 1, DD= 2                                                                                                                 
        line = line.replace("0/0", "0") # homo ancestral                                                                                                               
        line = line.replace("0/1", "1") # hetero 1                                                                                                                     
        line = line.replace("1/0", "1") # hetero 2                                                                                                                     
        line = line.replace("1/1", "2") # homo derived
                                                                                                                 
        # We remove the non needed vcf columns                                                                                                                         
        variant=line.split()
        del variant[0]
        del variant[4:8]

        # We change the positions in order to have snp_id, aa, dd, pos and genotype                                                                                    
        pos=variant[0]
        snp_id=variant[1]
        aa_allele=variant[2]
        derived_allele=variant[3]
        variant[0]=snp_id
        variant[1]=aa_allele
        variant[2]=derived_allele
        variant[3]=pos

        # We check that at least one of each genotypes is present for each variant site.                                                                               
        if '0' in variant[4:] and '1' in variant[4:] and '2' in variant[4:] :
            t_file.write('\t'.join(variant)+ "\n")

refalt_vcf.close()
t_file.close()
