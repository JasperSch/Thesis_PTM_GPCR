#introducing NME and ACE caps to the protein gaps

import os
import re
import pandas as pd
from pymol import cmd

def Gaps_fix():
    home = 'C:/Users/Jasper/Documents/School/2e_master_bioinformatica/Thesis'
    input_list = 'tables/gaps.csv'
    #pdb_dir = 'files_clean'
    pdb_dir = 'files_clean'
    #output_directory = 'files_capped20190126'
    output_directory = 'C:/Users/Jasper/Documents/School/2e_master_bioinformatica/Thesis/files_capped20190126'

    #df = pd.read_csv(os.path.join(home, input_list), sep=' ', names=['pdb_id', 'chain', 'nme_site', 'ace_site'])
    df = pd.read_csv(os.path.join(home, input_list), sep=' ', names=['pdb_id', 'chain', 'nme_site', 'ace_site'], nrows=5 , header=20)

    for j in range(len(df)):
        pdb_id = str(df.loc[j, 'pdb_id'])
        chain = str(df.loc[j, 'chain'])
        nme_site = list(str(df.loc[j, 'nme_site']).split(','))
        ace_site = list(str(df.loc[j, 'ace_site']).split(','))

        save_var = ''

        if j <= len(df) - 2:
            pdb_id_next = str(df.loc[j + 1, 'pdb_id'])
            if pdb_id != pdb_id_next:
                save_var = 'yes'
            else:
                save_var = 'no'
        else:
            save_var = 'yes'


        pdb_file = os.path.join(home, pdb_dir, pdb_id + '_clean.pdb')

        print('Capping ' + pdb_id + ' in progress...')

        try:
            cmd.load(pdb_file, 'protein')
            cmd.hide()


            for i in range(len(nme_site)):
                cmd.do('sele pk1, chain ' + chain + ' and resi ' + nme_site[i] + ' and n. C')
                cmd.do('editor.attach_amino_acid("pk1","nme")')
                print('Added NME to chain '+ chain +' amino acid ' + nme_site[i])

            for i in range(len(ace_site)):
                cmd.do('sele pk1, chain ' + chain + ' and resi ' + ace_site[i] + ' and n. N')
                cmd.do('editor.attach_amino_acid("pk1","ace")')
                print('Added ACE to chain '+ chain +' amino acid ' + ace_site[i])


            if save_var == 'yes':
                output_path = os.path.join(home, output_directory, pdb_id + '_capped.ent')
                cmd.save(output_path, 'protein')
                output_name = pdb_id + "_capped.ent"
                print('Saving file: ' + output_name)

                cmd.delete('all')

            elif save_var == 'no':
                pass


        except IOError as e:
            print("I/O_error({0}:_{1}".format(e.erno, e.sterror))

cmd.extend("Gaps_fix", Gaps_fix)
