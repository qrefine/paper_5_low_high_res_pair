from __future__ import absolute_import, division, print_function
import iotbx.pdb
import os
import mmtbx.alignment
from libtbx import easy_pickle

path_to_pdb = os.environ.get("PDB_MIRROR_PDB")

def file_from_code(code):
  file_path = None
  if(os.path.isdir(path_to_pdb)):
    file_path = os.path.join(
      path_to_pdb, code.lower()[1:3],"pdb"+code.lower()+".ent.gz")
  if(file_path is None): return None
  assert os.path.isfile(file_path)
  return file_path

def check_data(code):
  data_path = "/home/pdb/structure_factors/"
  file_path = os.path.join(
    data_path, code.lower()[1:3], "r"+code.lower()+"sf.ent.gz")
  if (not os.path.isfile(file_path)): return False
  return True

def run():
  results = easy_pickle.load("results.pkl")
  for result in results:
    if (result.method is "ELECTRON MICROSCOPY"): continue
    low_code = result.pdb_id
    if (not check_data(low_code)): continue
    low_file = file_from_code(low_code)

    pdb_inp = iotbx.pdb.input(low_file)
    hierarchy = pdb_inp.construct_hierarchy()
    atom_size = hierarchy.atoms().size()

    low_res = result.resolution
    chains = []
    codes = []
    ## result.model_rs means each chain of low_res pdbs pairs for sames model in model_rs
    for code in result.model_rs:
      codes.append(code)
    
    for chain_rs in result.rs:
      chains.append(chain_rs.chain_ref)

    chain_str = ".".join(chains)
    low_chain = result.rs[0].chain_ref
    for m in result.rs[0].match:
      high_code = m.pdb_code
      if (high_code in codes):
        high_file = file_from_code(high_code)
        high_res = m.resolution
        high_chain = m.chain_id

        ali = align(low_file,low_chain,high_file,high_chain)
        if ali < 0.99: continue
        print(",".join([str(atom_size), low_code, str(low_res), chain_str, high_code, high_chain, \
          str(high_res),str(ali)]))

def align(low_file,low_chain,high_file,high_chain):
  he = iotbx.pdb.input(file_name=low_file).construct_hierarchy()
  hx = iotbx.pdb.input(file_name=high_file).construct_hierarchy()

  sel = he.atom_selection_cache().selection("protein and chain %s"%low_chain)
  he = he.select(sel)
  sel = hx.atom_selection_cache().selection("protein and chain %s"%high_chain)
  hx = hx.select(sel)

  try:
    se = he.only_chain().as_padded_sequence()
    sx = hx.only_chain().as_padded_sequence()
  except Exception as e:
    # print ("only_chain() failed", str(e))
  #   continue
    return 0
  #
  ali = mmtbx.alignment.align(
    seq_a = se,
    seq_b = sx).extract_alignment()
  # print ("sequence_identity: %6.4f"%ali.calculate_sequence_identity(),\
  #   len(se), len(sx))
  return ali.calculate_sequence_identity()

if __name__ == '__main__':
  run()
