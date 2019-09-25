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

def run():
  results = easy_pickle.load("results.pkl")
  for result in results:
    low_code = result.pdb_id
    low_file = file_from_code(low_code)
    low_res = result.resolution
    cntr = 0
    chains = []
    for chain_rs in result.rs:
      chains.append(chain_rs.chain_ref)

    low_chain = result.rs[0].chain_ref
    high_code = result.rs[0].match[0].pdb_code
    high_file = file_from_code(high_code)
    high_res = result.rs[0].match[0].resolution
    high_chain = result.rs[0].match[0].chain_id

    ali = align(low_file,low_chain,high_file,high_chain)
    if ali < 0.99: continue
    print(low_code, low_res, chains, high_code, high_chain, high_res,ali)
      # if ali < 0.99: 
      #   cntr += 1
      #   print(low_code,low_chain,high_code,high_chain,ali)
      #   continue
      # print(low_code,low_chain,high_code,high_chain)
      # print(ali)
      # print(chain_rs)
      # print("_"*80)
    # if (cntr is 0):
    #   print(low_code)

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
    print ("only_chain() failed", str(e))
  #   continue
  #
  ali = mmtbx.alignment.align(
    seq_a = se,
    seq_b = sx).extract_alignment()
  # print ("sequence_identity: %6.4f"%ali.calculate_sequence_identity(),\
  #   len(se), len(sx))
  return ali.calculate_sequence_identity()

if __name__ == '__main__':
  run()
