from __future__ import division
import os, random
import libtbx.load_env
from libtbx.test_utils import show_diff
from libtbx import easy_run
from phenix.programs import homology

def exercise_00():
  """
  Exercise phenix.homology 
  For diff type files, get appro equal results
  """
  pdb_res = """[group_args
  chain_ref                      : A
  match                          : [group_args
  chain_id                       : A
  identity                       : 97.2602739726
  pdb_code                       : 3DS4
  resolution                     : 1.12], group_args
  chain_ref                      : B
  match                          : [group_args
  chain_id                       : A
  identity                       : 97.2602739726
  pdb_code                       : 3DS4
  resolution                     : 1.12], group_args
  chain_ref                      : C
  match                          : [group_args
  chain_id                       : A
  identity                       : 97.2602739726
  pdb_code                       : 3DS4
  resolution                     : 1.12], group_args
  chain_ref                      : D
  match                          : [group_args
  chain_id                       : A
  identity                       : 97.2602739726
  pdb_code                       : 3DS4
  resolution                     : 1.12]]"""

  seq_res="""[group_args
  chain_ref                      : A
  match                          : [group_args
  chain_id                       : A
  identity                       : 97.6744186047
  pdb_code                       : 3DS4
  resolution                     : 1.12], group_args
  chain_ref                      : B
  match                          : [group_args
  chain_id                       : A
  identity                       : 97.6744186047
  pdb_code                       : 3DS4
  resolution                     : 1.12], group_args
  chain_ref                      : C
  match                          : [group_args
  chain_id                       : A
  identity                       : 97.6744186047
  pdb_code                       : 3DS4
  resolution                     : 1.12], group_args
  chain_ref                      : D
  match                          : [group_args
  chain_id                       : A
  identity                       : 97.6744186047
  pdb_code                       : 3DS4
  resolution                     : 1.12]]"""

  pdb = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/phenix_homology/3dtj.pdb", 
    test=os.path.isfile)
  cif = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/phenix_homology/3dtj.cif", 
    test=os.path.isfile)
  seq = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/phenix_homology/3dtj.fa", 
    test=os.path.isfile)
  #
  params = homology.get_default_params()
  params.num_of_best_pdb = 1
  pdb_rs = homology.file_perfect_pair(pdb, params)
  cif_rs = homology.file_perfect_pair(cif, params)
  seq_rs = homology.file_perfect_pair(seq, params,is_input_sequence=True)
  assert not show_diff(str(pdb_rs), str(cif_rs))
  # print(pdb_rs)
  assert not show_diff(str(pdb_rs), pdb_res)
  assert not show_diff(str(seq_rs), seq_res)
  
def exercise_01():
  """
  Exercise phenix.homology
  Validate input
  """
  pdb1 = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/phenix_homology/3ds1.pdb", 
    test=os.path.isfile)
  pdb2 = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/phenix_homology/3dtj.pdb", 
    test=os.path.isfile)
  seq1 = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/phenix_homology/3ds1.fa", 
    test=os.path.isfile)
  seq2 = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/phenix_homology/3dtj.fa", 
    test=os.path.isfile)
  #
  base = "phenix.homology high_res=1.5 low_res=3.5 identity=95 %s %s"
  cmd1 = base%(pdb2, seq2)
  cmd2 = base%(pdb2, pdb1)
  for cmd in [cmd1, cmd2]:
    r = easy_run.fully_buffered(cmd)
    assert r.stderr_lines[0]== \
      'Sorry: One model or one sequence file is expected.'

def exercise_02():
  """
  Exercise phenix.homology
  """
  seq_res = """[group_args
  chain_id                       : A
  identity                       : 97.5
  pdb_code                       : 3DS4
  resolution                     : 1.12, group_args
  chain_id                       : A
  identity                       : 97.5
  pdb_code                       : 3DS2
  resolution                     : 1.2, group_args
  chain_id                       : A
  identity                       : 98.75
  pdb_code                       : 2XT1
  resolution                     : 1.32]"""
  sequence = "SPTSILDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTATLLVQNANPDCKTILKAL"+\
             "GPGATLEEMMTACQGVGGPG"
  params = homology.get_default_params()
  result = homology.perfect_pair(sequence, params)
  assert not show_diff(str(result), seq_res)

if (__name__ == "__main__"):
  exercise_00()
  exercise_01()
  exercise_02()
  print ("OK")
