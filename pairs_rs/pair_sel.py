from __future__ import absolute_import, division, print_function
from phenix.programs import homology
from libtbx import group_args
from libtbx import easy_pickle
import iotbx
import os

path_to_pdb = os.environ.get("PDB_MIRROR_PDB")

aa_codes = [
"resname ALA",
"resname ARG",
"resname ASN",
"resname ASP",
"resname CYS",
"resname GLN",
"resname GLU",
"resname GLY",
"resname HIS",
"resname ILE",
"resname LEU",
"resname LYS",
"resname MET",
"resname MSE",
"resname PHE",
"resname PRO",
"resname SER",
"resname THR",
"resname TRP",
"resname TYR",
"resname VAL"
]

def file_from_code(code):
  file_path = None
  if(os.path.isdir(path_to_pdb)):
    file_path = os.path.join(
      path_to_pdb, code.lower()[1:3],"pdb"+code.lower()+".ent.gz")
  if(file_path is None): return None
  assert os.path.isfile(file_path)
  return file_path

def pure_protein(hierarchy):
  asc = hierarchy.atom_selection_cache()
  ss = " or ".join(aa_codes)
  sel = asc.selection(ss)
  hierarchy = hierarchy.select(~sel)
  if(hierarchy.atoms().size()==0): return True
  return False

def check(match_pairs):
  m_ps = []
  for match_pair in match_pairs:
    m_p = []
    # print('{:6s} {:^4s}'.format("CHAIN", match_pair.chain_ref))
    for match in match_pair.match:
      m_p.append(match.pdb_code)
    #   print('{:^4s} {:^4s} {:^4.2f} {:^8.3f}'.format(match.pdb_code, \
    #     match.chain_id, match.resolution, match.identity))
    m_ps.append(m_p)
  return set.intersection(*[set(list) for list in m_ps])


def show(match_pairs):
  for match_pair in match_pairs:
    print('{:6s} {:^4s}'.format("CHAIN", match_pair.chain_ref))
    for match in match_pair.match:
      print('{:^4s} {:^4s} {:^4.2f} {:^8.3f}'.format(match.pdb_code, \
        match.chain_id, match.resolution, match.identity))

def run():
  params = homology.get_default_params()
  params.low_res = 3
  params.high_res = 2.5
  params.identity = 99
  params.chain_matching = False
  pdb_info = iotbx.bioinformatics.pdb_info.pdb_info_local().db_dict
  cntr = 0
  results = []
  for key,value in pdb_info.items():
    ## DEBUG
    # if key != "6F3A": continue
    # if key != "1M10": continue
    # if key != "3J05": continue
    ##
    if (not value[0]): continue 
    f1 = (float(value[0]) >= params.low_res)
    f3 = ("X-RAY DIFFRACTION" in value[4] or "NEUTRON DIFFRACTION" in value[4])
    f2 = (f3 and value[3]=="True") \
         or ("ELECTRON MICROSCOPY" in value[4])
    if(not (f1 and f2)):continue
    file_name = file_from_code(code=key.lower())
    if(file_name is None): continue

    pdb_inp = iotbx.pdb.input(file_name = low_file)
    model = mmtbx.model.manager(model_input = pdb_inp)
    hierarchy = model.get_hierarchy()
    atom_size = hierarchy.atoms().size()
    
    if atom_size > 4000: 
      # print("size problem")
      continue
    if pure_protein(hierarchy) is False:
      # print("not pure protein") 
      continue 
    try:
      rs = homology.file_perfect_pair(file_name, params)
      # print(rs)
      if (rs is None): continue
      model_rs = check(rs)
      if (len(model_rs) is 0): continue
      print(key.upper(), value[0], value[4])
      print(model_rs)
      show(rs)
      print ("*"*80)
      result = group_args(
        pdb_id     = key.upper(),
        resolution = value[0],
        method     = value[4],
        model_rs   = model_rs,
        rs         = rs)
      results.append(result)
    except KeyboardInterrupt: raise
    # except Exception
    except Exception, e:
      print ("FAILED:", file_name, str(e))
    cntr += 1
  print ("Processed:", cntr, "out of total:", len(pdb_info.keys()))
  easy_pickle.dump("results.pkl",results)



if __name__ == "__main__":
    run()
