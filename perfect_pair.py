from __future__ import division
import os
import mmtbx.model
import sys
import iotbx.pdb
from libtbx import group_args
from libtbx import easy_pickle
from iotbx.phil import process_command_line_with_files
import iotbx.bioinformatics.pdb_info
from iotbx.bioinformatics import local_blast
from iotbx.bioinformatics.structure import summarize_blast_output
from iotbx.pdb.fetch import fetch
from cctbx.array_family import flex

master_params_str = """\
model_name = None
  .type = path
  .multiple = False
  .help = Model file name
high_res = 1.5
  .type = float
  .help = High resolution pdb limit
low_res = 3.5
  .type = float
  .help = Low resolution pdb limit
identity = 100
  .type = int
  .help = sequence identity
match_for_piece = False
  .type = bool
  .help = piece of low resolution pdb match with high_res pdb
num_of_best_pdb = 3
  .type = int
  .help = num of high_res pdb for each chain
"""

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=True)

def get_inputs(args, log, master_params):
  cmdline = process_command_line_with_files(
    args         = args,
    master_phil  = master_params
  )
  params = cmdline.work.extract()
  return params

def get_perfect_pair(hierarchy, params):
  # Assume model is filtered to have protein etc
  pdb_info = iotbx.bioinformatics.pdb_info.pdb_info_local()
  pdb_id = os.path.basename(params.model_name).strip("pdb")[:4].upper()
  h = hierarchy
  results = {} 
  count = 0
  for chain in h.only_model().chains():
    sequence = chain.as_padded_sequence()
    l_blast = local_blast.pdbaa(seq=sequence)
    blast_xml_result = l_blast.run()
    blast_summary = summarize_blast_output("\n".join(blast_xml_result))
    pdb_ids_to_study = {}
    result = []
    for hit in blast_summary:
      #hit.show(out=sys.stdout)
      if hit.identity < params.identity:
        continue
      if(not params.match_for_piece):
        ali_identity = 100*hit.length/chain.residue_groups_size()
        if ali_identity < 95:
          continue
      pdb_ids_to_study[hit.pdb_id] = (hit.chain_id,hit.identity,ali_identity)
      for i in hit.all_ids:
        pdb_ids_to_study[str(i)] = (hit.chain_id,hit.identity,ali_identity)
    #
    info_lists = pdb_info.get_info_list(pdb_ids_to_study.keys())
    info_lists.sort(key=lambda tup: tup[1])
    if info_lists:
      for info_list in info_lists:
        best_pdb_id = info_list[0]
        best_pdb_chain = pdb_ids_to_study[info_list[0]][0]
        identity = pdb_ids_to_study[info_list[0]][1]
        if info_list[1] < params.high_res:
          result.append((best_pdb_id, best_pdb_chain,info_list[1],identity))
          result.sort(key=lambda tup: tup[3],reverse=True)
      if result:
        results[chain.id] = result[:params.num_of_best_pdb]
      else:
        count += 1
    else:
      count += 1
  return results,count

def file_from_code(code):
  work_dir_1 = "/home/pdb/pdb/"
  work_dir_2 = "/net/cci/pdb_mirror/mmcif/"
  file_path = None
  if(os.path.isdir(work_dir_1)):
    file_path = os.path.join(
      work_dir_1,code.lower()[1:3],"pdb"+code.lower()+".ent.gz")
  elif(work_dir_2):
    with open("".join([work_dir_2,"INDEX"]),"r") as fo:
      for l in fo.readlines():
        l = l.strip()
        if(l.count(code)>0):
          file_path = "".join([work_dir_2,l])
          break
  if(file_path is None): return None
  assert os.path.isfile(file_path)
  return file_path

def percent_of_single_atom_residues(hierarchy):
  sizes = flex.int()
  h = hierarchy
  for r in h.residue_groups():
    sizes.append(r.atoms().size())
  if(sizes.size()==0): return 0
  return sizes.count(1)*100./sizes.size()

def get_hierarchy(pdb_inp):
  hierarchy = pdb_inp.construct_hierarchy()
  asc = hierarchy.atom_selection_cache()
  ss = "protein and not (resname UNX or resname UNK or resname UNL)"
  sel = asc.selection(ss)
  hierarchy = hierarchy.select(sel)
  if(percent_of_single_atom_residues(hierarchy)>20):
    return None
  if(hierarchy.atoms().size()==0):
    return None
  if(len(hierarchy.models())>1):
    return None
  return hierarchy

def run(params):
  pdb_info = easy_pickle.load(file_name='pdb_info.pkl')
  cntr = 0
  results = []
  for key,value in pdb_info.items():
    ### DEBUG
    #if key != "6MYY": continue
    #if(not "ELECTRON MICROSCOPY" in value[4]): continue
    ### DEBUG
    fl = value[0] > params.low_res 
    f2 = ("X-RAY DIFFRACTION" in value[4] and value[3]==True) \
         or ("ELECTRON MICROSCOPY" in value[4])
    if(not (fl and f2)): continue
    file_name = file_from_code(code=key.lower())
    if(file_name is None): continue
    try:
      pdb_inp = iotbx.pdb.input(file_name = file_name)
      hierarchy = get_hierarchy(pdb_inp = pdb_inp)
      if(hierarchy is None): continue
      params.model_name = file_name
      rs,c = get_perfect_pair(hierarchy, params)
      if c==0 and rs:
        print key.upper(), value[0], hierarchy.atoms().size(), value[4]
        for k,v in rs.items():
          print k, v
        print "*"*80
        result = group_args(
          pdb_id     = key.upper(),
          resolution = value[0],
          size       = hierarchy.atoms().size(),
          method     = value[4],
          rs         = rs)
    except KeyboardInterrupt: raise
    except Exception, e:
      print "FAILED:", file_name, str(e)
    #
    cntr += 1
  results.append(result)
  print "Processed:", cntr, "out of total:", len(pdb_info.keys())
  easy_pickle.dump("pp_%d_%d_%d.pkl"\
    %(params.low_res*10,params.high_res*10,params.identity),results)

def run_one(params):
  if not os.path.isfile(params.model_name):
    data = fetch(id=params.model_name[:4])
    pdb_inp = iotbx.pdb.input(source_info=None, lines=data.readlines())
  else:
    pdb_inp = iotbx.pdb.input(params.model_name)
  hierarchy = get_hierarchy(pdb_inp = pdb_inp)
  assert (hierarchy is not None)
  rs,c = get_perfect_pair(hierarchy, params)
  if c==0 and rs:
    print params.model_name,model.size()
    for key,value in rs.items():
      print key,value

if __name__ == '__main__':
  args = sys.argv[1:]
  params = get_inputs(args=args, log=sys.stdout, master_params=master_params())
  if(params.model_name):
    run_one(params)
  else:
    run(params)
