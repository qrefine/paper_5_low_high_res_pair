from __future__ import division
import os
import mmtbx.model
import sys
import iotbx.pdb
from libtbx import group_args
from libtbx import easy_pickle
from libtbx.utils import Sorry
from iotbx.phil import process_command_line_with_files
from iotbx import bioinformatics
import iotbx.bioinformatics.pdb_info
from iotbx.bioinformatics import local_blast
from iotbx.bioinformatics.structure import summarize_blast_output
from iotbx.pdb.fetch import fetch
from cctbx.array_family import flex
import time

master_params_str = """\
model_name = None
  .type = path
  .multiple = False
  .help = Model file name
engine = *blastp blastall
  .type = choice(multi=False)
  .help = blast engine: blastp->blast+ blastall->legacy blast
high_res = None
  .type = float
  .help = High resolution pdb limit
low_res = 3.5
  .type = float
  .help = Low resolution pdb limit
identity = 95
  .type = int
  .help = sequence identity
piece_matching = False
  .type = bool
  .help = piece of low resolution pdb match with high_res pdb
num_of_best_pdb = 3
  .type = int
  .help = num of high_res pdb for each chain
chain_matching = True
  .type = bool
  
"""

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

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=True)

def get_inputs(args, log, master_params):
  cmdline = process_command_line_with_files(
    args         = args,
    master_phil  = master_params
  )
  params = cmdline.work.extract()
  return params

def pdb_perfect_pair(hierarchy, params):
  # Assume model is filtered to have protein etc
  pdb_id = os.path.basename(params.model_name).strip("pdb")[:4].upper()
  results = {} 
  count = 0
  for chain in hierarchy.only_model().chains():
    sequence = chain.as_padded_sequence()
    result = perfect_pair(sequence, params, pdb_id)
    if result:
      results[chain.id] = result
    else:
      count += 1
  return results, count

def sequence_perfect_pair(params):
  with open(params.model_name, "r") as f:
    fasta = f.read()
    (fastas, unknows) = bioinformatics.fasta_alignment_parse(fasta)
    for sequence in fastas.alignments:
      result = perfect_pair(sequence, params)
      print result

def perfect_pair(sequence, params, pdb_id=None):
  result=[]
  pdb_info = iotbx.bioinformatics.pdb_info.pdb_info_local()
  l_blast = local_blast.pdbaa(seq=sequence)
  blast_xml_result = l_blast.run(debug=True, binary=params.engine)
  try:
    blast_summary = summarize_blast_output("\n".join(blast_xml_result))
  except StopIteration:
    return result
  except Exception as e:
    if (("mismatched tag" in str(e))and params.engine=="blastall"):
      raise Sorry("setting engine=blastp and try again")
  pdb_ids_to_study = {}
  for hit in blast_summary:
    #hit.show(out=sys.stdout)
    hsp = hit.hsp
    ### The Gapped BLAST algorithm allows gaps (deletions and insertions) to 
    ### be introduced into the alignments
    # blastall: the surprising default value (None, None) instead of an integer.
    if hsp.gaps==(None, None): hsp.gaps=0
    identity = (hsp.identities-hsp.gaps)*100/(hsp.align_length-hsp.gaps)
    if(not params.piece_matching):
      ali_identity = len(hsp.query.replace('X',''))/len(sequence.replace('X',''))
      identity = identity * ali_identity
    if identity < params.identity:
      continue
    for i in hit.all_ids:
      if i[0] == pdb_id: continue
      if (i[0] not in pdb_ids_to_study) :
        pdb_ids_to_study[str(i[0])] = (str(i[1]),identity)
  info_lists = pdb_info.get_info_list(pdb_ids_to_study.keys())
  if info_lists:
    for info_list in info_lists:
      best_pdb_id = info_list[0]
      best_pdb_chain = pdb_ids_to_study[info_list[0]][0]
      identity = pdb_ids_to_study[info_list[0]][1]
      if params.high_res == None:
         result.append((best_pdb_id, best_pdb_chain,info_list[1],identity))
         result.sort(key=lambda tup:(tup[2],-tup[3]))
      elif params.high_res != None and info_list[1] <= params.high_res:
        result.append((best_pdb_id, best_pdb_chain,info_list[1],identity))
        result.sort(key=lambda tup:(-tup[3],tup[2]))
  result = result[:params.num_of_best_pdb]
  return result

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
  for r in hierarchy.residue_groups():
    sizes.append(r.atoms().size())
  if(sizes.size()==0): return 0
  return sizes.count(1)*100./sizes.size()

def get_hierarchy(pdb_inp):
  hierarchy = pdb_inp.construct_hierarchy()
  asc = hierarchy.atom_selection_cache()
  ss = " or ".join(aa_codes)
  ss = "(%s) and not (resname UNX or resname UNK or resname UNL)"%ss
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
      rs,c = pdb_perfect_pair(hierarchy, params)
      if (params.chain_matching and rs) or\
        (not params.chain_matching and c==0 and rs):
        print key.upper(), value[0], hierarchy.atoms().size(), value[4]
        for k,v in sorted(rs.items()):
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
  if (hierarchy is not None): 
    rs,c = pdb_perfect_pair(hierarchy, params)
    if (params.chain_matching and rs) or\
      (not params.chain_matching and c==0 and rs):
      print params.model_name, hierarchy.atoms().size()
      for key,value in sorted(rs.items()):
        print key,value
  else:
    raise Sorry("The structure not contain protein")

if __name__ == '__main__':
  t0 = time.time()
  args = sys.argv[1:]
  params = get_inputs(args=args, log=sys.stdout, master_params=master_params())
  if(params.model_name):
    if (os.path.isfile(params.model_name) \
      and (not iotbx.pdb.is_pdb_file(file_name=params.model_name)) \
      and (not iotbx.pdb.is_pdb_mmcif_file(file_name=params.model_name))):
      sequence_perfect_pair(params)
    else:
      run_one(params)
  else:
    run(params)
  print "Total time: %6.2f"%(time.time()-t0)
