from __future__ import absolute_import, division, print_function
from phenix.program_template import ProgramTemplate
from libtbx import group_args
import iotbx.pdb
import iotbx.pdb.fetch
import os
import sys
import iotbx.phil
from libtbx import easy_pickle
from libtbx.utils import Sorry
from iotbx import bioinformatics
import iotbx.bioinformatics.pdb_info
from iotbx.bioinformatics import local_blast
from iotbx.bioinformatics.structure import summarize_blast_output
from cctbx.array_family import flex

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

master_phil_str = """
run_over_database = False
  .type = bool
  .help = run over PDB dataset 
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

def get_default_params():
  return iotbx.phil.parse(
    input_string = master_phil_str,
    process_includes = True
  ).extract()

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

def percent_of_single_atom_residues(hierarchy):
  sizes = flex.int()
  for r in hierarchy.residue_groups():
    sizes.append(r.atoms().size())
  if(sizes.size()==0): return 0
  return sizes.count(1)*100./sizes.size()

def perfect_pair(sequence, params, pdb_id=None):
  result=[]
  pdb_info = iotbx.bioinformatics.pdb_info.pdb_info_local()
  l_blast = local_blast.pdbaa(seq=sequence)
  blast_xml_result = l_blast.run(debug=False, binary=params.engine)
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
      ali_identity=len(hsp.query.replace('X',''))/len(sequence.replace('X',''))
      identity = identity * ali_identity
    if identity < params.identity:
      continue
    for i in hit.all_ids:
      if i[0] == pdb_id: continue
      if (i[0] not in pdb_ids_to_study) :
        pdb_ids_to_study[str(i[0])] = (str(i[1]),identity)
  info_lists = pdb_info.get_info_list(pdb_ids_to_study.keys())
  info_lists.sort(key=lambda tup: tup[1])
  if info_lists:
    for info_list in info_lists:
      if not info_list[1]: continue
      best_pdb_id = info_list[0]
      best_pdb_chain = pdb_ids_to_study[info_list[0]][0]
      identity = pdb_ids_to_study[info_list[0]][1]
      resolution = float(info_list[1])
      g_a = group_args(
        pdb_code   = best_pdb_id,
        chain_id   = best_pdb_chain,
        resolution = resolution,
        identity   = identity)
      if params.high_res is None:
        result.append(g_a)
        result.sort(key=lambda tup:(tup.resolution,-tup.identity))
      elif resolution <= params.high_res:
        result.append(g_a)
        result.sort(key=lambda tup:(-tup.identity,tup.resolution))
  result = result[:params.num_of_best_pdb]
  return result

def file_perfect_pair(model_name, params, is_input_sequence=False):
  pdb_id = os.path.basename(model_name)[:4].upper()
  if (not iotbx.pdb.fetch.looks_like_pdb_id(pdb_id)): pdb_id = None
  seqs_dic = {}  
  if (is_input_sequence):
    with open(model_name, "r") as f:
      fasta = f.read()
      (seqs, unknows) = bioinformatics.seq_sequence_parse(fasta)
      for seq in seqs:
        sequence = seq.sequence
        chain_id = seq.name.split('|')[0][-1]
        seqs_dic[chain_id]=sequence
  else:
    pdb_inp = iotbx.pdb.input(model_name)
    hierarchy = get_hierarchy(pdb_inp = pdb_inp)
    if (hierarchy is not None):
      # Assume model is filtered to have protein etc
      for chain in hierarchy.only_model().chains():
        sequence = chain.as_padded_sequence()
        seqs_dic[chain.id] = sequence
    else: raise Sorry("The structure not contain protein")
  
  results = []
  count = 0
  for chain_id,seq in sorted(seqs_dic.items()):
    result = perfect_pair(seq, params, pdb_id)
    if result:
      results.append(group_args(
        chain_ref = chain_id,
        match     = result
      ))
    else: count += 1
  if (params.chain_matching and results) or\
      (not params.chain_matching and count==0 and results):
    return results
  else: return None
  

# =============================================================================

class Program(ProgramTemplate):

  description = """
Given low-resolution protein atomic model or corresponding sequence can search 
the PDB for a set of top highest-resolution models within specified sequence 
identify threshold.

  phenix.homology model.pdb  
  phenix.homology model.cif 
  phenix.homology fasta.fa
  phenix.homology run_over_database=True
"""

  datatypes = ['model', 'phil', 'sequence']

  master_phil_str = master_phil_str
  model_name = None
  is_input_sequence = False

  def validate(self):
    print('Validate inputs:', file=self.logger)    
    if self.params.run_over_database:
      self.path_to_pdb = os.environ.get("PDB_MIRROR_PDB")
      if(self.path_to_pdb is None):
        raise Sorry("Environmental variable PDB_MIRROR_PDB must be set.")
    else:
      model_names    = self.data_manager.get_model_names()
      sequence_names = self.data_manager.get_sequence_names()
      if((len(model_names) + len(sequence_names)) != 1):
         raise Sorry("One model or one sequence file is expected.")
      elif(len(model_names)>0):
        self.model_name = model_names[0]
        if(not (iotbx.pdb.is_pdb_file(      file_name=self.model_name) or
                iotbx.pdb.is_pdb_mmcif_file(file_name=self.model_name))):
          raise Sorry("Model file must be in PDB or mmCIF format.")
      else:
        self.model_name = sequence_names[0]
        self.is_input_sequence = True
  
  def _file_from_code(self, code):
    file_path = None
    if(os.path.isdir(self.path_to_pdb)):
      file_path = os.path.join(
        self.path_to_pdb, code.lower()[1:3],"pdb"+code.lower()+".ent.gz")
    if(file_path is None): return None
    assert os.path.isfile(file_path)
    return file_path

  def show(self, match_pairs):
    for match_pair in match_pairs:
      print("\nRef Chain ID: ", match_pair.chain_ref, file=self.logger)
      for match in match_pair.match:
        print('{:^4s} {:^4s} {:^4.2f} {:^8.3f}'.format(match.pdb_code, \
          match.chain_id, match.resolution, match.identity), file=self.logger)

  def run(self):
    if (self.params.run_over_database):
      pdb_info = iotbx.bioinformatics.pdb_info.pdb_info_local().db_dict
      cntr = 0
      results = []
      ## value read as string
      for key,value in pdb_info.items():
        ## DEBUG
        # if key != "6F3A": continue
        ##
        if (not value[0]): continue
        f1 = (float(value[0]) >= self.params.low_res)
        f2 = ("X-RAY DIFFRACTION" in value[4] and value[3]=="True") \
             or ("ELECTRON MICROSCOPY" in value[4])
        if(not (f1 and f2)): continue
        file_name = self._file_from_code(code=key.lower())
        if(file_name is None): continue
        try:
          rs = file_perfect_pair(file_name, self.params)
          if (rs is None): continue 
          print ("*"*80, file=self.logger)
          result = group_args(
            pdb_id     = key.upper(),
            resolution = value[0],
            method     = value[4],
            rs         = rs)
          # self.show(rs)
          results.append(result)
        except KeyboardInterrupt: raise
        except Exception, e:
          print ("FAILED:", file_name, str(e), file=self.logger)
        cntr += 1
      print ("Processed:", cntr, "out of total:", len(pdb_info.keys()),
        file=self.logger)
      easy_pickle.dump("results.pkl",results)

    else:
      rs = file_perfect_pair(self.model_name,self.params,self.is_input_sequence) 
      if rs: self.show(rs)
      else:
        print("\nNot find matching highest-resolution model!", file=self.logger)

    
