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

def get_perfect_pair(model, params, chain_ids=None):
  pdb_info = iotbx.bioinformatics.pdb_info.pdb_info_local()
  pdb_id = os.path.basename(params.model_name).strip("pdb")[:4].upper()
  h = model.get_hierarchy()
  results =  []
  for chain in h.only_model().chains():
    if not chain.is_protein():
      continue
    if chain_ids is not None and chain.id not in chain_ids():
      continue
    sequence = chain.as_padded_sequence()
    l_blast = local_blast.pdbaa(seq=sequence) 
    blast_xml_result = l_blast.run()
    blast_summary = summarize_blast_output("\n".join(blast_xml_result))
    pdb_ids_to_study = {}
    for hit in blast_summary:
      #hit.show(out=sys.stdout)
      if hit.identity < params.identity:
        continue
      pdb_ids_to_study[hit.pdb_id] = hit.chain_id 
    #
    info_list = pdb_info.get_info_list(pdb_ids_to_study.keys())
    info_list.sort(key=lambda tup: tup[1])
    if info_list :
      best_pdb_id = info_list[0][0]
      best_pdb_chain = pdb_ids_to_study[info_list[0][0]]
    #  print "Best pdb:", info_list[0], "chain:", pdb_ids_to_study[info_list[0][0]]
      if info_list[0][1] < params.high_res:
        result = group_args(
          low_res_pdb  = pdb_id+'_'+chain.id ,
          high_res_pdb = best_pdb_id +'_'+ best_pdb_chain)
        results.append(result) 
  return results

def run(params):
  work_dir = "/home/pdb/pdb/" 
  pdb_info = easy_pickle.load(file_name='pdb_dict.pickle')
  for key,value in pdb_info.items():
    if value[0] > params.low_res:
      file_name = "pdb"+key.lower()+".ent.gz"
      file_path = os.path.join(work_dir,key.lower()[1:3],file_name)
      if os.path.isfile(file_path):
        try:
          params.model_name = file_path
          model = mmtbx.model.manager(model_input=iotbx.pdb.input(file_path))
          rs = get_perfect_pair(model,params)
          if rs:
            for r in rs:
              print "%-10s  %-10s"%(r.low_res_pdb,r.high_res_pdb)
        except Exception, e:
          print "%s: %s"%(key, e)        

def run_one(params):
  if not os.path.isfile(params.model_name):
    data = fetch(id=params.model_name[:4])
    model = mmtbx.model.manager(
      model_input=iotbx.pdb.input(source_info=None, lines=data.readlines()))
  else:
    model = mmtbx.model.manager(model_input=iotbx.pdb.input(params.model_name))
  rs = get_perfect_pair(model, params)
  if rs:
    for r in rs:
      print "%-10s  %-10s"%(r.low_res_pdb,r.high_res_pdb)

if __name__ == '__main__':
  args = sys.argv[1:]
  params = get_inputs(args=args, log=sys.stdout, master_params=master_params())
  if  params.model_name:
    run_one(params)
  else:
    run(params)
