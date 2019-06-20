from __future__ import division
import os
import sys
from libtbx import easy_pickle
import iotbx.pdb
from libtbx import group_args

def get_experimental_pdb_info(file_name,have_experimental_data=False):
  prefix = os.path.basename(file_name)[3:7].upper()
  pdb_inp = iotbx.pdb.input(file_name = file_name)
  r_free = pdb_inp.get_r_rfree_sigma().r_free
  r_work = pdb_inp.get_r_rfree_sigma().r_work
  resolution = pdb_inp.resolution()
  data_type = pdb_inp.get_experiment_type()
  space_group = None
  try:
    crystal_symmetry = pdb_inp.crystal_symmetry()
    if crystal_symmetry.space_group():
      space_group = crystal_symmetry.space_group().type().lookup_symbol()
  except iotbx.pdb.records.FormatError:
    pass
  return (prefix,resolution,r_work,r_free,\
    have_experimental_data,data_type,space_group)

if __name__ == '__main__':
  if 1:
    rdict = {}
    path = "/home/pdb/pdb/"
    dpath = "/home/pdb/structure_factors/"    
    of = open("".join([path,"INDEX"]),"r")
    files = ["".join([path,f]).strip() for f in of.readlines()]
    of.close()
    of = open("".join([dpath,"INDEX"]),"r")
    dfiles = [os.path.basename("".join([path,f]).strip())[1:5] for f in of.readlines()]
    of.close()
    for f in files:
      pdb_code = os.path.basename(f)[3:7]
      if(pdb_code in dfiles):
        have_experimental_data = True
      try:
        tup = get_experimental_pdb_info(f,have_experimental_data)  
        rdict[tup[0]] = tup[1:]
      except Exception,e:
        print pdb_code,"Failed",str(e)
    easy_pickle.dump(file_name='pdb_info.pkl', obj=rdict)
  else:
#    info=get_experimental_pdb_info("/home/pdb/pdb/du/pdb6du8.ent.gz")
    info=get_experimental_pdb_info("/home/pdb/pdb/ny/pdb4ny6.ent.gz")
    print info
   # get_experimental_pdb_info("/home/pdb/pdb/us/pdb1us0.ent.gz")
