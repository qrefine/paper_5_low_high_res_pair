from __future__ import absolute_import, division, print_function
from libtbx import easy_run
import os,sys,libtbx.load_env

'''
This is a tool to run BLAST locally against selected databased such
as PDBaa ...etc.  The executables and databases are all distributed
with Phenix so no extra installation is required. PDBaa has been
implemted.  More classes, e.g. ScopE, PDBstructure ... will be
added later.
LWH 4/12/19

Useage:
>>>from iotbx.bioinformatics import local_blast
>>>myxml=local_blast.pdbaa(seq=myseq).run()

where
myseq is the query protein sequence string. You can use 'X' to fill gaps.
myxml is the stdout_lines object of the blast XML output.
'''
#setup lib dir
phenixpath=os.getenv('PHENIX')
ligand_lib_dir = libtbx.env.find_in_repositories(
  relative_path=os.path.join("chem_data","ligand_lib"),
  test=os.path.isdir)
cwd=os.getcwd()

#make sure blast exists.  Never had a problem but probably won't hurt.
def checkblast():
  phenix_blast_exe=''
  systype=sys.platform #linux2,darwin,win32
  sysname='Linux'
  phenix_blast_exe='blastp_%s'%(systype)
  if systype=='win32':
    phenix_blast_exe='blastp_%s.exe'%(systype)
    sysname='Windows'
  elif systype=='darwin':
    phenix_blast_exe='blastp_%s'%(systype)
    sysname='OSX'
  else:
    pass
  phenix_blast=os.path.join(ligand_lib_dir, phenix_blast_exe)
  blastexe=None
  if os.path.exists(phenix_blast):
    blastpath=phenix_blast
    blastexe=phenix_blast_exe
    #print('%s version is running...\n'%sysname)
  else:
    print('BLAST executable does not exist. please check your Phenix installation.')
    sys.exit(0)
  return blastpath



class pdbaa(object):
  def __init__(self, workdir=None, prefix=None, seq=None, output=None):
    self.workdir=workdir
    self.prefix=prefix
    self.seq=seq
    self.output=output

  def run(self, debug=False):
    blastpath=checkblast()
    curdir=os.getcwd()
    if self.workdir is None:
      self.workdir=curdir
    elif os.path.exists(self.workdir) is False:
      print("Input working directory does not exist. Try to work in current directory instead.")
      self.workdir=curdir
    fasta_file='myprotein.fasta'
    fasta_path=os.path.join(self.workdir,fasta_file)
    if self.prefix is None:
       self.prefix="myprotein"
    fastaline=">%s\n%s\n"%(self.prefix,self.seq)
    f=open(fasta_path,'w').writelines(fastaline)
    dbname="pdbaa.00"
    outfmt="-outfmt 5" #xml_out
    blastdb=os.path.join(ligand_lib_dir,dbname)

#    blastrun_seq=" -p blastp -i %s -a 8 -F F -W 3 -G 11 -E 2 \
#        -V F -e 1E-3 %s -d %s"%(fasta_path,outfmt, blastdb)
  
    blastrun_seq= "-query %s %s -evalue 1E-3 -num_threads 8 -db %s"%(fasta_path, outfmt, blastdb)

    cmds="%s %s"%(blastpath,blastrun_seq)
    #print(cmds)
    try:
      result = easy_run.fully_buffered(
        command=cmds,
        stdin_lines='')
    except KeyboardInterrupt :
      raise KeyboardInterrupt
    else :
      ## output contain wrong tag e.g. 6f3a
      ## should remove later
    #  if result.stdout_lines[-2] != "  </BlastOutput_iterations>":
    #    result.stdout_lines[-2] = "  </BlastOutput_iterations>"
      ##
      if debug:
        output='myprotein.xml'
        open(output, "w").write("\n".join(result.stdout_lines))
      return result.stdout_lines
