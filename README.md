# paper_5_low_high_res_pair

1. Get filess off PDB 

phenix.fetch_pdb 3dtj --mtz

2. Corrupt free-R flags. Re-set free-R flags.

phenix.cif_as_mtz 3dtj-sf.cif 

3. Run phenix.refine to remove memory (bias) from free-R set.

phenix.refine 3dtj.{pdb,mtz} main.number_of_macro_cycles=5 --overwrite ncs_search.enabled=true optimize_xyz=true main.nproc=16

4. Fetch high-res model and data

phenix.fetch_pdb 3ds1 --mtz

5. Rfree and Rfree-Rwork are bit high, so run some refinement:

phenix.refine 3ds1.mtz 3ds1.pdb ordered_solvent=true optimize_xyz=true optimize_adp=true --overwrite

6. Superpose high-res over low-res

phenix.superpose_pdbs 3dtj_refine_001.pdb 3ds1_refine_001.pdb

7. Run finalyze:

qr.finalyse 3dtj_refine_001.pdb

8. Run refinement

a. 3dtj start Rw: 0.2614 Rf: 0.3327 Rf-Rw: 0.0713 rmsd(b):  0.0039

# 3dtj xtb refinements 

1. mode=refine quantum.nproc=2 parallel.nproc=10 max_bond_rmsd=0.02 stpmax=0.2 gradient_only=true clustering=true use_convergence_test=true opt_log=1 restraints=qm engine_name=xtb

  a. maxnum_residues_in_cluster=15 -> 3dtj_xtb_refine_MaxRes15.log (49 clusters) - 
  Rw: 0.2659 Rf: 0.3186 Rf-Rw: 0.0527 rmsd(b):  0.0113 
  
  b. maxnum_residues_in_cluster=50 -> 3dtj_xtb_refine_MaxRes50.log (39 clusters) - run on cu38
  
  # Perfect pair select
  
  `phenix.homology 3dtj.pdb high_res=2 identity=95 n=3`
  
  `phenix.homology 3dtj.cif high_res=2 identity=95 n=3`
  
  `phenix.homology 3dtj.fa high_res=2 identity=95 n=3`
  
  `phenix.homology low_res=3 high_res=2 identity=95 run_over_database =True`
  
  params: 

  high_res(default: None): resolution limit for high_res pdb. If high _res setting to None, result sorted by high_res first, then by identity. Otherwise, will sorted by identity first, then by high_res
 
  low_res(default: 3.5): resolution limit for low_res pdb

  indentity(default: 95): Sequence Similarity Cutoff

  piece_matching(default: False): piece of low resolution pdb match with high_res pdb
  
  chain_matching(defalut: True): Output pdb chain that satisfies the condition
  
  n(default: 3): num of high_res pdb for each chain
  
  run_over_database(default: False): run over pdb database find pairs for model in question
  
  engine(default: blastp): engine choice for blastp/blastall
  
  
  pdb_dict.pickle: Get info (resolution, rwork, rfree) for all X-ray diffraction PDB from RCSB and dump into pickle file
  
  pdb_info.csv: Get info(resolution, rwork, rfree, have_experimental_data，data_type, space_group) for all PDB from RCSB and     store in .csv file

# Perfect pair result

pairs_rs/

script: pair_sel.py    * log file: pair_sel_99.log     * pairs sorted by model size:   sorted_pairs.csv

select rules:

   1. low_resolution > 3
   
   2. high_res < 2.5
   
   3. identity > 0.99
   
   4. structure atoms_size < 4000
   
   5. only contain common_amio_acid
   
   6. diffraction data are available
   
   7. only models that match fully( e.g. 2WBH_ A.B.C <> 2BU1_A )

