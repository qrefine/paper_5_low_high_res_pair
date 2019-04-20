# paper_5_low_high_res_pair

1) Get filess off PDB 

phenix.fetch_pdb 3dtj --mtz

2) Corrupt free-R flags. Re-set free-R flags.

phenix.cif_as_mtz 3dtj-sf.cif 

3) Run phenix.refine to remove memory (bias) from free-R set.

phenix.refine 3dtj.{pdb,mtz} main.number_of_macro_cycles=5 --overwrite ncs_search.enabled=true optimize_xyz=true main.nproc=16

4) Fetch high-res model and data

phenix.fetch_pdb 3ds1 --mtz

5) Rfree and Rfree-Rwork are bit high, so run some refinement:

phenix.refine 3ds1.mtz 3ds1.pdb ordered_solvent=true optimize_xyz=true optimize_adp=true --overwrite

6) Superpose high-res over low-res

phenix.superpose_pdbs 3dtj_refine_001.pdb 3ds1_refine_001.pdb

7) Run finalyze:

qr.finalyse 3dtj_refine_001.pdb

8) Run refinement
a. 3dtj start Rw: 0.2614 Rf: 0.3327 Rf-Rw: 0.0713 rmsd(b):  0.0039

# 3dtj xtb refinements 

1) mode=refine quantum.nproc=2 parallel.nproc=10 max_bond_rmsd=0.02 stpmax=0.2 gradient_only=true clustering=true use_convergence_test=true opt_log=1 restraints=qm engine_name=xtb

  a. maxnum_residues_in_cluster=15 -> 3dtj_xtb_refine_MaxRes15.log (49 clusters) - 
  Rw: 0.2659 Rf: 0.3186 Rf-Rw: 0.0527 rmsd(b):  0.0113 
  
  b. maxnum_residues_in_cluster=50 -> 3dtj_xtb_refine_MaxRes50.log (39 clusters) - run on cu38
