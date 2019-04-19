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

