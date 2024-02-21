Preparing the folding simulations has several steps:
1. Create the initial structure
    use pepgen to create a protein from the sequence
    parameterize with grappa
    run MD
    -> run create_peptide.sh with fasta file (could take an hour)
2. Pick a collapsed equilibrium structure and set up cluster simulations
    (most of the atoms in the simulations are water; if we have a small box, we simulate much faster)
    -> use vmd or other to pick a collapsed structure (e.g. by end-to-end distance) and save as gro; end-to-end distance should be lower than box-length in paper!
    -> copy structure (only protein, no solvent) and topology without solvent with grappa parameters (should be at protein/pep_grappa_initial.top) to input directory (name pep.gro and pep.top)
    -> make sure mdps/ directory exists and has files, empty mds/ and logfiles/ directory exist
3. Submit self-submitting job
    -> sbatch J <job name> submit_continue.sh <input directory> <Task ID (just some number)> <box length (from paper)>
    I also have submit_continue_test.sh which has a shorter run time for testing
4. Check on simulation regularly
    -> check out analysis script in bba/analysis or build your own