directories :
	mkdir WR
	mkdir Results/Test
	mkdir Results/MC
	mkdir Results/MC/Trials

build_large :
	cd C; make monte monte=mc_sim.cpp optim=simann.cpp simulate=simulate_crs.cpp model=model_crs.cpp; mv monte monte_sim; rm monte.o simulate.o model.o
	cd C; make monte monte=mc_epfq.cpp optim=simann.cpp simulate=simulate_crs.cpp model=model_crs.cpp; mv monte monte_epfq; rm monte.o simulate.o model.o
	cd C; make monte monte=mc_mom.cpp optim=simann.cpp simulate=simulate_crs.cpp model=model_crs.cpp; mv monte monte_mom; rm monte.o simulate.o model.o
	cd C; make monte monte=mc_epfq_diag.cpp optim=simann.cpp simulate=simulate_crs.cpp model=model_crs.cpp; mv monte monte_epfq_diag; rm monte.o simulate.o model.o
	cd C; make monte monte=mc_mom_diag.cpp optim=simann.cpp simulate=simulate_crs.cpp model=model_crs.cpp; mv monte monte_mom_diag; rm monte.o simulate.o model.o

build_small :
	cd C; make monte monte=mc_sim_small.cpp optim=simann.cpp simulate=simulate_crs.cpp model=model_crs.cpp; mv monte monte_sim_small; rm monte.o simulate.o model.o
	cd C; make monte monte=mc_epfq_small.cpp optim=simann.cpp simulate=simulate_crs.cpp model=model_crs.cpp; mv monte monte_epfq_small; rm monte.o simulate.o model.o
	cd C; make monte monte=mc_mom_small.cpp optim=simann.cpp simulate=simulate_crs.cpp model=model_crs.cpp; mv monte monte_mom_small; rm monte.o simulate.o model.o
	cd C; make monte monte=mc_epfq_small_diag.cpp optim=simann.cpp simulate=simulate_crs.cpp model=model_crs.cpp; mv monte monte_epfq_small_diag; rm monte.o simulate.o model.o
	cd C; make monte monte=mc_mom_small_diag.cpp optim=simann.cpp simulate=simulate_crs.cpp model=model_crs.cpp; mv monte monte_mom_small_diag; rm monte.o simulate.o model.o

build_mis :
	cd C; make monte monte=mc_sim_mis_debtc.cpp optim=simann.cpp simulate=simulate_crs.cpp model=model_crs.cpp; mv monte monte_sim_mis; rm monte.o simulate.o model.o
	cd C; make monte monte=mc_epfq_mis.cpp optim=simann.cpp simulate=simulate_crs.cpp model=model_crs.cpp; mv monte monte_epfq_mis; rm monte.o simulate.o model.o
	cd C; make monte monte=mc_mom_mis.cpp optim=simann.cpp simulate=simulate_crs.cpp model=model_crs.cpp; mv monte monte_mom_mis; rm monte.o simulate.o model.o

build_fixed :
	cd C; make monte monte=mc_sim_fixed.cpp optim=simann.cpp simulate=simulate_crs.cpp model=model_crs.cpp; mv monte monte_sim_fixed; rm monte.o simulate.o model.o
	cd C; make monte monte=mc_epfq_fixed.cpp optim=simann.cpp simulate=simulate_crs.cpp model=model_crs.cpp; mv monte monte_epfq_fixed; rm monte.o simulate.o model.o
	cd C; make monte monte=mc_epfc_fixed.cpp optim=simann.cpp simulate=simulate_crs.cpp model=model_crs.cpp; mv monte monte_epfc_fixed; rm monte.o simulate.o model.o

build_drs :
	cd C; make monte monte=mc_sim_drs.cpp optim=simann.cpp simulate=simulate_drs.cpp model=model_drs.cpp; mv monte monte_sim_drs; rm monte.o simulate.o model.o
	cd C; make monte monte=mc_epfq_drs.cpp optim=simann.cpp simulate=simulate_drs.cpp model=model_drs.cpp; mv monte monte_epfq_drs; rm monte.o simulate.o model.o
	cd C; make monte monte=mc_mom_drs.cpp optim=simann.cpp simulate=simulate_drs.cpp model=model_drs.cpp; mv monte monte_mom_drs; rm monte.o simulate.o model.o

build_robust :
	cd C; make monte monte=mc_sim.cpp optim=simann.cpp simulate=simulate_crs.cpp model=model_crs.cpp; mv monte monte_sim; rm monte.o simulate.o model.o
	cd C; make monte monte=mc_mom_cv.cpp optim=simann.cpp simulate=simulate_crs.cpp model=model_crs.cpp; mv monte monte_mom_cv; rm monte.o simulate.o model.o
	cd C; make monte monte=mc_epfq_trans.cpp optim=simann.cpp simulate=simulate_crs.cpp model=model_crs.cpp; mv monte monte_epfq_trans; rm monte.o simulate.o model.o
	cd C; make monte monte=mc_epfc.cpp optim=simann.cpp simulate=simulate_crs.cpp model=model_crs.cpp; mv monte monte_epfc; rm monte.o simulate.o model.o

run_sim_large :
	cd C; sbatch runhybrid_sim_bridges.sh

run_sim_small :
	cd C; sbatch runhybrid_sim_small_bridges.sh

run_sim_mis :
	cd C; sbatch runhybrid_sim_mis_bridges.sh

run_sim_fixed :
	cd C; sbatch runhybrid_sim_fixed_bridges.sh

run_sim_drs :
	cd C; sbatch runhybrid_sim_drs_bridges.sh

run_mc_large :
	cd Python; python invert_weight_matrices.py data_rr_
	cd C; sbatch runhybrid_epfq_bridges.sh
	cd C; sbatch runhybrid_mom_bridges.sh
	cd C; sbatch runhybrid_epfq_diag_bridges.sh
	cd C; sbatch runhybrid_mom_diag_bridges.sh

run_mc_correct :
	cd Python; python invert_weight_matrices.py data_rr_
	cd C; sbatch runhybrid_epfq_bridges.sh
	cd C; sbatch runhybrid_mom_bridges.sh

run_mc_small :
	cd Python; python invert_weight_matrices.py data_rr_small_
	cd C; sbatch runhybrid_epfq_small_bridges.sh
	cd C; sbatch runhybrid_mom_small_bridges.sh
	cd C; sbatch runhybrid_epfq_small_diag_bridges.sh
	cd C; sbatch runhybrid_mom_small_diag_bridges.sh

run_mc_mis :
	cd Python; python invert_weight_matrices.py data_rr_mis_
	cd C; sbatch runhybrid_epfq_mis_bridges.sh
	cd C; sbatch runhybrid_mom_mis_bridges.sh

run_mc_fixed :
	cd Python; python invert_weight_matrices.py data_rr_fixed_
	cd C; sbatch runhybrid_epfq_fixed_bridges.sh
	cd C; sbatch runhybrid_epfc_fixed_bridges.sh

run_mc_drs :
	cd Python; python invert_weight_matrices.py data_rr_drs_
	cd C; sbatch runhybrid_epfq_drs_bridges.sh
	cd C; sbatch runhybrid_mom_drs_bridges.sh

run_mc_robust :
	cd Python; python invert_weight_matrices.py data_rr_
	cd C; sbatch runhybrid_epfq_trans_bridges.sh
	cd C; sbatch runhybrid_epfc_bridges.sh
	cd C; sbatch runhybrid_mom_cv_bridges.sh

table_large :
	cd Python; python generate_table_large.py

table_small :
	cd Python; python generate_table_small.py

table_robust :
	cd Python; python generate_table_robust.py

table_size :
	cd Python; python generate_table_size.py

table_drs :
	cd Python; python generate_table_drs.py

table_fixed :
	cd Python; python generate_table_fixed.py

table_mis :
	cd Python; python generate_table_mis.py
