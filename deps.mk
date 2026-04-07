PdV.o: clover.o ideal_gas.o kernels/PdV_kernel.o report.o revert.o update_halo.o
accelerate.o: clover.o kernels/accelerate_kernel.o
advec_cell_driver.o: clover.o kernels/advec_cell_kernel.o
advec_mom_driver.o: clover.o kernels/advec_mom_kernel.o
advection.o: advec_cell_driver.o advec_mom_driver.o clover.o update_halo.o
build_field.o: clover.o
calc_dt.o: clover.o kernels/calc_dt_kernel.o
clover.o: data.o definitions.o kernels/pack_kernel.o
clover_leaf.o: clover.o
definitions.o: data.o
field_summary.o: clover.o ideal_gas.o kernels/field_summary_kernel.o
flux_calc.o: clover.o kernels/flux_calc_kernel.o
generate_chunk.o: clover.o kernels/generate_chunk_kernel.o
hydro.o: PdV.o accelerate.o advection.o clover.o flux_calc.o reset_field.o timestep.o viscosity.o
ideal_gas.o: clover.o kernels/ideal_gas_kernel.o
initialise.o: clover.o parse.o report.o
initialise_chunk.o: clover.o kernels/initialise_chunk_kernel.o
kernels/update_halo_kernel.o: data.o
kernels/update_tile_halo_kernel.o: data.o
parse.o: clover.o data.o report.o
read_input.o: clover.o parse.o report.o
report.o: clover.o data.o
reset_field.o: clover.o kernels/reset_field_kernel.o
revert.o: clover.o kernels/revert_kernel.o
start.o: clover.o ideal_gas.o parse.o update_halo.o
timestep.o: calc_dt.o clover.o definitions.o ideal_gas.o report.o update_halo.o viscosity.o
update_halo.o: clover.o kernels/update_halo_kernel.o update_tile_halo.o
update_tile_halo.o: clover.o kernels/update_tile_halo_kernel.o
viscosity.o: clover.o kernels/viscosity_kernel.o
visit.o: clover.o ideal_gas.o update_halo.o viscosity.o
