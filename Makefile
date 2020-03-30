# !!!IMPORTANT!!!
# make should be run in the appropriate python (conda) environment!
# CONDA_ACTIVATE=conda activate clim
# CONDA_DEACTIVATE=conda deactivate
EXEC_NB=python execute_notebook.py

CODE_DIR=code
FIG_DIR=plots
TAB_DIR=tables
SADIR=../modelling/um/results/sa/

PROC_DATA_MAIN=\
    $(SADIR)/trap1e_grcs/_processed/trap1e_grcs.nc \
    $(SADIR)/trap1e_llcs_all_rain/_processed/trap1e_llcs_all_rain.nc \
    $(SADIR)/trap1e_acoff/_processed/trap1e_acoff.nc \
    $(SADIR)/proxb_grcs/_processed/proxb_grcs.nc \
    $(SADIR)/proxb_llcs_all_rain/_processed/proxb_llcs_all_rain.nc \
    $(SADIR)/proxb_acoff/_processed/proxb_acoff.nc \

# $(CODE_DIR)/pp_main_gl_output.py -p trap1e -r grcs -s 1440
# $(CODE_DIR)/proc_global_um_output_ens.py -p proxb -s 720 -r NONE
# $(CODE_DIR)/pp_ns_data.py -p trap1e -r grcs
# $(CODE_DIR)/aggr_ens_output.py
# $(CODE_DIR)/aggr_ns_data.py
# $(CODE_DIR)/ns_mean_vprof.py

FIGURES=\
    $(FIG_DIR)/trap1e_grcs_110d_view3d__enlarged.png \
    $(FIG_DIR)/trap1e_proxb__grcs_llcs_all_rain_acoff__tsfc_winds__plev250hpa.pdf \
    $(FIG_DIR)/trap1e_proxb__grcs_llcs_all_rain_acoff__vprof_daynight_pm01deg__temp_shum.pdf \
    $(FIG_DIR)/trap1e_proxb__grcs_llcs_all_rain_acoff__ghgt_winds_eddy__250hpa__temp__700hpa.pdf \
    $(FIG_DIR)/trap1e_proxb__grcs_llcs_all_rain_acoff__nondim_rossby_rhines__hgt0-15km.pdf \
    $(FIG_DIR)/trap1e_proxb__grcs_llcs_all_rain_acoff__cloud_types__w_precip.pdf \
    $(FIG_DIR)/trap1e_proxb__grcs_llcs_all_rain_acoff__mse_div_merid_mean.pdf \
    $(FIG_DIR)/trap1e_proxb__grcs_llcs_all_rain_acoff__vcross_lon90E__hgt0-16__zonal_fluxes.pdf \
    $(FIG_DIR)/trap1e__grcs__toa_olr_map__102d.pdf \
    $(FIG_DIR)/proxb__grcs__toa_olr_map__100d.pdf \
    $(FIG_DIR)/trap1e_proxb__grcs__toa_olr_hist.pdf \
    $(FIG_DIR)/proxb__grcs__precip_w_map.pdf \
    $(FIG_DIR)/trap1e_proxb__grcs__vprof_qcl_qcf.pdf \
    $(FIG_DIR)/trap1e_proxb__grcs__vprof_t_incr_lh.pdf \
    $(FIG_DIR)/trap1e_proxb__grcs__scatter_w_linreg__hdiv_5_20km_ss__t_sfc_diff_dn.pdf
TABLES=\
    $(TAB_DIR)/mean_sfc_temp_table.tex \
    $(TAB_DIR)/cre_table.tex \

$(FIG_DIR)/trap1e_grcs_110d_view3d__enlarged.png: $(CODE_DIR)/Fig01-Model-Set-Up-3D.ipynb
$(FIG_DIR)/trap1e_proxb__grcs_llcs_all_rain_acoff__tsfc_winds__plev250hpa.pdf: $(CODE_DIR)/Fig02-Map-Sfc-Temp-Winds.ipynb
$(FIG_DIR)/trap1e_proxb__grcs_llcs_all_rain_acoff__vprof_daynight_pm01deg__temp_shum: $(CODE_DIR)/Fig03-Vert-Prof-Temp-Hum.ipynb
$(FIG_DIR)/trap1e_proxb__grcs_llcs_all_rain_acoff__ghgt_winds_eddy__250hpa__temp__700hpa.pdf: $(CODE_DIR)/Fig04-Map-Gpt-Hgt-Temp.ipynb
$(FIG_DIR)/trap1e_proxb__grcs_llcs_all_rain_acoff__nondim_rossby_rhines__hgt0-15km.pdf: $(CODE_DIR)/Fig05-Rossby-Rhines.ipynb
$(FIG_DIR)/trap1e_proxb__grcs_llcs_all_rain_acoff__cloud_types__w_precip.pdf: $(CODE_DIR)/Fig06-Map-Cloud-Precip.ipynb
$(FIG_DIR)/trap1e_proxb__grcs_llcs_all_rain_acoff__mse_div_merid_mean.pdf: $(CODE_DIR)/Fig07-MSE-Flux-Divergence.ipynb
$(FIG_DIR)/trap1e_proxb__grcs_llcs_all_rain_acoff__vcross_lon90E__hgt0-16__zonal_fluxes.pdf: $(CODE_DIR)/Fig08-Vert-Cross-Heat-Flux.ipynb
$(FIG_DIR)/trap1e__grcs__toa_olr_map__102d.pdf: $(CODE_DIR)/Fig09-Map-TOA-OLR-Trap1e.ipynb
$(FIG_DIR)/proxb__grcs__toa_olr_map__100d.pdf: $(CODE_DIR)/Fig10-Map-TOA-OLR-Proxb.ipynb
$(FIG_DIR)/trap1e_proxb__grcs__toa_olr_hist.pdf: $(CODE_DIR)/Fig11-Hist-TOA-OLR.ipynb
$(FIG_DIR)/proxb__grcs__precip_w_map.pdf: $(CODE_DIR)/Fig12-Map-Precip-Vert-Wind.ipynb
$(FIG_DIR)/trap1e_proxb__grcs__vprof_qcl_qcf.pdf: $(CODE_DIR)/Fig13-Vert-Prof-Cloud-Condensate.ipynb
$(FIG_DIR)/trap1e_proxb__grcs__vprof_t_incr_lh.pdf: $(CODE_DIR)/Fig14-Vert-Prof-Latent-Heating.ipynb
$(FIG_DIR)/trap1e_proxb__grcs__scatter_w_linreg__hdiv_5_20km_ss__t_sfc_diff_dn.pdf: $(CODE_DIR)/Fig15-Day-Night-Impact-Lin-Reg.ipynb

$(TAB_DIR)/mean_sfc_temp_table.tex: $(CODE_DIR)/Tab03-Mean-Sfc-Temp.ipynb
$(TAB_DIR)/cre_table.tex: $(CODE_DIR)/Tab04-CRE.ipynb

all: $(FIGURES)

$(FIGURES):
	@echo "making figure: $@"
	$(EXEC_NB) $<


.PHONY: clean help
clean:
	@echo "Cleaning..."
	rm -f $(FIG_DIR)/*

help:
	@echo ""
	@echo "Usage:"
	@echo "    make all: run Jupyter Notebooks to create all figures"
	@echo "    make clean: Delete files in plots/ folder"
	@echo "    make help: Print this message and exit"
	@echo ""
