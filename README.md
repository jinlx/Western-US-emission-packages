# WECAN_emission_packages

### Goal of this repo
This is using to reproduce the figures in the BB emission paper below. In this paper, we used **Unix shell, Python, IDL, and R** for data analysis and processing. The main script will be provided here.

- **Lixu Jin**, Wade Permar, Vanessa Selimovic, Damien Ketcherside, Robert J. Yokelson, Rebecca S. Hornbrook, Eric C. Apel, I-Ting Ku, Jeffrey L. Collett, Amy P. Sullivan, Daniel A. Jaffe, Jeffrey R. Pierce, Alan Fried, Matthew M. Coggon, Georgios I. Gkatzelis, Carsten Warneke, Emily V. Fischer, and Lu Hu (submitted) <br>
[Constraining emissions of volatile organic compounds from western US wildfires with WE-CAN and FIREX-AQ airborne observations](TBD for the web). <br>
  *Atmospheric Chemistry and Physics*.

***Citation is encouraged if your work is related to my BB emission paper! Also, please do not hestiate to reach out if there is anything I could help!***

### Description of each figure
- [Figure 1](https://github.com/jinlx/Western-US-emission-packages/blob/main/Data_Analysis/Fig1_Budget_Map_Bpch_Based.pro)
This script is used for plotting spatial distributions of total VOCs from different sources, used as figure 1 in BB emission paper

- [Figures 2 and S3](https://github.com/jinlx/Western-US-emission-packages/blob/main/Data_Analysis/Fig2_Emission_Flux_Ers_Comp.ipynb)
This script is used for analyzing emission fluxes data and emission-flux-based regional ER. The emission flux data is based on GCv12.5.0 bpch files and it's currently obsolete and future anlaysis is recommended to use hmco netcdf output. I recommend to save out data in a csv and analyze data by using panda in the future rather than hand-writing here. This script seems a little bit messy as I included anlayze codes here. Still, the **bootstrapping_CI** and **plotting tech** here are still useful.  Just run through it!

- [Figures 3-5 and S4-S5](https://github.com/jinlx/Western-US-emission-packages/blob/main/Data_Analysis/Fig345_comp_mod_obs_vpro_wecan_clean.pro)
This script is used to plot vertical profile of WE-CAN observation and corresponding GEOS-Chem outputs for Figs 3-5. 
1) To reproduce Fig.3, you need to compile the script and use the command below. <br />
comp_mod_obs_vpro_wecan_clean,plume=0,/ind_co,/each,/nested,/save,/test,/prs,/med,/wus,/whole,/obs_major,/codeployed,/errorbar
2) To reproduce Fig.4, you need to compile the script and use the command below. <br />
comp_mod_obs_vpro_wecan_clean,plume=0,/all,/nested,/save,/test,/prs,/med,/wus,/whole,/obs_major,/codeployed,/errorbar
3) To reproduce Fig.5, you need to compile the script and use the command below. <br />
comp_mod_obs_vpro_wecan_clean,plume=1,/filter2_nobb,/filter3_nobb,/all,/nested,/save,/test,/prs,/med,/wus,/whole,/obs_major,/codeployed,/errorbar
4) To reproduce Fig. S4, you need to compile the script and use the command below. <br />
comp_mod_obs_vpro_wecan_clean,plume=0,/all,/nested,/save,/test,/prs,/med,/wus,/whole,/obs_major,/codeployed,/errorbar
4) To reproduce Fig. S5, you need to compile the script and use the command below. <br />
comp_mod_obs_vpro_wecan_clean,plume=1,/filter2_nobb,/filter3_nobb,/ind_co,/each,/nested,/save,/test,/prs,/med,/wus,/whole,/obs_major,/codeployed,/errorbar


- [Figure 6](https://github.com/jinlx/Western-US-emission-packages/blob/main/Data_Analysis/Fig6_ERs_wecan_not_clean_yet_plots.pro)
This script is used to plot scatter plot of WE-CAN observation and corresponding GEOS-Chem outputs for Fig. 6 and provide stat data (i.e., ER/slope and R) for Fig. 7.
To reproduce Fig. 6, you need to compile the script and use the command below. <br />
ERs_wecan_not_clean_yet_plots,plume=0,/emipass,/all,/test,/save,inventories='gfas',/nested,/xeqco

- [Figures 7 and S7](https://github.com/jinlx/Western-US-emission-packages/blob/main/Data_Analysis/Fig7_ER_compilation_campaigns.ipynb)
This script is used for plotting barplot of VOC ER, used as figs. 7 and S7 in WE-CAN emission paper. Just run through it!

- [Figure 8](https://github.com/jinlx/Western-US-emission-packages/blob/main/Data_Analysis/Fig8_grd_ts.ipynb)
This script is used for plotting CO time series from modelling output and EPA ground measurement. Just run through it!

- [Figures 9 and S9](https://github.com/jinlx/Western-US-emission-packages/blob/main/Data_Analysis/Fig9_comp_mod_obs_vpro_firexaq_clean.pro)
This script is used to plot vertical profile of FIREX-AQ observation and corresponding GEOS-Chem outputs. To reproduce figs, you need to compile the script and use the command below. <br />
comp_mod_obs_vpro_firexaq_clean,plume=0,/all,/nested,/save,/test,/prs,/med,/wus,/whole,/obs_major,/codeployed,/errorbar

- [Figuer S1a](https://github.com/jinlx/Western-US-emission-packages/blob/main/Data_Analysis/FigS1a_EnRs_comp_multi_instru_wecan_clean.pro)
This script is used to plot scatter plot for instrument comparison in the WE-CAN for Fig S1a. <br />
You need to compile the script and use the command below. <br />
EnRs_comp_multi_instru_wecan_clean,plume=0,/emipass,/all,/test,/save,inventories='gfas',/nested,/xeqco

- [Figuer S1b](https://github.com/jinlx/Western-US-emission-packages/blob/main/Data_Analysis/FigS1b_EnRs_comp_multi_instru_firexaq_clean.pro)
This script is used to plot scatter plot for instrument comparison in the FIREX-AQ for Fig S1b. <br />
You need to compile the script and use the command below. <br />
EnRs_comp_multi_instru_firexaq_clean,plume=0,/emipass,/all,/test,/save,inventories='gfas',/nested,/xeqco

- [Figuer S2](https://github.com/jinlx/Western-US-emission-packages/blob/main/Data_Analysis/FigS2_pdf_comp_wecan_firexaq.pro)
This script is used to do the distribution plot  
scatter plot for instrument comparison in the FIREX-AQ for Fig S1b. You need to compile the script and use the command below. <br />
pdf_comp_wecan_firexaq,/save

- [Figure S6](https://github.com/jinlx/Western-US-emission-packages/blob/main/Data_Analysis/FigS6_comp_mod_obs_vpro_wecan_celan_inj.pro)
This script is used to do the vertical profile of WE-CAN observations and corresponding injection height sensitivity experiments. <br />
To reproduce the figure, you need to compile the script and use the command below. <br />
comp_mod_obs_vpro_wecan_celan_inj,plume=0,/all,/test,/save,/prs,/med,/errorbar,/wus,/whole,/codeployed,/obs_major

- [Figure S8](https://github.com/jinlx/Western-US-emission-packages/blob/main/Data_Analysis/FigS8_EnRs_OH_evaluation.pro)
This script is used to examine the OH level in the model and in the observations. <br /> 
To reproduce the figure, you need to compile the script and use the command below. <br />
EnRs_OH_evaluation,plume=0,/test,/save



