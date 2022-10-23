;; modified by lixu, 10/17/2019, used for compare the profiles from wecan and gc 
;; original script mk_flight_v2, lhu
;; modicication: parameters, set missing value for observation data, cancel filter.
;; update observation WE-CAN data: R4, 09/25/2020, add pressure
;; comp_mod_obs_vpro_all_v6, plume=1,/filter2,/all,/test,/cloud1,/save,/med,/nested,/prs
;; add 'addup' from version 5

;; add nonBB into calculation
;; this script is used for ND40 output vertical profiles. Observation here is comparing with simulations within GFED4, FINN, GFAS, QFED and non-BB. 

;; add CO into analysis;
;; add up cleanp option and nocleanup option.


;; v11, the time stamp is directly from Kate's help
;; trying to divide into different environments, 0208/2021
;; 1) BB-impacted: CO>85ppb&HCN>275ppt&CH3CN>200ppt
;;    - young (<1 day): 2-methylfuran > 0.7 ppt (95th percentile of nonsmoke background observations)
;;    - median (1-3 days): 2-methylfuran was not elevated but acrolein was >7.4 ppt
;;    - old (>3 days): neither 2-methylfuran nor acrolein was elevated, but acrylonitrile was >2.9 ppt; none of these age tracers is elevated but the smoke tracers(CO, HCN and CH3CN) are elevated

;; 2) urban-impacted: 2,2,4-trimethylpentane > 20 ppt & tetrachloroethene > 2 ppt & HFC-134a > 125 ppt or HCFC-22 > 275 ppt 

;; 3) bkg: transects where CO is below the 5th percentile

;; 4) transects: plume transect where CO is above the 5th percentile for that crossing
;;    - cores: plume where CO > 75th percentile
;;    - edges: plume transect where 25th < CO < 75th percentile
;;    - wings: plume transect where 5th < CO < 25th percentile

;; v11: introduce timestampes from Kate Odell

;; v13 add nonbb 

;; v14 add normalized keyword for checking injection height
;; and errorbar keyword
;; In this version, we scale up GFAS with a factor of 2.5 and 7 in fbf
;; We may need to change file directory after finished the test.

;; v15: add O3, NO2, PAN in individual VOCs, may need to give them a better position later
;; v15_v2: add OVOCS and NMHCs, methonal, formic acid and acetic acids are waiting to be added
;; fix no, no2
;; detailed pedges from sree
;; v16, try ch3cn
@bin_vert
    pro comp_mod_obs_vpro_wecan_celan_inj, $
        plume=plume, $
        filter2_nobb=filter2_nobb,filter3_nobb=filter3_nobb,$
        all=all,each=each, $
        ind_co=ind_co,ind_benz=ind_benz,ind_o3=ind_o3,ind_pan=ind_pan,ind_no2=ind_no2, $;; it should be used with each keyword_set
        test=test,save=save,$
        nested=nested,fbf=fbf,$
        metrics=metrics, $
        prs=prs, $
        med=med,average=average,$
        normalized=normalized, $
        errorbar=errorbar, $
        wus=wus,sus=sus, $ ; meaningless but being consistent with FIREX-AQ
        whole=whole,individual=individual,$
        no_obs_include=no_obs_include, codeployed= codeployed, obs_major=obs_major
        

myct,/WhGrYlRd
;///////////////////////////////////////////////////////////////////////////////////////////
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;++                                      STEP 1                                           ++
;++                      Set working dir, data dir, and parameters                        ++
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;///////////////////////////////////////////////////////////////////////////////////////////
  ; Get defaults (bmy, 6/7/11)
  X_OMARGIN   = !X.OMARGIN
  Y_OMARGIN   = !Y.OMARGIN
  X_MARGIN    = !X.MARGIN
  Y_MARGIN    = !Y.MARGIN
  P_CHARTHICK = !P.CHARTHICK
  P_THICK     = !P.THICK
  X_THICK     = !X.THICK
  Y_THICK     = !Y.THICK

  ; Plot parameters
  !X.OMARGIN=[6,1]
  !Y.OMARGIN=[2.5, 0.1]
  !X.MARGIN=[2,4]
  !Y.MARGIN=[1,1.5]
  !P.CHARTHICK=1.5
  !P.THICK=1.5
  !X.THICK=1.5
  !Y.THICK=1.5
  
  PEdges = [0,0.5,1,1.5,2.,2.5,3.,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12]
  ;if keyword_set(prs) then PEdges = [2000.,1000, 950, 900, 850, 800, 750, 700, 650, 600, 550, 500, 450, 400]
  ;if keyword_set(prs) then PEdges = 1000-25*findgen(25)  
  if keyword_set(prs) then PEdges = 1000-30*findgen(21)  
  ;if keyword_set(prs) then PEdges = 1000-40*findgen(16)  
  
  avo= 6.022e23  ;; molec/mol

  ; PLOT SETTING
  if keyword_set(test) then fi = './test_vp'  
    if keyword_set(save) then $
    open_device, /ps, /color, /landscape, $
               fi = './ps/' + fi + '.ps'           
   DEVICE, XSIZE=16, YSIZE=12, /INCHES

;///////////////////////////////////////////////////////////////////////////////////////////
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;++                                      STEP 2                                           ++
;++                       Read merge file and plane log files                             ++
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;///////////////////////////////////////////////////////////////////////////////////////////          
    dates_all=  ['20180724','20180726','20180730','20180731','20180802',$
                 '20180803','20180806','20180808','20180809','20180813','20180815','20180816','20180820',$
                 '20180823','20180826','20180828','20180906','20180910','20180913']
               
;; emission passes time stamp from Wade
    dates_emipass = ['20180809','20180809','20180815','20180815','20180726','20180726','20180809','20180806',$
                        '20180806','20180806',$
                        '20180813','20180815',$
                        '20180820','20180820'$
                        ,'20180813','20180803','20180813',$
                        '20180815','20180815','20180815','20180815','20180815',$
                        '20180910','20180910',$
                        '20180731','20180731',$
                        '20180731','20180913','20180913',$
                        '20180730','20180730',$
                        '20180813']
                        
                        
    start_time_emipass = [76935,82840,78905,88160,78735,81365,87380,81530,$
                    87550,87790,$
                    79710,89340,$
                    86655,87280,$
                    82270,91940,78195,$
                    74660,75230,75980,77275, 77870,$
                    73395,80765,$
                    74365,75255,$
                    86600,76607,76695,$
                    82060,82290,$
                    82140]
    
    end_time_emipass = [77040,83125,78985,88285,78850,81585,87480,81630,$
                     87775,87860,$
                     79860,89430,$
                     87015,87645,$
                     82370,92010,78320,$
                     74770,75565,76195,77545,78180,$
                     73455,80795,$
                     74500,75330,$
                     86740,76680,76740,$
                     82120,82365,$
                     82215]
                         
    avo= 6.022e23  ;; molec/mol
    
    for dd = 0, n_elements(dates_all) do begin
        ;;only consider whole flight
        if keyword_set(whole) then if dd ne 0 then break
        ;;consider each flight
        if keyword_set(individual) then if dd eq 0 then continue
    
;; ==============
;; WECAN DATA
;; ==============
        ;co_mod = [0]
        co_obs = [0]
        ;o3_mod = [0]
        o3_obs = [0]
        ;pan_mod = [0]
        pan_obs = [0]
        ;hcho_mod = [0]
        hcho_obs = [0]
        ;acet_mod = [0]
        acet_obs = [0]
        ;benz_mod = [0]
        benz_obs = [0]
        ;ch3oh_mod = [0]
        ch3oh_obs = [0]
        ald2_obs=[0]

        ;no_mod = [0]
        no_obs = [0]
        ;no2_mod = [0]
        no2_obs = [0]
        ;oh_mod = [0]
        oh_obs = [0]
        ;so2_mod = [0]
        so2_obs = [0]

        lat_obs  = [0]  
        lon_obs  = [0]
        prs_obs = [0]
        alt_obs = [0]

        ; used for time series
        utc_obs = [0]
        jday_obs = [0]
        lstime_obs = [0]


        noy_obs   = [0]
        ch3cn_obs = [0]
        

    ;; PTR and TOGA
    ;; PTR    
        isop_obs_ptr = [0]
        acet_obs_ptr = [0]
        propanal_obs_ptr = [0]
        mek_obs_ptr  = [0]
        butanal_obs_ptr = [0]
        macr_mvk_obs_ptr = [0]
        ALD2_obs_ptr = [0]
        ch2o_obs_ptr = [0]
        ch3cn_obs_ptr = [0]
        benz_obs_ptr = [0]
        tolu_obs_ptr = [0]

        c8h10_obs_ptr = [0]
        Xylenes_obs_ptr =  [0]
        
        monoterpenes_obs_ptr = [0]   
        
        dms_obs_ptr = [0]

        moh_obs_ptr = [0]
     
        ;; organic acid
        hcooh_obs_ptr  = [0]
        acta_obs_ptr   = [0]

     ; TOGA
        isop_obs_toga = [0]
        acet_obs_toga = [0]
        mek_obs_toga  = [0]
        macr_mvk_obs_toga     = [0]
        ALD2_obs_toga = [0]
        ch2o_obs_toga = [0]
        benz_obs_toga = [0]
        tolu_obs_toga = [0]

        mbo_obs_toga  = [0]
        propanal_obs_toga =[0]

        c2Butenal_obs_toga = [0]
        t2Butenal_obs_toga = [0]

        butanal_obs_toga = [0]

        etbenzene_obs_toga   = [0]

        tricyclene_obs_toga  = [0]
        apinene_obs_toga     = [0]
        camphene_obs_toga    = [0]
        bpinenemyrcene_obs_toga   = [0]
        limonened3carene_obs_toga = [0]

        dms_obs_toga = [0]

        moh_obs_toga = [0]
        hcn_obs_TOGA    = [0]

        c8h10_obs_toga = [0]
        Xylenes_obs_toga = [0]
        
        MPXYLENE_obs_toga = [0]
        OXYLENE_obs_toga  = [0]
        EtBenzene_obs_toga  = [0]
        
        ;;AWAS
        c2h2_obs_awas = [0]
        c2h4_obs_awas = [0]
        
        c3h6_obs_awas = [0]
        trans2butene_obs_awas = [0]
        x1Butene_obs_awas = [0]
        cisx2butene_obs_awas = [0]
        x1Hexene_obs_awas = [0]
        
        prpe_obs_awas = [0]
        
        cyclopentane_obs_awas = [0]
        cyclohexane_obs_awas = [0]
        methylcyclohexane_obs_awas = [0]
                
        c2h6_obs_awas = [0]
        isobutane_obs_awas = [0]
        nbutane_obs_awas = [0]
        isopentane_obs_awas = [0]
        nPentane_obs_awas = [0]
        hexane_obs_awas = [0]
        x24Dimethylpentane_obs_awas = [0]
        x22Dimethylbutane_obs_awas = [0]
        x3Methylpentane_obs_awas = [0]
        x23Dimethylpentane_obs_awas = [0]
        x2Methylhexane_obs_awas = [0]
        x3Methylhexane_obs_awas = [0]
        nHeptan_obs_awas = [0]
        x234Trimethylpentane_obs_awas = [0]
        x2Methylheptane_obs_awas = [0]
        x3Methylheptane_obs_awas = [0]
        nOctane_obs_awas = [0]
        nNonane_obs_awas = [0]
        undecane_obs_awas = [0]
        ndecane_obs_awas = [0]
        
        ALK4_obs_awas = [0]

        
        c5h8_obs_awas = [0]

        
        ;;for alk4 and prpe
        ;;alk4
        alk4_obs_toga = [0]
        propane_obs_toga = [0]
        butane_obs_toga = [0]
        pentane_obs_toga = [0]
        mepentane_obs_toga = [0]
        hexane_obs_toga = [0]
        trimepentane_obs_toga = [0]
        heptane_obs_toga = [0]
        octane_obs_toga = [0]
        
        ;;prpe
        prpe_obs_toga = [0]
        butene_obs_toga = [0]
        
        
        
        ;; CIMS
        ;; organic acid
        hcooh_obs_cims  = [0]
        acta_obs_cims   = [0]
        
    ;; for relative humidity
        rhum = [0]
        
    ;; for temperature
        temperature_obs = [0]

    ;; empty column for time index
        ind_enter_time = [0]
        ind_exit_time = [0]

;; ==============
;; GEOS-Chem DATA
;; ==============
    ;; GFAS_BOP_TOP
        co_gc_gfas_bop_top   = [0]
        O3_gc_gfas_bop_top   = [0]
        hcho_gc_gfas_bop_top = [0]
        pan_gc_gfas_bop_top  = [0]
        acet_gc_gfas_bop_top = [0]
        benz_gc_gfas_bop_top = [0]
        ch3oh_gc_gfas_bop_top= [0]
        ald2_gc_gfas_bop_top = [0]

        no_gc_gfas_bop_top   = [0]
        no2_gc_gfas_bop_top  = [0]
        so2_gc_gfas_bop_top  = [0]
        oh_gc_gfas_bop_top   = [0]

        date_gc_gfas_bop_top = [0]
        utc_gc_gfas_bop_top  = [0]
        doy_gc_gfas_bop_top  = [0]
        lat_gc_gfas_bop_top  = [0]
        lon_gc_gfas_bop_top  = [0]
        alt_gc_gfas_bop_top  = [0]
        prs_gc_gfas_bop_top  = [0]

     ;; add vocs
        c3h8_gc_gfas_bop_top = [0]
        c3h6_gc_gfas_bop_top = [0]
        c2h5oh_gc_gfas_bop_top = [0]
        c5h8_gc_gfas_bop_top   = [0]
        c7h8_gc_gfas_bop_top   = [0]
        dms_gc_gfas_bop_top    = [0]
        mek_gc_gfas_bop_top    = [0]
        c8h10_gc_gfas_bop_top  = [0]

        acta_gc_gfas_bop_top   = [0]
        macr_mvk_gc_gfas_bop_top = [0]  
        hcooh_gc_gfas_bop_top  = [0]

        mtpa_gc_gfas_bop_top   = [0]
        limo_gc_gfas_bop_top   = [0]
        mtpo_gc_gfas_bop_top   = [0]
        
        c2h4_gc_gfas_bop_top   = [0]
        c2h6_gc_gfas_bop_top   = [0]
    
    ;;add lumped species
        alk4_gc_gfas_bop_top   = [0]
        prpe_gc_gfas_bop_top   = [0]

        
    ;; GFAS_MAMI
        co_gc_gfas_mami   = [0]
        O3_gc_gfas_mami   = [0]
        hcho_gc_gfas_mami = [0]
        pan_gc_gfas_mami  = [0]
        acet_gc_gfas_mami = [0]
        benz_gc_gfas_mami = [0]
        ch3oh_gc_gfas_mami= [0]
        ald2_gc_gfas_mami = [0]

        no_gc_gfas_mami   = [0]
        no2_gc_gfas_mami  = [0]
        so2_gc_gfas_mami  = [0]
        oh_gc_gfas_mami   = [0]

        date_gc_gfas_mami = [0]
        utc_gc_gfas_mami  = [0]
        doy_gc_gfas_mami  = [0]
        lat_gc_gfas_mami  = [0]
        lon_gc_gfas_mami  = [0]
        alt_gc_gfas_mami  = [0]
        prs_gc_gfas_mami  = [0]

     ;; add vocs
        c3h8_gc_gfas_mami = [0]
        c3h6_gc_gfas_mami = [0]
        c2h5oh_gc_gfas_mami = [0]
        c5h8_gc_gfas_mami   = [0]
        c7h8_gc_gfas_mami   = [0]
        dms_gc_gfas_mami    = [0]
        mek_gc_gfas_mami    = [0]
        c8h10_gc_gfas_mami  = [0]

        acta_gc_gfas_mami   = [0]
        macr_mvk_gc_gfas_mami = [0]  
        hcooh_gc_gfas_mami  = [0]

        mtpa_gc_gfas_mami   = [0]
        limo_gc_gfas_mami   = [0]
        mtpo_gc_gfas_mami   = [0]
        
        c2h4_gc_gfas_mami   = [0]
        c2h6_gc_gfas_mami   = [0]
    
    ;;add lumped species
        alk4_gc_gfas_mami   = [0]
        prpe_gc_gfas_mami   = [0]
    ;; GFAS_NOBB
        co_gc_gfas_nobb   = [0]
        O3_gc_gfas_nobb   = [0]
        hcho_gc_gfas_nobb = [0]
        pan_gc_gfas_nobb  = [0]
        acet_gc_gfas_nobb = [0]
        benz_gc_gfas_nobb = [0]
        ch3oh_gc_gfas_nobb= [0]
        ald2_gc_gfas_nobb = [0]

        no_gc_gfas_nobb   = [0]
        no2_gc_gfas_nobb  = [0]
        so2_gc_gfas_nobb  = [0]
        oh_gc_gfas_nobb   = [0]

        date_gc_gfas_nobb = [0]
        utc_gc_gfas_nobb  = [0]
        doy_gc_gfas_nobb  = [0]
        lat_gc_gfas_nobb  = [0]
        lon_gc_gfas_nobb  = [0]
        alt_gc_gfas_nobb  = [0]
        prs_gc_gfas_nobb  = [0]

     ;; add vocs
        c3h8_gc_gfas_nobb = [0]
        c3h6_gc_gfas_nobb = [0]
        c2h5oh_gc_gfas_nobb = [0]
        c5h8_gc_gfas_nobb   = [0]
        c7h8_gc_gfas_nobb   = [0]
        dms_gc_gfas_nobb    = [0]
        mek_gc_gfas_nobb    = [0]
        c8h10_gc_gfas_nobb  = [0]

        acta_gc_gfas_nobb   = [0]
        macr_mvk_gc_gfas_nobb = [0]  
        hcooh_gc_gfas_nobb  = [0]

        mtpa_gc_gfas_nobb   = [0]
        limo_gc_gfas_nobb   = [0]
        mtpo_gc_gfas_nobb   = [0]
        
        c2h4_gc_gfas_nobb   = [0]
        c2h6_gc_gfas_nobb   = [0]
    
    ;;add lumped species
        alk4_gc_gfas_nobb   = [0]
        prpe_gc_gfas_nobb   = [0]
    ;; GFAS_PBL65_FT35
        co_gc_gfas_pbl65_ft35   = [0]
        O3_gc_gfas_pbl65_ft35   = [0]
        hcho_gc_gfas_pbl65_ft35 = [0]
        pan_gc_gfas_pbl65_ft35  = [0]
        acet_gc_gfas_pbl65_ft35 = [0]
        benz_gc_gfas_pbl65_ft35 = [0]
        ch3oh_gc_gfas_pbl65_ft35= [0]
        ald2_gc_gfas_pbl65_ft35 = [0]

        no_gc_gfas_pbl65_ft35   = [0]
        no2_gc_gfas_pbl65_ft35  = [0]
        so2_gc_gfas_pbl65_ft35  = [0]
        oh_gc_gfas_pbl65_ft35   = [0]

        date_gc_gfas_pbl65_ft35 = [0]
        utc_gc_gfas_pbl65_ft35  = [0]
        doy_gc_gfas_pbl65_ft35  = [0]
        lat_gc_gfas_pbl65_ft35  = [0]
        lon_gc_gfas_pbl65_ft35  = [0]
        alt_gc_gfas_pbl65_ft35  = [0]
        prs_gc_gfas_pbl65_ft35  = [0]

     ;; add vocs
        c3h8_gc_gfas_pbl65_ft35 = [0]
        c3h6_gc_gfas_pbl65_ft35 = [0]
        c2h5oh_gc_gfas_pbl65_ft35 = [0]
        c5h8_gc_gfas_pbl65_ft35   = [0]
        c7h8_gc_gfas_pbl65_ft35   = [0]
        dms_gc_gfas_pbl65_ft35    = [0]
        mek_gc_gfas_pbl65_ft35    = [0]
        c8h10_gc_gfas_pbl65_ft35  = [0]

        acta_gc_gfas_pbl65_ft35   = [0]
        macr_mvk_gc_gfas_pbl65_ft35 = [0]  
        hcooh_gc_gfas_pbl65_ft35  = [0]

        mtpa_gc_gfas_pbl65_ft35   = [0]
        limo_gc_gfas_pbl65_ft35   = [0]
        mtpo_gc_gfas_pbl65_ft35   = [0]
        
        c2h4_gc_gfas_pbl65_ft35   = [0]
        c2h6_gc_gfas_pbl65_ft35   = [0]
    
    ;;add lumped species
        alk4_gc_gfas_pbl65_ft35   = [0]
        prpe_gc_gfas_pbl65_ft35   = [0]

    ;; GFAS_SF_MAMI
        co_gc_gfas_sf_mami   = [0]
        O3_gc_gfas_sf_mami   = [0]
        hcho_gc_gfas_sf_mami = [0]
        pan_gc_gfas_sf_mami  = [0]
        acet_gc_gfas_sf_mami = [0]
        benz_gc_gfas_sf_mami = [0]
        ch3oh_gc_gfas_sf_mami= [0]
        ald2_gc_gfas_sf_mami = [0]

        no_gc_gfas_sf_mami   = [0]
        no2_gc_gfas_sf_mami  = [0]
        so2_gc_gfas_sf_mami  = [0]
        oh_gc_gfas_sf_mami   = [0]

        date_gc_gfas_sf_mami = [0]
        utc_gc_gfas_sf_mami  = [0]
        doy_gc_gfas_sf_mami  = [0]
        lat_gc_gfas_sf_mami  = [0]
        lon_gc_gfas_sf_mami  = [0]
        alt_gc_gfas_sf_mami  = [0]
        prs_gc_gfas_sf_mami  = [0]

     ;; add vocs
        c3h8_gc_gfas_sf_mami = [0]
        c3h6_gc_gfas_sf_mami = [0]
        c2h5oh_gc_gfas_sf_mami = [0]
        c5h8_gc_gfas_sf_mami   = [0]
        c7h8_gc_gfas_sf_mami   = [0]
        dms_gc_gfas_sf_mami    = [0]
        mek_gc_gfas_sf_mami    = [0]
        c8h10_gc_gfas_sf_mami  = [0]

        acta_gc_gfas_sf_mami   = [0]
        macr_mvk_gc_gfas_sf_mami = [0]  
        hcooh_gc_gfas_sf_mami  = [0]

        mtpa_gc_gfas_sf_mami   = [0]
        limo_gc_gfas_sf_mami   = [0]
        mtpo_gc_gfas_sf_mami   = [0]
        
        c2h4_gc_gfas_sf_mami   = [0]
        c2h6_gc_gfas_sf_mami   = [0]
    
    ;;add lumped species
        alk4_gc_gfas_sf_mami   = [0]
        prpe_gc_gfas_sf_mami   = [0]

    ;; GFAS_SURFACE
        co_gc_gfas_surface   = [0]
        O3_gc_gfas_surface   = [0]
        hcho_gc_gfas_surface = [0]
        pan_gc_gfas_surface  = [0]
        acet_gc_gfas_surface = [0]
        benz_gc_gfas_surface = [0]
        ch3oh_gc_gfas_surface= [0]
        ald2_gc_gfas_surface = [0]

        no_gc_gfas_surface   = [0]
        no2_gc_gfas_surface  = [0]
        so2_gc_gfas_surface  = [0]
        oh_gc_gfas_surface   = [0]

        date_gc_gfas_surface = [0]
        utc_gc_gfas_surface  = [0]
        doy_gc_gfas_surface  = [0]
        lat_gc_gfas_surface  = [0]
        lon_gc_gfas_surface  = [0]
        alt_gc_gfas_surface  = [0]
        prs_gc_gfas_surface  = [0]

     ;; add vocs
        c3h8_gc_gfas_surface = [0]
        c3h6_gc_gfas_surface = [0]
        c2h5oh_gc_gfas_surface = [0]
        c5h8_gc_gfas_surface   = [0]
        c7h8_gc_gfas_surface   = [0]
        dms_gc_gfas_surface    = [0]
        mek_gc_gfas_surface    = [0]
        c8h10_gc_gfas_surface  = [0]

        acta_gc_gfas_surface   = [0]
        macr_mvk_gc_gfas_surface = [0]  
        hcooh_gc_gfas_surface  = [0]

        mtpa_gc_gfas_surface   = [0]
        limo_gc_gfas_surface   = [0]
        mtpo_gc_gfas_surface   = [0]
        
        c2h4_gc_gfas_surface   = [0]
        c2h6_gc_gfas_surface   = [0]
    
    ;;add lumped species
        alk4_gc_gfas_surface   = [0]
        prpe_gc_gfas_surface   = [0]


    ;; empty column for time index
        ind_enter_time = [0]
        ind_exit_time = [0]

    ; smart loop for each fligiht and total flights from lu
        if dd eq 0 then dates = dates_all else dates=dates_all[dd-1]
        dates = STRTRIM(string(dates),1)
        for n=0,n_elements(dates)-1 do begin
    ;; GC file for all data
            str_date=dates[n]
            print,'processing '+str_date


;; ======================
;; wecan observation
;; ======================
;; -------------------------------------------------------------------------------------
            ;c130fi = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
            ;    'mrg2sav/R4_merges/wecan-mrg60-c130_merge_'+dates[n]+'_R4.sav'

            ;c130fi = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
            ;    'mrg2sav/R4_merges_avg/wecan-mrg1m-c130_merge_'+dates[n]+'_R4.sav'
            
            c130fi = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                'mrg2sav/WECAN/R4_merges_avg/wecan-mrg5m-c130_merge_'+dates[n]+'_R4.sav'
            restore, c130fi
            tmp_co_obs = c130.CO_PICARRO
            tmp_o3_obs =  c130.O3
            tmp_pan_obs = (c130.pan)/1e3 ;; ppt to ppb
            tmp_hcho_obs = c130.FORMALDEHYDE_MIXINGRATIO_PTR
            tmp_acet_obs = c130.ACETONE_MIXINGRATIO_PTR*0.78
            tmp_benz_obs = c130.BENZENE_MIXINGRATIO_PTR
            tmp_ch3oh_obs = c130.METHANOL_MIXINGRATIO_PTR
            tmp_ald2_obs = c130.ACETALDEHYDE_MIXINGRATIO_PTR
            tmp_no_obs = c130.no/1e3 ;; ppt to ppb
            tmp_no2_obs = c130.no2/1e3 ;; ppt to ppb
            tmp_so2_obs = c130.SO2_UWCIMS
            
            tmp_lat_obs  = c130.latitude
            filter_lat_obs = c130.latitude
            tmp_lon_obs  = c130.longitude
            tmp_prs_obs = c130.pressure  ;; Pa
            tmp_alt_obs = c130.altp

            ;noy_obs   = [noy_obs,dc8.NOy_ESRL]
            tmp_ch3cn_obs = c130.Acetonitrile_MixingRatio_PTR ; we also have mixing ratio of ch3cn, need to check the unit lixu
            tmp_jday_obs = c130.jday
            tmp_utc_obs = c130.utc_mid 
            tmp_lstime_obs = c130.LOCAL_SUN_TIME

    ;; PTR and TOGA
    ;; PTR    
            tmp_isop_obs_ptr = c130.Isoprene_MixingRatio_PTR
            tmp_acet_obs_ptr = c130.Acetone_MixingRatio_PTR*0.78
            tmp_propanal_obs_ptr = c130.Acetone_MixingRatio_PTR*0.22
            tmp_mek_obs_ptr  = c130.MEK_MixingRatio_PTR*0.8
            tmp_butanal_obs_ptr  = c130.MEK_MixingRatio_PTR*0.2
            tmp_macr_mvk_obs_ptr     = c130.MACR_MVK_MixingRatio_PTR
            tmp_ALD2_obs_ptr = c130.Acetaldehyde_MixingRatio_PTR
            tmp_ch3cn_obs_ptr = c130.Acetonitrile_MixingRatio_PTR ; we also have mixing ratio of ch3cn, need to check the unit lixu

            tmp_ch2o_obs_ptr = c130.Formaldehyde_MixingRatio_PTR
            tmp_benz_obs_ptr = c130.Benzene_MixingRatio_PTR
            tmp_tolu_obs_ptr = c130.Toluene_MixingRatio_PTR

            tmp_c8h10_obs_ptr   = c130.C8_aromatics_MixingRatio_PTR 
            tmp_Xylenes_obs_ptr = c130.C8_aromatics_MixingRatio_PTR*0.65

            tmp_monoterpenes_obs_ptr = c130.Monoterpenes_MixingRatio_PTR  
            
            tmp_moh_obs_ptr = c130.Methanol_MixingRatio_PTR
            
            tmp_acta_obs_ptr   = c130.AceticAcid_MixingRatio_PTR
            tmp_hcooh_obs_ptr  = c130.FormicAcid_PTR

    ; TOGA
            tmp_isop_obs_toga = c130.Isoprene_TOGA/1e3
            tmp_acet_obs_toga = c130.Acetone_TOGA/1e3
            tmp_mek_obs_toga  = c130.MEK_TOGA/1e3
            tmp_macr_mvk_obs_toga     = (c130.MVK_TOGA+c130.MAC_TOGA)/1e3
            tmp_ALD2_obs_toga = c130.ch3cho_toga/1e3
            tmp_ch2o_obs_toga = c130.CH2O_TOGA/1e3
            tmp_benz_obs_toga = c130.Benzene_TOGA/1e3
            tmp_tolu_obs_toga = c130.Toluene_TOGA/1e3

            tmp_mbo_obs_toga  =c130.MBO_TOGA/1e3
            tmp_propanal_obs_toga =c130.Propanal_TOGA/1e3

            tmp_moh_obs_toga = c130.CH3OH_TOGA/1e3
;;need to changed later!!!!!!!!!!!!!!!!!!!!!!!
            tmp_c2Butenal_obs_toga = 0
            tmp_t2Butenal_obs_toga = 0

            tmp_butanal_obs_toga = c130.Butanal_TOGA/1e3

            tmp_etbenzene_obs_toga   = c130.EtBenzene_TOGA/1e3

            tmp_tricyclene_obs_toga  = c130.Tricyclene_TOGA/1e3
            tmp_apinene_obs_toga     = c130.aPinene_TOGA/1e3
            tmp_camphene_obs_toga  = c130.Camphene_TOGA/1e3
            tmp_bpinenemyrcene_obs_toga   = c130.bPineneMyrcene_TOGA/1e3
            tmp_limonened3carene_obs_toga = c130.LimoneneD3Carene_TOGA/1e3 
            tmp_hcn_obs_toga = c130.HCN_TOGA/1000 ; pptv to ppbv

            tmp_c8h10_obs_toga     = (c130.MPXYLENE_TOGA + c130.OXYLENE_TOGA +  c130.EtBenzene_TOGA)/1e3 ; 
            tmp_Xylenes_obs_toga   = (c130.MPXYLENE_TOGA + c130.OXYLENE_TOGA)/1e3 ; 
            tmp_MPXYLENE_obs_toga  = (c130.MPXYLENE_TOGA)/1e3 ; 
            tmp_OXYLENE_obs_toga   = (c130.OXYLENE_TOGA)/1e3 ; 
            tmp_EtBenzene_obs_toga = (c130.EtBenzene_TOGA)/1e3 ; 

        ;;AWAS
            tmp_c2h2_obs_awas = c130.Ethyne_AWAS
            tmp_c2h4_obs_awas   = c130.Ethene_AWAS
            tmp_c2h6_obs_awas   = c130.Ethane_AWAS
            
            ;; PRPE
            tmp_c3h6_obs_awas = c130.Propane_AWAS
            tmp_trans2butene_obs_awas = c130.transx2butene_AWAS
            tmp_x1Butene_obs_awas = c130.x1Butene_AWAS
            tmp_cisx2butene_obs_awas = c130.cisx2butene_awas
            tmp_x1Hexene_obs_awas = c130.x1Hexene_AWAS
            
            tmp_cyclopentane_obs_awas = c130.Cyclopentane_AWAS
            tmp_cyclohexane_obs_awas = c130.Cyclohexane_CHB_AWAS
            tmp_methylcyclohexane_obs_awas = c130.Methylcyclohexane_CHB_AWAS
            
            tmp_prpe_obs_awas = tmp_c3h6_obs_awas + tmp_trans2butene_obs_awas + tmp_x1Butene_obs_awas + tmp_cisx2butene_obs_awas + tmp_x1Hexene_obs_awas + $
                                tmp_cyclopentane_obs_awas + tmp_cyclohexane_obs_awas + tmp_methylcyclohexane_obs_awas
            ;; ALK4
            tmp_isobutane_obs_awas = c130.Isobutane_AWAS
            tmp_nbutane_obs_awas = c130.N_butane_AWAS
            tmp_isopentane_obs_awas = c130.Isopentane_AWAS
            tmp_nPentane_obs_awas = c130.nPentane_AWAS
            tmp_hexane_obs_awas = c130.Hexane_AWAS
            tmp_x24Dimethylpentane_obs_awas = c130.x24Dimethylpentane_AWAS
            tmp_x22Dimethylbutane_obs_awas = c130.x22Dimethylbutane_AWAS
            tmp_x3Methylpentane_obs_awas = c130.x3Methylpentane_AWAS
            tmp_x23Dimethylpentane_obs_awas = c130.x23Dimethylpentane_AWAS
            tmp_x2Methylhexane_obs_awas = c130.x2Methylhexane_AWAS
            tmp_x3Methylhexane_obs_awas = c130.x3Methylhexane_AWAS
            tmp_nHeptan_obs_awas       = c130.nHeptan_CHB_AWAS
            tmp_x234Trimethylpentane_obs_awas = c130.x234Trimethylpentane_AWAS
            tmp_x2Methylheptane_obs_awas = c130.x2Methylheptane_AWAS
            tmp_x3Methylheptane_obs_awas = c130.x3Methylheptane_AWAS
            tmp_nOctane_obs_awas = c130.nOctane_AWAS
            tmp_nNonane_obs_awas= c130.nNonane_AWAS
            tmp_undecane_obs_awas = c130.Undecane_AWAS
            tmp_ndecane_obs_awas = c130.ndecane_AWAS
            
            tmp_alk4_obs_awas = tmp_isobutane_obs_awas + tmp_nbutane_obs_awas + tmp_isopentane_obs_awas + tmp_nPentane_obs_awas + $
                                    tmp_hexane_obs_awas + tmp_x24Dimethylpentane_obs_awas + tmp_x22Dimethylbutane_obs_awas + tmp_x3Methylpentane_obs_awas + $
                                    tmp_x23Dimethylpentane_obs_awas + tmp_x2Methylhexane_obs_awas + tmp_x3Methylhexane_obs_awas + tmp_nHeptan_obs_awas + $
                                    tmp_x234Trimethylpentane_obs_awas + tmp_x2Methylheptane_obs_awas + tmp_x3Methylheptane_obs_awas + $
                                    tmp_nOctane_obs_awas + tmp_nNonane_obs_awas + tmp_undecane_obs_awas + tmp_ndecane_obs_awas
                                    
            ;;for alk4 and prpe from TOGA, may need double check later!!!!!!!!
            ;;alk4
            tmp_propane_obs_toga = c130.Propane_TOGA/1e3
            tmp_butane_obs_toga = (c130.iButane_TOGA + c130.nButane_TOGA)/1e3
            tmp_pentane_obs_toga = (c130.iPentane_TOGA + c130.nPentane_TOGA)/1e3
            tmp_mepentane_obs_toga = (c130.x2MePentane_TOGA + c130.x3MePentane_TOGA)/1e3
            tmp_hexane_obs_toga = c130.nHexane_TOGA/1e3
            tmp_trimepentane_obs_toga = c130.x224TrimePentane_TOGA/1e3
            tmp_heptane_obs_toga = c130.nHeptane_TOGA/1e3
            tmp_octane_obs_toga = c130.nOctane_TOGA/1e3
            
            tmp_alk4_obs_toga = tmp_propane_obs_toga + tmp_butane_obs_toga + tmp_pentane_obs_toga + tmp_mepentane_obs_toga + $
                                    tmp_hexane_obs_toga + tmp_trimepentane_obs_toga + tmp_heptane_obs_toga + tmp_octane_obs_toga

            ;;prpe
            tmp_prpe_obs_toga = c130.iButene1Butene_TOGA
            tmp_butene_obs_toga = c130.iButene1Butene_TOGA

            ;; CIMS
            tmp_hcooh_obs_cims  = c130.CH2O2_UWCIMS/1e3
            tmp_acta_obs_cims  = c130.C2H4O2_UWCIMS/1e3

            
    ;; for cloud
            tmp_rhum = c130.rhum        
    ;; for temperature
            tmp_temperature_obs = c130.TEMPERATURE

    ;; prepare data for 20180826, lat is less than zero and it is not simulated by the obs,fixed
            i = where(filter_lat_obs le 0,ct) 
            if ct gt 0 then begin
                tmp_lat_obs[i] = !VALUES.F_NAN
                remove_ind=where(finite(tmp_lat_obs),ct)
                if remove_ind[0] ne -1 then begin
                    tmp_co_obs = tmp_co_obs[remove_ind]
                    tmp_o3_obs =  tmp_o3_obs[remove_ind]
                    tmp_pan_obs = tmp_pan_obs[remove_ind]
                    tmp_hcho_obs = tmp_hcho_obs[remove_ind]
                    tmp_acet_obs = tmp_acet_obs[remove_ind]
                    tmp_benz_obs = tmp_benz_obs[remove_ind]
                    tmp_ch3oh_obs = tmp_ch3oh_obs[remove_ind]
                    tmp_ald2_obs = tmp_ald2_obs[remove_ind]
                    tmp_no_obs = tmp_no_obs[remove_ind]
                    tmp_no2_obs =tmp_no2_obs[remove_ind]
                    tmp_so2_obs = tmp_so2_obs[remove_ind]
                    tmp_lat_obs  = tmp_lat_obs[remove_ind]
                    tmp_lon_obs  = tmp_lon_obs[remove_ind]
                    tmp_prs_obs = tmp_prs_obs[remove_ind]  ;; Pa
                    tmp_alt_obs = tmp_alt_obs[remove_ind]

                    ;noy_obs   = [noy_obs,dc8.NOy_ESRL]
                    tmp_ch3cn_obs = tmp_ch3cn_obs[remove_ind]; we also have mixing ratio of ch3cn, need to check the unit lixu
                    tmp_jday_obs = tmp_jday_obs[remove_ind]
                    tmp_utc_obs = tmp_utc_obs[remove_ind]
                    tmp_lstime_obs = tmp_lstime_obs[remove_ind]


            ;; PTR and TOGA
            ;; PTR    
                    tmp_isop_obs_ptr = tmp_isop_obs_ptr[remove_ind]
                    tmp_acet_obs_ptr = tmp_acet_obs_ptr[remove_ind]
                    tmp_propanal_obs_ptr = tmp_propanal_obs_ptr[remove_ind]
               
                    tmp_mek_obs_ptr  = tmp_mek_obs_ptr[remove_ind]
                    tmp_butanal_obs_ptr  = tmp_butanal_obs_ptr[remove_ind]
                    
                    tmp_macr_mvk_obs_ptr     = tmp_macr_mvk_obs_ptr[remove_ind]
                    tmp_ALD2_obs_ptr = tmp_ALD2_obs_ptr[remove_ind]
                    tmp_ch3cn_obs_ptr = tmp_ch3cn_obs_ptr[remove_ind]; we also have mixing ratio of ch3cn, need to check the unit lixu

                    tmp_ch2o_obs_ptr = tmp_ch2o_obs_ptr[remove_ind]
                    tmp_benz_obs_ptr = tmp_benz_obs_ptr[remove_ind]
                    tmp_tolu_obs_ptr = tmp_tolu_obs_ptr[remove_ind]

                    tmp_c8h10_obs_ptr   = tmp_c8h10_obs_ptr[remove_ind]
                    tmp_Xylenes_obs_ptr = tmp_Xylenes_obs_ptr[remove_ind]

                    tmp_monoterpenes_obs_ptr = tmp_monoterpenes_obs_ptr[remove_ind]
                    
                    tmp_moh_obs_ptr = tmp_moh_obs_ptr[remove_ind]
                    
                    tmp_hcooh_obs_ptr  = tmp_hcooh_obs_ptr[remove_ind]
                    tmp_acta_obs_ptr   = tmp_acta_obs_ptr[remove_ind]
                    
            ; TOGA
                    tmp_isop_obs_toga = tmp_isop_obs_toga[remove_ind]
                    tmp_acet_obs_toga =tmp_acet_obs_toga[remove_ind]
                    tmp_mek_obs_toga  = tmp_mek_obs_toga[remove_ind]
                    tmp_macr_mvk_obs_toga     = tmp_macr_mvk_obs_toga[remove_ind]
                    tmp_ALD2_obs_toga =tmp_ALD2_obs_toga[remove_ind]
                    tmp_ch2o_obs_toga = tmp_ch2o_obs_toga[remove_ind]
                    tmp_benz_obs_toga = tmp_benz_obs_toga[remove_ind]
                    tmp_tolu_obs_toga = tmp_tolu_obs_toga[remove_ind]

                    tmp_mbo_obs_toga  =tmp_mbo_obs_toga[remove_ind]
                    tmp_propanal_obs_toga =tmp_propanal_obs_toga[remove_ind]

                    tmp_c2Butenal_obs_toga =tmp_c2Butenal_obs_toga[remove_ind]
                    tmp_t2Butenal_obs_toga = tmp_t2Butenal_obs_toga[remove_ind]
                    tmp_butanal_obs_toga = tmp_butanal_obs_toga[remove_ind]

                    tmp_etbenzene_obs_toga   = tmp_etbenzene_obs_toga[remove_ind]
                    tmp_tricyclene_obs_toga  = tmp_tricyclene_obs_toga[remove_ind]
                    tmp_apinene_obs_toga     =tmp_apinene_obs_toga[remove_ind]
                    tmp_camphene_obs_toga  = tmp_camphene_obs_toga[remove_ind]
                    tmp_bpinenemyrcene_obs_toga   = tmp_bpinenemyrcene_obs_toga[remove_ind]
                    tmp_limonened3carene_obs_toga = tmp_limonened3carene_obs_toga[remove_ind]
                    tmp_hcn_obs_toga = tmp_hcn_obs_toga[remove_ind]

                    tmp_c8h10_obs_toga     = tmp_c8h10_obs_toga[remove_ind]
                    tmp_Xylenes_obs_toga   = tmp_Xylenes_obs_toga[remove_ind]
                    tmp_MPXYLENE_obs_toga  = tmp_MPXYLENE_obs_toga[remove_ind]
                    tmp_OXYLENE_obs_toga   = tmp_OXYLENE_obs_toga[remove_ind]
                    tmp_EtBenzene_obs_toga = tmp_EtBenzene_obs_toga[remove_ind]
                    
                    tmp_moh_obs_toga = tmp_moh_obs_toga[remove_ind]
                    
                    tmp_rhum = tmp_rhum[remove_ind]
                    tmp_temperature_obs = tmp_temperature_obs[remove_ind]


                ;;for alk4 and prpe, toga, may need to double check itlater!!!!!!
                ;;alk4
                    tmp_propane_obs_toga = tmp_propane_obs_toga[remove_ind]
                    tmp_butane_obs_toga = tmp_butane_obs_toga[remove_ind]
                    tmp_pentane_obs_toga = tmp_pentane_obs_toga[remove_ind]
                    tmp_mepentane_obs_toga = tmp_mepentane_obs_toga[remove_ind]
                    tmp_hexane_obs_toga = tmp_hexane_obs_toga[remove_ind]
                    tmp_trimepentane_obs_toga = tmp_trimepentane_obs_toga[remove_ind]
                    tmp_heptane_obs_toga = tmp_heptane_obs_toga[remove_ind]
                    tmp_octane_obs_toga = tmp_octane_obs_toga[remove_ind]

                    tmp_alk4_obs_toga = tmp_alk4_obs_toga[remove_ind]

                    ;;prpe
                    tmp_prpe_obs_toga = tmp_prpe_obs_toga[remove_ind]
                    tmp_butene_obs_toga = tmp_butene_obs_toga[remove_ind]
                    
                    
                ;;AWAS
                    tmp_c2h2_obs_awas = tmp_c2h2_obs_awas[remove_ind]
                    tmp_c2h4_obs_awas   = tmp_c2h4_obs_awas[remove_ind]
                    tmp_c2h6_obs_awas   = tmp_c2h6_obs_awas[remove_ind]
                    
                    ;; PRPE
                    tmp_c3h6_obs_awas = tmp_c3h6_obs_awas[remove_ind]
                    tmp_trans2butene_obs_awas =tmp_trans2butene_obs_awas[remove_ind]
                    tmp_x1Butene_obs_awas = tmp_x1Butene_obs_awas[remove_ind]
                    tmp_cisx2butene_obs_awas = tmp_cisx2butene_obs_awas[remove_ind]
                    tmp_x1Hexene_obs_awas = tmp_x1Hexene_obs_awas[remove_ind]

                    tmp_cyclopentane_obs_awas = tmp_cyclopentane_obs_awas[remove_ind]
                    tmp_cyclohexane_obs_awas = tmp_cyclohexane_obs_awas[remove_ind]
                    tmp_methylcyclohexane_obs_awas = tmp_methylcyclohexane_obs_awas[remove_ind]

                    tmp_prpe_obs_awas = tmp_prpe_obs_awas[remove_ind]
                    
                    ;; ALK4
                    tmp_isobutane_obs_awas = tmp_isobutane_obs_awas[remove_ind]
                    tmp_nbutane_obs_awas = tmp_nbutane_obs_awas[remove_ind]
                    tmp_isopentane_obs_awas = tmp_isopentane_obs_awas[remove_ind]
                    tmp_nPentane_obs_awas = tmp_nPentane_obs_awas[remove_ind]
                    tmp_hexane_obs_awas = tmp_hexane_obs_awas[remove_ind]
                    tmp_x24Dimethylpentane_obs_awas = tmp_x24Dimethylpentane_obs_awas[remove_ind]
                    tmp_x22Dimethylbutane_obs_awas = tmp_x22Dimethylbutane_obs_awas[remove_ind]
                    tmp_x3Methylpentane_obs_awas = tmp_x3Methylpentane_obs_awas[remove_ind]
                    tmp_x23Dimethylpentane_obs_awas = tmp_x23Dimethylpentane_obs_awas[remove_ind]
                    tmp_x2Methylhexane_obs_awas = tmp_x2Methylhexane_obs_awas[remove_ind]
                    tmp_x3Methylhexane_obs_awas = tmp_x3Methylhexane_obs_awas[remove_ind]
                    tmp_nHeptan_obs_awas       = tmp_nHeptan_obs_awas[remove_ind]
                    tmp_x234Trimethylpentane_obs_awas = tmp_x234Trimethylpentane_obs_awas[remove_ind]
                    tmp_x2Methylheptane_obs_awas = tmp_x2Methylheptane_obs_awas[remove_ind]
                    tmp_x3Methylheptane_obs_awas = tmp_x3Methylheptane_obs_awas[remove_ind]
                    tmp_nOctane_obs_awas = tmp_nOctane_obs_awas[remove_ind]
                    tmp_nNonane_obs_awas=tmp_nNonane_obs_awas[remove_ind]
                    tmp_undecane_obs_awas = tmp_undecane_obs_awas[remove_ind]
                    tmp_ndecane_obs_awas =tmp_ndecane_obs_awas[remove_ind]      
                    tmp_alk4_obs_awas = tmp_alk4_obs_awas[remove_ind]

                    ;; CIMS
                    tmp_hcooh_obs_cims  = tmp_hcooh_obs_cims[remove_ind]
                    tmp_acta_obs_cims   = tmp_acta_obs_cims[remove_ind]
                    
                endif
            endif
            
            help, tmp_co_obs
            co_obs = [co_obs,tmp_co_obs]
            ;pan_obs = [pan_obs,tmp_pan_obs]
            o3_obs = [o3_obs, tmp_o3_obs]
            pan_obs = [pan_obs, tmp_pan_obs]   
            hcho_obs = [hcho_obs,tmp_hcho_obs]
            acet_obs = [acet_obs,tmp_acet_obs]
            benz_obs = [benz_obs,tmp_benz_obs]
            ch3oh_obs= [ch3oh_obs,tmp_ch3oh_obs]
            ald2_obs = [ald2_obs,tmp_ald2_obs]

            no_obs = [no_obs,tmp_no_obs]
            no2_obs = [no2_obs,tmp_no2_obs]
            ; oh_obs = [oh_obs,tmp.spcoh] ;; dummy
            so2_obs = [so2_obs,tmp_so2_obs]

            lat_obs  = [lat_obs,tmp_lat_obs]
            lon_obs  = [lon_obs,tmp_lon_obs]
            prs_obs = [prs_obs, tmp_prs_obs]  ;; Pa
            alt_obs = [alt_obs, tmp_alt_obs]

            ;noy_obs   = [noy_obs,dc8.NOy_ESRL]
            ch3cn_obs = [ch3cn_obs,tmp_ch3cn_obs] ; we also have mixing ratio of ch3cn, need to check the unit lixu
            jday_obs = [jday_obs,tmp_jday_obs]
            ;utc_obs = [utc_obs,utc_obs] 
            lstime_obs = [lstime_obs,tmp_lstime_obs]


    ;; PTR and TOGA
    ;; PTR    
            isop_obs_ptr = [isop_obs_ptr,tmp_isop_obs_ptr]
            acet_obs_ptr = [acet_obs_ptr,tmp_acet_obs_ptr]
            propanal_obs_ptr = [propanal_obs_ptr,tmp_propanal_obs_ptr]
            mek_obs_ptr  = [mek_obs_ptr,tmp_mek_obs_ptr]
            butanal_obs_ptr  = [butanal_obs_ptr,tmp_butanal_obs_ptr]
            macr_mvk_obs_ptr     = [macr_mvk_obs_ptr,tmp_macr_mvk_obs_ptr]
            ALD2_obs_ptr = [ALD2_obs_ptr,tmp_ALD2_obs_ptr]
            ch3cn_obs_ptr = [ch3cn_obs_ptr,tmp_ch3cn_obs_ptr] ; we also have mixing ratio of ch3cn, need to check the unit lixu

            ch2o_obs_ptr = [ch2o_obs_ptr,tmp_ch2o_obs_ptr]
            benz_obs_ptr = [benz_obs_ptr,tmp_benz_obs_ptr]
            tolu_obs_ptr = [tolu_obs_ptr,tmp_tolu_obs_ptr]

            c8h10_obs_ptr   = [c8h10_obs_ptr,tmp_c8h10_obs_ptr]
            Xylenes_obs_ptr = [Xylenes_obs_ptr,tmp_Xylenes_obs_ptr]
            
            monoterpenes_obs_ptr = [monoterpenes_obs_ptr,tmp_monoterpenes_obs_ptr]
            
            moh_obs_ptr = [moh_obs_ptr,tmp_moh_obs_ptr]

            hcooh_obs_ptr  = [hcooh_obs_ptr,tmp_hcooh_obs_ptr]
            acta_obs_ptr   = [acta_obs_ptr,tmp_acta_obs_ptr]            
    ; TOGA
            isop_obs_toga = [isop_obs_toga,tmp_isop_obs_toga]
            acet_obs_toga = [acet_obs_toga,tmp_acet_obs_toga]
            mek_obs_toga  = [mek_obs_toga,tmp_mek_obs_toga]
            macr_mvk_obs_toga     = [macr_mvk_obs_toga,tmp_macr_mvk_obs_toga]
            ALD2_obs_toga = [ALD2_obs_toga,tmp_ALD2_obs_toga]
            ch2o_obs_toga = [ch2o_obs_toga,tmp_ch2o_obs_toga]
            benz_obs_toga = [benz_obs_toga,tmp_benz_obs_toga]
            tolu_obs_toga = [tolu_obs_toga,tmp_tolu_obs_toga]

            mbo_obs_toga  = [mbo_obs_toga,tmp_mbo_obs_toga]
            propanal_obs_toga =[propanal_obs_toga,tmp_propanal_obs_toga]

            c2Butenal_obs_toga = [c2Butenal_obs_toga, tmp_c2Butenal_obs_toga]
            t2Butenal_obs_toga = [t2Butenal_obs_toga,tmp_t2Butenal_obs_toga]

            butanal_obs_toga = [butanal_obs_toga, tmp_butanal_obs_toga]

            etbenzene_obs_toga   = [etbenzene_obs_toga, tmp_etbenzene_obs_toga]

            tricyclene_obs_toga  = [tricyclene_obs_toga,tmp_tricyclene_obs_toga]
            apinene_obs_toga     = [apinene_obs_toga,tmp_apinene_obs_toga]
            camphene_obs_toga  = [camphene_obs_toga, tmp_camphene_obs_toga]
            bpinenemyrcene_obs_toga   = [bpinenemyrcene_obs_toga, tmp_bpinenemyrcene_obs_toga]
            limonened3carene_obs_toga = [limonened3carene_obs_toga, tmp_limonened3carene_obs_toga]     
            hcn_obs_toga = [hcn_obs_toga,tmp_hcn_obs_toga] ; pptv to ppbv

            c8h10_obs_toga     = [c8h10_obs_toga,tmp_c8h10_obs_toga]
            Xylenes_obs_toga   = [Xylenes_obs_toga,tmp_Xylenes_obs_toga]
            MPXYLENE_obs_toga  = [MPXYLENE_obs_toga,tmp_MPXYLENE_obs_toga]
            OXYLENE_obs_toga   = [OXYLENE_obs_toga,tmp_OXYLENE_obs_toga]
            EtBenzene_obs_toga = [EtBenzene_obs_toga,tmp_EtBenzene_obs_toga]

            moh_obs_toga = [moh_obs_toga,tmp_moh_obs_toga]

            
            ;;for alk4 and prpe
            ;;alk4
            propane_obs_toga = [propane_obs_toga,tmp_propane_obs_toga]
            butane_obs_toga = [butane_obs_toga,tmp_butane_obs_toga]
            pentane_obs_toga = [pentane_obs_toga,tmp_pentane_obs_toga]
            mepentane_obs_toga = [mepentane_obs_toga,tmp_mepentane_obs_toga]
            hexane_obs_toga = [hexane_obs_toga,tmp_hexane_obs_toga]
            trimepentane_obs_toga = [trimepentane_obs_toga,tmp_trimepentane_obs_toga]
            heptane_obs_toga = [heptane_obs_toga,tmp_heptane_obs_toga]
            octane_obs_toga = [octane_obs_toga,tmp_octane_obs_toga]

            alk4_obs_toga = [alk4_obs_toga,tmp_alk4_obs_toga]
            ;;prpe
            prpe_obs_toga = [prpe_obs_toga,tmp_prpe_obs_toga]
            butene_obs_toga = [butene_obs_toga,tmp_butene_obs_toga] 

    ;;for cloud
            rhum = [rhum,tmp_rhum]
    ;;for temperature
            temperature_obs = [temperature_obs,tmp_temperature_obs]
            
    ;;AWAS
            ;c2h2_obs_awas = [c2h2_obs_awas, tmp_c2h2_obs_awas]
            c2h4_obs_awas   = [c2h4_obs_awas, tmp_c2h4_obs_awas]
            c2h6_obs_awas   = [c2h6_obs_awas, tmp_c2h6_obs_awas]

            ;; PRPE
            c3h6_obs_awas = [c3h6_obs_awas,tmp_c3h6_obs_awas]
            trans2butene_obs_awas = [trans2butene_obs_awas, tmp_trans2butene_obs_awas]
            x1Butene_obs_awas = [x1Butene_obs_awas, tmp_x1Butene_obs_awas]
            cisx2butene_obs_awas = [cisx2butene_obs_awas, tmp_cisx2butene_obs_awas]
            x1Hexene_obs_awas = [x1Hexene_obs_awas, tmp_x1Hexene_obs_awas]

            cyclopentane_obs_awas = [cyclopentane_obs_awas, tmp_cyclopentane_obs_awas]
            cyclohexane_obs_awas = [cyclohexane_obs_awas, tmp_cyclohexane_obs_awas]
            methylcyclohexane_obs_awas = [methylcyclohexane_obs_awas,tmp_methylcyclohexane_obs_awas]
            prpe_obs_awas = [prpe_obs_awas,tmp_prpe_obs_awas]

            ;; ALK4
            isobutane_obs_awas = [isobutane_obs_awas,tmp_isobutane_obs_awas]
            nbutane_obs_awas = [nbutane_obs_awas,tmp_nbutane_obs_awas]
            isopentane_obs_awas = [isopentane_obs_awas,tmp_isopentane_obs_awas]
            nPentane_obs_awas = [nPentane_obs_awas,tmp_nPentane_obs_awas]
            hexane_obs_awas = [hexane_obs_awas,tmp_hexane_obs_awas]
            x24Dimethylpentane_obs_awas = [x24Dimethylpentane_obs_awas,tmp_x24Dimethylpentane_obs_awas]
            x22Dimethylbutane_obs_awas = [x22Dimethylbutane_obs_awas, tmp_x22Dimethylbutane_obs_awas]
            x3Methylpentane_obs_awas = [x3Methylpentane_obs_awas, tmp_x3Methylpentane_obs_awas]
            x23Dimethylpentane_obs_awas = [x23Dimethylpentane_obs_awas, tmp_x23Dimethylpentane_obs_awas]
            x2Methylhexane_obs_awas = [x2Methylhexane_obs_awas, tmp_x2Methylhexane_obs_awas]
            x3Methylhexane_obs_awas = [x3Methylhexane_obs_awas, tmp_x3Methylhexane_obs_awas]
            nHeptan_obs_awas       = [nHeptan_obs_awas, tmp_nHeptan_obs_awas]
            x234Trimethylpentane_obs_awas = [x234Trimethylpentane_obs_awas, tmp_x234Trimethylpentane_obs_awas]
            x2Methylheptane_obs_awas = [x2Methylheptane_obs_awas, tmp_x2Methylheptane_obs_awas]
            x3Methylheptane_obs_awas = [x3Methylheptane_obs_awas, tmp_x3Methylheptane_obs_awas]
            nOctane_obs_awas = [nOctane_obs_awas, tmp_nOctane_obs_awas]
            nNonane_obs_awas= [nNonane_obs_awas, tmp_nNonane_obs_awas]
            undecane_obs_awas = [undecane_obs_awas, tmp_undecane_obs_awas]
            ndecane_obs_awas = [ndecane_obs_awas, tmp_ndecane_obs_awas]
            alk4_obs_awas = [alk4_obs_awas,tmp_alk4_obs_awas]

        ;; CIMS
            hcooh_obs_cims  = [hcooh_obs_cims,tmp_hcooh_obs_cims]
            acta_obs_cims   = [acta_obs_cims,tmp_acta_obs_cims] 
            
            undefine,c130
;;--------------------------------------------------------------------------------------------------------------
;; ============================
;; GEOS-Chem 0.25x0.3125: GFAS + BOTTOME OF THE PLUME TO THE TOP
;; ============================
            ;if keyword_set(nested) then begin
            ;gcfi_gfas_bop_top   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
            ;    'planelog2sav/sen_inj_v2/output_gfas' + '_bop_top' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'
                
            gcfi_gfas_bop_top   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                'planelog2sav/sen_inj_v2/output_gfas' + '_bop_top' + '/mrg5m_wecan_c130_'+dates[n]+'.sav'
            ;endif
            restore, gcfi_gfas_bop_top  
            tmp_utc_gc_gfas_bop_top = gc.utc             
            tmp_co_gc_gfas_bop_top   = gc.co*1e9
            tmp_O3_gc_gfas_bop_top   = gc.o3*1e9
            tmp_pan_gc_gfas_bop_top  = gc.pan*1e9
            tmp_hcho_gc_gfas_bop_top = gc.ch2o*1e9
            tmp_acet_gc_gfas_bop_top = gc.acet*1e9;/3; 3C acetone in ppb
            tmp_benz_gc_gfas_bop_top = gc.benz*1e9;/6; 6C benzene in ppb
    ;        ch3oh_gc= [ch3oh_gc, gc.MOH*1e9]
            tmp_ald2_gc_gfas_bop_top = gc.ald2*1e9;/2; 2C ch3cho in ppb

            tmp_no_gc_gfas_bop_top   = gc.no*1e9
            tmp_no2_gc_gfas_bop_top  = gc.no2*1e9
            tmp_so2_gc_gfas_bop_top  = gc.so2*1e9

            tmp_na = avo * gc.pres * 100./(8.31 * gc.temp) * 1e-6;; air density molec/cm3
            tmp_oh_gc_gfas_bop_top  = gc.oh*tmp_na ;; v/v --> molec/cm3 

            tmp_date_gc_gfas_bop_top = gc.date
            ;tmp_utc_gc_gfas_bop_top  = gc.utc
            tmp_doy_gc_gfas_bop_top  = gc.doy
            tmp_lat_gc_gfas_bop_top  = gc.lat
            tmp_lon_gc_gfas_bop_top  = gc.lon
            tmp_alt_gc_gfas_bop_top  = gc.alt
            tmp_prs_gc_gfas_bop_top  = gc.pres
     ;; add vocs
            tmp_c3h8_gc_gfas_bop_top   = gc.c3h8*1e9;/3; 3C c3h8 in ppb
            tmp_c3h6_gc_gfas_bop_top   = gc.prpe*1e9;/3; 3C PRPE in ppb
            tmp_c2h5oh_gc_gfas_bop_top = gc.eoh*1e9;/2; 2C eoh in ppb
            tmp_c5h8_gc_gfas_bop_top   = gc.isop*1e9;/5; 5C isoprene in ppb
            tmp_c7h8_gc_gfas_bop_top   = gc.tolu*1e9;/7; 7C toluene in ppb
            tmp_dms_gc_gfas_bop_top    = gc.dms*1e9
            tmp_mek_gc_gfas_bop_top    = gc.mek*1e9;/4 ;4C  MEK in ppb
            tmp_c8h10_gc_gfas_bop_top  = gc.xyle*1e9;/8 ; 8C xylenes in ppb

            tmp_acta_gc_gfas_bop_top   = gc.acta*1e9
            tmp_macr_mvk_gc_gfas_bop_top = (gc.MVK+gc.MACR)*1e9  
            tmp_hcooh_gc_gfas_bop_top  = gc.hcooh*1e9

            tmp_mtpa_gc_gfas_bop_top  = gc.mtpa*1e9
            tmp_limo_gc_gfas_bop_top  = gc.limo*1e9
            tmp_mtpo_gc_gfas_bop_top  = gc.mtpo*1e9
        ;;add lumped species
            tmp_alk4_gc_gfas_bop_top   = gc.alk4*1e9;/4.3 ; 4.3C
            tmp_prpe_gc_gfas_bop_top   = gc.prpe*1e9;/3   ; 3C 

 

            co_gc_gfas_bop_top = [co_gc_gfas_bop_top,tmp_co_gc_gfas_bop_top]
            O3_gc_gfas_bop_top   = [o3_gc_gfas_bop_top,tmp_o3_gc_gfas_bop_top]
            pan_gc_gfas_bop_top  = [pan_gc_gfas_bop_top,tmp_pan_gc_gfas_bop_top]
            hcho_gc_gfas_bop_top = [hcho_gc_gfas_bop_top,tmp_hcho_gc_gfas_bop_top]
            acet_gc_gfas_bop_top = [acet_gc_gfas_bop_top,tmp_acet_gc_gfas_bop_top]
            benz_gc_gfas_bop_top = [benz_gc_gfas_bop_top,tmp_benz_gc_gfas_bop_top]

            ald2_gc_gfas_bop_top = [ald2_gc_gfas_bop_top,tmp_ald2_gc_gfas_bop_top]


            no_gc_gfas_bop_top   = [no_gc_gfas_bop_top,tmp_no_gc_gfas_bop_top]
            no2_gc_gfas_bop_top  = [no2_gc_gfas_bop_top,tmp_no2_gc_gfas_bop_top]
            so2_gc_gfas_bop_top  = [so2_gc_gfas_bop_top,tmp_so2_gc_gfas_bop_top]

            oh_gc_gfas_bop_top  = [oh_gc_gfas_bop_top,tmp_oh_gc_gfas_bop_top] ;; v/v --> molec/cm3 

            date_gc_gfas_bop_top = [date_gc_gfas_bop_top,tmp_date_gc_gfas_bop_top]
            ;utc_gc_gfas_bop_top  = [utc_gc_gfas_bop_top,tmp_utc_gc_gfas_bop_top]
            doy_gc_gfas_bop_top  = [doy_gc_gfas_bop_top,tmp_doy_gc_gfas_bop_top]
            lat_gc_gfas_bop_top  = [lat_gc_gfas_bop_top,tmp_lat_gc_gfas_bop_top]
            lon_gc_gfas_bop_top  = [lon_gc_gfas_bop_top,tmp_lon_gc_gfas_bop_top]
            alt_gc_gfas_bop_top  = [alt_gc_gfas_bop_top,tmp_alt_gc_gfas_bop_top]
            prs_gc_gfas_bop_top  = [prs_gc_gfas_bop_top,tmp_prs_gc_gfas_bop_top]
     ;; add vocs
            c3h8_gc_gfas_bop_top   = [c3h8_gc_gfas_bop_top,tmp_c3h8_gc_gfas_bop_top]
            c3h6_gc_gfas_bop_top   = [c3h6_gc_gfas_bop_top,tmp_c3h6_gc_gfas_bop_top]
            c2h5oh_gc_gfas_bop_top = [c2h5oh_gc_gfas_bop_top,tmp_c2h5oh_gc_gfas_bop_top]
            c5h8_gc_gfas_bop_top   = [c5h8_gc_gfas_bop_top,tmp_c5h8_gc_gfas_bop_top]
            c7h8_gc_gfas_bop_top   = [c7h8_gc_gfas_bop_top,tmp_c7h8_gc_gfas_bop_top]
            dms_gc_gfas_bop_top    = [dms_gc_gfas_bop_top,tmp_dms_gc_gfas_bop_top]
            mek_gc_gfas_bop_top    = [mek_gc_gfas_bop_top,tmp_mek_gc_gfas_bop_top]
            c8h10_gc_gfas_bop_top  = [c8h10_gc_gfas_bop_top,tmp_c8h10_gc_gfas_bop_top]

            acta_gc_gfas_bop_top   = [acta_gc_gfas_bop_top,tmp_acta_gc_gfas_bop_top]
            macr_mvk_gc_gfas_bop_top = [macr_mvk_gc_gfas_bop_top,tmp_macr_mvk_gc_gfas_bop_top]  
            hcooh_gc_gfas_bop_top  = [hcooh_gc_gfas_bop_top,tmp_hcooh_gc_gfas_bop_top]

            mtpa_gc_gfas_bop_top  = [mtpa_gc_gfas_bop_top,tmp_mtpa_gc_gfas_bop_top]
            limo_gc_gfas_bop_top  = [limo_gc_gfas_bop_top,tmp_limo_gc_gfas_bop_top]
            mtpo_gc_gfas_bop_top  = [mtpo_gc_gfas_bop_top,tmp_mtpo_gc_gfas_bop_top]
        ;;add lumped species
            alk4_gc_gfas_bop_top  = [alk4_gc_gfas_bop_top,tmp_alk4_gc_gfas_bop_top]
            prpe_gc_gfas_bop_top  = [prpe_gc_gfas_bop_top,tmp_prpe_gc_gfas_bop_top]
    
            undefine,gc


;; ============================
;; GEOS-Chem 0.25x0.3125: GFAS + MAMI
;; ============================
            ;gcfi_gfas_mami   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
            ;    'planelog2sav/sen_inj_v2/output_'+'gfas_mami_v3' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'
            gcfi_gfas_mami   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                'planelog2sav/sen_inj_v2/output_'+'gfas_mami_v3' + '/mrg5m_wecan_c130_'+dates[n]+'.sav'
            restore, gcfi_gfas_mami
            
            tmp_utc_gc_gfas_mami = gc.utc 
            tmp_co_gc_gfas_mami   = gc.co*1e9
            tmp_O3_gc_gfas_mami   = gc.o3*1e9
            tmp_pan_gc_gfas_mami  = gc.pan*1e9
            tmp_hcho_gc_gfas_mami = gc.ch2o*1e9
            tmp_acet_gc_gfas_mami = gc.acet*1e9;/3; 3 acetone in ppb
            tmp_benz_gc_gfas_mami = gc.benz*1e9;/6; 6C benzene in ppb
    ;        ch3oh_gc= [ch3oh_gc, gc.MOH*1e9]
            tmp_ald2_gc_gfas_mami = gc.ald2*1e9;/2; 2C ch3cho in ppb
            
            tmp_no_gc_gfas_mami   = gc.no*1e9
            tmp_no2_gc_gfas_mami  = gc.no2*1e9
            tmp_so2_gc_gfas_mami  = gc.so2*1e9

            tmp_na = avo * gc.pres * 100./(8.31 * gc.temp) * 1e-6;; air density molec/cm3
            tmp_oh_gc_gfas_mami  = gc.oh*tmp_na ;; v/v --> molec/cm3 

            tmp_date_gc_gfas_mami = gc.date
            ;tmp_utc_gc_gfas_mami  = gc.utc
            tmp_doy_gc_gfas_mami  = gc.doy
            tmp_lat_gc_gfas_mami  = gc.lat
            tmp_lon_gc_gfas_mami  = gc.lon
            tmp_alt_gc_gfas_mami  = gc.alt
            tmp_prs_gc_gfas_mami  = gc.pres
     ;; add vocs
            tmp_c3h8_gc_gfas_mami   = gc.c3h8*1e9;/3; 3C c3h8 in ppb
            tmp_c3h6_gc_gfas_mami   = gc.prpe*1e9;/3; 3C PRPE in ppb
            tmp_c2h5oh_gc_gfas_mami = gc.eoh*1e9;/2; 2C eoh in ppb
            tmp_c5h8_gc_gfas_mami   = gc.isop*1e9;/5; 5C isoprene in ppb
            tmp_c7h8_gc_gfas_mami   = gc.tolu*1e9;/7; 7C toluene in ppb
            tmp_dms_gc_gfas_mami    = gc.dms*1e9
            tmp_mek_gc_gfas_mami    = gc.mek*1e9;/4 ;4C  MEK in ppb
            tmp_c8h10_gc_gfas_mami  = gc.xyle*1e9;/8 ; 8C xylenes in ppb

            tmp_acta_gc_gfas_mami   = gc.acta*1e9
            tmp_macr_mvk_gc_gfas_mami = (gc.MVK+gc.MACR)*1e9  
            tmp_hcooh_gc_gfas_mami  = gc.hcooh*1e9

            tmp_mtpa_gc_gfas_mami  = gc.mtpa*1e9
            tmp_limo_gc_gfas_mami  = gc.limo*1e9
            tmp_mtpo_gc_gfas_mami  = gc.mtpo*1e9
        ;;add lumped species
            tmp_alk4_gc_gfas_mami   = gc.alk4*1e9
            tmp_prpe_gc_gfas_mami   = gc.prpe*1e9

            co_gc_gfas_mami = [co_gc_gfas_mami,tmp_co_gc_gfas_mami]
            O3_gc_gfas_mami   = [o3_gc_gfas_mami,tmp_o3_gc_gfas_mami]
            pan_gc_gfas_mami  = [pan_gc_gfas_mami,tmp_pan_gc_gfas_mami]
            hcho_gc_gfas_mami = [hcho_gc_gfas_mami,tmp_hcho_gc_gfas_mami]
            acet_gc_gfas_mami = [acet_gc_gfas_mami,tmp_acet_gc_gfas_mami]
            benz_gc_gfas_mami = [benz_gc_gfas_mami,tmp_benz_gc_gfas_mami]

            ald2_gc_gfas_mami = [ald2_gc_gfas_mami,tmp_ald2_gc_gfas_mami]


            no_gc_gfas_mami   = [no_gc_gfas_mami,tmp_no_gc_gfas_mami]
            no2_gc_gfas_mami  = [no2_gc_gfas_mami,tmp_no2_gc_gfas_mami]
            so2_gc_gfas_mami  = [so2_gc_gfas_mami,tmp_so2_gc_gfas_mami]

            oh_gc_gfas_mami  = [oh_gc_gfas_mami,tmp_oh_gc_gfas_mami] ;; v/v --> molec/cm3 

            date_gc_gfas_mami = [date_gc_gfas_mami,tmp_date_gc_gfas_mami]
            ;utc_gc_gfas_mami  = [utc_gc_gfas_mami,tmp_utc_gc_gfas_mami]
            doy_gc_gfas_mami  = [doy_gc_gfas_mami,tmp_doy_gc_gfas_mami]
            lat_gc_gfas_mami  = [lat_gc_gfas_mami,tmp_lat_gc_gfas_mami]
            lon_gc_gfas_mami  = [lon_gc_gfas_mami,tmp_lon_gc_gfas_mami]
            alt_gc_gfas_mami  = [alt_gc_gfas_mami,tmp_alt_gc_gfas_mami]
            prs_gc_gfas_mami  = [prs_gc_gfas_mami,tmp_prs_gc_gfas_mami]
     ;; add vocs
            c3h8_gc_gfas_mami   = [c3h8_gc_gfas_mami,tmp_c3h8_gc_gfas_mami]
            c3h6_gc_gfas_mami   = [c3h6_gc_gfas_mami,tmp_c3h6_gc_gfas_mami]
            c2h5oh_gc_gfas_mami = [c2h5oh_gc_gfas_mami,tmp_c2h5oh_gc_gfas_mami]
            c5h8_gc_gfas_mami   = [c5h8_gc_gfas_mami,tmp_c5h8_gc_gfas_mami]
            c7h8_gc_gfas_mami   = [c7h8_gc_gfas_mami,tmp_c7h8_gc_gfas_mami]
            dms_gc_gfas_mami    = [dms_gc_gfas_mami,tmp_dms_gc_gfas_mami]
            mek_gc_gfas_mami    = [mek_gc_gfas_mami,tmp_mek_gc_gfas_mami]
            c8h10_gc_gfas_mami  = [c8h10_gc_gfas_mami,tmp_c8h10_gc_gfas_mami]

            acta_gc_gfas_mami   = [acta_gc_gfas_mami,tmp_acta_gc_gfas_mami]
            macr_mvk_gc_gfas_mami = [macr_mvk_gc_gfas_mami,tmp_macr_mvk_gc_gfas_mami]  
            hcooh_gc_gfas_mami  = [hcooh_gc_gfas_mami,tmp_hcooh_gc_gfas_mami]

            mtpa_gc_gfas_mami  = [mtpa_gc_gfas_mami,tmp_mtpa_gc_gfas_mami]
            limo_gc_gfas_mami  = [limo_gc_gfas_mami,tmp_limo_gc_gfas_mami]
            mtpo_gc_gfas_mami  = [mtpo_gc_gfas_mami,tmp_mtpo_gc_gfas_mami]
        ;;add lumped species
            alk4_gc_gfas_mami  = [alk4_gc_gfas_mami,tmp_alk4_gc_gfas_mami]
            prpe_gc_gfas_mami  = [prpe_gc_gfas_mami,tmp_prpe_gc_gfas_mami]
    
            undefine,gc

;; ============================
;; GEOS-Chem 0.25x0.3125: GFAS + NOBB
;; ============================
            ;gcfi_gfas_nobb   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
            ;    'planelog2sav/sen_inj_v2/output_gfas_'+'nobb' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'
            gcfi_gfas_nobb   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                'planelog2sav/sen_inj_v2/output_gfas_'+'nobb' + '/mrg5m_wecan_c130_'+dates[n]+'.sav'
            restore, gcfi_gfas_nobb       
            tmp_utc_gc_gfas_nobb = gc.utc 
            tmp_co_gc_gfas_nobb   = gc.co*1e9
            tmp_O3_gc_gfas_nobb   = gc.o3*1e9
            tmp_pan_gc_gfas_nobb  = gc.pan*1e9
            tmp_hcho_gc_gfas_nobb = gc.ch2o*1e9
            tmp_acet_gc_gfas_nobb = gc.acet*1e9;/3; 3 acetone in ppb
            tmp_benz_gc_gfas_nobb = gc.benz*1e9;/6; 6C benzene in ppb
    ;        ch3oh_gc= [ch3oh_gc, gc.MOH*1e9]
            tmp_ald2_gc_gfas_nobb = gc.ald2*1e9;/2; 2C ch3cho in ppb

            tmp_no_gc_gfas_nobb   = gc.no*1e9
            tmp_no2_gc_gfas_nobb  = gc.no2*1e9
            tmp_so2_gc_gfas_nobb  = gc.so2*1e9

            tmp_na = avo * gc.pres * 100./(8.31 * gc.temp) * 1e-6;; air density molec/cm3
            tmp_oh_gc_gfas_nobb  = gc.oh*tmp_na ;; v/v --> molec/cm3 

            tmp_date_gc_gfas_nobb = gc.date
            ;tmp_utc_gc_gfas_nobb  = gc.utc
            tmp_doy_gc_gfas_nobb  = gc.doy
            tmp_lat_gc_gfas_nobb  = gc.lat
            tmp_lon_gc_gfas_nobb  = gc.lon
            tmp_alt_gc_gfas_nobb  = gc.alt
            tmp_prs_gc_gfas_nobb  = gc.pres
     ;; add vocs
            tmp_c3h8_gc_gfas_nobb   = gc.c3h8*1e9;/3; 3C c3h8 in ppb
            tmp_c3h6_gc_gfas_nobb   = gc.prpe*1e9;/3; 3C PRPE in ppb
            tmp_c2h5oh_gc_gfas_nobb = gc.eoh*1e9;/2; 2C eoh in ppb
            tmp_c5h8_gc_gfas_nobb   = gc.isop*1e9;/5; 5C isoprene in ppb
            tmp_c7h8_gc_gfas_nobb   = gc.tolu*1e9;/7; 7C toluene in ppb
            tmp_dms_gc_gfas_nobb    = gc.dms*1e9
            tmp_mek_gc_gfas_nobb    = gc.mek*1e9;/4 ;4C  MEK in ppb
            tmp_c8h10_gc_gfas_nobb  = gc.xyle*1e9;/8 ; 8C xylenes in ppb


            tmp_acta_gc_gfas_nobb   = gc.acta*1e9
            tmp_macr_mvk_gc_gfas_nobb = (gc.MVK+gc.MACR)*1e9  
            tmp_hcooh_gc_gfas_nobb  = gc.hcooh*1e9

            tmp_mtpa_gc_gfas_nobb  = gc.mtpa*1e9
            tmp_limo_gc_gfas_nobb  = gc.limo*1e9
            tmp_mtpo_gc_gfas_nobb  = gc.mtpo*1e9
        ;;add lumped species
            tmp_alk4_gc_gfas_nobb   = gc.alk4*1e9
            tmp_prpe_gc_gfas_nobb   = gc.prpe*1e9

            co_gc_gfas_nobb = [co_gc_gfas_nobb,tmp_co_gc_gfas_nobb]
            O3_gc_gfas_nobb   = [o3_gc_gfas_nobb,tmp_o3_gc_gfas_nobb]
            pan_gc_gfas_nobb  = [pan_gc_gfas_nobb,tmp_pan_gc_gfas_nobb]
            hcho_gc_gfas_nobb = [hcho_gc_gfas_nobb,tmp_hcho_gc_gfas_nobb]
            acet_gc_gfas_nobb = [acet_gc_gfas_nobb,tmp_acet_gc_gfas_nobb]
            benz_gc_gfas_nobb = [benz_gc_gfas_nobb,tmp_benz_gc_gfas_nobb]

            ald2_gc_gfas_nobb = [ald2_gc_gfas_nobb,tmp_ald2_gc_gfas_nobb]


            no_gc_gfas_nobb   = [no_gc_gfas_nobb,tmp_no_gc_gfas_nobb]
            no2_gc_gfas_nobb  = [no2_gc_gfas_nobb,tmp_no2_gc_gfas_nobb]
            so2_gc_gfas_nobb  = [so2_gc_gfas_nobb,tmp_so2_gc_gfas_nobb]

            oh_gc_gfas_nobb  = [oh_gc_gfas_nobb,tmp_oh_gc_gfas_nobb] ;; v/v --> molec/cm3 

            date_gc_gfas_nobb = [date_gc_gfas_nobb,tmp_date_gc_gfas_nobb]
            ;utc_gc_gfas_nobb  = [utc_gc_gfas_nobb,tmp_utc_gc_gfas_nobb]
            doy_gc_gfas_nobb  = [doy_gc_gfas_nobb,tmp_doy_gc_gfas_nobb]
            lat_gc_gfas_nobb  = [lat_gc_gfas_nobb,tmp_lat_gc_gfas_nobb]
            lon_gc_gfas_nobb  = [lon_gc_gfas_nobb,tmp_lon_gc_gfas_nobb]
            alt_gc_gfas_nobb  = [alt_gc_gfas_nobb,tmp_alt_gc_gfas_nobb]
            prs_gc_gfas_nobb  = [prs_gc_gfas_nobb,tmp_prs_gc_gfas_nobb]
     ;; add vocs
            c3h8_gc_gfas_nobb   = [c3h8_gc_gfas_nobb,tmp_c3h8_gc_gfas_nobb]
            c3h6_gc_gfas_nobb   = [c3h6_gc_gfas_nobb,tmp_c3h6_gc_gfas_nobb]
            c2h5oh_gc_gfas_nobb = [c2h5oh_gc_gfas_nobb,tmp_c2h5oh_gc_gfas_nobb]
            c5h8_gc_gfas_nobb   = [c5h8_gc_gfas_nobb,tmp_c5h8_gc_gfas_nobb]
            c7h8_gc_gfas_nobb   = [c7h8_gc_gfas_nobb,tmp_c7h8_gc_gfas_nobb]
            dms_gc_gfas_nobb    = [dms_gc_gfas_nobb,tmp_dms_gc_gfas_nobb]
            mek_gc_gfas_nobb    = [mek_gc_gfas_nobb,tmp_mek_gc_gfas_nobb]
            c8h10_gc_gfas_nobb  = [c8h10_gc_gfas_nobb,tmp_c8h10_gc_gfas_nobb]

            acta_gc_gfas_nobb   = [acta_gc_gfas_nobb,tmp_acta_gc_gfas_nobb]
            macr_mvk_gc_gfas_nobb = [macr_mvk_gc_gfas_nobb,tmp_macr_mvk_gc_gfas_nobb]  
            hcooh_gc_gfas_nobb  = [hcooh_gc_gfas_nobb,tmp_hcooh_gc_gfas_nobb]

            mtpa_gc_gfas_nobb  = [mtpa_gc_gfas_nobb,tmp_mtpa_gc_gfas_nobb]
            limo_gc_gfas_nobb  = [limo_gc_gfas_nobb,tmp_limo_gc_gfas_nobb]
            mtpo_gc_gfas_nobb  = [mtpo_gc_gfas_nobb,tmp_mtpo_gc_gfas_nobb]
        ;;add lumped species
            alk4_gc_gfas_nobb  = [alk4_gc_gfas_nobb,tmp_alk4_gc_gfas_nobb]
            prpe_gc_gfas_nobb  = [prpe_gc_gfas_nobb,tmp_prpe_gc_gfas_nobb]
    
            undefine,gc

;; ============================
;; GEOS-Chem 0.25x0.3125: GFAS + PBL65_FT35
;; ============================
            ;gcfi_gfas_pbl65_ft35   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
            ;    'planelog2sav/sen_inj_v2/output_'+'gfas_65PBL_35FT' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'
            gcfi_gfas_pbl65_ft35   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                'planelog2sav/sen_inj_v2/output_'+'gfas_65PBL_35FT' + '/mrg5m_wecan_c130_'+dates[n]+'.sav'
            restore, gcfi_gfas_pbl65_ft35       
            tmp_utc_gc_gfas_pbl65_ft35 = gc.utc 
            tmp_co_gc_gfas_pbl65_ft35   = gc.co*1e9
            tmp_O3_gc_gfas_pbl65_ft35   = gc.o3*1e9
            tmp_pan_gc_gfas_pbl65_ft35  = gc.pan*1e9
            tmp_hcho_gc_gfas_pbl65_ft35 = gc.ch2o*1e9
            tmp_acet_gc_gfas_pbl65_ft35 = gc.acet*1e9;/3; 3 acetone in ppb
            tmp_benz_gc_gfas_pbl65_ft35 = gc.benz*1e9;/6; 6C benzene in ppb
    ;        ch3oh_gc= [ch3oh_gc, gc.MOH*1e9]
            tmp_ald2_gc_gfas_pbl65_ft35 = gc.ald2*1e9;/2; 2C ch3cho in ppb

            tmp_no_gc_gfas_pbl65_ft35   = gc.no*1e9
            tmp_no2_gc_gfas_pbl65_ft35  = gc.no2*1e9
            tmp_so2_gc_gfas_pbl65_ft35  = gc.so2*1e9

            tmp_na = avo * gc.pres * 100./(8.31 * gc.temp) * 1e-6;; air density molec/cm3
            tmp_oh_gc_gfas_pbl65_ft35  = gc.oh*tmp_na ;; v/v --> molec/cm3 

            tmp_date_gc_gfas_pbl65_ft35 = gc.date
            ;tmp_utc_gc_gfas_pbl65_ft35  = gc.utc
            tmp_doy_gc_gfas_pbl65_ft35  = gc.doy
            tmp_lat_gc_gfas_pbl65_ft35  = gc.lat
            tmp_lon_gc_gfas_pbl65_ft35  = gc.lon
            tmp_alt_gc_gfas_pbl65_ft35  = gc.alt
            tmp_prs_gc_gfas_pbl65_ft35  = gc.pres
     ;; add vocs
            tmp_c3h8_gc_gfas_pbl65_ft35   = gc.c3h8*1e9;/3; 3C c3h8 in ppb
            tmp_c3h6_gc_gfas_pbl65_ft35   = gc.prpe*1e9;/3; 3C PRPE in ppb
            tmp_c2h5oh_gc_gfas_pbl65_ft35 = gc.eoh*1e9;/2; 2C eoh in ppb
            tmp_c5h8_gc_gfas_pbl65_ft35   = gc.isop*1e9;/5; 5C isoprene in ppb
            tmp_c7h8_gc_gfas_pbl65_ft35   = gc.tolu*1e9;/7; 7C toluene in ppb
            tmp_dms_gc_gfas_pbl65_ft35    = gc.dms*1e9
            tmp_mek_gc_gfas_pbl65_ft35    = gc.mek*1e9;/4 ;4C  MEK in ppb
            tmp_c8h10_gc_gfas_pbl65_ft35  = gc.xyle*1e9;/8 ; 8C xylenes in ppb


            tmp_acta_gc_gfas_pbl65_ft35   = gc.acta*1e9
            tmp_macr_mvk_gc_gfas_pbl65_ft35 = (gc.MVK+gc.MACR)*1e9  
            tmp_hcooh_gc_gfas_pbl65_ft35  = gc.hcooh*1e9

            tmp_mtpa_gc_gfas_pbl65_ft35  = gc.mtpa*1e9
            tmp_limo_gc_gfas_pbl65_ft35  = gc.limo*1e9
            tmp_mtpo_gc_gfas_pbl65_ft35  = gc.mtpo*1e9
        ;;add lumped species
            tmp_alk4_gc_gfas_pbl65_ft35   = gc.alk4*1e9
            tmp_prpe_gc_gfas_pbl65_ft35   = gc.prpe*1e9

            co_gc_gfas_pbl65_ft35 = [co_gc_gfas_pbl65_ft35,tmp_co_gc_gfas_pbl65_ft35]
            O3_gc_gfas_pbl65_ft35   = [o3_gc_gfas_pbl65_ft35,tmp_o3_gc_gfas_pbl65_ft35]
            pan_gc_gfas_pbl65_ft35  = [pan_gc_gfas_pbl65_ft35,tmp_pan_gc_gfas_pbl65_ft35]
            hcho_gc_gfas_pbl65_ft35 = [hcho_gc_gfas_pbl65_ft35,tmp_hcho_gc_gfas_pbl65_ft35]
            acet_gc_gfas_pbl65_ft35 = [acet_gc_gfas_pbl65_ft35,tmp_acet_gc_gfas_pbl65_ft35]
            benz_gc_gfas_pbl65_ft35 = [benz_gc_gfas_pbl65_ft35,tmp_benz_gc_gfas_pbl65_ft35]

            ald2_gc_gfas_pbl65_ft35 = [ald2_gc_gfas_pbl65_ft35,tmp_ald2_gc_gfas_pbl65_ft35]


            no_gc_gfas_pbl65_ft35   = [no_gc_gfas_pbl65_ft35,tmp_no_gc_gfas_pbl65_ft35]
            no2_gc_gfas_pbl65_ft35  = [no2_gc_gfas_pbl65_ft35,tmp_no2_gc_gfas_pbl65_ft35]
            so2_gc_gfas_pbl65_ft35  = [so2_gc_gfas_pbl65_ft35,tmp_so2_gc_gfas_pbl65_ft35]

            oh_gc_gfas_pbl65_ft35  = [oh_gc_gfas_pbl65_ft35,tmp_oh_gc_gfas_pbl65_ft35] ;; v/v --> molec/cm3 

            date_gc_gfas_pbl65_ft35 = [date_gc_gfas_pbl65_ft35,tmp_date_gc_gfas_pbl65_ft35]
            ;utc_gc_gfas_pbl65_ft35  = [utc_gc_gfas_pbl65_ft35,tmp_utc_gc_gfas_pbl65_ft35]
            doy_gc_gfas_pbl65_ft35  = [doy_gc_gfas_pbl65_ft35,tmp_doy_gc_gfas_pbl65_ft35]
            lat_gc_gfas_pbl65_ft35  = [lat_gc_gfas_pbl65_ft35,tmp_lat_gc_gfas_pbl65_ft35]
            lon_gc_gfas_pbl65_ft35  = [lon_gc_gfas_pbl65_ft35,tmp_lon_gc_gfas_pbl65_ft35]
            alt_gc_gfas_pbl65_ft35  = [alt_gc_gfas_pbl65_ft35,tmp_alt_gc_gfas_pbl65_ft35]
            prs_gc_gfas_pbl65_ft35  = [prs_gc_gfas_pbl65_ft35,tmp_prs_gc_gfas_pbl65_ft35]
     ;; add vocs
            c3h8_gc_gfas_pbl65_ft35   = [c3h8_gc_gfas_pbl65_ft35,tmp_c3h8_gc_gfas_pbl65_ft35]
            c3h6_gc_gfas_pbl65_ft35   = [c3h6_gc_gfas_pbl65_ft35,tmp_c3h6_gc_gfas_pbl65_ft35]
            c2h5oh_gc_gfas_pbl65_ft35 = [c2h5oh_gc_gfas_pbl65_ft35,tmp_c2h5oh_gc_gfas_pbl65_ft35]
            c5h8_gc_gfas_pbl65_ft35   = [c5h8_gc_gfas_pbl65_ft35,tmp_c5h8_gc_gfas_pbl65_ft35]
            c7h8_gc_gfas_pbl65_ft35   = [c7h8_gc_gfas_pbl65_ft35,tmp_c7h8_gc_gfas_pbl65_ft35]
            dms_gc_gfas_pbl65_ft35    = [dms_gc_gfas_pbl65_ft35,tmp_dms_gc_gfas_pbl65_ft35]
            mek_gc_gfas_pbl65_ft35    = [mek_gc_gfas_pbl65_ft35,tmp_mek_gc_gfas_pbl65_ft35]
            c8h10_gc_gfas_pbl65_ft35  = [c8h10_gc_gfas_pbl65_ft35,tmp_c8h10_gc_gfas_pbl65_ft35]

            acta_gc_gfas_pbl65_ft35   = [acta_gc_gfas_pbl65_ft35,tmp_acta_gc_gfas_pbl65_ft35]
            macr_mvk_gc_gfas_pbl65_ft35 = [macr_mvk_gc_gfas_pbl65_ft35,tmp_macr_mvk_gc_gfas_pbl65_ft35]  
            hcooh_gc_gfas_pbl65_ft35  = [hcooh_gc_gfas_pbl65_ft35,tmp_hcooh_gc_gfas_pbl65_ft35]

            mtpa_gc_gfas_pbl65_ft35  = [mtpa_gc_gfas_pbl65_ft35,tmp_mtpa_gc_gfas_pbl65_ft35]
            limo_gc_gfas_pbl65_ft35  = [limo_gc_gfas_pbl65_ft35,tmp_limo_gc_gfas_pbl65_ft35]
            mtpo_gc_gfas_pbl65_ft35  = [mtpo_gc_gfas_pbl65_ft35,tmp_mtpo_gc_gfas_pbl65_ft35]
        ;;add lumped species
            alk4_gc_gfas_pbl65_ft35  = [alk4_gc_gfas_pbl65_ft35,tmp_alk4_gc_gfas_pbl65_ft35]
            prpe_gc_gfas_pbl65_ft35  = [prpe_gc_gfas_pbl65_ft35,tmp_prpe_gc_gfas_pbl65_ft35]
    
            undefine,gc

;; ============================
;; GEOS-Chem 0.25x0.3125: GFAS + SF_MAMI
;; ============================
            ;gcfi_gfas_sf_mami   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
            ;    'planelog2sav/sen_inj_v2/output_'+'gfas_sf_mami' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'
            gcfi_gfas_sf_mami   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                'planelog2sav/sen_inj_v2/output_'+'gfas_sf_mami' + '/mrg5m_wecan_c130_'+dates[n]+'.sav'
            restore, gcfi_gfas_sf_mami       
            tmp_utc_gc_gfas_sf_mami = gc.utc 
            tmp_co_gc_gfas_sf_mami   = gc.co*1e9
            tmp_O3_gc_gfas_sf_mami   = gc.o3*1e9
            tmp_pan_gc_gfas_sf_mami  = gc.pan*1e9
            tmp_hcho_gc_gfas_sf_mami = gc.ch2o*1e9
            tmp_acet_gc_gfas_sf_mami = gc.acet*1e9;/3; 3 acetone in ppb
            tmp_benz_gc_gfas_sf_mami = gc.benz*1e9;/6; 6C benzene in ppb
    ;        ch3oh_gc= [ch3oh_gc, gc.MOH*1e9]
            tmp_ald2_gc_gfas_sf_mami = gc.ald2*1e9;/2; 2C ch3cho in ppb

            tmp_no_gc_gfas_sf_mami   = gc.no*1e9
            tmp_no2_gc_gfas_sf_mami  = gc.no2*1e9
            tmp_so2_gc_gfas_sf_mami  = gc.so2*1e9

            tmp_na = avo * gc.pres * 100./(8.31 * gc.temp) * 1e-6;; air density molec/cm3
            tmp_oh_gc_gfas_sf_mami  = gc.oh*tmp_na ;; v/v --> molec/cm3 

            tmp_date_gc_gfas_sf_mami = gc.date
            ;tmp_utc_gc_gfas_sf_mami  = gc.utc
            tmp_doy_gc_gfas_sf_mami  = gc.doy
            tmp_lat_gc_gfas_sf_mami  = gc.lat
            tmp_lon_gc_gfas_sf_mami  = gc.lon
            tmp_alt_gc_gfas_sf_mami  = gc.alt
            tmp_prs_gc_gfas_sf_mami  = gc.pres
     ;; add vocs
            tmp_c3h8_gc_gfas_sf_mami   = gc.c3h8*1e9;/3; 3C c3h8 in ppb
            tmp_c3h6_gc_gfas_sf_mami   = gc.prpe*1e9;/3; 3C PRPE in ppb
            tmp_c2h5oh_gc_gfas_sf_mami = gc.eoh*1e9;/2; 2C eoh in ppb
            tmp_c5h8_gc_gfas_sf_mami   = gc.isop*1e9;/5; 5C isoprene in ppb
            tmp_c7h8_gc_gfas_sf_mami   = gc.tolu*1e9;/7; 7C toluene in ppb
            tmp_dms_gc_gfas_sf_mami    = gc.dms*1e9
            tmp_mek_gc_gfas_sf_mami    = gc.mek*1e9;/4 ;4C  MEK in ppb
            tmp_c8h10_gc_gfas_sf_mami  = gc.xyle*1e9;/8 ; 8C xylenes in ppb


            tmp_acta_gc_gfas_sf_mami   = gc.acta*1e9
            tmp_macr_mvk_gc_gfas_sf_mami = (gc.MVK+gc.MACR)*1e9  
            tmp_hcooh_gc_gfas_sf_mami  = gc.hcooh*1e9

            tmp_mtpa_gc_gfas_sf_mami  = gc.mtpa*1e9
            tmp_limo_gc_gfas_sf_mami  = gc.limo*1e9
            tmp_mtpo_gc_gfas_sf_mami  = gc.mtpo*1e9
        ;;add lumped species
            tmp_alk4_gc_gfas_sf_mami   = gc.alk4*1e9
            tmp_prpe_gc_gfas_sf_mami   = gc.prpe*1e9

            co_gc_gfas_sf_mami = [co_gc_gfas_sf_mami,tmp_co_gc_gfas_sf_mami]
            O3_gc_gfas_sf_mami   = [o3_gc_gfas_sf_mami,tmp_o3_gc_gfas_sf_mami]
            pan_gc_gfas_sf_mami  = [pan_gc_gfas_sf_mami,tmp_pan_gc_gfas_sf_mami]
            hcho_gc_gfas_sf_mami = [hcho_gc_gfas_sf_mami,tmp_hcho_gc_gfas_sf_mami]
            acet_gc_gfas_sf_mami = [acet_gc_gfas_sf_mami,tmp_acet_gc_gfas_sf_mami]
            benz_gc_gfas_sf_mami = [benz_gc_gfas_sf_mami,tmp_benz_gc_gfas_sf_mami]

            ald2_gc_gfas_sf_mami = [ald2_gc_gfas_sf_mami,tmp_ald2_gc_gfas_sf_mami]


            no_gc_gfas_sf_mami   = [no_gc_gfas_sf_mami,tmp_no_gc_gfas_sf_mami]
            no2_gc_gfas_sf_mami  = [no2_gc_gfas_sf_mami,tmp_no2_gc_gfas_sf_mami]
            so2_gc_gfas_sf_mami  = [so2_gc_gfas_sf_mami,tmp_so2_gc_gfas_sf_mami]

            oh_gc_gfas_sf_mami  = [oh_gc_gfas_sf_mami,tmp_oh_gc_gfas_sf_mami] ;; v/v --> molec/cm3 

            date_gc_gfas_sf_mami = [date_gc_gfas_sf_mami,tmp_date_gc_gfas_sf_mami]
            ;utc_gc_gfas_sf_mami  = [utc_gc_gfas_sf_mami,tmp_utc_gc_gfas_sf_mami]
            doy_gc_gfas_sf_mami  = [doy_gc_gfas_sf_mami,tmp_doy_gc_gfas_sf_mami]
            lat_gc_gfas_sf_mami  = [lat_gc_gfas_sf_mami,tmp_lat_gc_gfas_sf_mami]
            lon_gc_gfas_sf_mami  = [lon_gc_gfas_sf_mami,tmp_lon_gc_gfas_sf_mami]
            alt_gc_gfas_sf_mami  = [alt_gc_gfas_sf_mami,tmp_alt_gc_gfas_sf_mami]
            prs_gc_gfas_sf_mami  = [prs_gc_gfas_sf_mami,tmp_prs_gc_gfas_sf_mami]
     ;; add vocs
            c3h8_gc_gfas_sf_mami   = [c3h8_gc_gfas_sf_mami,tmp_c3h8_gc_gfas_sf_mami]
            c3h6_gc_gfas_sf_mami   = [c3h6_gc_gfas_sf_mami,tmp_c3h6_gc_gfas_sf_mami]
            c2h5oh_gc_gfas_sf_mami = [c2h5oh_gc_gfas_sf_mami,tmp_c2h5oh_gc_gfas_sf_mami]
            c5h8_gc_gfas_sf_mami   = [c5h8_gc_gfas_sf_mami,tmp_c5h8_gc_gfas_sf_mami]
            c7h8_gc_gfas_sf_mami   = [c7h8_gc_gfas_sf_mami,tmp_c7h8_gc_gfas_sf_mami]
            dms_gc_gfas_sf_mami    = [dms_gc_gfas_sf_mami,tmp_dms_gc_gfas_sf_mami]
            mek_gc_gfas_sf_mami    = [mek_gc_gfas_sf_mami,tmp_mek_gc_gfas_sf_mami]
            c8h10_gc_gfas_sf_mami  = [c8h10_gc_gfas_sf_mami,tmp_c8h10_gc_gfas_sf_mami]

            acta_gc_gfas_sf_mami   = [acta_gc_gfas_sf_mami,tmp_acta_gc_gfas_sf_mami]
            macr_mvk_gc_gfas_sf_mami = [macr_mvk_gc_gfas_sf_mami,tmp_macr_mvk_gc_gfas_sf_mami]  
            hcooh_gc_gfas_sf_mami  = [hcooh_gc_gfas_sf_mami,tmp_hcooh_gc_gfas_sf_mami]

            mtpa_gc_gfas_sf_mami  = [mtpa_gc_gfas_sf_mami,tmp_mtpa_gc_gfas_sf_mami]
            limo_gc_gfas_sf_mami  = [limo_gc_gfas_sf_mami,tmp_limo_gc_gfas_sf_mami]
            mtpo_gc_gfas_sf_mami  = [mtpo_gc_gfas_sf_mami,tmp_mtpo_gc_gfas_sf_mami]
        ;;add lumped species
            alk4_gc_gfas_sf_mami  = [alk4_gc_gfas_sf_mami,tmp_alk4_gc_gfas_sf_mami]
            prpe_gc_gfas_sf_mami  = [prpe_gc_gfas_sf_mami,tmp_prpe_gc_gfas_sf_mami]
    
            undefine,gc
       
       
;; ============================
;; GEOS-Chem 0.25x0.3125: GFAS + SURFACE
;; ============================
            ;gcfi_gfas_surface   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
            ;    'planelog2sav/sen_inj_v2/output_'+'gfas_surface' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'
            gcfi_gfas_surface   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                'planelog2sav/sen_inj_v2/output_'+'gfas_surface' + '/mrg5m_wecan_c130_'+dates[n]+'.sav'
            restore, gcfi_gfas_surface       
            tmp_utc_gc_gfas_surface = gc.utc 
            tmp_co_gc_gfas_surface   = gc.co*1e9
            tmp_O3_gc_gfas_surface   = gc.o3*1e9
            tmp_pan_gc_gfas_surface  = gc.pan*1e9
            tmp_hcho_gc_gfas_surface = gc.ch2o*1e9
            tmp_acet_gc_gfas_surface = gc.acet*1e9;/3; 3 acetone in ppb
            tmp_benz_gc_gfas_surface = gc.benz*1e9;/6; 6C benzene in ppb
    ;        ch3oh_gc= [ch3oh_gc, gc.MOH*1e9]
            tmp_ald2_gc_gfas_surface = gc.ald2*1e9;/2; 2C ch3cho in ppb

            tmp_no_gc_gfas_surface   = gc.no*1e9
            tmp_no2_gc_gfas_surface  = gc.no2*1e9
            tmp_so2_gc_gfas_surface  = gc.so2*1e9

            tmp_na = avo * gc.pres * 100./(8.31 * gc.temp) * 1e-6;; air density molec/cm3
            tmp_oh_gc_gfas_surface  = gc.oh*tmp_na ;; v/v --> molec/cm3 

            tmp_date_gc_gfas_surface = gc.date
            ;tmp_utc_gc_gfas_surface  = gc.utc
            tmp_doy_gc_gfas_surface  = gc.doy
            tmp_lat_gc_gfas_surface  = gc.lat
            tmp_lon_gc_gfas_surface  = gc.lon
            tmp_alt_gc_gfas_surface  = gc.alt
            tmp_prs_gc_gfas_surface  = gc.pres
     ;; add vocs
            tmp_c3h8_gc_gfas_surface   = gc.c3h8*1e9;/3; 3C c3h8 in ppb
            tmp_c3h6_gc_gfas_surface   = gc.prpe*1e9;/3; 3C PRPE in ppb
            tmp_c2h5oh_gc_gfas_surface = gc.eoh*1e9;/2; 2C eoh in ppb
            tmp_c5h8_gc_gfas_surface   = gc.isop*1e9;/5; 5C isoprene in ppb
            tmp_c7h8_gc_gfas_surface   = gc.tolu*1e9;/7; 7C toluene in ppb
            tmp_dms_gc_gfas_surface    = gc.dms*1e9
            tmp_mek_gc_gfas_surface    = gc.mek*1e9;/4 ;4C  MEK in ppb
            tmp_c8h10_gc_gfas_surface  = gc.xyle*1e9;/8 ; 8C xylenes in ppb


            tmp_acta_gc_gfas_surface   = gc.acta*1e9
            tmp_macr_mvk_gc_gfas_surface = (gc.MVK+gc.MACR)*1e9  
            tmp_hcooh_gc_gfas_surface  = gc.hcooh*1e9

            tmp_mtpa_gc_gfas_surface  = gc.mtpa*1e9
            tmp_limo_gc_gfas_surface  = gc.limo*1e9
            tmp_mtpo_gc_gfas_surface  = gc.mtpo*1e9
        ;;add lumped species
            tmp_alk4_gc_gfas_surface   = gc.alk4*1e9
            tmp_prpe_gc_gfas_surface   = gc.prpe*1e9

    ;; concert hhmm(utc) into ss(utc)      
            ;print,tmp_utc_gc_gfas_surface
            hh = floor(tmp_utc_gc_gfas_surface/100)
            ind = where(hh gt 12, ct)
            if ct gt 0 then hh[ind] = hh[ind]
            ind = where(hh lt 12, ct)
            if ct gt 0 then hh[ind] = hh[ind] + 24
            mm = tmp_utc_gc_gfas_surface- floor(tmp_utc_gc_gfas_surface/100)*100
            tmp_utc_gc_gfas_surface = float(hh)*60*60+float(mm)*60


            co_gc_gfas_surface = [co_gc_gfas_surface,tmp_co_gc_gfas_surface]
            O3_gc_gfas_surface   = [o3_gc_gfas_surface,tmp_o3_gc_gfas_surface]
            pan_gc_gfas_surface  = [pan_gc_gfas_surface,tmp_pan_gc_gfas_surface]
            hcho_gc_gfas_surface = [hcho_gc_gfas_surface,tmp_hcho_gc_gfas_surface]
            acet_gc_gfas_surface = [acet_gc_gfas_surface,tmp_acet_gc_gfas_surface]
            benz_gc_gfas_surface = [benz_gc_gfas_surface,tmp_benz_gc_gfas_surface]

            ald2_gc_gfas_surface = [ald2_gc_gfas_surface,tmp_ald2_gc_gfas_surface]


            no_gc_gfas_surface   = [no_gc_gfas_surface,tmp_no_gc_gfas_surface]
            no2_gc_gfas_surface  = [no2_gc_gfas_surface,tmp_no2_gc_gfas_surface]
            so2_gc_gfas_surface  = [so2_gc_gfas_surface,tmp_so2_gc_gfas_surface]

            oh_gc_gfas_surface  = [oh_gc_gfas_surface,tmp_oh_gc_gfas_surface] ;; v/v --> molec/cm3 

            date_gc_gfas_surface = [date_gc_gfas_surface,tmp_date_gc_gfas_surface]
            ;utc_gc_gfas_surface  = [utc_gc_gfas_surface,tmp_utc_gc_gfas_surface]
            doy_gc_gfas_surface  = [doy_gc_gfas_surface,tmp_doy_gc_gfas_surface]
            lat_gc_gfas_surface  = [lat_gc_gfas_surface,tmp_lat_gc_gfas_surface]
            lon_gc_gfas_surface  = [lon_gc_gfas_surface,tmp_lon_gc_gfas_surface]
            alt_gc_gfas_surface  = [alt_gc_gfas_surface,tmp_alt_gc_gfas_surface]
            prs_gc_gfas_surface  = [prs_gc_gfas_surface,tmp_prs_gc_gfas_surface]
     ;; add vocs
            c3h8_gc_gfas_surface   = [c3h8_gc_gfas_surface,tmp_c3h8_gc_gfas_surface]
            c3h6_gc_gfas_surface   = [c3h6_gc_gfas_surface,tmp_c3h6_gc_gfas_surface]
            c2h5oh_gc_gfas_surface = [c2h5oh_gc_gfas_surface,tmp_c2h5oh_gc_gfas_surface]
            c5h8_gc_gfas_surface   = [c5h8_gc_gfas_surface,tmp_c5h8_gc_gfas_surface]
            c7h8_gc_gfas_surface   = [c7h8_gc_gfas_surface,tmp_c7h8_gc_gfas_surface]
            dms_gc_gfas_surface    = [dms_gc_gfas_surface,tmp_dms_gc_gfas_surface]
            mek_gc_gfas_surface    = [mek_gc_gfas_surface,tmp_mek_gc_gfas_surface]
            c8h10_gc_gfas_surface  = [c8h10_gc_gfas_surface,tmp_c8h10_gc_gfas_surface]

            acta_gc_gfas_surface   = [acta_gc_gfas_surface,tmp_acta_gc_gfas_surface]
            macr_mvk_gc_gfas_surface = [macr_mvk_gc_gfas_surface,tmp_macr_mvk_gc_gfas_surface]  
            hcooh_gc_gfas_surface  = [hcooh_gc_gfas_surface,tmp_hcooh_gc_gfas_surface]

            mtpa_gc_gfas_surface  = [mtpa_gc_gfas_surface,tmp_mtpa_gc_gfas_surface]
            limo_gc_gfas_surface  = [limo_gc_gfas_surface,tmp_limo_gc_gfas_surface]
            mtpo_gc_gfas_surface  = [mtpo_gc_gfas_surface,tmp_mtpo_gc_gfas_surface]
        ;;add lumped species
            alk4_gc_gfas_surface  = [alk4_gc_gfas_surface,tmp_alk4_gc_gfas_surface]
            prpe_gc_gfas_surface  = [prpe_gc_gfas_surface,tmp_prpe_gc_gfas_surface]
    
            undefine,gc   
        endfor ;; loop for dates, read files
        
;; ==========================================================
;; Remove placeholder
;; ==========================================================
        ;obs
        co_obs = co_obs[1:*] 
        o3_obs = o3_obs[1:*]
        pan_obs= pan_obs[1:*]
        hcho_obs=hcho_obs[1:*]
        acet_obs = acet_obs[1:*]
        benz_obs = benz_obs[1:*]
        ch3oh_obs= ch3oh_obs[1:*]
        ald2_obs = ald2_obs[1:*]


        no_obs = no_obs[1:*]
        no2_obs= no2_obs[1:*]
        so2_obs = so2_obs[1:*]
        prs_obs = prs_obs[1:*] ;; --> hPa

        ;noy_obs   = noy_obs[1:*]
        ch3cn_obs = ch3cn_obs[1:*]

        lat_obs  = lat_obs[1:*]  
        lon_obs  = lon_obs[1:*]
        alt_obs = alt_obs[1:*]
        jday_obs = jday_obs[1:*]
        ;utc_obs = utc_obs[1:*]
        lstime_obs = lstime_obs[1:*]

    ;; PTR and TOGA
    ;; PTR    
        isop_obs_ptr = isop_obs_ptr[1:*]
        acet_obs_ptr = acet_obs_ptr[1:*]
        propanal_obs_ptr = propanal_obs_ptr[1:*]
        mek_obs_ptr  = mek_obs_ptr[1:*]
        butanal_obs_ptr  = butanal_obs_ptr[1:*]
        macr_mvk_obs_ptr     = macr_mvk_obs_ptr[1:*]
        ALD2_obs_ptr = ALD2_obs_ptr[1:*]
        ch3cn_obs_ptr = ch3cn_obs_ptr[1:*]; we also have mixing ratio of ch3cn, need to check the unit lixu

        ch2o_obs_ptr = ch2o_obs_ptr[1:*]
        benz_obs_ptr = benz_obs_ptr[1:*]
        tolu_obs_ptr = tolu_obs_ptr[1:*]

        c8h10_obs_ptr     = c8h10_obs_ptr[1:*]
        Xylenes_obs_ptr   = Xylenes_obs_ptr[1:*]

        monoterpenes_obs_ptr = monoterpenes_obs_ptr[1:*]

        moh_obs_ptr = moh_obs_ptr[1:*]

        hcooh_obs_ptr  = hcooh_obs_ptr[1:*]
        acta_obs_ptr   = acta_obs_ptr[1:*]
        
    ; TOGA
        isop_obs_toga = isop_obs_toga[1:*]
        acet_obs_toga = acet_obs_toga[1:*]
        mek_obs_toga  = mek_obs_toga[1:*]
        macr_mvk_obs_toga     = macr_mvk_obs_toga[1:*]
        ALD2_obs_toga = ALD2_obs_toga[1:*]
        ch2o_obs_toga = ch2o_obs_toga[1:*]
        benz_obs_toga = benz_obs_toga[1:*]
        tolu_obs_toga = tolu_obs_toga[1:*]

        mbo_obs_toga  = mbo_obs_toga[1:*]
        propanal_obs_toga = propanal_obs_toga[1:*]

        c2Butenal_obs_toga = c2Butenal_obs_toga[1:*]
        t2Butenal_obs_toga = t2Butenal_obs_toga[1:*]

        butanal_obs_toga = butanal_obs_toga[1:*]

        etbenzene_obs_toga   = etbenzene_obs_toga[1:*]

        tricyclene_obs_toga  = tricyclene_obs_toga[1:*]
        apinene_obs_toga     = apinene_obs_toga[1:*]
        camphene_obs_toga    = camphene_obs_toga[1:*]
        bpinenemyrcene_obs_toga   = bpinenemyrcene_obs_toga[1:*]
        limonened3carene_obs_toga = limonened3carene_obs_toga[1:*]
        hcn_obs_toga = hcn_obs_toga[1:*]

        C8H10_obs_toga     = C8H10_obs_toga[1:*] 
        Xylenes_obs_toga   = Xylenes_obs_toga[1:*]
        MPXYLENE_obs_toga  = MPXYLENE_obs_toga[1:*] 
        OXYLENE_obs_toga   = OXYLENE_obs_toga[1:*] 
        EtBenzene_obs_toga = EtBenzene_obs_toga[1:*] 
        
        moh_obs_toga = moh_obs_toga[1:*]
    ;; CIMS
        hcooh_obs_cims  = hcooh_obs_cims[1:*]
        acta_obs_cims   = acta_obs_cims[1:*]
    

    ;;for cloud
        rhum = rhum[1:*]
        temperature_obs = temperature_obs[1:*]
    
    ;;add lumped species
        alk4_obs_toga = alk4_obs_toga[1:*]
        prpe_obs_toga = prpe_obs_toga[1:*]
        alk4_obs_awas = alk4_obs_awas[1:*]
        prpe_obs_awas = prpe_obs_awas[1:*]
        
        c2h6_obs_awas = c2h6_obs_awas[1:*]
;; ============================
;; GEOS-Chem 0.25x0.3125
;; ============================
    ;; GFAS_BOP_TOP
        co_gc_gfas_bop_top   = co_gc_gfas_bop_top[1:*]
        o3_gc_gfas_bop_top   = o3_gc_gfas_bop_top[1:*]
        pan_gc_gfas_bop_top  = pan_gc_gfas_bop_top[1:*]
        hcho_gc_gfas_bop_top = hcho_gc_gfas_bop_top[1:*]
        acet_gc_gfas_bop_top = acet_gc_gfas_bop_top[1:*]
        benz_gc_gfas_bop_top = benz_gc_gfas_bop_top[1:*]
    ;    ch3oh_gc= ch3oh_gc[1:*]
        ald2_gc_gfas_bop_top = ald2_gc_gfas_bop_top[1:*]

        no_gc_gfas_bop_top   = no_gc_gfas_bop_top[1:*]
        no2_gc_gfas_bop_top  = no2_gc_gfas_bop_top[1:*]
        so2_gc_gfas_bop_top  = so2_gc_gfas_bop_top[1:*]
        oh_gc_gfas_bop_top   = oh_gc_gfas_bop_top[1:*] 

        date_gc_gfas_bop_top = date_gc_gfas_bop_top[1:*]
        ;utc_gc_gfas_bop_top  = utc_gc_gfas_bop_top[1:*]
        doy_gc_gfas_bop_top  = doy_gc_gfas_bop_top[1:*]
        lat_gc_gfas_bop_top  = lat_gc_gfas_bop_top[1:*]
        lon_gc_gfas_bop_top  = lon_gc_gfas_bop_top[1:*]
        alt_gc_gfas_bop_top  = alt_gc_gfas_bop_top[1:*]
        prs_gc_gfas_bop_top  = prs_gc_gfas_bop_top[1:*]

     ;; add vocs
        c3h8_gc_gfas_bop_top   = c3h8_gc_gfas_bop_top[1:*]
        c3h6_gc_gfas_bop_top   = c3h6_gc_gfas_bop_top[1:*]
        c2h5oh_gc_gfas_bop_top = c2h5oh_gc_gfas_bop_top[1:*]
        c5h8_gc_gfas_bop_top   = c5h8_gc_gfas_bop_top[1:*]
        c7h8_gc_gfas_bop_top   = c7h8_gc_gfas_bop_top[1:*]
        dms_gc_gfas_bop_top    = dms_gc_gfas_bop_top[1:*]
        mek_gc_gfas_bop_top    = mek_gc_gfas_bop_top[1:*]
        c8h10_gc_gfas_bop_top  = c8h10_gc_gfas_bop_top[1:*]

        acta_gc_gfas_bop_top   = acta_gc_gfas_bop_top[1:*]
        macr_mvk_gc_gfas_bop_top = macr_mvk_gc_gfas_bop_top[1:*]
        hcooh_gc_gfas_bop_top  = hcooh_gc_gfas_bop_top[1:*]

        mtpa_gc_gfas_bop_top  = mtpa_gc_gfas_bop_top[1:*]
        limo_gc_gfas_bop_top  = limo_gc_gfas_bop_top[1:*]
        mtpo_gc_gfas_bop_top  = mtpo_gc_gfas_bop_top[1:*]
        
        ;;add lumped species
        alk4_gc_gfas_bop_top = alk4_gc_gfas_bop_top[1:*]
        prpe_gc_gfas_bop_top = prpe_gc_gfas_bop_top[1:*]
        


    ;; GFAS_MAMI
        co_gc_gfas_mami   = co_gc_gfas_mami[1:*]
        o3_gc_gfas_mami   = o3_gc_gfas_mami[1:*]
        pan_gc_gfas_mami  = pan_gc_gfas_mami[1:*]
        hcho_gc_gfas_mami = hcho_gc_gfas_mami[1:*]
        acet_gc_gfas_mami = acet_gc_gfas_mami[1:*]
        benz_gc_gfas_mami = benz_gc_gfas_mami[1:*]
    ;    ch3oh_gc= ch3oh_gc[1:*]
        ald2_gc_gfas_mami = ald2_gc_gfas_mami[1:*]

        no_gc_gfas_mami   = no_gc_gfas_mami[1:*]
        no2_gc_gfas_mami  = no2_gc_gfas_mami[1:*]
        so2_gc_gfas_mami  = so2_gc_gfas_mami[1:*]
        oh_gc_gfas_mami   = oh_gc_gfas_mami[1:*] 

        date_gc_gfas_mami = date_gc_gfas_mami[1:*]
        ;utc_gc_gfas_mami  = utc_gc_gfas_mami[1:*]
        doy_gc_gfas_mami  = doy_gc_gfas_mami[1:*]
        lat_gc_gfas_mami  = lat_gc_gfas_mami[1:*]
        lon_gc_gfas_mami  = lon_gc_gfas_mami[1:*]
        alt_gc_gfas_mami  = alt_gc_gfas_mami[1:*]
        prs_gc_gfas_mami  = prs_gc_gfas_mami[1:*]

     ;; add vocs
        c3h8_gc_gfas_mami   = c3h8_gc_gfas_mami[1:*]
        c3h6_gc_gfas_mami   = c3h6_gc_gfas_mami[1:*]
        c2h5oh_gc_gfas_mami = c2h5oh_gc_gfas_mami[1:*]
        c5h8_gc_gfas_mami   = c5h8_gc_gfas_mami[1:*]
        c7h8_gc_gfas_mami   = c7h8_gc_gfas_mami[1:*]
        dms_gc_gfas_mami    = dms_gc_gfas_mami[1:*]
        mek_gc_gfas_mami    = mek_gc_gfas_mami[1:*]
        c8h10_gc_gfas_mami  = c8h10_gc_gfas_mami[1:*]

        acta_gc_gfas_mami   = acta_gc_gfas_mami[1:*]
        macr_mvk_gc_gfas_mami = macr_mvk_gc_gfas_mami[1:*]
        hcooh_gc_gfas_mami  = hcooh_gc_gfas_mami[1:*]

        mtpa_gc_gfas_mami  = mtpa_gc_gfas_mami[1:*]
        limo_gc_gfas_mami  = limo_gc_gfas_mami[1:*]
        mtpo_gc_gfas_mami  = mtpo_gc_gfas_mami[1:*]
        
        ;;add lumped species
        alk4_gc_gfas_mami = alk4_gc_gfas_mami[1:*]
        prpe_gc_gfas_mami = prpe_gc_gfas_mami[1:*]
        
    ;; GFAS_NOBB
        co_gc_gfas_nobb   = co_gc_gfas_nobb[1:*]
        o3_gc_gfas_nobb   = o3_gc_gfas_nobb[1:*]
        pan_gc_gfas_nobb  = pan_gc_gfas_nobb[1:*]
        hcho_gc_gfas_nobb = hcho_gc_gfas_nobb[1:*]
        acet_gc_gfas_nobb = acet_gc_gfas_nobb[1:*]
        benz_gc_gfas_nobb = benz_gc_gfas_nobb[1:*]
    ;    ch3oh_gc= ch3oh_gc[1:*]
        ald2_gc_gfas_nobb = ald2_gc_gfas_nobb[1:*]

        no_gc_gfas_nobb   = no_gc_gfas_nobb[1:*]
        no2_gc_gfas_nobb  = no2_gc_gfas_nobb[1:*]
        so2_gc_gfas_nobb  = so2_gc_gfas_nobb[1:*]
        oh_gc_gfas_nobb   = oh_gc_gfas_nobb[1:*] 

        date_gc_gfas_nobb = date_gc_gfas_nobb[1:*]
        ;utc_gc_gfas_nobb  = utc_gc_gfas_nobb[1:*]
        doy_gc_gfas_nobb  = doy_gc_gfas_nobb[1:*]
        lat_gc_gfas_nobb  = lat_gc_gfas_nobb[1:*]
        lon_gc_gfas_nobb  = lon_gc_gfas_nobb[1:*]
        alt_gc_gfas_nobb  = alt_gc_gfas_nobb[1:*]
        prs_gc_gfas_nobb  = prs_gc_gfas_nobb[1:*]

     ;; add vocs
        c3h8_gc_gfas_nobb   = c3h8_gc_gfas_nobb[1:*]
        c3h6_gc_gfas_nobb   = c3h6_gc_gfas_nobb[1:*]
        c2h5oh_gc_gfas_nobb = c2h5oh_gc_gfas_nobb[1:*]
        c5h8_gc_gfas_nobb   = c5h8_gc_gfas_nobb[1:*]
        c7h8_gc_gfas_nobb   = c7h8_gc_gfas_nobb[1:*]
        dms_gc_gfas_nobb    = dms_gc_gfas_nobb[1:*]
        mek_gc_gfas_nobb    = mek_gc_gfas_nobb[1:*]
        c8h10_gc_gfas_nobb  = c8h10_gc_gfas_nobb[1:*]

        acta_gc_gfas_nobb   = acta_gc_gfas_nobb[1:*]
        macr_mvk_gc_gfas_nobb = macr_mvk_gc_gfas_nobb[1:*]
        hcooh_gc_gfas_nobb  = hcooh_gc_gfas_nobb[1:*]

        mtpa_gc_gfas_nobb  = mtpa_gc_gfas_nobb[1:*]
        limo_gc_gfas_nobb  = limo_gc_gfas_nobb[1:*]
        mtpo_gc_gfas_nobb  = mtpo_gc_gfas_nobb[1:*]
        
        ;;add lumped species
        alk4_gc_gfas_nobb = alk4_gc_gfas_nobb[1:*]
        prpe_gc_gfas_nobb = prpe_gc_gfas_nobb[1:*]
        

    ;; GFAS_PBL65_FT35
        co_gc_gfas_pbl65_ft35   = co_gc_gfas_pbl65_ft35[1:*]
        o3_gc_gfas_pbl65_ft35   = o3_gc_gfas_pbl65_ft35[1:*]
        pan_gc_gfas_pbl65_ft35  = pan_gc_gfas_pbl65_ft35[1:*]
        hcho_gc_gfas_pbl65_ft35 = hcho_gc_gfas_pbl65_ft35[1:*]
        acet_gc_gfas_pbl65_ft35 = acet_gc_gfas_pbl65_ft35[1:*]
        benz_gc_gfas_pbl65_ft35 = benz_gc_gfas_pbl65_ft35[1:*]
    ;    ch3oh_gc= ch3oh_gc[1:*]
        ald2_gc_gfas_pbl65_ft35 = ald2_gc_gfas_pbl65_ft35[1:*]

        no_gc_gfas_pbl65_ft35   = no_gc_gfas_pbl65_ft35[1:*]
        no2_gc_gfas_pbl65_ft35  = no2_gc_gfas_pbl65_ft35[1:*]
        so2_gc_gfas_pbl65_ft35  = so2_gc_gfas_pbl65_ft35[1:*]
        oh_gc_gfas_pbl65_ft35   = oh_gc_gfas_pbl65_ft35[1:*] 

        date_gc_gfas_pbl65_ft35 = date_gc_gfas_pbl65_ft35[1:*]
        ;utc_gc_gfas_pbl65_ft35  = utc_gc_gfas_pbl65_ft35[1:*]
        doy_gc_gfas_pbl65_ft35  = doy_gc_gfas_pbl65_ft35[1:*]
        lat_gc_gfas_pbl65_ft35  = lat_gc_gfas_pbl65_ft35[1:*]
        lon_gc_gfas_pbl65_ft35  = lon_gc_gfas_pbl65_ft35[1:*]
        alt_gc_gfas_pbl65_ft35  = alt_gc_gfas_pbl65_ft35[1:*]
        prs_gc_gfas_pbl65_ft35  = prs_gc_gfas_pbl65_ft35[1:*]

     ;; add vocs
        c3h8_gc_gfas_pbl65_ft35   = c3h8_gc_gfas_pbl65_ft35[1:*]
        c3h6_gc_gfas_pbl65_ft35   = c3h6_gc_gfas_pbl65_ft35[1:*]
        c2h5oh_gc_gfas_pbl65_ft35 = c2h5oh_gc_gfas_pbl65_ft35[1:*]
        c5h8_gc_gfas_pbl65_ft35   = c5h8_gc_gfas_pbl65_ft35[1:*]
        c7h8_gc_gfas_pbl65_ft35   = c7h8_gc_gfas_pbl65_ft35[1:*]
        dms_gc_gfas_pbl65_ft35    = dms_gc_gfas_pbl65_ft35[1:*]
        mek_gc_gfas_pbl65_ft35    = mek_gc_gfas_pbl65_ft35[1:*]
        c8h10_gc_gfas_pbl65_ft35  = c8h10_gc_gfas_pbl65_ft35[1:*]

        acta_gc_gfas_pbl65_ft35   = acta_gc_gfas_pbl65_ft35[1:*]
        macr_mvk_gc_gfas_pbl65_ft35 = macr_mvk_gc_gfas_pbl65_ft35[1:*]
        hcooh_gc_gfas_pbl65_ft35  = hcooh_gc_gfas_pbl65_ft35[1:*]

        mtpa_gc_gfas_pbl65_ft35  = mtpa_gc_gfas_pbl65_ft35[1:*]
        limo_gc_gfas_pbl65_ft35  = limo_gc_gfas_pbl65_ft35[1:*]
        mtpo_gc_gfas_pbl65_ft35  = mtpo_gc_gfas_pbl65_ft35[1:*]
        
        ;;add lumped species
        alk4_gc_gfas_pbl65_ft35 = alk4_gc_gfas_pbl65_ft35[1:*]
        prpe_gc_gfas_pbl65_ft35 = prpe_gc_gfas_pbl65_ft35[1:*]
        

    ;; GFAS_SF_MAMI
        co_gc_gfas_sf_mami   = co_gc_gfas_sf_mami[1:*]
        o3_gc_gfas_sf_mami   = o3_gc_gfas_sf_mami[1:*]
        pan_gc_gfas_sf_mami  = pan_gc_gfas_sf_mami[1:*]
        hcho_gc_gfas_sf_mami = hcho_gc_gfas_sf_mami[1:*]
        acet_gc_gfas_sf_mami = acet_gc_gfas_sf_mami[1:*]
        benz_gc_gfas_sf_mami = benz_gc_gfas_sf_mami[1:*]
    ;    ch3oh_gc= ch3oh_gc[1:*]
        ald2_gc_gfas_sf_mami = ald2_gc_gfas_sf_mami[1:*]

        no_gc_gfas_sf_mami   = no_gc_gfas_sf_mami[1:*]
        no2_gc_gfas_sf_mami  = no2_gc_gfas_sf_mami[1:*]
        so2_gc_gfas_sf_mami  = so2_gc_gfas_sf_mami[1:*]
        oh_gc_gfas_sf_mami   = oh_gc_gfas_sf_mami[1:*] 

        date_gc_gfas_sf_mami = date_gc_gfas_sf_mami[1:*]
        ;utc_gc_gfas_sf_mami  = utc_gc_gfas_sf_mami[1:*]
        doy_gc_gfas_sf_mami  = doy_gc_gfas_sf_mami[1:*]
        lat_gc_gfas_sf_mami  = lat_gc_gfas_sf_mami[1:*]
        lon_gc_gfas_sf_mami  = lon_gc_gfas_sf_mami[1:*]
        alt_gc_gfas_sf_mami  = alt_gc_gfas_sf_mami[1:*]
        prs_gc_gfas_sf_mami  = prs_gc_gfas_sf_mami[1:*]

     ;; add vocs
        c3h8_gc_gfas_sf_mami   = c3h8_gc_gfas_sf_mami[1:*]
        c3h6_gc_gfas_sf_mami   = c3h6_gc_gfas_sf_mami[1:*]
        c2h5oh_gc_gfas_sf_mami = c2h5oh_gc_gfas_sf_mami[1:*]
        c5h8_gc_gfas_sf_mami   = c5h8_gc_gfas_sf_mami[1:*]
        c7h8_gc_gfas_sf_mami   = c7h8_gc_gfas_sf_mami[1:*]
        dms_gc_gfas_sf_mami    = dms_gc_gfas_sf_mami[1:*]
        mek_gc_gfas_sf_mami    = mek_gc_gfas_sf_mami[1:*]
        c8h10_gc_gfas_sf_mami  = c8h10_gc_gfas_sf_mami[1:*]

        acta_gc_gfas_sf_mami   = acta_gc_gfas_sf_mami[1:*]
        macr_mvk_gc_gfas_sf_mami = macr_mvk_gc_gfas_sf_mami[1:*]
        hcooh_gc_gfas_sf_mami  = hcooh_gc_gfas_sf_mami[1:*]

        mtpa_gc_gfas_sf_mami  = mtpa_gc_gfas_sf_mami[1:*]
        limo_gc_gfas_sf_mami  = limo_gc_gfas_sf_mami[1:*]
        mtpo_gc_gfas_sf_mami  = mtpo_gc_gfas_sf_mami[1:*]
        
        ;;add lumped species
        alk4_gc_gfas_sf_mami = alk4_gc_gfas_sf_mami[1:*]
        prpe_gc_gfas_sf_mami = prpe_gc_gfas_sf_mami[1:*]
        
    ;; GFAS_SURFACE
        co_gc_gfas_surface   = co_gc_gfas_surface[1:*]
        o3_gc_gfas_surface   = o3_gc_gfas_surface[1:*]
        pan_gc_gfas_surface  = pan_gc_gfas_surface[1:*]
        hcho_gc_gfas_surface = hcho_gc_gfas_surface[1:*]
        acet_gc_gfas_surface = acet_gc_gfas_surface[1:*]
        benz_gc_gfas_surface = benz_gc_gfas_surface[1:*]
    ;    ch3oh_gc= ch3oh_gc[1:*]
        ald2_gc_gfas_surface = ald2_gc_gfas_surface[1:*]

        no_gc_gfas_surface   = no_gc_gfas_surface[1:*]
        no2_gc_gfas_surface  = no2_gc_gfas_surface[1:*]
        so2_gc_gfas_surface  = so2_gc_gfas_surface[1:*]
        oh_gc_gfas_surface   = oh_gc_gfas_surface[1:*] 

        date_gc_gfas_surface = date_gc_gfas_surface[1:*]
        ;utc_gc_gfas_surface  = utc_gc_gfas_surface[1:*]
        doy_gc_gfas_surface  = doy_gc_gfas_surface[1:*]
        lat_gc_gfas_surface  = lat_gc_gfas_surface[1:*]
        lon_gc_gfas_surface  = lon_gc_gfas_surface[1:*]
        alt_gc_gfas_surface  = alt_gc_gfas_surface[1:*]
        prs_gc_gfas_surface  = prs_gc_gfas_surface[1:*]

     ;; add vocs
        c3h8_gc_gfas_surface   = c3h8_gc_gfas_surface[1:*]
        c3h6_gc_gfas_surface   = c3h6_gc_gfas_surface[1:*]
        c2h5oh_gc_gfas_surface = c2h5oh_gc_gfas_surface[1:*]
        c5h8_gc_gfas_surface   = c5h8_gc_gfas_surface[1:*]
        c7h8_gc_gfas_surface   = c7h8_gc_gfas_surface[1:*]
        dms_gc_gfas_surface    = dms_gc_gfas_surface[1:*]
        mek_gc_gfas_surface    = mek_gc_gfas_surface[1:*]
        c8h10_gc_gfas_surface  = c8h10_gc_gfas_surface[1:*]

        acta_gc_gfas_surface   = acta_gc_gfas_surface[1:*]
        macr_mvk_gc_gfas_surface = macr_mvk_gc_gfas_surface[1:*]
        hcooh_gc_gfas_surface  = hcooh_gc_gfas_surface[1:*]

        mtpa_gc_gfas_surface  = mtpa_gc_gfas_surface[1:*]
        limo_gc_gfas_surface  = limo_gc_gfas_surface[1:*]
        mtpo_gc_gfas_surface  = mtpo_gc_gfas_surface[1:*]
        
        ;;add lumped species
        alk4_gc_gfas_surface = alk4_gc_gfas_surface[1:*]
        prpe_gc_gfas_surface = prpe_gc_gfas_surface[1:*]
                
        if n_elements(co_obs) ne n_elements(co_gc_gfas_bop_top) or $
            n_elements(co_obs) ne n_elements(co_gc_gfas_mami) or $
            n_elements(co_obs) ne n_elements(co_gc_gfas_nobb) or $
            n_elements(co_obs) ne n_elements(co_gc_gfas_PBL65_FT35) or $
            n_elements(co_obs) ne n_elements(co_gc_gfas_sf_mami) or $
            n_elements(co_obs) ne n_elements(co_gc_gfas_surface) then stop 

;; test for percentile
;data_test = ch3cn_obs
;ind = where(data_test gt 0, ct)
;if ct gt 0 then data_test = data_test[ind]
;print, cgPercentiles(data_test, Percentiles=[0.1,0.25, 0.5, 0.75,0.9])
;continue

;; ==============================
;; set up detection limit for PTR
;; ==============================
;        LoD = 50.0/1000
;        ind = where(hcho_obs_ptr le LoD, ct)
;        if ct gt 0 then begin
;            ;hcho_obs_ptr[ind] = 4.0*LoD
;            ;hcho_obs_toga[ind] = 4.0*LoD
;            hcho_obs_ptr[ind] = hcho_obs_toga[ind]
;        endif
;        
;        LoD = 10.0/1000
;        ind = where(ald2_obs_ptr le LoD, ct)
;        if ct gt 0 then begin
;            ;ald2_obs_ptr[ind] = 4.0*LoD
;            ;ald2_obs_toga[ind] = 4.0*LoD
;            ald2_obs_ptr[ind] = ald2_obs_toga[ind]
;        endif
;        
;        LoD = 10.0/1000
;        ind = where(acet_obs_ptr le LoD, ct)
;        if ct gt 0 then begin
;            ;acet_obs_ptr[ind] = 4.0*LoD
;            ;acet_obs_toga[ind] = 4.0*LoD
;            acet_obs_ptr[ind] = acet_obs_toga[ind]
;        endif
;
;        LoD = 10.0/1000
;        ind = where(mek_obs_ptr le LoD, ct)
;        if ct gt 0 then begin
;            ;mek_obs_ptr[ind] = 4.0*LoD
;            ;mek_obs_toga[ind] = 4.0*LoD
;            mek_obs_ptr[ind] = mek_obs_toga[ind]
;        endif
;        
;        LoD = 10.0/1000  ; change 15.0 to 10.0
;        ind = where(benz_obs_ptr le LoD, ct)
;        if ct gt 0 then begin
;            ;benz_obs_ptr[ind] = 4.0*LoD
;            ;benz_obs_toga[ind] = 4.0*LoD
;            benz_obs_ptr[ind] = benz_obs_toga[ind]
;        endif
        
        
;        LoD = 10.0/1000 ; change 20.0 to 10.0
;        ind = where(tolu_obs_ptr le LoD, ct)
;        if ct gt 0 then begin
;            ;tolu_obs_ptr[ind] = 4.0*LoD
;            ;tolu_obs_toga[ind] = 4.0*LoD
;            tolu_obs_ptr[ind] = tolu_obs_toga[ind]
;        endif
;        
;        LoD = 10.0/1000 
;        ind = where(c8h10_obs_ptr le LoD, ct)
;        if ct gt 0 then begin
;            ;c8h10_obs_ptr[ind] = 4.0*LoD
;            ;c8h10_obs_toga[ind] = 4.0*LoD
;            c8h10_obs_ptr[ind] = c8h10_obs_toga[ind]
;        endif

        LoD = 50.0/1000
        ind = where(ch2o_obs_ptr le LoD and ch2o_obs_ptr gt 0, ct) ; try make PTR lt 0
        if ct gt 0 then begin
            ;ch2o_obs_ptr[ind] = !VALUES.F_NAN
            ch2o_obs_ptr[ind] = ch2o_obs_toga[ind]
        endif

        LoD = 10.0/1000
        ind = where(ald2_obs_ptr le LoD and ald2_obs_ptr gt 0, ct)
        if ct gt 0 then begin
            ;ald2_obs_ptr[ind] = !VALUES.F_NAN
            ald2_obs_ptr[ind] = ald2_obs_toga[ind]
        endif

        LoD = 20.0/1000
        ind = where(acet_obs_ptr le LoD and acet_obs_ptr gt 0, ct)
        if ct gt 0 then begin
            ;acet_obs_ptr[ind] = !VALUES.F_NAN
            acet_obs_ptr[ind] = acet_obs_toga[ind]
         endif

        LoD = 10.0/1000
        ind = where(mek_obs_ptr le LoD and mek_obs_ptr gt 0, ct)
            if ct gt 0 then begin
            ;mek_obs_ptr[ind] = !VALUES.F_NAN
            mek_obs_ptr[ind] = mek_obs_toga[ind]
        endif

        LoD = 30.0/1000 ; change 15.0 to 10.0
        ind = where(benz_obs_ptr le LoD and benz_obs_ptr gt 0, ct)
        if ct gt 0 then begin
            ;benz_obs_ptr[ind] = !VALUES.F_NAN
            benz_obs_ptr[ind] = benz_obs_toga[ind]
        endif

        LoD = 50.0/1000 ; change 20.0 to 10.0
        ind = where(tolu_obs_ptr le LoD and tolu_obs_ptr le 0, ct)
        if ct gt 0 then begin
            ;tolu_obs_ptr[ind] = !VALUES.F_NAN
            tolu_obs_ptr[ind] = tolu_obs_toga[ind]
        endif

        LoD = 50.0/1000
        ind = where(Xylenes_obs_ptr le LoD and Xylenes_obs_ptr gt 0, ct)
        if ct gt 0 then begin
            ;c8h10_obs_ptr[ind] = !VALUES.F_NAN
            Xylenes_obs_ptr[ind] = Xylenes_obs_toga[ind]
        endif

        LoD = 50.0/1000
        ind = where(c8h10_obs_ptr le LoD and c8h10_obs_ptr gt 0, ct)
        if ct gt 0 then begin
            ;c8h10_obs_ptr[ind] = !VALUES.F_NAN
            c8h10_obs_ptr[ind] = c8h10_obs_toga[ind]
        endif
;; ========================
;; Prepare for the filters 
;; ========================
    ;; jupyter time
        doy_obs = jday_obs + utc_obs / 3600. / 24. 

    ;; lixu, 12/18/2019, 0208/2021   
    ;; plume filter setting, use in the plotting part
        no2_thresh = 4. ; ppb
        co_thresh  = 150. ; ppb
    ;; urban influence
        x224TrimePentane_obs_toga = 20/1000
        C2Cl4_obs_toga = 2/1000
        HFC134a_obs_toga = 125/1000
        HCFC22_obs_toga = 275/1000
    
    ;; biomass burning threshhold 
        hcn_thresh   = 275./1000 ;pptv to ppbv
        ch3cn_thresh = 131/1000. ;;ppt to ppb, 131
        co_thresh =  85
    ;; strat filter
        strat_thresh=1.25 ; O3/CO ratio
        
    ;; CO threshold of fresh and aged smoke 
        ;; young
        methylfuran_thresh = 0.7/1000 
        ;; median
        acrolein_thresh = 7.4/1000
        ;; old
        acrylonitrile_thresh = 2.9/1000
        
    ;; get hhmm of each day
        utc_obs = 24.*(doy_obs mod 1.)
        utc_gc_gfas  = 24.*(doy_gc_gfas_bop_top  mod 1.)

        
    ;; for vertical profile
        yrange = [0,7]
    ;; lxu, 06202021
        acn_co_tresh = 2.01 ;ppb/ppm
;; ========keep the same for regression, profiles and time series=================
        if keyword_set(all) then rows = 3
        if keyword_set(all) then cols = 3

        if keyword_set(each) then rows = 1
        if keyword_set(each) then cols = 1
    
        if keyword_set(onebyone) then rows = 1
        if keyword_set(onebyone) then cols = 1
        
        
        
        multipanel, rows = rows, cols = cols

        if keyword_set(all) then s_end = 9
        if keyword_set(obs_only) then s_end = 8
        if keyword_set(each) then s_end = 0


        for s = 0, s_end do begin
    ;;set for individual species, e.g., CO, benzene
            if keyword_set(each) then begin
                if keyword_set(ind_co) then begin
                    ovoc_ptr = co_obs
                    ovoc_other = co_obs
                    gvoc_bop_top = co_gc_gfas_bop_top
                    gvoc_mami = co_gc_gfas_mami
                    gvoc_nobb = co_gc_gfas_nobb
                    gvoc_pbl65_ft35 = co_gc_gfas_pbl65_ft35
                    gvoc_sf_mami = co_gc_gfas_sf_mami
                    gvoc_surface = co_gc_gfas_surface
                    title = 'CO'
                    xrange = [0,600]
                    yrange = [0,7]

                endif
                if keyword_set(ind_benz) then begin
                    ovoc = benz_obs_ptr
                    gvoc_bop_top = benz_gc_gfas_bop_top
                    gvoc_mami = benz_gc_gfas_mami
                    gvoc_nobb = benz_gc_gfas_nobb
                    gvoc_pbl65_ft35 = benz_gc_gfas_pbl65_ft35
                    gvoc_sf_mami = benz_gc_gfas_sf_mami
                    gvoc_surface = benz_gc_gfas_surface
                    title = 'Benzene'
                    xrange = [0,0.6]

                endif        
                if keyword_set(ind_o3) then begin
                    ovoc = o3_obs
                    gvoc_bop_top = o3_gc_gfas_bop_top
                    gvoc_mami = o3_gc_gfas_mami
                    gvoc_nobb = o3_gc_gfas_nobb
                    gvoc_pbl65_ft35 = o3_gc_gfas_pbl65_ft35
                    gvoc_sf_mami = o3_gc_gfas_sf_mami
                    gvoc_surface = o3_gc_gfas_surface
                    title = 'O3'
                    xrange = [0,100]

                endif  
                if keyword_set(ind_pan) then begin
                    ovoc = pan_obs
                    gvoc_bop_top = pan_gc_gfas_bop_top
                    gvoc_mami = pan_gc_gfas_mami
                    gvoc_nobb = pan_gc_gfas_nobb
                    gvoc_pbl65_ft35 = pan_gc_gfas_pbl65_ft35
                    gvoc_sf_mami = pan_gc_gfas_sf_mami
                    gvoc_surface = pan_gc_gfas_surface
                    title = 'PAN'
                    xrange = [0,1]

                endif  
                if keyword_set(ind_no2) then begin
                    ovoc = no2_obs
                    gvoc_bop_top = no2_gc_gfas_bop_top
                    gvoc_mami = no2_gc_gfas_mami
                    gvoc_nobb = no2_gc_gfas_nobb
                    gvoc_pbl65_ft35 = no2_gc_gfas_pbl65_ft35
                    gvoc_sf_mami = no2_gc_gfas_sf_mami
                    gvoc_surface = no2_gc_gfas_surface
                    title = 'NO2'
                    xrange = [0,1]

                endif  
                if keyword_set(ind_dms) then begin
                    ovoc = dms_obs
                    gvoc_bop_top = dms_gc_gfas_bop_top
                    gvoc_mami = dms_gc_gfas_mami
                    gvoc_nobb = dms_gc_gfas_nobb
                    gvoc_pbl65_ft35 = dms_gc_gfas_pbl65_ft35
                    gvoc_sf_mami = dms_gc_gfas_sf_mami
                    gvoc_surface = dms_gc_gfas_surface
                    title = 'DMS'
                    xrange = [0,1]
                endif  
            endif

    ;; ignore add up keyword here
    ;;all
            if keyword_set(all) then begin
                case s of                     
                    ;;c2h6(AWAS, discrete data)

                    ;;hcho
                    0:begin
                        ovoc_ptr = ch2o_obs_ptr 
                        ovoc_other = ch2o_obs_toga 
                        gvoc_bop_top = hcho_gc_gfas_bop_top
                        gvoc_mami = hcho_gc_gfas_mami
                        gvoc_nobb = hcho_gc_gfas_nobb
                        gvoc_pbl65_ft35 = hcho_gc_gfas_pbl65_ft35
                        gvoc_sf_mami = hcho_gc_gfas_sf_mami
                        gvoc_surface = hcho_gc_gfas_surface
                        title = 'Formaldehyde'
                        xrange = [0,6]
                        
                        if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) then xrange = [0,2]    
                        
                        if keyword_set(no_obs_include) then xrange = [0, 2]
                    end
                    ;;ald2
                    1:begin
                        ovoc_ptr = ald2_obs_ptr
                        ovoc_other = ald2_obs_toga
                        
                        gvoc_bop_top = ald2_gc_gfas_bop_top/2
                        gvoc_mami = ald2_gc_gfas_mami/2
                        gvoc_nobb = ald2_gc_gfas_nobb/2
                        gvoc_pbl65_ft35 = ald2_gc_gfas_pbl65_ft35/2
                        gvoc_sf_mami = ald2_gc_gfas_sf_mami/2
                        gvoc_surface = ald2_gc_gfas_surface/2
                        title = 'Acetaldehyde'
                        xrange = [0,3]

                        if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) then xrange = [0,0.6]    
                        if keyword_set(OHrplot) then xrange = [0,0.8]
                        if keyword_set(no_obs_include) then xrange = [0, 0.5]

                    end
                    ;;acet
                    2:begin
                        ovoc_ptr = acet_obs_ptr 
                        ovoc_other = acet_obs_toga
                        
                        gvoc_bop_top = acet_gc_gfas_bop_top/3
                        gvoc_mami = acet_gc_gfas_mami/3
                        gvoc_nobb = acet_gc_gfas_nobb/3
                        gvoc_pbl65_ft35 = acet_gc_gfas_pbl65_ft35/3
                        gvoc_sf_mami = acet_gc_gfas_sf_mami/3
                        gvoc_surface = acet_gc_gfas_surface/3
                        title = 'Acetone'
                        xrange = [0,4]

                        if keyword_set(no_obs_include) then xrange = [0, 2]

                    end
                    ;;c3h8
                    3:begin
                        ovoc_ptr = propane_obs_toga
                        ovoc_other = propane_obs_toga
                        gvoc_bop_top = c3h8_gc_gfas_bop_top/3
                        gvoc_mami = c3h8_gc_gfas_mami/3
                        gvoc_nobb = c3h8_gc_gfas_nobb/3
                        gvoc_pbl65_ft35 = c3h8_gc_gfas_pbl65_ft35/3
                        gvoc_sf_mami = c3h8_gc_gfas_sf_mami/3
                        gvoc_surface = c3h8_gc_gfas_surface/3
                        title = 'Propane'
                        xrange = [0,1]

                        yrange = [0,7]
                        
                        if keyword_set(no_obs_include) then xrange = [0, 0.3]

                        
                    end
                    ;;mek
                    4:begin
                        ovoc_ptr = mek_obs_ptr
                        ovoc_other = mek_obs_toga
                        
                        gvoc_bop_top = mek_gc_gfas_bop_top/5
                        gvoc_mami = mek_gc_gfas_mami/5
                        gvoc_nobb = mek_gc_gfas_nobb/5
                        gvoc_pbl65_ft35 = mek_gc_gfas_pbl65_ft35/5
                        gvoc_sf_mami = mek_gc_gfas_sf_mami/5
                        gvoc_surface = mek_gc_gfas_surface/5
                        title = 'MEK'
                        xrange = [0,0.4]

                        if keyword_set(filter1) and plume eq 1 then xrange = [0,0.3]

                    end
                    ;;eoh
                    5:begin
                        ovoc_ptr = co_obs
                        ovoc_other = co_obs
                        
                        gvoc_bop_top = c2h5oh_gc_gfas_bop_top/2
                        gvoc_mami = c2h5oh_gc_gfas_mami/2
                        gvoc_nobb = c2h5oh_gc_gfas_nobb/2
                        gvoc_pbl65_ft35 = c2h5oh_gc_gfas_pbl65_ft35/2
                        gvoc_sf_mami = c2h5oh_gc_gfas_sf_mami/2
                        gvoc_surface = c2h5oh_gc_gfas_surface/2
                        title = 'Ethanol'
                        xrange = [0,0.6]

                    end
                    ;;benz
                    6:begin
                        ovoc_ptr = benz_obs_ptr
                        ovoc_other = benz_obs_toga
                        
                        gvoc_bop_top = benz_gc_gfas_bop_top/6
                        gvoc_mami = benz_gc_gfas_mami/6
                        gvoc_nobb = benz_gc_gfas_nobb/6
                        gvoc_pbl65_ft35 = benz_gc_gfas_pbl65_ft35/6
                        gvoc_sf_mami = benz_gc_gfas_sf_mami/6
                        gvoc_surface = benz_gc_gfas_surface/6
                        title = 'Benzene'
                        xrange = [0,0.6]


                        ;if keyword_set(no_obs) then xrange = [0, 0.3]

                    end
                    ;;tolu
                    7:begin
                        ovoc_ptr = tolu_obs_ptr
                        ovoc_other = tolu_obs_toga
                        
                        gvoc_bop_top = c7h8_gc_gfas_bop_top/7
                        gvoc_mami = c7h8_gc_gfas_mami/7
                        gvoc_nobb = c7h8_gc_gfas_nobb/7
                        gvoc_pbl65_ft35 = c7h8_gc_gfas_pbl65_ft35/7
                        gvoc_sf_mami = c7h8_gc_gfas_sf_mami/7
                        gvoc_surface = c7h8_gc_gfas_surface/7
                        title = 'Toluene'
                        xrange = [0,0.25]

                        if keyword_set(no_obs_include) then xrange = [0, 0.1]

                    end  
                    ;;xyle
                    8:begin
                        ovoc_ptr = Xylenes_obs_ptr ;; PTR would be different 
                        ovoc_other = Xylenes_obs_toga ;; PTR would be different 
                        
                        gvoc_bop_top = c8h10_gc_gfas_bop_top/8
                        gvoc_mami = c8h10_gc_gfas_mami/8
                        gvoc_nobb = c8h10_gc_gfas_nobb/8
                        gvoc_pbl65_ft35 = c8h10_gc_gfas_pbl65_ft35/8
                        gvoc_sf_mami = c8h10_gc_gfas_sf_mami/8
                        gvoc_surface = c8h10_gc_gfas_surface/8
                        title = 'Xylene'
                        xrange = [0, 0.025]

                        if keyword_set(no_obs_include) then xrange = [0, 0.02]

                    end

                    ;;CO
                    9:begin
                        ovoc_ptr = co_obs
                        ovoc_other = co_obs
                        
                        
                        gvoc_bop_top = co_gc_gfas_bop_top
                        gvoc_mami = co_gc_gfas_mami
                        gvoc_nobb = co_gc_gfas_nobb
                        gvoc_pbl65_ft35 = co_gc_gfas_pbl65_ft35
                        gvoc_sf_mami = co_gc_gfas_sf_mami
                        gvoc_surface = co_gc_gfas_surface
                        title = 'CO'
                        ;xrange = [0,1.5]
                    end
                endcase
            endif
            
            print, 'processing ' + title
            ;stop


;///////////////////////////////////////////////////////////////////////////////////////////
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;++                                      STEP 3                                           ++
;++                     Setting filters:NA/bb emission/fresh/aged                         ++
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;///////////////////////////////////////////////////////////////////////////////////////////
    ;for vertical profile
            zvar_obs  = alt_obs
            zvar_gc_gfas   = alt_gc_gfas_bop_top

            ;test,lxu
            if keyword_set(prs) then begin
                    zvar_obs  = prs_obs
                    zvar_gc_gfas   = prs_obs
            endif
        
    ;for time series, reset the time array every loop (diff)
            xtime_obs  = doy_obs
            xtime_gc   = doy_gc_gfas_bop_top

    ;for voc .vs. co
            oco_tmp=co_obs
            gco_bop_top = co_gc_gfas_bop_top
            gco_mami = co_gc_gfas_mami
            gco_nobb = co_gc_gfas_nobb
            gco_pbl65_ft35 = co_gc_gfas_pbl65_ft35
            gco_sf_mami = co_gc_gfas_sf_mami
            gco_surface = co_gc_gfas_surface
            
    ;;add for OVOC VS. ISOP
            oisop_ptr = isop_obs_ptr
            oisop_toga = isop_obs_toga
            gisop_bop_top = c5h8_gc_gfas_bop_top
            gisop_mami = c5h8_gc_gfas_mami
            gisop_nobb = c5h8_gc_gfas_nobb
            gisop_pbl65_ft35 = c5h8_gc_gfas_pbl65_ft35
            gisop_sf_mami = c5h8_gc_gfas_sf_mami
            gisop_surface = c5h8_gc_gfas_surface

    ;for filters
            no2_filter=no2_obs
            o3_filter=o3_obs
            co_filter=co_obs
            ch3cn_filter=ch3cn_obs_ptr
            hcn_filter = hcn_obs_toga
            acn_co_filter = ch3cn_obs_ptr/co_obs*1000 ; convert it into ppb/ppm


    ;for new criteria
            isop_filter = isop_obs_ptr
            isop_filter_toga = isop_obs_toga
            mek_filter = mek_obs_ptr 
            benz_filter = benz_obs_ptr
            tolu_filter = tolu_obs_ptr
            xyle_filter = c8h10_obs_ptr
            mono_filter = tricyclene_obs_toga + aPinene_obs_TOGA + Camphene_obs_TOGA + bPineneMyrcene_obs_toga

    ;;tracks
            lat_gc = lat_gc_gfas_bop_top
            lon_gc = lon_gc_gfas_bop_top
    
    ;;setting tmp value
            tmp_isop_ptr = oisop_ptr
            tmp_isop_toga = oisop_toga
            tmp_lat_obs = lat_obs
            tmp_lon_obs = lon_obs
            tmp_lat_gc  = lat_gc
            tmp_lon_gc  = lon_gc
            ;cloud
            tmp_rhum = rhum

;///////////////////////////////////////////////////////////////////////////////////////////
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;++                                      STEP 3                                           ++
;++                     Setting filters:NA/bb emission/fresh/aged                         ++
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;///////////////////////////////////////////////////////////////////////////////////////////

    ;; For removing missing values!!!!!!!
    ;;PTR and TOGA
    
            if keyword_set(codeployed) then begin
                ind = WHERE(ovoc_ptr LT 0 or ovoc_other LT 0 or zvar_obs LT 0, count)
                if count gt 0 then begin
                    ovoc_ptr[ind] = !VALUES.F_NAN
                    ovoc_other[ind] = !VALUES.F_NAN
                endif
            endif else begin
                ind = WHERE(ovoc_ptr LT 0 or zvar_obs LT 0, count)
                if count gt 0 then begin
                    ovoc_ptr[ind] = !VALUES.F_NAN
                endif
            endelse

;;=================
;;plume filters
;;=================
            if plume eq 1 then begin
                ;;filter2:use 25 percentile of acetonitrile
                if keyword_set(filter2_nobb) then begin
                    ind = where(ch3cn_filter ge ch3cn_thresh,ct)
                    if ct gt 0 then begin
                        ovoc_ptr[ind] = !VALUES.F_NAN
                        ovoc_other[ind] = !VALUES.F_NAN
                    endif
                    ;print, 'this is the data points to be deleted for flitler2', ct
                endif

                if keyword_set(filter3_nobb) then begin
                    ind = where(acn_co_filter ge acn_co_tresh,ct)
                    if ct gt 0 then begin
                        ovoc_ptr[ind] = !VALUES.F_NAN
                        ovoc_other[ind] = !VALUES.F_NAN
                    endif
                    ;print, 'this is the data points to be deleted for flitler3', ct
                endif
            endif

            if keyword_set(codeployed) then remove_ind = where(finite(ovoc_ptr) and finite(ovoc_other)) else remove_ind = where(finite(ovoc_ptr))
;///////////////////////////////////////////////////////////////////////////////////////////
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;++                                      STEP 4                                           ++
;++                               Removing missing value                                  ++
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;///////////////////////////////////////////////////////////////////////////////////////////
            if remove_ind[0] ne -1 then begin

                ovoc_ptr = ovoc_ptr[remove_ind]
                ovoc_other = ovoc_other[remove_ind]
                gvoc_bop_top = gvoc_bop_top[remove_ind]
                gvoc_mami = gvoc_mami[remove_ind]
                gvoc_nobb = gvoc_nobb[remove_ind]
                gvoc_pbl65_ft35 = gvoc_pbl65_ft35[remove_ind]
                gvoc_sf_mami = gvoc_sf_mami[remove_ind]
                gvoc_surface = gvoc_surface[remove_ind]

                ;; addup
                
                ;for vertical profile
                zvar_obs=zvar_obs[remove_ind]
                zvar_gc_gfas=zvar_gc_gfas[remove_ind]

                ;for time series
                xtime_obs=xtime_obs[remove_ind] 
                xtime_gc = xtime_gc[remove_ind]
                
                ;for filter
                no2_filter=no2_filter[remove_ind]
                o3_filter=o3_filter[remove_ind]
                co_filter = co_filter[remove_ind]
                ch3cn_filter=ch3cn_filter[remove_ind]
                hcn_filter = hcn_filter[remove_ind]
                
                ;for voc vs. co
                oco_tmp=oco_tmp[remove_ind]
                gco_bop_top = co_gc_gfas_bop_top[remove_ind]
                gco_mami = co_gc_gfas_mami[remove_ind]
                gco_nobb = co_gc_gfas_nobb[remove_ind]
                gco_pbl65_ft35 = co_gc_gfas_pbl65_ft35[remove_ind]
                gco_sf_mami = co_gc_gfas_sf_mami[remove_ind]
                gco_surface = co_gc_gfas_surface[remove_ind]
            
                ;for voc vs. isop
                tmp_isop_ptr = tmp_isop_ptr[remove_ind]
                tmp_isop_toga = tmp_isop_toga[remove_ind]
                gisop_bop_top = c5h8_gc_gfas_bop_top[remove_ind]
                gisop_mami = c5h8_gc_gfas_mami[remove_ind]
                gisop_nobb = c5h8_gc_gfas_nobb[remove_ind]
                gisop_pbl65_ft35 = c5h8_gc_gfas_pbl65_ft35[remove_ind]
                gisop_sf_mami = c5h8_gc_gfas_sf_mami[remove_ind]
                gisop_surface = c5h8_gc_gfas_surface[remove_ind]

                ;for new criteria
                isop_filter = isop_filter[remove_ind]
                mek_filter = mek_filter[remove_ind]
                benz_filter = benz_filter[remove_ind]
                tolu_filter = tolu_filter[remove_ind]
                xyle_filter = xyle_filter[remove_ind]
                mono_filter = mono_filter[remove_ind]
                ;;tracks
                tmp_lat_obs = tmp_lat_obs[remove_ind]
                tmp_lon_obs = tmp_lon_obs[remove_ind]
                tmp_lat_gc = tmp_lat_gc[remove_ind]
                tmp_lon_gc = tmp_lon_gc[remove_ind]
            endif
            ;;set NA 
            if remove_ind[0] eq -1 then gvoc_bop_top = !values.f_nan
            if remove_ind[0] eq -1 then gvoc_mami = !values.f_nan
            if remove_ind[0] eq -1 then gvoc_nobb = !values.f_nan
            if remove_ind[0] eq -1 then gvoc_pbl65_ft35 = !values.f_nan
            if remove_ind[0] eq -1 then gvoc_sf_mami = !values.f_nan
            if remove_ind[0] eq -1 then gvoc_surface = !values.f_nan


            if remove_ind[0] eq -1 then ovoc_ptr = !values.f_nan
            if remove_ind[0] eq -1 then ovoc_other = !values.f_nan


;///////////////////////////////////////////////////////////////////////////////////////////
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;++                                      STEP 5                                           ++
;++                                  DO THE PLOTTING                                      ++
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;///////////////////////////////////////////////////////////////////////////////////////////

    ;; Vertical profiles
    ;; pressure
    ;; just use one altitude
            zvar_gc    = zvar_gc_gfas
            
            oVert_ptr  = bin_vert(Data=ovoc_ptr, ZData=zvar_obs, $
                                ZEdges=PEdges)
            oVert_other  = bin_vert(Data=ovoc_other, ZData=zvar_obs, $
                                ZEdges=PEdges)

            gVert_bop_top  = bin_vert(Data=gvoc_bop_top, ZData=zvar_gc, $
                                ZEdges=PEdges)
            gVert_mami  = bin_vert(Data=gvoc_mami, ZData=zvar_gc, $
                                ZEdges=PEdges)
            gVert_nobb  = bin_vert(Data=gvoc_nobb, ZData=zvar_gc, $
                                ZEdges=PEdges)
            gVert_pbl65_ft35  = bin_vert(Data=gvoc_pbl65_ft35, ZData=zvar_gc, $
                                ZEdges=PEdges)
            gVert_sf_mami  = bin_vert(Data=gvoc_sf_mami, ZData=zvar_gc, $
                                ZEdges=PEdges)             
                                
            gVert_surface  = bin_vert(Data=gvoc_surface, ZData=zvar_gc, $
                                ZEdges=PEdges)   
    
            if keyword_set(prs) then begin
                oVert_ptr = bin_vert(Data=ovoc_ptr, ZData = zvar_obs, $
                                    ZEdges=PEdges,/Press)
                                        
                oVert_other = bin_vert(Data=ovoc_other, ZData = zvar_obs, $
                                    ZEdges=PEdges,/Press)
                                    
                gVert_bop_top  = bin_vert(Data=gvoc_bop_top, ZData=zvar_gc, $
                                    ZEdges=PEdges,/Press)
                gVert_mami  = bin_vert(Data=gvoc_mami, ZData=zvar_gc, $
                                    ZEdges=PEdges,/Press)
                gVert_nobb  = bin_vert(Data=gvoc_nobb, ZData=zvar_gc, $
                                    ZEdges=PEdges,/Press)
                gVert_pbl65_ft35  = bin_vert(Data=gvoc_pbl65_ft35, ZData=zvar_gc, $
                                    ZEdges=PEdges,/Press)
                gVert_sf_mami  = bin_vert(Data=gvoc_sf_mami, ZData=zvar_gc, $
                                    ZEdges=PEdges,/Press)             
                gVert_surface  = bin_vert(Data=gvoc_surface, ZData=zvar_gc, $
                                    ZEdges=PEdges,/Press)   
            endif
            
            ;; charsize
            if keyword_set(primary) then charsize=2.2
            if keyword_set(secondary) then charsize = 1.2
            if keyword_set(all) or keyword_set(ovocs) then charsize = 1
            if keyword_set(chemistry) or keyword_set(nmhcs) then charsize = 2

            ;; err_thick
            if keyword_set(all) or keyword_set(ovocs) then err_thick = 4
            if keyword_set(each) then err_thick = 6

    ;; Plotting part
            ;if dd eq 0 then title = 'All SEAC4RS data' ;else title = 'Single flight: '+dates         
            ;plot,oVert_ptr.DataMean,oVert_ptr.zMean,col=1,xrange=xrange,yrange=[0,7],$;title=title,
            ;xtitle=xtitle,$;,ytitle='Altitude[km]',$     
            ;thick=8,charsize=charsize;,ystyle=4,$
            ;YTICKFORMAT="(A1)",XTICKFORMAT="(A1)",title=title
            
            if keyword_set(average) then begin
                obs_ptr = oVert_ptr.DataMean
                obs_other = oVert_other.DataMean
                                
                bop_top = gVert_bop_top.DataMean
                mami    = gVert_mami.DataMean
                nobb    = gVert_nobb.DataMean
                pbl65_ft35  = gVert_pbl65_ft35.DataMean
                sf_mami  = gVert_sf_mami.DataMean
                sf       = gVert_surface.DataMean                
            endif
            if keyword_set(med) then begin
                obs_ptr = oVert_ptr.DataMed
                obs_other = oVert_other.DataMed
                
                bop_top = gVert_bop_top.DataMed
                mami    = gVert_mami.DataMed
                nobb    = gVert_nobb.DataMed
                pbl65_ft35  = gVert_pbl65_ft35.DataMed
                sf_mami  = gVert_sf_mami.DataMed
                sf       = gVert_surface.DataMed       
            endif  
            

            
            
            ;; OUTPUTS
            ;  alt
            z_obs_ptr = oVert_ptr.zMean
            z_obs_other = oVert_other.zMean
            
            z_bop_top = gVert_bop_top.zMean
            z_mami    = gVert_mami.zMean
            z_nobb    = gVert_nobb.zMean
            z_pbl65_ft35  = gVert_pbl65_ft35.zMean
            z_sf_mami  = gVert_sf_mami.zMean
            z_sf       = gVert_surface.zMean       

            ;; PERCENTILE and NUMBER OF DATAPOINTS
            ;; PTR
            obs_ptr_q10 = oVert_ptr.DataQ10
            obs_ptr_q90 = oVert_ptr.DataQ90
            obs_ptr_q25 = oVert_ptr.DataQ25
            obs_ptr_q75 = oVert_ptr.DataQ75
            obs_ptr_numpts = oVert_ptr.NumPts
            
            ;; others
            obs_other_q10 = oVert_other.DataQ10
            obs_other_q90 = oVert_other.DataQ90
            obs_other_q25 = oVert_other.DataQ25
            obs_other_q75 = oVert_other.DataQ75
            obs_other_numpts = oVert_other.NumPts     
            
            ;; bop_top
            bop_top_q10 = gVert_bop_top.DataQ10
            bop_top_q90 = gVert_bop_top.DataQ90
            bop_top_q25 = gVert_bop_top.DataQ25
            bop_top_q75 = gVert_bop_top.DataQ75
            bop_top_numpts = gVert_bop_top.NumPts

            ;; mami
            mami_q10 = gVert_mami.DataQ10
            mami_q90 = gVert_mami.DataQ90
            mami_q25 = gVert_mami.DataQ25
            mami_q75 = gVert_mami.DataQ75
            mami_numpts = gVert_mami.NumPts

            ;; nobb
            nobb_q10 = gVert_nobb.DataQ10
            nobb_q90 = gVert_nobb.DataQ90
            nobb_q25 = gVert_nobb.DataQ25
            nobb_q75 = gVert_nobb.DataQ75
            nobb_numpts = gVert_nobb.NumPts
            
            ;; pbl65_ft35
            pbl65_ft35_q10 = gVert_pbl65_ft35.DataQ10
            pbl65_ft35_q90 = gVert_pbl65_ft35.DataQ90
            pbl65_ft35_q25 = gVert_pbl65_ft35.DataQ25
            pbl65_ft35_q75 = gVert_pbl65_ft35.DataQ75
            pbl65_ft35_numpts = gVert_pbl65_ft35.NumPts
            
            ;; sf_mami
            sf_mami_q10 = gVert_sf_mami.DataQ10
            sf_mami_q90 = gVert_sf_mami.DataQ90
            sf_mami_q25 = gVert_sf_mami.DataQ25
            sf_mami_q75 = gVert_sf_mami.DataQ75
            sf_mami_numpts = gVert_sf_mami.NumPts

            ;; sf
            sf_q10 = gVert_surface.DataQ10
            sf_q90 = gVert_surface.DataQ90
            sf_q25 = gVert_surface.DataQ25
            sf_q75 = gVert_surface.DataQ75
            sf_numpts = gVert_surface.NumPts


            ;; delete data point less than 20
            ;; PTR
            ind_num = where(obs_ptr_numpts lt 10, ct)
            if ct gt 0 then obs_ptr[ind_num] = !VALUES.F_NAN
            
            ;; other
            ind_num = where(obs_other_numpts lt 5, ct)
            if ct gt 0 then obs_other[ind_num] = !VALUES.F_NAN
            
            ind =  where(finite(obs_ptr),ct)
            if ct gt 0 then begin
                obs_ptr = obs_ptr[ind]
                obs_other = obs_other[ind]
                
                bop_top = bop_top[ind]
                mami  = mami[ind]
                nobb  = nobb[ind]
                pbl65_ft35  = pbl65_ft35[ind]
                sf_mami  = sf_mami[ind]
                sf = sf[ind]

                z_obs_ptr = z_obs_ptr[ind]
                z_obs_other = z_obs_other[ind]
                                
                z_bop_top = z_bop_top[ind]
                z_mami    = z_mami[ind]
                z_nobb    = z_nobb[ind]
                z_pbl65_ft35  = z_pbl65_ft35[ind]
                z_sf_mami  = z_sf_mami[ind]
                z_sf       = z_sf[ind]
                
                obs_ptr_q10 = obs_ptr_q10[ind]
                obs_ptr_q90 = obs_ptr_q90[ind]
                obs_ptr_q25 = obs_ptr_q25[ind]
                obs_ptr_q75 = obs_ptr_q75[ind]
                obs_ptr_numpts = obs_ptr_numpts[ind]
                
                obs_other_q10 = obs_other_q10[ind]
                obs_other_q90 = obs_other_q90[ind]
                obs_other_q25 = obs_other_q25[ind]
                obs_other_q75 = obs_other_q75[ind]
                obs_other_numpts = obs_other_numpts[ind]
                
                bop_top_q10 = bop_top_q10[ind]
                bop_top_q90 = bop_top_q90[ind]
                bop_top_q25 = bop_top_q25[ind]
                bop_top_q75 = bop_top_q75[ind]
                bop_top_numpts = bop_top_numpts[ind]
                
                mami_q10 = mami_q10[ind]
                mami_q90 = mami_q90[ind]
                mami_q25 = mami_q25[ind]
                mami_q75 = mami_q75[ind]
                mami_numpts = mami_numpts[ind]
                
                nobb_q10 = nobb_q10[ind]
                nobb_q90 = nobb_q90[ind]
                nobb_q25 = nobb_q25[ind]
                nobb_q75 = nobb_q75[ind]
                nobb_numpts = nobb_numpts[ind]
                
                pbl65_ft35_q10 = pbl65_ft35_q10[ind]
                pbl65_ft35_q90 = pbl65_ft35_q90[ind]
                pbl65_ft35_q25 = pbl65_ft35_q25[ind]
                pbl65_ft35_q75 = pbl65_ft35_q75[ind]
                pbl65_ft35_numpts = pbl65_ft35_numpts[ind]
                
                sf_mami_q10 = sf_mami_q10[ind]
                sf_mami_q90 = sf_mami_q90[ind]
                sf_mami_q25 = sf_mami_q25[ind]
                sf_mami_q75 = sf_mami_q75[ind]
                sf_mami_numpts = sf_mami_numpts[ind]
                
                sf_q10 = sf_q10[ind]
                sf_q90 = sf_q90[ind]
                sf_q25 = sf_q25[ind]
                sf_q75 = sf_q75[ind]
                sf_numpts = sf_numpts[ind]

            endif
            
            ;; plotting setting
            
            if keyword_set(all) or keyword_set(ovocs) then begin
                charsize = 2.5
                charthick = 4

                thick1 = 20
                thick2 = 12
            endif

        
            if keyword_set(each) then begin
                charsize = 3
                charthick = 8
                thick1 = 20
                thick2 = 15
            endif


            
            yrange = [1000,400]
            

            
            if keyword_set(no_obs_include) then begin
                obs_ptr = bop_top
                z_obs_ptr = z_bop_top
                obs_ptr_q25 = bop_top_q25
                obs_ptr_q75 = bop_top_q75
            endif
            
            if keyword_set(all) or keyword_set(ovocs) then begin
                if s eq 0 or s eq 3 or s eq 6 then begin                
                    if keyword_set(errorbar) then begin
                        ;; coyote library by using error bar to represent percentile lines
                        ;; sources: http://www.idlcoyote.com/idldoc/cg/cgplot.html
                        cgplot,obs_ptr,z_obs_ptr,col=1,xrange=xrange,yrange=yrange,$;, title=title,,$; title=title,$
                        ;xtitle=xtitle,$;,ytitle='Altitude[km]',$     
                        thick=thick1,charsize=charsize,charthick=charthick, $
                        ERR_XLow=(obs_ptr-obs_ptr_q25), ERR_XHigh=(obs_ptr_q75-obs_ptr), $
                        ERR_Color='black', $
                        ERR_THICK = ERR_THICK, $
                        /ERR_CLIP
                        ;YTICKFORMAT="(A1)", XTICKFORMAT="(A1)"
                        ;,ystyle=4,$

                    endif else begin
                        ;; default setting by using other lines to represent 25th/75th percentile
                        plot,obs_ptr,z_obs_ptr,col=1,xrange=xrange,yrange=yrange,$;, title=title,,$; title=title,$
                        ;xtitle=xtitle,$;,ytitle='Altitude[km]',$     
                        thick=thick1,charsize=charsize,charthick=charthick
                        ;YTICKFORMAT="(A1)", XTICKFORMAT="(A1)"
                        ;,ystyle=4,$
                    endelse
                endif else begin
                    if keyword_set(errorbar) then begin
                        ;; coyote library by using error bar to represent percentile lines
                        ;; sources: http://www.idlcoyote.com/idldoc/cg/cgplot.html
                        cgplot,obs_ptr,z_obs_ptr,col=1,xrange=xrange,yrange=yrange,$;, title=title,,$; title=title,$
                        ;xtitle=xtitle,$;,ytitle='Altitude[km]',$     
                        thick=thick1,charsize=charsize,charthick=charthick, $
                        YTICKFORMAT="(A1)",$
                        ERR_XLow=(obs_ptr-obs_ptr_q25), ERR_XHigh=(obs_ptr_q75-obs_ptr), $
                        ERR_Color='black', $
                        ERR_THICK = ERR_THICK, $
                        /ERR_CLIP
                        ;,ystyle=4,$
                    endif else begin
                        ; default setting by using other lines to represent 25th/75th percentile
                        plot,obs_ptr,z_obs_ptr,col=1,xrange=xrange,yrange=yrange,$;, title=title,,$; title=title,$
                        ;xtitle=xtitle,$;,ytitle='Altitude[km]',$     
                        thick=thick1,charsize=charsize,charthick=charthick, $
                        YTICKFORMAT="(A1)"
                        ;,ystyle=4,$
                    endelse
                endelse
            endif
            
        
            if keyword_set(each) then begin
                if keyword_set(errorbar) and not keyword_set(no_obs_include) then begin
                ;if keyword_set(errorbar) then begin
                    ;; coyote library by using error bar to represent percentile lines
                    ;; sources: http://www.idlcoyote.com/idldoc/cg/cgplot.html
                    cgplot,obs_ptr,z_obs_ptr,col=1,xrange=xrange,yrange=yrange,$;, title=title,,$; title=title,$
                    ;xtitle=xtitle,$;,ytitle='Altitude[km]',$     
                    thick=thick1,charsize=charsize,charthick=charthick, $
                    ERR_XLow=(obs_ptr-obs_ptr_q25), ERR_XHigh=(obs_ptr_q75-obs_ptr), $
                    ERR_Color='black', $
                    ERR_THICK = ERR_THICK, $
                    /ERR_CLIP
                    ;YTICKFORMAT="(A1)", XTICKFORMAT="(A1)"
                    ;,ystyle=4,$

                endif else begin
                    ;; default setting by using other lines to represent 25th/75th percentile
                    plot,obs_ptr,z_obs_ptr,col=1,xrange=xrange,yrange=yrange,$;, title=title,,$; title=title,$
                    ;xtitle=xtitle,$;,ytitle='Altitude[km]',$     
                    thick=thick1,charsize=charsize,charthick=charthick
                    ;YTICKFORMAT="(A1)", XTICKFORMAT="(A1)"
                    ;,ystyle=4,$
                endelse
            endif
            
                
;; only set sepcies implemented in the model
            oplot,bop_top,z_bop_top,col=2,thick=thick1,line=0
            oplot,mami,z_mami,col=3,thick=thick1,line=0 
            ;oplot,nobb,z_nobb,col=4,thick=thick1,line=0
            oplot,pbl65_ft35,z_pbl65_ft35,col=5,thick=thick1,line=0
            oplot,sf_mami,z_sf_mami,col=6,thick=thick1,line=0
            oplot,sf,z_sf,col=7,thick=thick1,line=0
            
            
            test = 0 
            if test eq 1 then begin
                cgplot,bop_top,z_bop_top,col=2,xrange=xrange,yrange=yrange,$;, title=title,,$; title=title,$
                        ;xtitle=xtitle,$;,ytitle='Altitude[km]',$     
                        thick=thick1,charsize=charsize,charthick=charthick, $
                        ERR_XLow=(bop_top-bop_top_q25), ERR_XHigh=(bop_top_q75-bop_top), $
                        ERR_Color=2, $
                        ERR_THICK = ERR_THICK, $
                        /ERR_CLIP,/OVERPLOT

                cgplot,mami,z_mami,col=3,xrange=xrange,yrange=yrange,$;, title=title,,$; title=title,$
                        ;xtitle=xtitle,$;,ytitle='Altitude[km]',$     
                        thick=thick1,charsize=charsize,charthick=charthick, $
                        ERR_XLow=(mami-mami_q25), ERR_XHigh=(mami_q75-mami), $
                        ERR_Color=3, $
                        ERR_THICK = ERR_THICK, $
                        /ERR_CLIP,/OVERPLOT

                ;cgplot,nobb,z_nobb,col=4,xrange=xrange,yrange=yrange,$;, title=title,,$; title=title,$
                ;        ;xtitle=xtitle,$;,ytitle='Altitude[km]',$     
                ;        thick=thick1,charsize=charsize,charthick=charthick, $
                ;        ERR_XLow=(nobb-nobb_q25), ERR_XHigh=(nobb_q75-nobb), $
                ;        ERR_Color=4, $
                ;        ERR_THICK = ERR_THICK, $
                ;        /ERR_CLIP,/OVERPLOT

                cgplot,pbl65_ft35,z_pbl65_ft35,col=5,xrange=xrange,yrange=yrange,$;, title=title,,$; title=title,$
                        ;xtitle=xtitle,$;,ytitle='Altitude[km]',$     
                        thick=thick1,charsize=charsize,charthick=charthick, $
                        ERR_XLow=(pbl65_ft35-pbl65_ft35_q25), ERR_XHigh=(pbl65_ft35_q75-pbl65_ft35), $
                        ERR_Color=5, $
                        ERR_THICK = ERR_THICK, $
                        /ERR_CLIP,/OVERPLOT

                cgplot,sf_mami,z_sf_mami,col=6,xrange=xrange,yrange=yrange,$;, title=title,,$; title=title,$
                        ;xtitle=xtitle,$;,ytitle='Altitude[km]',$     
                        thick=thick1,charsize=charsize,charthick=charthick, $
                        ERR_XLow=(sf_mami-sf_mami_q25), ERR_XHigh=(sf_mami_q75-sf_mami), $
                        ERR_Color=6, $
                        ERR_THICK = ERR_THICK, $
                        /ERR_CLIP,/OVERPLOT

                cgplot,sf,z_sf,col=7,xrange=xrange,yrange=yrange,$;, title=title,,$; title=title,$
                        ;xtitle=xtitle,$;,ytitle='Altitude[km]',$     
                        thick=thick1,charsize=charsize,charthick=charthick, $
                        ERR_XLow=(sf-sf_q25), ERR_XHigh=(sf_q75-sf), $
                        ERR_Color=7, $
                        ERR_THICK = ERR_THICK, $
                        /ERR_CLIP,/OVERPLOT
            endif





            if keyword_set(each) then begin
                charsize =2
                charthick =5
            endif 
            if keyword_set(all) or keyword_set(ovocs) then begin
                charsize =1.1
                charthick =3
            endif 
            ;; test
            ;; change the postition later.
            ;; xyouts for NumPoints for each layers            
            for num = 0, n_elements(obs_ptr_numpts)-1 do begin
                if z_obs_ptr[num] ge yrange[0] then continue            
                xyouts, 1.015*max(xrange),z_obs_ptr[num], STRTRIM(string(obs_ptr_numpts[num]), 1), /data, col=1, charsize=charsize, $
                    charthick=charthick, alignment=0.      
            endfor
        
;///////////////////////////////////////////////////////////////////////////////////////////
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;++                                      STEP 6                                           ++
;++                                   PLOTTING xyouts                                     ++
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;///////////////////////////////////////////////////////////////////////////////////////////
        ;; VOCs names, add title to each plot
        ;; add bar
            if keyword_set(each) then begin
                ;xpos1 = 0.88
                xpos1 = 0.7
                ypos1 = [0.9,0.85,0.80,0.75,0.70,0.65,0.60,0.55]
                ;ypos1 = [0.9,0.87,0.84,0.81,0.78,0.75,0.72]
                charsize= 3
                charthick = 7
            endif
            if keyword_set(each) then begin
                xyouts, xpos1,ypos1[1], 'Observation', /normal, col=1, charsize=charsize, $
                    charthick=charthick, alignment=0. 
                xyouts, xpos1,ypos1[2], 'Bop2top', /normal, col=2, charsize=charsize, $
                    charthick=charthick, alignment=0
                xyouts, xpos1,ypos1[3], 'MAMI', /normal, col=3, charsize=charsize, $
                    charthick=charthick, alignment=0. 
                ;xyouts, xpos1,ypos1[4], 'NOBB', /normal, col=4, charsize=charsize, $
                ;    charthick=charthick, alignment=0.
                xyouts, xpos1,ypos1[4], 'PBL65_FT35', /normal, col=5, charsize=charsize, $
                    charthick=charthick, alignment=0. 
                xyouts, xpos1,ypos1[5], 'Sf2mami', /normal, col=6, charsize=charsize, $
                    charthick=charthick, alignment=0. 
                xyouts, xpos1,ypos1[6], 'Surface', /normal, col=7, charsize=charsize, $
                    charthick=charthick, alignment=0. 
            endif
        endfor ;; s loop
    endfor ;; dd loop
    close_device
    ;DEVICE, /CLOSE
    print,'done!'
end
