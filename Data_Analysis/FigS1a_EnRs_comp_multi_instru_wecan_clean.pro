;; modified by lixu, 10/17/2019, used for compare the profiles from wecan and gc 
;; original script mk_flight_v2, lhu
;; modicication: parameters, set missing value for observation data, cancel filter. 
;@tvmap
; This is used to answer coauthor's comments

@org_boot
pro EnRs_comp_multi_instru_wecan_clean,$
        plume=plume,$ 
        emipass = emipass, filter2_nobb=filter2_nobb, filter3_nobb = filter3_nobb,$
        all=all,$
        test=test, save=save,$
        nested=nested,fbf=fbf,$
        mod_only = mod_only, $
        raw_60s = raw_60s, avg_5min=avg_5min, TOGA_merge=TOGA_merge

;myct,/WhGrYlRd
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

  ;; need to check it later!!!!!!!!
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
  avo= 6.022e23  ;; molec/mol

  ; PLOT SETTING
  ;;v2 means add for those species have both two measurement
    fi = 'emipass_ERs'
    if keyword_set(nested) then fi = fi + '_nested'
    if keyword_set(fbf) then fi = fi + '_fbf'
    if keyword_set(test) then fi = 'test_emipass'  
    if keyword_set(save) then $
    open_device, /ps, /color, /landscape, $
               fi = './ps/' + fi + '.ps'
    if n_elements(region) eq 0 then region=''
    case strlowcase(region) of 
        'world'    : limit = [0,-160,90,160]
        'southeast': limit=[25,-100,40,-75]
        'se'       : limit=[25,-100,40,-75]
        'west'     : limit=[36,-127,49.5,-105]
        ;lixu, 10/10/2019
        'na'       : limit=[10,-140,70,-70]    
        else       : limit=[25,-127,52,-65]
    endcase               

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
    
    if keyword_set(TOGA_merge) then     dates_all=  ['20180724','20180726','20180730','20180731','20180802',$
                                                     '20180803','20180806','20180808','20180809','20180813','20180815','20180816','20180820',$
                                                     '20180823','20180826','20180828']
    
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
    if keyword_set(emipass) then begin
        dates_all = dates_emipass
        enter_time = start_time_emipass
        exit_time = end_time_emipass
        help, dates_all
    endif
   
    avo= 6.022e23  ;; molec/mol
    
    for dd = 0, n_elements(dates_all) do begin
        ;;only consider whole flight
        if dd ne 0 then break
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

    ;; add vocs
        c3h8_obs = [0]
        c3h6_obs = [0]
        c2h5oh_obs = [0]
        c5h8_obs   = [0]
        c7h8_obs   = [0]
        dms_obs    = [0]
        mek_obs    = [0]
        hcn_obs    = [0]
        macr_mvk_obs = [0]

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
            if keyword_set(avg_5min) then c130fi = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                                                    'mrg2sav/WECAN/R4_merges_avg/wecan-mrg5m-c130_merge_'+dates[n]+'_R4.sav'
                
            if keyword_set(raw_60s) then c130fi = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                                                    'mrg2sav/WECAN/R4_merges_raw/wecan-mrg60-c130_merge_'+dates[n]+'_R4.sav'

            if keyword_set(TOGA_merge) then c130fi = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                                                    'mrg2sav/WECAN/TOGA_merges_raw/wecan-mrgTOGA-c130_merge_'+dates[n]+'_R4.sav'
                                                    
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

    ;; add VOCs
            tmp_c3h8_obs = c130.PROPANE_TOGA/1000;change pptv into ppbv
            ;c3h6_obs = [c3h6_obs,c130.PROPANE_TOGA]
            tmp_c2h5oh_obs = c130.c2h5oh_toga/1000 ; convert pptv into ppbv
            tmp_c5h8_obs   = c130.ISOPRENE_MIXINGRATIO_PTR
            tmp_c7h8_obs   = c130.TOLUENE_MIXINGRATIO_PTR
            tmp_dms_obs    = c130.DMS_MIXINGRATIO_PTR
            tmp_mek_obs    = c130.MEK_MIXINGRATIO_PTR*0.8
            tmp_hcn_obs = c130.HCN_TOGA/1000 ; pptv to ppbv

            tmp_macr_mvk_obs = c130.MACR_MVK_MixingRatio_PTR 

    ;; PTR and TOGA
    ;; PTR    
            tmp_isop_obs_ptr = c130.Isoprene_MixingRatio_PTR
            tmp_acet_obs_ptr = c130.Acetone_MixingRatio_PTR*0.78
            tmp_propanal_obs_ptr = c130.Acetone_MixingRatio_PTR*0.22
            tmp_mek_obs_ptr  = c130.MEK_MixingRatio_PTR*0.8
            tmp_butanal_obs_ptr  = c130.MEK_MixingRatio_PTR*0.2
            tmp_macr_mvk_obs_ptr     = c130.MACR_MVK_MixingRatio_PTR
            tmp_ALD2_obs_ptr = c130.Acetaldehyde_MixingRatio_PTR
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
            ;; add VOCs
                    tmp_c3h8_obs = tmp_c3h8_obs[remove_ind];change pptv into ppbv
                    
                    ;c3h6_obs = [c3h6_obs,c130.PROPANE_TOGA]
                    tmp_c2h5oh_obs = tmp_c2h5oh_obs[remove_ind] ; convert pptv into ppbv
                    tmp_c5h8_obs   = tmp_c5h8_obs[remove_ind]
                    tmp_c7h8_obs   = tmp_c7h8_obs[remove_ind]
                    tmp_dms_obs    = tmp_dms_obs[remove_ind]
                    tmp_mek_obs    = tmp_mek_obs[remove_ind]
                    tmp_hcn_obs =tmp_hcn_obs[remove_ind]

                    tmp_macr_mvk_obs = tmp_macr_mvk_obs[remove_ind]


            ;; PTR and TOGA
            ;; PTR    
                    tmp_isop_obs_ptr = tmp_isop_obs_ptr[remove_ind]
                    tmp_acet_obs_ptr = tmp_acet_obs_ptr[remove_ind]
                    tmp_propanal_obs_ptr = tmp_propanal_obs_ptr[remove_ind]
               
                    tmp_mek_obs_ptr  = tmp_mek_obs_ptr[remove_ind]
                    tmp_butanal_obs_ptr  = tmp_butanal_obs_ptr[remove_ind]
                    
                    tmp_macr_mvk_obs_ptr     = tmp_macr_mvk_obs_ptr[remove_ind]
                    tmp_ALD2_obs_ptr = tmp_ALD2_obs_ptr[remove_ind]
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


            if keyword_set(emipass) or keyword_set(young) or $
                keyword_set(intermed) or keyword_set(GT4) or keyword_set(LT4) or $
                keyword_set(nonsmoke) or keyword_set(aged) then begin
                ind_enter_time = where(abs(tmp_utc_obs - enter_time[n]) le 30,ct1)
                ind_exit_time = where(abs(tmp_utc_obs - exit_time[n]) le 30, ct2)
                
                print, ind_enter_time[0] , ind_exit_time[-1]
                if ct1 gt 0 and ct2 gt 0 then begin
                    tmp_co_obs = tmp_co_obs[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_o3_obs =  tmp_o3_obs[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_pan_obs = tmp_pan_obs[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_hcho_obs = tmp_hcho_obs[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_acet_obs = tmp_acet_obs[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_benz_obs = tmp_benz_obs[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_ch3oh_obs = tmp_ch3oh_obs[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_ald2_obs = tmp_ald2_obs[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_no_obs = tmp_no_obs[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_no2_obs =tmp_no2_obs[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_so2_obs = tmp_so2_obs[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_lat_obs  = tmp_lat_obs[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_lon_obs  = tmp_lon_obs[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_prs_obs = tmp_prs_obs[ind_enter_time[0]:ind_exit_time[-1]]  ;; Pa
                    tmp_alt_obs = tmp_alt_obs[ind_enter_time[0]:ind_exit_time[-1]]

                    ;noy_obs   = [noy_obs,dc8.NOy_ESRL]
                    tmp_ch3cn_obs = tmp_ch3cn_obs[ind_enter_time[0]:ind_exit_time[-1]]; we also have mixing ratio of ch3cn, need to check the unit lixu
                    tmp_jday_obs = tmp_jday_obs[ind_enter_time[0]:ind_exit_time[-1]]
                    ;tmp_utc_obs = tmp_utc_obs[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_lstime_obs = tmp_lstime_obs[ind_enter_time[0]:ind_exit_time[-1]]
            ;; add VOCs
                    tmp_c3h8_obs = tmp_c3h8_obs[ind_enter_time[0]:ind_exit_time[-1]];change pptv into ppbv
                    ;c3h6_obs = [c3h6_obs,c130.PROPANE_TOGA]
                    tmp_c2h5oh_obs = tmp_c2h5oh_obs[ind_enter_time[0]:ind_exit_time[-1]] ; convert pptv into ppbv
                    tmp_c5h8_obs   = tmp_c5h8_obs[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c7h8_obs   = tmp_c7h8_obs[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_dms_obs    = tmp_dms_obs[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_mek_obs    = tmp_mek_obs[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_hcn_obs =tmp_hcn_obs[ind_enter_time[0]:ind_exit_time[-1]]

                    tmp_macr_mvk_obs = tmp_macr_mvk_obs[ind_enter_time[0]:ind_exit_time[-1]]
                    


            ;; PTR and TOGA
            ;; PTR    
                    tmp_isop_obs_ptr = tmp_isop_obs_ptr[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_acet_obs_ptr = tmp_acet_obs_ptr[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_propanal_obs_ptr = tmp_propanal_obs_ptr[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_mek_obs_ptr  = tmp_mek_obs_ptr[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_butanal_obs_ptr = tmp_butanal_obs_ptr[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_macr_mvk_obs_ptr     = tmp_macr_mvk_obs_ptr[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_ALD2_obs_ptr = tmp_ALD2_obs_ptr[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_ch2o_obs_ptr = tmp_ch2o_obs_ptr[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_benz_obs_ptr = tmp_benz_obs_ptr[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_tolu_obs_ptr = tmp_tolu_obs_ptr[ind_enter_time[0]:ind_exit_time[-1]]

                    tmp_c8h10_obs_ptr        = tmp_c8h10_obs_ptr[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_Xylenes_obs_ptr      = tmp_Xylenes_obs_ptr[ind_enter_time[0]:ind_exit_time[-1]]
                    
                    tmp_monoterpenes_obs_ptr = tmp_monoterpenes_obs_ptr[ind_enter_time[0]:ind_exit_time[-1]]
            
                    tmp_moh_obs_ptr = tmp_moh_obs_ptr[ind_enter_time[0]:ind_exit_time[-1]]
                    
                    tmp_hcooh_obs_ptr  = tmp_hcooh_obs_ptr[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_acta_obs_ptr   = tmp_acta_obs_ptr[ind_enter_time[0]:ind_exit_time[-1]]
            ; TOGA
                    tmp_isop_obs_toga = tmp_isop_obs_toga[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_acet_obs_toga =tmp_acet_obs_toga[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_mek_obs_toga  = tmp_mek_obs_toga[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_macr_mvk_obs_toga     = tmp_macr_mvk_obs_toga[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_ALD2_obs_toga =tmp_ALD2_obs_toga[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_ch2o_obs_toga = tmp_ch2o_obs_toga[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_benz_obs_toga = tmp_benz_obs_toga[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_tolu_obs_toga = tmp_tolu_obs_toga[ind_enter_time[0]:ind_exit_time[-1]]

                    tmp_mbo_obs_toga  =tmp_mbo_obs_toga[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_propanal_obs_toga =tmp_propanal_obs_toga[ind_enter_time[0]:ind_exit_time[-1]]

                    ;tmp_c2Butenal_obs_toga =tmp_c2Butenal_obs_toga[ind_enter_time[0]:ind_exit_time[-1]]
                    ;tmp_t2Butenal_obs_toga = tmp_t2Butenal_obs_toga[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_butanal_obs_toga = tmp_butanal_obs_toga[ind_enter_time[0]:ind_exit_time[-1]]

                    tmp_etbenzene_obs_toga   = tmp_etbenzene_obs_toga[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_tricyclene_obs_toga  = tmp_tricyclene_obs_toga[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_apinene_obs_toga     =tmp_apinene_obs_toga[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_camphene_obs_toga  = tmp_camphene_obs_toga[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_bpinenemyrcene_obs_toga   = tmp_bpinenemyrcene_obs_toga[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_limonened3carene_obs_toga = tmp_limonened3carene_obs_toga[ind_enter_time[0]:ind_exit_time[-1]]

                    tmp_c8h10_obs_toga     = tmp_c8h10_obs_toga[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_Xylenes_obs_toga   = tmp_Xylenes_obs_toga[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_MPXYLENE_obs_toga  = tmp_MPXYLENE_obs_toga[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_OXYLENE_obs_toga   = tmp_OXYLENE_obs_toga[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_EtBenzene_obs_toga = tmp_EtBenzene_obs_toga[ind_enter_time[0]:ind_exit_time[-1]]
                                        
                    tmp_moh_obs_toga = tmp_moh_obs_toga[ind_enter_time[0]:ind_exit_time[-1]]
                    
                    tmp_rhum = tmp_rhum[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_temperature_obs = tmp_temperature_obs[ind_enter_time[0]:ind_exit_time[-1]]

                ;;for alk4 and prpe
                ;;alk4
                    tmp_propane_obs_toga = tmp_propane_obs_toga[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_butane_obs_toga = tmp_butane_obs_toga[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_pentane_obs_toga = tmp_pentane_obs_toga[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_mepentane_obs_toga = tmp_mepentane_obs_toga[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_hexane_obs_toga = tmp_hexane_obs_toga[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_trimepentane_obs_toga = tmp_trimepentane_obs_toga[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_heptane_obs_toga = tmp_heptane_obs_toga[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_octane_obs_toga = tmp_octane_obs_toga[ind_enter_time[0]:ind_exit_time[-1]]

                    tmp_alk4_obs_toga = tmp_alk4_obs_toga[ind_enter_time[0]:ind_exit_time[-1]]

                    ;;prpe
                    tmp_prpe_obs_toga = tmp_prpe_obs_toga[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_butene_obs_toga = tmp_butene_obs_toga[ind_enter_time[0]:ind_exit_time[-1]]

                ;;AWAS
                    tmp_c2h2_obs_awas = tmp_c2h2_obs_awas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c2h4_obs_awas   = tmp_c2h4_obs_awas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c2h6_obs_awas   = tmp_c2h6_obs_awas[ind_enter_time[0]:ind_exit_time[-1]]

                    ;; PRPE
                    tmp_c3h6_obs_awas = tmp_c3h6_obs_awas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_trans2butene_obs_awas = tmp_trans2butene_obs_awas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_x1Butene_obs_awas = tmp_x1Butene_obs_awas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_cisx2butene_obs_awas = tmp_cisx2butene_obs_awas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_x1Hexene_obs_awas = tmp_x1Hexene_obs_awas[ind_enter_time[0]:ind_exit_time[-1]]
                    
                    tmp_prpe_obs_awas = tmp_prpe_obs_awas[ind_enter_time[0]:ind_exit_time[-1]]
                    ;; ALK4
                    tmp_cyclopentane_obs_awas = tmp_cyclopentane_obs_awas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_cyclohexane_obs_awas = tmp_cyclohexane_obs_awas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_methylcyclohexane_obs_awas = tmp_methylcyclohexane_obs_awas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_isobutane_obs_awas = tmp_isobutane_obs_awas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_nbutane_obs_awas = tmp_nbutane_obs_awas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_isopentane_obs_awas = tmp_isopentane_obs_awas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_nPentane_obs_awas = tmp_nPentane_obs_awas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_hexane_obs_awas = tmp_hexane_obs_awas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_x24Dimethylpentane_obs_awas = tmp_x24Dimethylpentane_obs_awas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_x22Dimethylbutane_obs_awas = tmp_x22Dimethylbutane_obs_awas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_x3Methylpentane_obs_awas = tmp_x3Methylpentane_obs_awas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_x23Dimethylpentane_obs_awas = tmp_x23Dimethylpentane_obs_awas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_x2Methylhexane_obs_awas = tmp_x2Methylhexane_obs_awas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_x3Methylhexane_obs_awas = tmp_x3Methylhexane_obs_awas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_nHeptan_obs_awas       = tmp_nHeptan_obs_awas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_x234Trimethylpentane_obs_awas = tmp_x234Trimethylpentane_obs_awas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_x2Methylheptane_obs_awas = tmp_x2Methylheptane_obs_awas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_x3Methylheptane_obs_awas = tmp_x3Methylheptane_obs_awas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_nOctane_obs_awas = tmp_nOctane_obs_awas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_nNonane_obs_awas= tmp_nNonane_obs_awas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_undecane_obs_awas = tmp_undecane_obs_awas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_ndecane_obs_awas = tmp_ndecane_obs_awas[ind_enter_time[0]:ind_exit_time[-1]]
                    
                    tmp_alk4_obs_awas = tmp_alk4_obs_awas[ind_enter_time[0]:ind_exit_time[-1]]
                    
                    ;; CIMS
                    tmp_hcooh_obs_cims  = tmp_hcooh_obs_cims[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_acta_obs_cims   = tmp_acta_obs_cims[ind_enter_time[0]:ind_exit_time[-1]]
                endif
                if ct1 le 0 or ct2 le 0 then continue
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

    ;; add VOCs
            c3h8_obs = [c3h8_obs,tmp_c3h8_obs];change pptv into ppbv
        
            ;c3h6_obs = [c3h6_obs,c130.PROPANE_TOGA]
            c2h5oh_obs = [c2h5oh_obs,tmp_c2h5oh_obs] ; convert pptv into ppbv
            c5h8_obs   = [c5h8_obs,tmp_c5h8_obs]
            c7h8_obs   = [c7h8_obs,tmp_c7h8_obs]
            dms_obs    = [dms_obs,tmp_dms_obs]
            mek_obs    = [mek_obs,tmp_mek_obs]
            hcn_obs = [hcn_obs,tmp_hcn_obs] ; pptv to ppbv

            macr_mvk_obs = [macr_mvk_obs,tmp_macr_mvk_obs] 


    ;; PTR and TOGA
    ;; PTR    
            isop_obs_ptr = [isop_obs_ptr,tmp_isop_obs_ptr]
            acet_obs_ptr = [acet_obs_ptr,tmp_acet_obs_ptr]
            propanal_obs_ptr = [propanal_obs_ptr,tmp_propanal_obs_ptr]
            mek_obs_ptr  = [mek_obs_ptr,tmp_mek_obs_ptr]
            butanal_obs_ptr  = [butanal_obs_ptr,tmp_butanal_obs_ptr]
            macr_mvk_obs_ptr     = [macr_mvk_obs_ptr,tmp_macr_mvk_obs_ptr]
            ALD2_obs_ptr = [ALD2_obs_ptr,tmp_ALD2_obs_ptr]
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

    ;; add VOCs
        c3h8_obs = c3h8_obs[1:*]
        ;c3h6_obs = c3h6_obs[1:*]
        c2h5oh_obs = c2h5oh_obs[1:*]
        ch3oh_obs = ch3oh_obs[1:*]
        c5h8_obs   = c5h8_obs[1:*]
        c7h8_obs   = c7h8_obs[1:*]
        dms_obs    = dms_obs[1:*]
        mek_obs    = mek_obs[1:*]
        hcn_obs    = hcn_obs[1:*]
        macr_mvk_obs = macr_mvk_obs[1:*]


    ;; PTR and TOGA
    ;; PTR    
        isop_obs_ptr = isop_obs_ptr[1:*]
        acet_obs_ptr = acet_obs_ptr[1:*]
        propanal_obs_ptr = propanal_obs_ptr[1:*]
        mek_obs_ptr  = mek_obs_ptr[1:*]
        butanal_obs_ptr  = butanal_obs_ptr[1:*]
        macr_mvk_obs_ptr     = macr_mvk_obs_ptr[1:*]
        ALD2_obs_ptr = ALD2_obs_ptr[1:*]
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

            
;; ==============================
;; set up detection limit for PTR
;; not for missing values
;; ==============================   
        ;if keyword_set(lod_ptr) then begin
            LoD = 50.0/1000
            ind = where(ch2o_obs_ptr le LoD and ch2o_obs_ptr gt 0, ct)
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
            ind = where(tolu_obs_ptr le LoD and tolu_obs_ptr gt 0, ct)
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
        ;endif
;; ========================
;; Prepare for the filters 
;; ========================
    ;; lixu, 12/18/2019, 0208/2021   
    ;; plume filter setting, use in the plotting part
        no2_thresh = 4. ; ppb
        co_thresh  = 150. ; ppb
    ;; biomass burning threshhold 
        hcn_thresh   = 275./1000 ;pptv to ppbv
        ch3cn_thresh = 155/1000. ;;ppt to ppb, 131
        co_thresh =  100
    ;; strat filter
        strat_thresh=1.25 ; O3/CO ratio
        
    ;; CO threshold of fresh and aged smoke 
        ;; young
        methylfuran_thresh = 0.7/1000 
        ;; median
        acrolein_thresh = 7.4/1000
        ;; old
        acrylonitrile_thresh = 2.9/1000
    ;; for vertical profile
    ;    yrange = [0,7]
    ;; lxu, 06202021
        acn_co_tresh = 2.01 ;ppb/ppm
;; ========keep the same for regression, profiles and time series=================
        if keyword_set(all) then rows = 3
        if keyword_set(all) then cols = 3

        multipanel, rows = rows, cols = cols
        n_species = rows*cols

        ; initialize save-out data
        spec_obs_total   = []
        spec_obs_std_total   = []
        corr_obs_total   = []
        title_all = []
        
        for s = 0, n_species - 1 do begin      
            if keyword_set(all) then begin
                case s of
                    ;;hcho
                    0:begin
                        ovoc_ptr = ch2o_obs_ptr
                        ovoc_other = ch2o_obs_toga
                        title = 'Formaldehyde'
                        if plume eq 0 then begin
                            xrange = [0, 120]
                            yrange = [0, 120]
                        endif
                        if plume eq 1 then begin
                            xrange = [0, 5]
                            yrange = [0, 5]
                        endif
                    end
                    ;;ald2
                    1:begin
                        ovoc_ptr = ald2_obs_ptr
                        ovoc_other = ald2_obs_toga
                        title = 'Acetaldehyde'
                        if plume eq 0 then begin
                            xrange = [0, 50]
                            yrange = [0, 50]
                        endif
                        if plume eq 1 then begin
                            xrange = [0, 1.5]
                            yrange = [0, 1.5]
                        endif
                        
                    end
                    ;;acet
                    2:begin
                        ovoc_ptr = acet_obs_ptr
                        ovoc_other = acet_obs_toga
                        title = 'Acetone'
                        if plume eq 0 then begin
                            xrange = [0, 12]
                            yrange = [0, 12]
                        endif
                        if plume eq 1 then begin
                            xrange = [0, 3]
                            yrange = [0, 3]
                        endif
                        
                    end
                    
                    ;; Formic acid
                    3:begin
                        ovoc_ptr = hcooh_obs_ptr
                        ovoc_other = hcooh_obs_cims
                        title = 'Formic acid' 
                        if plume eq 0 then begin
                            xrange = [0, 50]
                            yrange = [0, 50]
                        endif
                        if plume eq 1 then begin
                            xrange = [0, 10]
                            yrange = [0, 10]
                        endif
                        
                    end
                    
                    ;;MEK
                    4:begin
                        ovoc_ptr = mek_obs_ptr
                        ovoc_other = mek_obs_toga
                        title = 'MEK'
                        if plume eq 0 then begin
                            xrange = [0, 3.0]
                            yrange = [0, 3.0]
                        endif
                        if plume eq 1 then begin
                            xrange = [0, 0.15]
                            yrange = [0, 0.15]
                        endif
                        
                    end
                    ;;benzene
                    5:begin
                        ovoc_ptr = benz_obs_ptr
                        ovoc_other = benz_obs_toga
                        title = 'Benzene'
                        if plume eq 0 then begin
                            xrange = [0, 5.0]
                            yrange = [0, 5.0]
                        endif
                        if plume eq 1 then begin
                            xrange = [0, 1]
                            yrange = [0, 1]
                        endif
                    end
                    
                    ;;benzene
                    6:begin
                        ovoc_ptr = benz_obs_ptr
                        ovoc_other = benz_obs_toga
                        title = 'Benzene'
                        if plume eq 0 then begin
                            xrange = [0, 8.0]
                            yrange = [0, 8.0]
                        endif
                        if plume eq 1 then begin
                            xrange = [0, 1]
                            yrange = [0, 1]
                        endif
                    end  
                    ;;tolu
                    7:begin
                        ovoc_ptr = tolu_obs_ptr
                        ovoc_other = tolu_obs_toga
                        title = 'Toluene'
                        if plume eq 0 then begin
                            xrange = [0, 5.0]
                            yrange = [0, 5.0]
                        endif
                        if plume eq 1 then begin
                            xrange = [0, 0.3]
                            yrange = [0, 0.3]
                        endif
                    end  
                    ;;xyle
                    8:begin
                        ovoc_ptr = Xylenes_obs_ptr
                        ovoc_other = Xylenes_obs_toga
                        ;ovoc_ptr = c8h10_obs_ptr
                        ;ovoc_other = c8h10_obs_toga

                        
                        title = 'Xylene'
                        if plume eq 0 then begin
                            xrange = [0, 2.0]
                            yrange = [0, 2.0]
                        endif
                        if plume eq 1 then begin
                            xrange = [0, 0.05]
                            yrange = [0, 0.05]
                        endif
                    end
                endcase
            endif
;///////////////////////////////////////////////////////////////////////////////////////////
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;++                                      STEP 3                                           ++
;++                     Setting filters:NA/bb emission/fresh/aged                         ++
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;///////////////////////////////////////////////////////////////////////////////////////////
    ;for voc .vs. co
            oco_tmp=co_obs
    ;for filters
            no2_filter=no2_obs
            o3_filter=o3_obs
            co_filter=co_obs
            ch3cn_filter=ch3cn_obs
            hcn_filter = hcn_obs
            acn_co_filter = ch3cn_obs/co_obs*1000 ; convert it into ppb/ppm
      
    ;cloud
            tmp_rhum = rhum  
;;=================
;;plume filters
;;=================
            if plume eq 1 then begin
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
        ;;remove missing value in either PTR or TOGA
            i = where(ovoc_ptr le 0 or ovoc_other le 0, ct)
            if ct gt 0 then begin
                ovoc_ptr[i] = !VALUES.F_NAN
                ovoc_other[i] = !VALUES.F_NAN
            end
            remove_ind = where(finite(ovoc_ptr) and finite(ovoc_other))

;///////////////////////////////////////////////////////////////////////////////////////////
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;++                                      STEP 4                                           ++
;++                               Removing missing value                                  ++
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;///////////////////////////////////////////////////////////////////////////////////////////

            print, n_elements(ovoc_ptr)
            control_ind = 1
            if not keyword_set(lumped) then control_ind = 1
            if control_ind eq 1 then begin
                com_ind = remove_ind
                if com_ind[0] ne -1 then begin
                    ovoc_ptr = ovoc_ptr[com_ind]
                    ovoc_other = ovoc_other[com_ind] 

                    ;for filters
                    no2_filter=no2_filter[com_ind]
                    o3_filter=o3_filter[com_ind]
                    co_filter = co_filter[com_ind]
                    ch3cn_filter=ch3cn_filter[com_ind]
                    hcn_filter = hcn_filter[com_ind]

                    ;for voc vs. co
                    oco_tmp=oco_tmp[com_ind]
                    
                    ;;tracks
                    tmp_lat_obs = tmp_lat_obs[com_ind]
                    tmp_lon_obs = tmp_lon_obs[com_ind]
                endif          
            endif
    
;///////////////////////////////////////////////////////////////////////////////////////////
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;++                                      STEP 5                                           ++
;++                                  DO THE PLOTTING                                      ++
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;///////////////////////////////////////////////////////////////////////////////////////////
        ;; unit scale factor
            y1=ovoc_ptr
            x1=ovoc_other                         
            print,'this is ',title
            print,'slope, slope_se, intercept, R'
            
;; this is used for PTR vs other emission ratio.
test = 1
if test eq 1 then begin
unitscale =1
    ;;PTR vs Other
            stat1 = org_boot(X=x1,Y=y1,ntrials=1000l) ;M (slope), MSE, MCI_boot, B, BSE, BCI_Boot, R, RCI_Boot, N            
            print,':OBS PTR vs Other',' slope',stat1.M*unitscale, $
                ' SE of slope', stat1.MSE*unitscale, $
                ' CI of slope', stat1.MCI_boot*unitscale, $
                ' intercept', stat1.B, $
                ' SE of intercept', stat1.BSE, $
                ' CI of intercept', stat1.BCI_Boot, $
                ' correlations', stat1.R, $
                ' CI of correlations', stat1.RCI_Boot, $
                ' # of dp', stat1.N
                
;; assign values
    slope_obs = (stat1.MCI_boot[0]+stat1.MCI_boot[1])/2.
    slope_se_obs = (stat1.MCI_boot[1]-stat1.MCI_boot[0])/2.
    intercept_obs = (stat1.BCI_boot[0]+stat1.BCI_boot[1])/2.    
    corr_obs = (stat1.RCI_boot[0]+stat1.RCI_boot[1])/2.


;; Plotting part
            ;if control_ind eq 1 then xrange=[0,10]
            ;ovoc_ptr .VS. oco_picarro
            SCATTERPLOT, x1, y1,_Extra = _Extra, $
                charsize=2.5, charthick = 4, thick = 10,xrange=xrange,yrange=yrange,/log;,xrange=[0,200];,yrange=yrange,xrange=xrange;,$
                ;,XTICKFORMAT="(A1)",YTICKFORMAT="(A1)", 
            oplot,!x.crange,intercept_obs+!x.crange*slope_obs,col=1,line=0, thick=10

;; add math symbol, waiting for test!!
;; http://www.idlcoyote.com/ps_tips/greeksym.php
        thisletter = "262B
        greekLetter = '!9' + String(thisLetter) + '!X'
;; setting 

        if n_species eq 4 then begin
            charsize = 2
            charthick = 6
        endif

        if n_species eq 9 then begin
            charsize = 1.5
            charthick = 5
            second_sf = 0.83
            third_sf = 0.76
        endif
        
        
        if n_species eq 1  then begin
            charsize = 3.5
            charthick = 8
            second_sf = 0.85
            third_sf = 0.80
        endif
        
        thesymbol =   cgSymbol('+-')
        
        format_R = '(F4.2)'
        
        
        ;; Observation
        if slope_obs*unitscale ge 10 then format_slope_obs = '(F5.2)'
        if slope_obs*unitscale lt 10 then format_slope_obs = '(F4.2)'
        
        ;if keyword_set(all) then $
            xyouts,!x.crange(0) + (!x.crange(1)-!x.crange(0))/25,!y.crange(1)*0.90,$
               /data,col=1,'Slope='+string(strtrim(slope_obs*unitscale, 1), format=format_slope_obs) +thesymbol + string(strtrim(slope_se_obs*unitscale, 1), format=format_R) +', R='+$
               string(strtrim(corr_obs, 1),format=format_R),charsize = charsize, charthick = charthick
endif

    endfor ;; s loop

endfor ;; dd loop
close_device
;DEVICE, /CLOSE
print,'done!'
end
