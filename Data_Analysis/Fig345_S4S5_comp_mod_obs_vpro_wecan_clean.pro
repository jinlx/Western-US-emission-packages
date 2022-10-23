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
;; This script is from v16 but adding OH reactivity field for each VOCs, wokring...
;.run GCRateConst
; adding LoD, 0914/2022


; comp_mod_obs_vpro_wecan_clean,plume=0,/all,/nested,/save,/prs,/med,/errorbar,/wus,/whole,/obs_major,/test
; comp_mod_obs_vpro_wecan_clean,plume=0,/each,/ind_co,/nested,/save,/prs,/med,/errorbar,/wus,/whole,/obs_major,/test
; comp_mod_obs_vpro_wecan_clean,plume=1,/filter2_nobb,/filter3_nobb,/all,/nested,/save,/prs,/med,/errorbar,/wus,/whole,/obs_major,/test

@bin_vert
@GCRateConst_co
@GCRateConst_acet
@GCRateConst_c3h8
@GCRateConst
pro comp_mod_obs_vpro_wecan_clean,$
        plume=plume,$
        filter2_nobb=filter2_nobb,filter3_nobb=filter3_nobb,$
        all=all, obs_only=obs_only,each=each,onebyone=onebyone,addup=addup,$
        ind_co=ind_co,ind_benz=ind_benz,ind_o3=ind_o3,ind_pan=ind_pan,ind_no2=ind_no2, $;; it should be used with 'each' keyword_set
        n_criteria=n_criteria,$
        nested=nested,fbf=fbf,save=save, test=test, $
        prs=prs, $
        med=med,average=average,$
        errorbar=errorbar,$
        wus=wus,sus=sus, $ ; meaningless but being consistent with FIREX-AQ
        whole=whole,individual=individual,$
        codeployed=codeployed, obs_major=obs_major
        

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
  ;;v2 means add for those species have both two measurement
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
    ;; Temperature in model
        temperature_mod = [0]
    ;; GFED4
        co_gc_gfed4   = [0]
        O3_gc_gfed4   = [0]
        hcho_gc_gfed4 = [0]
        pan_gc_gfed4  = [0]
        acet_gc_gfed4 = [0]
        benz_gc_gfed4 = [0]
        ch3oh_gc_gfed4= [0]
        ald2_gc_gfed4 = [0]

        no_gc_gfed4   = [0]
        no2_gc_gfed4  = [0]
        so2_gc_gfed4  = [0]
        oh_gc_gfed4   = [0]

        date_gc_gfed4 = [0]
        utc_gc_gfed4  = [0]
        doy_gc_gfed4  = [0]
        lat_gc_gfed4  = [0]
        lon_gc_gfed4  = [0]
        alt_gc_gfed4  = [0]
        prs_gc_gfed4  = [0]

     ;; add vocs
        c3h8_gc_gfed4 = [0]
        c3h6_gc_gfed4 = [0]
        c2h5oh_gc_gfed4 = [0]
        c5h8_gc_gfed4   = [0]
        c7h8_gc_gfed4   = [0]
        dms_gc_gfed4    = [0]
        mek_gc_gfed4    = [0]
        c8h10_gc_gfed4  = [0]

        acta_gc_gfed4   = [0]
        macr_mvk_gc_gfed4 = [0]
        hcooh_gc_gfed4  = [0]

        mtpa_gc_gfed4   = [0]
        limo_gc_gfed4   = [0]
        mtpo_gc_gfed4   = [0]

        c2h4_gc_gfed4   = [0]
        c2h6_gc_gfed4   = [0]

    ;;add lumped species
        alk4_gc_gfed4   = [0]
        prpe_gc_gfed4   = [0]
        rcho_gc_gfed4   = [0]


    ;; FINN
        co_gc_finn   = [0]
        O3_gc_finn   = [0]
        hcho_gc_finn = [0]
        pan_gc_finn  = [0]
        acet_gc_finn = [0]
        benz_gc_finn = [0]
        ch3oh_gc_finn= [0]
        ald2_gc_finn = [0]

        no_gc_finn   = [0]
        no2_gc_finn  = [0]
        so2_gc_finn  = [0]
        oh_gc_finn   = [0]

        date_gc_finn = [0]
        utc_gc_finn  = [0]
        doy_gc_finn  = [0]
        lat_gc_finn  = [0]
        lon_gc_finn  = [0]
        alt_gc_finn  = [0]
        prs_gc_finn  = [0]


     ;; add vocs
        c3h8_gc_finn = [0]
        c3h6_gc_finn = [0]
        c2h5oh_gc_finn = [0]
        c5h8_gc_finn   = [0]
        c7h8_gc_finn   = [0]
        dms_gc_finn    = [0]
        mek_gc_finn    = [0]
        c8h10_gc_finn  = [0]  
        acta_gc_finn   = [0]
        macr_mvk_gc_finn = [0] 
        hcooh_gc_finn  = [0]

        mtpa_gc_finn  = [0]
        limo_gc_finn  = [0]
        mtpo_gc_finn  = [0]

        c2h4_gc_finn   = [0]
        c2h6_gc_finn   = [0]
        
    ;; add lumped species
        alk4_gc_finn   = [0]
        prpe_gc_finn   = [0]
        rcho_gc_finn   = [0]

    ;; GFAS
        co_gc_gfas   = [0]
        O3_gc_gfas   = [0]
        hcho_gc_gfas = [0]
        pan_gc_gfas  = [0]
        acet_gc_gfas = [0]
        benz_gc_gfas = [0]
        ch3oh_gc_gfas= [0]
        ald2_gc_gfas = [0]

        no_gc_gfas   = [0]
        no2_gc_gfas  = [0]
        so2_gc_gfas  = [0]
        oh_gc_gfas   = [0]

        date_gc_gfas = [0]
        utc_gc_gfas  = [0]
        doy_gc_gfas  = [0]
        lat_gc_gfas  = [0]
        lon_gc_gfas  = [0]
        alt_gc_gfas  = [0]
        prs_gc_gfas  = [0]
        temperature_gc_gfas = [0]
     ;; add vocs
        c3h8_gc_gfas = [0]
        c3h6_gc_gfas = [0]
        c2h5oh_gc_gfas = [0]
        c5h8_gc_gfas   = [0]
        c7h8_gc_gfas   = [0]
        dms_gc_gfas    = [0]
        mek_gc_gfas    = [0]
        c8h10_gc_gfas  = [0]

        acta_gc_gfas   = [0]
        macr_mvk_gc_gfas = [0]  
        hcooh_gc_gfas  = [0]

        mtpa_gc_gfas   = [0]
        limo_gc_gfas   = [0]
        mtpo_gc_gfas   = [0]
        
        c2h4_gc_gfas   = [0]
        c2h6_gc_gfas   = [0]
    
    ;;add lumped species
        alk4_gc_gfas   = [0]
        prpe_gc_gfas   = [0]
        rcho_gc_gfas   = [0]

    ;; QFED
        co_gc_qfed   = [0]
        O3_gc_qfed   = [0]
        hcho_gc_qfed = [0]
        pan_gc_qfed  = [0]
        acet_gc_qfed = [0]
        benz_gc_qfed = [0]
        ch3oh_gc_qfed= [0]
        ald2_gc_qfed = [0]

        no_gc_qfed   = [0]
        no2_gc_qfed  = [0]
        so2_gc_qfed  = [0]
        oh_gc_qfed   = [0]

        date_gc_qfed = [0]
        utc_gc_qfed  = [0]
        doy_gc_qfed  = [0]
        lat_gc_qfed  = [0]
        lon_gc_qfed  = [0]
        alt_gc_qfed  = [0]
        prs_gc_qfed  = [0]

     ;; add vocs
        c3h8_gc_qfed = [0]
        c3h6_gc_qfed = [0]
        c2h5oh_gc_qfed = [0]
        c5h8_gc_qfed   = [0]
        c7h8_gc_qfed   = [0]
        dms_gc_qfed    = [0]
        mek_gc_qfed    = [0]
        c8h10_gc_qfed  = [0]  
        acta_gc_qfed   = [0]
        macr_mvk_gc_qfed = [0] 
        hcooh_gc_qfed  = [0]

        mtpa_gc_qfed  = [0]
        limo_gc_qfed  = [0]
        mtpo_gc_qfed   = [0]

        c2h4_gc_qfed   = [0]
        c2h6_gc_qfed   = [0]
        
    ;;add lumped species
        alk4_gc_qfed   = [0]
        prpe_gc_qfed   = [0]
        rcho_gc_qfed   = [0]
        
    ;; THREEGFAS
        co_gc_threegfas   = [0]
        O3_gc_threegfas   = [0]
        hcho_gc_threegfas = [0]
        pan_gc_threegfas  = [0]
        acet_gc_threegfas = [0]
        benz_gc_threegfas = [0]
        ch3oh_gc_threegfas= [0]
        ald2_gc_threegfas = [0]

        no_gc_threegfas   = [0]
        no2_gc_threegfas  = [0]
        so2_gc_threegfas  = [0]
        oh_gc_threegfas   = [0]

        date_gc_threegfas = [0]
        utc_gc_threegfas  = [0]
        doy_gc_threegfas  = [0]
        lat_gc_threegfas  = [0]
        lon_gc_threegfas  = [0]
        alt_gc_threegfas  = [0]
        prs_gc_threegfas  = [0]

     ;; add vocs
        c3h8_gc_threegfas = [0]
        c3h6_gc_threegfas = [0]
        c2h5oh_gc_threegfas = [0]
        c5h8_gc_threegfas   = [0]
        c7h8_gc_threegfas   = [0]
        dms_gc_threegfas    = [0]
        mek_gc_threegfas    = [0]
        c8h10_gc_threegfas  = [0]  
        acta_gc_threegfas   = [0]
        macr_mvk_gc_threegfas = [0] 
        hcooh_gc_threegfas  = [0]

        mtpa_gc_threegfas  = [0]
        limo_gc_threegfas  = [0]
        mtpo_gc_threegfas   = [0]

        c2h4_gc_threegfas   = [0]
        c2h6_gc_threegfas   = [0]
        
    ;;add lumped species
        alk4_gc_threegfas   = [0]
        prpe_gc_threegfas   = [0]
        rcho_gc_threegfas   = [0]


    ;; NOBB
        co_gc_nobb   = [0]
        O3_gc_nobb   = [0]
        hcho_gc_nobb = [0]
        pan_gc_nobb  = [0]
        acet_gc_nobb = [0]
        benz_gc_nobb = [0]
        ch3oh_gc_nobb= [0]
        ald2_gc_nobb = [0]

        no_gc_nobb   = [0]
        no2_gc_nobb  = [0]
        so2_gc_nobb  = [0]
        oh_gc_nobb   = [0]

        date_gc_nobb = [0]
        utc_gc_nobb  = [0]
        doy_gc_nobb  = [0]
        lat_gc_nobb  = [0]
        lon_gc_nobb  = [0]
        alt_gc_nobb  = [0]
        prs_gc_nobb  = [0]


     ;; add vocs
        c3h8_gc_nobb = [0]
        c3h6_gc_nobb = [0]
        c2h5oh_gc_nobb = [0]
        c5h8_gc_nobb   = [0]
        c7h8_gc_nobb   = [0]
        dms_gc_nobb    = [0]
        mek_gc_nobb    = [0]
        c8h10_gc_nobb  = [0]  
        acta_gc_nobb   = [0]
        macr_mvk_gc_nobb = [0] 
        hcooh_gc_nobb  = [0]

        mtpa_gc_nobb  = [0]
        limo_gc_nobb  = [0]
        mtpo_gc_nobb  = [0]

        c2h4_gc_nobb   = [0]
        c2h6_gc_nobb   = [0]
        
    ;; add lumped species
        alk4_gc_nobb   = [0]
        prpe_gc_nobb   = [0]
        rcho_gc_nobb   = [0]

    ;; empty column for time index
        ind_enter_time = [0]
        ind_exit_time = [0]

    ; smart loop for each fligiht and total flights from lu
        if dd eq 0 then dates = dates_all else dates=dates_all[dd-1]
        dates = STRTRIM(string(dates),1)
        for n=0,n_elements(dates)-1 do begin
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
;; GEOS-Chem 0.25x0.3125: GFED4
;; ============================
            if keyword_set(nested) then begin
                ;gcfi_gfed4   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                ;    'planelog2sav/output_'+'gfed4' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'
                ;gcfi_gfed4   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                ;    'planelog2sav/output_'+'gfed4_tmp' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'
                gcfi_gfed4   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                    'planelog2sav/output_'+'gfed4_tmp' + '/mrg5m_wecan_c130_'+dates[n]+'.sav'
            endif
            if keyword_set(fbf) then begin
                ;gcfi_gfed4   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                ;    'planelog2sav/output_'+'gfed4_4x5' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'
                gcfi_gfed4   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                    'planelog2sav/'+'output_WECAN_C2raw' + '/mrg30m_wecan_c130_'+dates[n]+'.sav'
            endif
            restore, gcfi_gfed4
            
            tmp_utc_gc_gfed4 = gc.utc 
            tmp_co_gc_gfed4   = gc.co*1e9
            tmp_O3_gc_gfed4   = gc.o3*1e9
            tmp_pan_gc_gfed4  = gc.pan*1e9
            tmp_hcho_gc_gfed4 = gc.ch2o*1e9
            tmp_acet_gc_gfed4 = gc.acet*1e9/3; 3 acetone in ppb
            tmp_benz_gc_gfed4 = gc.benz*1e9/6; 6C benzene in ppb
    ;        ch3oh_gc= [ch3oh_gc, gc.MOH*1e9]
            tmp_ald2_gc_gfed4 = gc.ald2*1e9/2; 2C ch3cho in ppb
            
            tmp_no_gc_gfed4   = gc.no*1e9
            tmp_no2_gc_gfed4  = gc.no2*1e9
            tmp_so2_gc_gfed4  = gc.so2*1e9

            tmp_na = avo * gc.pres * 100./(8.31 * gc.temp) * 1e-6;; air density molec/cm3
            tmp_oh_gc_gfed4  = gc.oh*tmp_na ;; v/v --> molec/cm3 

            tmp_date_gc_gfed4 = gc.date
            ;tmp_utc_gc_gfed4  = gc.utc
            tmp_doy_gc_gfed4  = gc.doy
            tmp_lat_gc_gfed4  = gc.lat
            tmp_lon_gc_gfed4  = gc.lon
            tmp_alt_gc_gfed4  = gc.alt
            tmp_prs_gc_gfed4  = gc.pres

     ;; add vocs
            tmp_c3h8_gc_gfed4   = gc.c3h8*1e9/3; 3C c3h8 in ppb
            tmp_c3h6_gc_gfed4   = gc.prpe*1e9/3; 3C PRPE in ppb
            tmp_c2h5oh_gc_gfed4 = gc.eoh*1e9/2; 2C eoh in ppb
            tmp_c5h8_gc_gfed4   = gc.isop*1e9/5; 5C isoprene in ppb
            tmp_c7h8_gc_gfed4   = gc.tolu*1e9/7; 7C toluene in ppb
            tmp_dms_gc_gfed4    = gc.dms*1e9
            tmp_mek_gc_gfed4    = gc.mek*1e9/4 ;4C  MEK in ppb
            tmp_c8h10_gc_gfed4  = gc.xyle*1e9/8 ; 8C xylenes in ppb

            tmp_acta_gc_gfed4   = gc.acta*1e9
            tmp_macr_mvk_gc_gfed4 = (gc.MVK+gc.MACR)*1e9  
            tmp_hcooh_gc_gfed4  = gc.hcooh*1e9

            tmp_mtpa_gc_gfed4  = gc.mtpa*1e9
            tmp_limo_gc_gfed4  = gc.limo*1e9
            tmp_mtpo_gc_gfed4  = gc.mtpo*1e9
        ;;add lumped species
            tmp_alk4_gc_gfed4   = gc.alk4*1e9/4.3 ; 4.3C
            tmp_prpe_gc_gfed4   = gc.prpe*1e9/3   ; 3C 
            tmp_rcho_gc_gfed4   = gc.rcho*1e9 

            ;tmp_c2h4_gc_gfed4   = gc.c2h4*1e9
            tmp_c2h6_gc_gfed4   = gc.c2h6*1e9/2 ;; 2C

            co_gc_gfed4 = [co_gc_gfed4,tmp_co_gc_gfed4]
            O3_gc_gfed4   = [o3_gc_gfed4,tmp_o3_gc_gfed4]
            pan_gc_gfed4  = [pan_gc_gfed4,tmp_pan_gc_gfed4]
            hcho_gc_gfed4 = [hcho_gc_gfed4,tmp_hcho_gc_gfed4]
            acet_gc_gfed4 = [acet_gc_gfed4,tmp_acet_gc_gfed4]
            benz_gc_gfed4 = [benz_gc_gfed4,tmp_benz_gc_gfed4]

            ald2_gc_gfed4 = [ald2_gc_gfed4,tmp_ald2_gc_gfed4]


            no_gc_gfed4   = [no_gc_gfed4,tmp_no_gc_gfed4]
            no2_gc_gfed4  = [no2_gc_gfed4,tmp_no2_gc_gfed4]
            so2_gc_gfed4  = [so2_gc_gfed4,tmp_so2_gc_gfed4]

            oh_gc_gfed4  = [oh_gc_gfed4,tmp_oh_gc_gfed4] ;; v/v --> molec/cm3 

            date_gc_gfed4 = [date_gc_gfed4,tmp_date_gc_gfed4]
            ;utc_gc_gfed4  = [utc_gc_gfed4,tmp_utc_gc_gfed4]
            doy_gc_gfed4  = [doy_gc_gfed4,tmp_doy_gc_gfed4]
            lat_gc_gfed4  = [lat_gc_gfed4,tmp_lat_gc_gfed4]
            lon_gc_gfed4  = [lon_gc_gfed4,tmp_lon_gc_gfed4]
            alt_gc_gfed4  = [alt_gc_gfed4,tmp_alt_gc_gfed4]
            prs_gc_gfed4  = [prs_gc_gfed4,tmp_prs_gc_gfed4]
     ;; add vocs
            c3h8_gc_gfed4   = [c3h8_gc_gfed4,tmp_c3h8_gc_gfed4]
            c3h6_gc_gfed4   = [c3h6_gc_gfed4,tmp_c3h6_gc_gfed4]
            c2h5oh_gc_gfed4 = [c2h5oh_gc_gfed4,tmp_c2h5oh_gc_gfed4]
            c5h8_gc_gfed4   = [c5h8_gc_gfed4,tmp_c5h8_gc_gfed4]
            c7h8_gc_gfed4   = [c7h8_gc_gfed4,tmp_c7h8_gc_gfed4]
            dms_gc_gfed4    = [dms_gc_gfed4,tmp_dms_gc_gfed4]
            mek_gc_gfed4    = [mek_gc_gfed4,tmp_mek_gc_gfed4]
            c8h10_gc_gfed4  = [c8h10_gc_gfed4,tmp_c8h10_gc_gfed4]

            acta_gc_gfed4   = [acta_gc_gfed4,tmp_acta_gc_gfed4]
            macr_mvk_gc_gfed4 = [macr_mvk_gc_gfed4,tmp_macr_mvk_gc_gfed4]  
            hcooh_gc_gfed4  = [hcooh_gc_gfed4,tmp_hcooh_gc_gfed4]

            mtpa_gc_gfed4  = [mtpa_gc_gfed4,tmp_mtpa_gc_gfed4]
            limo_gc_gfed4  = [limo_gc_gfed4,tmp_limo_gc_gfed4]
            mtpo_gc_gfed4  = [mtpo_gc_gfed4,tmp_mtpo_gc_gfed4]
        ;;add lumped species
            alk4_gc_gfed4  = [alk4_gc_gfed4,tmp_alk4_gc_gfed4]
            prpe_gc_gfed4  = [prpe_gc_gfed4,tmp_prpe_gc_gfed4]
            rcho_gc_gfed4  = [rcho_gc_gfed4,tmp_rcho_gc_gfed4]

            ;c2h4_gc_gfed4   = [c2h4_gc_gfed4, tmp_c2h4_gc_gfed4]
            c2h6_gc_gfed4   = [c2h6_gc_gfed4, tmp_c2h6_gc_gfed4]
            undefine,gc
            
;; ============================
;; GEOS-Chem 0.25x0.3125: FINN
;; ============================
            if keyword_set(nested) then begin
                ;gcfi_finn   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                ;    'planelog2sav/output_'+'finn' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'
                ;gcfi_finn   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                ;    'planelog2sav/output_'+'finn_tmp' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'
                gcfi_finn   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                    'planelog2sav/output_'+'gfed4_tmp' + '/mrg5m_wecan_c130_'+dates[n]+'.sav'
            endif
            if keyword_set(fbf) then begin
                gcfi_finn   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                    'planelog2sav/'+'output_WECAN_C2raw' + '/mrg30m_wecan_c130_'+dates[n]+'.sav'
            endif
            
            restore, gcfi_finn
            tmp_utc_gc_finn = gc.utc 
            tmp_co_gc_finn   = gc.co*1e9
            tmp_O3_gc_finn   = gc.o3*1e9
            tmp_pan_gc_finn  = gc.pan*1e9
            tmp_hcho_gc_finn = gc.ch2o*1e9
            tmp_acet_gc_finn = gc.acet*1e9/3; 3 acetone in ppb
            tmp_benz_gc_finn = gc.benz*1e9/6; 6C benzene in ppb
    ;        ch3oh_gc= [ch3oh_gc, gc.MOH*1e9]
            tmp_ald2_gc_finn = gc.ald2*1e9/2; 2C ch3cho in ppb

            tmp_no_gc_finn   = gc.no*1e9
            tmp_no2_gc_finn  = gc.no2*1e9
            tmp_so2_gc_finn  = gc.so2*1e9

            tmp_na = avo * gc.pres * 100./(8.31 * gc.temp) * 1e-6;; air density molec/cm3
            tmp_oh_gc_finn  = gc.oh*tmp_na ;; v/v --> molec/cm3 

            tmp_date_gc_finn = gc.date
            ;tmp_utc_gc_finn  = gc.utc
            tmp_doy_gc_finn  = gc.doy
            tmp_lat_gc_finn  = gc.lat
            tmp_lon_gc_finn  = gc.lon
            tmp_alt_gc_finn  = gc.alt
            tmp_prs_gc_finn  = gc.pres
     ;; add vocs
            tmp_c3h8_gc_finn   = gc.c3h8*1e9/3; 3C c3h8 in ppb
            tmp_c3h6_gc_finn   = gc.prpe*1e9/3; 3C PRPE in ppb
            tmp_c2h5oh_gc_finn = gc.eoh*1e9/2; 2C eoh in ppb
            tmp_c5h8_gc_finn   = gc.isop*1e9/5; 5C isoprene in ppb
            tmp_c7h8_gc_finn   = gc.tolu*1e9/7; 7C toluene in ppb
            tmp_dms_gc_finn    = gc.dms*1e9
            tmp_mek_gc_finn    = gc.mek*1e9/4 ;4C  MEK in ppb
            tmp_c8h10_gc_finn  = gc.xyle*1e9/8 ; 8C xylenes in ppb
            
            tmp_acta_gc_finn   = gc.acta*1e9
            tmp_macr_mvk_gc_finn = (gc.MVK+gc.MACR)*1e9  
            tmp_hcooh_gc_finn  = gc.hcooh*1e9

            tmp_mtpa_gc_finn  = gc.mtpa*1e9
            tmp_limo_gc_finn  = gc.limo*1e9
            tmp_mtpo_gc_finn  = gc.mtpo*1e9
        ;;add lumped species
            tmp_alk4_gc_finn   = gc.alk4*1e9/4.3 ; 4.3C
            tmp_prpe_gc_finn   = gc.prpe*1e9/3   ; 3C 
            tmp_rcho_gc_finn   = gc.rcho*1e9

            ;tmp_c2h4_gc_finn   = gc.c2h4*1e9
            tmp_c2h6_gc_finn   = gc.c2h6*1e9/2 ;;2C

;; prepare data for 20180826, lat is less than zero and it is not simulated by the obs,fixed
            i = where(filter_lat_obs le 0,ct) 
            ;; test, since we did something wrong in FINN
            if ct gt 0 and keyword_set(nested) then begin
                tmp_lat_gc_finn[i] = !VALUES.F_NAN
                remove_ind=where(finite(tmp_lat_gc_finn),ct)
                if remove_ind[0] ne -1 then begin
                    tmp_co_gc_finn   = tmp_co_gc_finn[remove_ind]
                    tmp_O3_gc_finn   = tmp_O3_gc_finn[remove_ind]
                    tmp_pan_gc_finn  = tmp_pan_gc_finn[remove_ind]
                    tmp_hcho_gc_finn = tmp_hcho_gc_finn[remove_ind]
                    tmp_acet_gc_finn = tmp_acet_gc_finn[remove_ind]
                    tmp_benz_gc_finn = tmp_benz_gc_finn[remove_ind]
            ;        ch3oh_gc= [ch3oh_gc, gc.MOH*1e9]
                    tmp_ald2_gc_finn = tmp_ald2_gc_finn[remove_ind]

                    tmp_no_gc_finn   = tmp_no_gc_finn[remove_ind]
                    tmp_no2_gc_finn  = tmp_no2_gc_finn[remove_ind]
                    tmp_so2_gc_finn  = tmp_so2_gc_finn[remove_ind]

                    ;tmp_na = tmp_na[remove_ind]
                    tmp_oh_gc_finn  = tmp_oh_gc_finn[remove_ind]

                    tmp_date_gc_finn = tmp_date_gc_finn[remove_ind]
                    tmp_utc_gc_finn  = tmp_utc_gc_finn[remove_ind]
                    tmp_doy_gc_finn  = tmp_doy_gc_finn[remove_ind]
                    tmp_lat_gc_finn  = tmp_lat_gc_finn[remove_ind]
                    tmp_lon_gc_finn  = tmp_lon_gc_finn[remove_ind]
                    tmp_alt_gc_finn  = tmp_alt_gc_finn[remove_ind]
                    tmp_prs_gc_finn  = tmp_prs_gc_finn[remove_ind]
             ;; add vocs
                    tmp_c3h8_gc_finn   = tmp_c3h8_gc_finn[remove_ind]
                    tmp_c3h6_gc_finn   = tmp_c3h6_gc_finn[remove_ind]
                    tmp_c2h5oh_gc_finn = tmp_c2h5oh_gc_finn[remove_ind]
                    tmp_c5h8_gc_finn   = tmp_c5h8_gc_finn[remove_ind]
                    tmp_c7h8_gc_finn   = tmp_c7h8_gc_finn[remove_ind]
                    tmp_dms_gc_finn    = tmp_dms_gc_finn[remove_ind]
                    tmp_mek_gc_finn    = tmp_mek_gc_finn[remove_ind]
                    tmp_c8h10_gc_finn  = tmp_c8h10_gc_finn[remove_ind]

                    tmp_acta_gc_finn   = tmp_acta_gc_finn[remove_ind]
                    tmp_macr_mvk_gc_finn = tmp_macr_mvk_gc_finn[remove_ind]
                    tmp_hcooh_gc_finn  = tmp_hcooh_gc_finn[remove_ind]

                    tmp_mtpa_gc_finn  = tmp_mtpa_gc_finn[remove_ind]
                    tmp_limo_gc_finn  = tmp_limo_gc_finn[remove_ind]
                    tmp_mtpo_gc_finn  = tmp_mtpo_gc_finn[remove_ind]

                ;;add lumped species
                    tmp_alk4_gc_finn   = tmp_alk4_gc_finn[remove_ind]
                    tmp_prpe_gc_finn   = tmp_prpe_gc_finn[remove_ind]
                    tmp_rcho_gc_finn   = tmp_rcho_gc_finn[remove_ind]

                    ;tmp_c2h4_gc_finn   = tmp_c2h4_gc_finn[remove_ind]
                    tmp_c2h6_gc_finn   = tmp_c2h6_gc_finn[remove_ind]
                endif
            endif
            
            co_gc_finn = [co_gc_finn,tmp_co_gc_finn]
            O3_gc_finn   = [o3_gc_finn,tmp_o3_gc_finn]
            pan_gc_finn  = [pan_gc_finn,tmp_pan_gc_finn]
            hcho_gc_finn = [hcho_gc_finn,tmp_hcho_gc_finn]
            acet_gc_finn = [acet_gc_finn,tmp_acet_gc_finn]
            benz_gc_finn = [benz_gc_finn,tmp_benz_gc_finn]

            ald2_gc_finn = [ald2_gc_finn,tmp_ald2_gc_finn]


            no_gc_finn   = [no_gc_finn,tmp_no_gc_finn]
            no2_gc_finn  = [no2_gc_finn,tmp_no2_gc_finn]
            so2_gc_finn  = [so2_gc_finn,tmp_so2_gc_finn]

            oh_gc_finn  = [oh_gc_finn,tmp_oh_gc_finn] ;; v/v --> molec/cm3 

            date_gc_finn = [date_gc_finn,tmp_date_gc_finn]
            ;utc_gc_finn  = [utc_gc_finn,tmp_utc_gc_finn]
            doy_gc_finn  = [doy_gc_finn,tmp_doy_gc_finn]
            lat_gc_finn  = [lat_gc_finn,tmp_lat_gc_finn]
            lon_gc_finn  = [lon_gc_finn,tmp_lon_gc_finn]
            alt_gc_finn  = [alt_gc_finn,tmp_alt_gc_finn]
            prs_gc_finn  = [prs_gc_finn,tmp_prs_gc_finn]
     ;; add vocs
            c3h8_gc_finn   = [c3h8_gc_finn,tmp_c3h8_gc_finn]
            c3h6_gc_finn   = [c3h6_gc_finn,tmp_c3h6_gc_finn]
            c2h5oh_gc_finn = [c2h5oh_gc_finn,tmp_c2h5oh_gc_finn]
            c5h8_gc_finn   = [c5h8_gc_finn,tmp_c5h8_gc_finn]
            c7h8_gc_finn   = [c7h8_gc_finn,tmp_c7h8_gc_finn]
            dms_gc_finn    = [dms_gc_finn,tmp_dms_gc_finn]
            mek_gc_finn    = [mek_gc_finn,tmp_mek_gc_finn]
            c8h10_gc_finn  = [c8h10_gc_finn,tmp_c8h10_gc_finn]

            acta_gc_finn   = [acta_gc_finn,tmp_acta_gc_finn]
            macr_mvk_gc_finn = [macr_mvk_gc_finn,tmp_macr_mvk_gc_finn]  
            hcooh_gc_finn  = [hcooh_gc_finn,tmp_hcooh_gc_finn]

            mtpa_gc_finn  = [mtpa_gc_finn,tmp_mtpa_gc_finn]
            limo_gc_finn  = [limo_gc_finn,tmp_limo_gc_finn]
            mtpo_gc_finn  = [mtpo_gc_finn,tmp_mtpo_gc_finn]
        ;;add lumped species
            alk4_gc_finn  = [alk4_gc_finn,tmp_alk4_gc_finn]
            prpe_gc_finn  = [prpe_gc_finn,tmp_prpe_gc_finn]
            rcho_gc_finn  = [rcho_gc_finn,tmp_rcho_gc_finn]

            ;c2h4_gc_finn = [c2h4_gc_finn,tmp_c2h4_gc_finn]
            c2h6_gc_finn = [c2h6_gc_finn,tmp_c2h6_gc_finn]
    
            undefine,gc
;; ============================
;; GEOS-Chem 0.25x0.3125: GFAS
;; ============================
            if keyword_set(nested) then begin
                ;gcfi_gfas   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                ;    'planelog2sav/output_'+'gfas' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'
                ;gcfi_gfas   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                ;    'planelog2sav/output_'+'gfas_tmp' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'         
                gcfi_gfas   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                    'planelog2sav/output_'+'gfas_tmp' + '/mrg5m_wecan_c130_'+dates[n]+'.sav'   
            endif
            if keyword_set(fbf) then begin
                ;gcfi_gfas   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                ;    'planelog2sav/output_'+'gfas_4x5' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'
                gcfi_gfas   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                    'planelog2sav/'+'output_WECAN_C2test' + '/mrg30m_wecan_c130_'+dates[n]+'.sav'
            endif
            restore, gcfi_gfas
            tmp_utc_gc_gfas = gc.utc             
            tmp_co_gc_gfas   = gc.co*1e9
            tmp_O3_gc_gfas   = gc.o3*1e9
            tmp_pan_gc_gfas  = gc.pan*1e9
            tmp_hcho_gc_gfas = gc.ch2o*1e9
            tmp_acet_gc_gfas = gc.acet*1e9/3; 3C acetone in ppb
            tmp_benz_gc_gfas = gc.benz*1e9/6; 6C benzene in ppb
    ;        ch3oh_gc= [ch3oh_gc, gc.MOH*1e9]
            tmp_ald2_gc_gfas = gc.ald2*1e9/2; 2C ch3cho in ppb

            tmp_no_gc_gfas   = gc.no*1e9
            tmp_no2_gc_gfas  = gc.no2*1e9
            tmp_so2_gc_gfas  = gc.so2*1e9

            tmp_na = avo * gc.pres * 100./(8.31 * gc.temp) * 1e-6;; air density molec/cm3
            tmp_oh_gc_gfas  = gc.oh*tmp_na ;; v/v --> molec/cm3 

            tmp_date_gc_gfas = gc.date
            ;tmp_utc_gc_gfas  = gc.utc
            tmp_doy_gc_gfas  = gc.doy
            tmp_lat_gc_gfas  = gc.lat
            tmp_lon_gc_gfas  = gc.lon
            tmp_alt_gc_gfas  = gc.alt
            tmp_prs_gc_gfas  = gc.pres
            ;;temperature
            tmp_temperature_gc_gfas = gc.temp
     ;; add vocs
            tmp_c3h8_gc_gfas   = gc.c3h8*1e9/3; 3C c3h8 in ppb
            tmp_c3h6_gc_gfas   = gc.prpe*1e9/3; 3C PRPE in ppb
            tmp_c2h5oh_gc_gfas = gc.eoh*1e9/2; 2C eoh in ppb
            tmp_c5h8_gc_gfas   = gc.isop*1e9/5; 5C isoprene in ppb
            tmp_c7h8_gc_gfas   = gc.tolu*1e9/7; 7C toluene in ppb
            tmp_dms_gc_gfas    = gc.dms*1e9
            tmp_mek_gc_gfas    = gc.mek*1e9/4 ;4C  MEK in ppb
            tmp_c8h10_gc_gfas  = gc.xyle*1e9/8 ; 8C xylenes in ppb

            tmp_acta_gc_gfas   = gc.acta*1e9
            tmp_macr_mvk_gc_gfas = (gc.MVK+gc.MACR)*1e9  
            tmp_hcooh_gc_gfas  = gc.hcooh*1e9

            tmp_mtpa_gc_gfas  = gc.mtpa*1e9
            tmp_limo_gc_gfas  = gc.limo*1e9
            tmp_mtpo_gc_gfas  = gc.mtpo*1e9
        ;;add lumped species
            tmp_alk4_gc_gfas   = gc.alk4*1e9/4.3 ; 4.3C
            tmp_prpe_gc_gfas   = gc.prpe*1e9/3   ; 3C 
            tmp_rcho_gc_gfas   = gc.rcho*1e9

            ;tmp_c2h4_gc_gfas   = gc.c2h4*1e9
            tmp_c2h6_gc_gfas   = gc.c2h6*1e9/2



            co_gc_gfas = [co_gc_gfas,tmp_co_gc_gfas]
            O3_gc_gfas   = [o3_gc_gfas,tmp_o3_gc_gfas]
            pan_gc_gfas  = [pan_gc_gfas,tmp_pan_gc_gfas]
            hcho_gc_gfas = [hcho_gc_gfas,tmp_hcho_gc_gfas]
            acet_gc_gfas = [acet_gc_gfas,tmp_acet_gc_gfas]
            benz_gc_gfas = [benz_gc_gfas,tmp_benz_gc_gfas]

            ald2_gc_gfas = [ald2_gc_gfas,tmp_ald2_gc_gfas]


            no_gc_gfas   = [no_gc_gfas,tmp_no_gc_gfas]
            no2_gc_gfas  = [no2_gc_gfas,tmp_no2_gc_gfas]
            so2_gc_gfas  = [so2_gc_gfas,tmp_so2_gc_gfas]

            oh_gc_gfas  = [oh_gc_gfas,tmp_oh_gc_gfas] ;; v/v --> molec/cm3 

            date_gc_gfas = [date_gc_gfas,tmp_date_gc_gfas]
            ;utc_gc_gfas  = [utc_gc_gfas,tmp_utc_gc_gfas]
            doy_gc_gfas  = [doy_gc_gfas,tmp_doy_gc_gfas]
            lat_gc_gfas  = [lat_gc_gfas,tmp_lat_gc_gfas]
            lon_gc_gfas  = [lon_gc_gfas,tmp_lon_gc_gfas]
            alt_gc_gfas  = [alt_gc_gfas,tmp_alt_gc_gfas]
            prs_gc_gfas  = [prs_gc_gfas,tmp_prs_gc_gfas]
            ;;temperature        
            temperature_gc_gfas = [temperature_gc_gfas,tmp_temperature_gc_gfas]
     ;; add vocs
            c3h8_gc_gfas   = [c3h8_gc_gfas,tmp_c3h8_gc_gfas]
            c3h6_gc_gfas   = [c3h6_gc_gfas,tmp_c3h6_gc_gfas]
            c2h5oh_gc_gfas = [c2h5oh_gc_gfas,tmp_c2h5oh_gc_gfas]
            c5h8_gc_gfas   = [c5h8_gc_gfas,tmp_c5h8_gc_gfas]
            c7h8_gc_gfas   = [c7h8_gc_gfas,tmp_c7h8_gc_gfas]
            dms_gc_gfas    = [dms_gc_gfas,tmp_dms_gc_gfas]
            mek_gc_gfas    = [mek_gc_gfas,tmp_mek_gc_gfas]
            c8h10_gc_gfas  = [c8h10_gc_gfas,tmp_c8h10_gc_gfas]

            acta_gc_gfas   = [acta_gc_gfas,tmp_acta_gc_gfas]
            macr_mvk_gc_gfas = [macr_mvk_gc_gfas,tmp_macr_mvk_gc_gfas]  
            hcooh_gc_gfas  = [hcooh_gc_gfas,tmp_hcooh_gc_gfas]

            mtpa_gc_gfas  = [mtpa_gc_gfas,tmp_mtpa_gc_gfas]
            limo_gc_gfas  = [limo_gc_gfas,tmp_limo_gc_gfas]
            mtpo_gc_gfas  = [mtpo_gc_gfas,tmp_mtpo_gc_gfas]
        ;;add lumped species
            alk4_gc_gfas  = [alk4_gc_gfas,tmp_alk4_gc_gfas]
            prpe_gc_gfas  = [prpe_gc_gfas,tmp_prpe_gc_gfas]
            rcho_gc_gfas  = [rcho_gc_gfas,tmp_rcho_gc_gfas]

            ;c2h4_gc_gfas   = [c2h4_gc_gfas, tmp_c2h4_gc_gfas]
            c2h6_gc_gfas   = [c2h6_gc_gfas, tmp_c2h6_gc_gfas]
    
            undefine,gc

;; ============================
;; GEOS-Chem 0.25x0.3125: QFED
;; ============================
            if keyword_set(nested) then begin
                ;gcfi_qfed   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                ;    'planelog2sav/output_'+'qfed' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'
                ;gcfi_qfed   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                ;    'planelog2sav/output_'+'qfed_tmp' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'
                gcfi_qfed   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                    'planelog2sav/output_'+'qfed_tmp' + '/mrg5m_wecan_c130_'+dates[n]+'.sav'
            endif
            if keyword_set(fbf) then begin
                gcfi_qfed   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                    'planelog2sav/'+'output_WECAN_C2test' + '/mrg30m_wecan_c130_'+dates[n]+'.sav'
            endif
            restore, gcfi_qfed       
            tmp_utc_gc_qfed = gc.utc 
            tmp_co_gc_qfed   = gc.co*1e9
            tmp_O3_gc_qfed   = gc.o3*1e9
            tmp_pan_gc_qfed  = gc.pan*1e9
            tmp_hcho_gc_qfed = gc.ch2o*1e9
            tmp_acet_gc_qfed = gc.acet*1e9/3; 3 acetone in ppb
            tmp_benz_gc_qfed = gc.benz*1e9/6; 6C benzene in ppb
    ;        ch3oh_gc= [ch3oh_gc, gc.MOH*1e9]
            tmp_ald2_gc_qfed = gc.ald2*1e9/2; 2C ch3cho in ppb

            tmp_no_gc_qfed   = gc.no*1e9
            tmp_no2_gc_qfed  = gc.no2*1e9
            tmp_so2_gc_qfed  = gc.so2*1e9

            tmp_na = avo * gc.pres * 100./(8.31 * gc.temp) * 1e-6;; air density molec/cm3
            tmp_oh_gc_qfed  = gc.oh*tmp_na ;; v/v --> molec/cm3 

            tmp_date_gc_qfed = gc.date
            ;tmp_utc_gc_qfed  = gc.utc
            tmp_doy_gc_qfed  = gc.doy
            tmp_lat_gc_qfed  = gc.lat
            tmp_lon_gc_qfed  = gc.lon
            tmp_alt_gc_qfed  = gc.alt
            tmp_prs_gc_qfed  = gc.pres
     ;; add vocs
            tmp_c3h8_gc_qfed   = gc.c3h8*1e9/3; 3C c3h8 in ppb
            tmp_c3h6_gc_qfed   = gc.prpe*1e9/3; 3C PRPE in ppb
            tmp_c2h5oh_gc_qfed = gc.eoh*1e9/2; 2C eoh in ppb
            tmp_c5h8_gc_qfed   = gc.isop*1e9/5; 5C isoprene in ppb
            tmp_c7h8_gc_qfed   = gc.tolu*1e9/7; 7C toluene in ppb
            tmp_dms_gc_qfed    = gc.dms*1e9
            tmp_mek_gc_qfed    = gc.mek*1e9/4 ;4C  MEK in ppb
            tmp_c8h10_gc_qfed  = gc.xyle*1e9/8 ; 8C xylenes in ppb


            tmp_acta_gc_qfed   = gc.acta*1e9
            tmp_macr_mvk_gc_qfed = (gc.MVK+gc.MACR)*1e9  
            tmp_hcooh_gc_qfed  = gc.hcooh*1e9

            tmp_mtpa_gc_qfed  = gc.mtpa*1e9
            tmp_limo_gc_qfed  = gc.limo*1e9
            tmp_mtpo_gc_qfed  = gc.mtpo*1e9
        ;;add lumped species
            tmp_alk4_gc_qfed  = gc.alk4*1e9/4.3 ; 4.3C
            tmp_prpe_gc_qfed   = gc.prpe*1e9/3   ; 3C 
            tmp_rcho_gc_qfed   = gc.rcho*1e9

            ;tmp_c2h4_gc_qfed   = gc.c2h4*1e9
            tmp_c2h6_gc_qfed   = gc.c2h6*1e9/2 ;; 2C


            co_gc_qfed = [co_gc_qfed,tmp_co_gc_qfed]
            O3_gc_qfed   = [o3_gc_qfed,tmp_o3_gc_qfed]
            pan_gc_qfed  = [pan_gc_qfed,tmp_pan_gc_qfed]
            hcho_gc_qfed = [hcho_gc_qfed,tmp_hcho_gc_qfed]
            acet_gc_qfed = [acet_gc_qfed,tmp_acet_gc_qfed]
            benz_gc_qfed = [benz_gc_qfed,tmp_benz_gc_qfed]

            ald2_gc_qfed = [ald2_gc_qfed,tmp_ald2_gc_qfed]


            no_gc_qfed   = [no_gc_qfed,tmp_no_gc_qfed]
            no2_gc_qfed  = [no2_gc_qfed,tmp_no2_gc_qfed]
            so2_gc_qfed  = [so2_gc_qfed,tmp_so2_gc_qfed]

            oh_gc_qfed  = [oh_gc_qfed,tmp_oh_gc_qfed] ;; v/v --> molec/cm3 

            date_gc_qfed = [date_gc_qfed,tmp_date_gc_qfed]
            ;utc_gc_qfed  = [utc_gc_qfed,tmp_utc_gc_qfed]
            doy_gc_qfed  = [doy_gc_qfed,tmp_doy_gc_qfed]
            lat_gc_qfed  = [lat_gc_qfed,tmp_lat_gc_qfed]
            lon_gc_qfed  = [lon_gc_qfed,tmp_lon_gc_qfed]
            alt_gc_qfed  = [alt_gc_qfed,tmp_alt_gc_qfed]
            prs_gc_qfed  = [prs_gc_qfed,tmp_prs_gc_qfed]
     ;; add vocs
            c3h8_gc_qfed   = [c3h8_gc_qfed,tmp_c3h8_gc_qfed]
            c3h6_gc_qfed   = [c3h6_gc_qfed,tmp_c3h6_gc_qfed]
            c2h5oh_gc_qfed = [c2h5oh_gc_qfed,tmp_c2h5oh_gc_qfed]
            c5h8_gc_qfed   = [c5h8_gc_qfed,tmp_c5h8_gc_qfed]
            c7h8_gc_qfed   = [c7h8_gc_qfed,tmp_c7h8_gc_qfed]
            dms_gc_qfed    = [dms_gc_qfed,tmp_dms_gc_qfed]
            mek_gc_qfed    = [mek_gc_qfed,tmp_mek_gc_qfed]
            c8h10_gc_qfed  = [c8h10_gc_qfed,tmp_c8h10_gc_qfed]

            acta_gc_qfed   = [acta_gc_qfed,tmp_acta_gc_qfed]
            macr_mvk_gc_qfed = [macr_mvk_gc_qfed,tmp_macr_mvk_gc_qfed]  
            hcooh_gc_qfed  = [hcooh_gc_qfed,tmp_hcooh_gc_qfed]

            mtpa_gc_qfed  = [mtpa_gc_qfed,tmp_mtpa_gc_qfed]
            limo_gc_qfed  = [limo_gc_qfed,tmp_limo_gc_qfed]
            mtpo_gc_qfed  = [mtpo_gc_qfed,tmp_mtpo_gc_qfed]
        ;;add lumped species
            alk4_gc_qfed  = [alk4_gc_qfed,tmp_alk4_gc_qfed]
            prpe_gc_qfed  = [prpe_gc_qfed,tmp_prpe_gc_qfed]
            rcho_gc_qfed  = [rcho_gc_qfed,tmp_rcho_gc_qfed]

            ;c2h4_gc_qfed  = [c2h4_gc_qfed,tmp_c2h4_gc_qfed]
            c2h6_gc_qfed  = [c2h6_gc_qfed,tmp_c2h6_gc_qfed]
            
            undefine,gc
            
;; ============================
;; GEOS-Chem 0.25x0.3125: THREEGFAS
;; ============================
            if keyword_set(nested) then begin
                ;gcfi_threegfas   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                ;    'planelog2sav/output_'+'threegfas' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'
                gcfi_threegfas   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                    'planelog2sav/output_'+'3Xgfas_tmp' + '/mrg5m_wecan_c130_'+dates[n]+'.sav'
            endif
            if keyword_set(fbf) then begin
                gcfi_threegfas   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                    'planelog2sav/'+'output_WECAN_C2test' + '/mrg30m_wecan_c130_'+dates[n]+'.sav'
            endif
            restore, gcfi_threegfas       
            tmp_utc_gc_threegfas = gc.utc 
            tmp_co_gc_threegfas   = gc.co*1e9
            tmp_O3_gc_threegfas   = gc.o3*1e9
            tmp_pan_gc_threegfas  = gc.pan*1e9
            tmp_hcho_gc_threegfas = gc.ch2o*1e9
            tmp_acet_gc_threegfas = gc.acet*1e9/3; 3 acetone in ppb
            tmp_benz_gc_threegfas = gc.benz*1e9/6; 6C benzene in ppb
    ;        ch3oh_gc= [ch3oh_gc, gc.MOH*1e9]
            tmp_ald2_gc_threegfas = gc.ald2*1e9/2; 2C ch3cho in ppb

            tmp_no_gc_threegfas   = gc.no*1e9
            tmp_no2_gc_threegfas  = gc.no2*1e9
            tmp_so2_gc_threegfas  = gc.so2*1e9

            tmp_na = avo * gc.pres * 100./(8.31 * gc.temp) * 1e-6;; air density molec/cm3
            tmp_oh_gc_threegfas  = gc.oh*tmp_na ;; v/v --> molec/cm3 

            tmp_date_gc_threegfas = gc.date
            ;tmp_utc_gc_threegfas  = gc.utc
            tmp_doy_gc_threegfas  = gc.doy
            tmp_lat_gc_threegfas  = gc.lat
            tmp_lon_gc_threegfas  = gc.lon
            tmp_alt_gc_threegfas  = gc.alt
            tmp_prs_gc_threegfas  = gc.pres
     ;; add vocs
            tmp_c3h8_gc_threegfas   = gc.c3h8*1e9/3; 3C c3h8 in ppb
            tmp_c3h6_gc_threegfas   = gc.prpe*1e9/3; 3C PRPE in ppb
            tmp_c2h5oh_gc_threegfas = gc.eoh*1e9/2; 2C eoh in ppb
            tmp_c5h8_gc_threegfas   = gc.isop*1e9/5; 5C isoprene in ppb
            tmp_c7h8_gc_threegfas   = gc.tolu*1e9/7; 7C toluene in ppb
            tmp_dms_gc_threegfas    = gc.dms*1e9
            tmp_mek_gc_threegfas    = gc.mek*1e9/4 ;4C  MEK in ppb
            tmp_c8h10_gc_threegfas  = gc.xyle*1e9/8 ; 8C xylenes in ppb


            tmp_acta_gc_threegfas   = gc.acta*1e9
            tmp_macr_mvk_gc_threegfas = (gc.MVK+gc.MACR)*1e9  
            tmp_hcooh_gc_threegfas  = gc.hcooh*1e9

            tmp_mtpa_gc_threegfas  = gc.mtpa*1e9
            tmp_limo_gc_threegfas  = gc.limo*1e9
            tmp_mtpo_gc_threegfas  = gc.mtpo*1e9
        ;;add lumped species
            tmp_alk4_gc_threegfas  = gc.alk4*1e9/4.3 ; 4.3C
            tmp_prpe_gc_threegfas   = gc.prpe*1e9/3   ; 3C 
            tmp_rcho_gc_threegfas   = gc.rcho*1e9

            ;tmp_c2h4_gc_threegfas   = gc.c2h4*1e9
            tmp_c2h6_gc_threegfas   = gc.c2h6*1e9/2 ;; 2C


            co_gc_threegfas = [co_gc_threegfas,tmp_co_gc_threegfas]
            O3_gc_threegfas   = [o3_gc_threegfas,tmp_o3_gc_threegfas]
            pan_gc_threegfas  = [pan_gc_threegfas,tmp_pan_gc_threegfas]
            hcho_gc_threegfas = [hcho_gc_threegfas,tmp_hcho_gc_threegfas]
            acet_gc_threegfas = [acet_gc_threegfas,tmp_acet_gc_threegfas]
            benz_gc_threegfas = [benz_gc_threegfas,tmp_benz_gc_threegfas]

            ald2_gc_threegfas = [ald2_gc_threegfas,tmp_ald2_gc_threegfas]


            no_gc_threegfas   = [no_gc_threegfas,tmp_no_gc_threegfas]
            no2_gc_threegfas  = [no2_gc_threegfas,tmp_no2_gc_threegfas]
            so2_gc_threegfas  = [so2_gc_threegfas,tmp_so2_gc_threegfas]

            oh_gc_threegfas  = [oh_gc_threegfas,tmp_oh_gc_threegfas] ;; v/v --> molec/cm3 

            date_gc_threegfas = [date_gc_threegfas,tmp_date_gc_threegfas]
            ;utc_gc_threegfas  = [utc_gc_threegfas,tmp_utc_gc_threegfas]
            doy_gc_threegfas  = [doy_gc_threegfas,tmp_doy_gc_threegfas]
            lat_gc_threegfas  = [lat_gc_threegfas,tmp_lat_gc_threegfas]
            lon_gc_threegfas  = [lon_gc_threegfas,tmp_lon_gc_threegfas]
            alt_gc_threegfas  = [alt_gc_threegfas,tmp_alt_gc_threegfas]
            prs_gc_threegfas  = [prs_gc_threegfas,tmp_prs_gc_threegfas]
     ;; add vocs
            c3h8_gc_threegfas   = [c3h8_gc_threegfas,tmp_c3h8_gc_threegfas]
            c3h6_gc_threegfas   = [c3h6_gc_threegfas,tmp_c3h6_gc_threegfas]
            c2h5oh_gc_threegfas = [c2h5oh_gc_threegfas,tmp_c2h5oh_gc_threegfas]
            c5h8_gc_threegfas   = [c5h8_gc_threegfas,tmp_c5h8_gc_threegfas]
            c7h8_gc_threegfas   = [c7h8_gc_threegfas,tmp_c7h8_gc_threegfas]
            dms_gc_threegfas    = [dms_gc_threegfas,tmp_dms_gc_threegfas]
            mek_gc_threegfas    = [mek_gc_threegfas,tmp_mek_gc_threegfas]
            c8h10_gc_threegfas  = [c8h10_gc_threegfas,tmp_c8h10_gc_threegfas]

            acta_gc_threegfas   = [acta_gc_threegfas,tmp_acta_gc_threegfas]
            macr_mvk_gc_threegfas = [macr_mvk_gc_threegfas,tmp_macr_mvk_gc_threegfas]  
            hcooh_gc_threegfas  = [hcooh_gc_threegfas,tmp_hcooh_gc_threegfas]

            mtpa_gc_threegfas  = [mtpa_gc_threegfas,tmp_mtpa_gc_threegfas]
            limo_gc_threegfas  = [limo_gc_threegfas,tmp_limo_gc_threegfas]
            mtpo_gc_threegfas  = [mtpo_gc_threegfas,tmp_mtpo_gc_threegfas]
        ;;add lumped species
            alk4_gc_threegfas  = [alk4_gc_threegfas,tmp_alk4_gc_threegfas]
            prpe_gc_threegfas  = [prpe_gc_threegfas,tmp_prpe_gc_threegfas]
            rcho_gc_threegfas  = [rcho_gc_threegfas,tmp_rcho_gc_threegfas]

            ;c2h4_gc_threegfas  = [c2h4_gc_threegfas,tmp_c2h4_gc_threegfas]
            c2h6_gc_threegfas  = [c2h6_gc_threegfas,tmp_c2h6_gc_threegfas]
            
            undefine,gc

;; ============================
;; GEOS-Chem 0.25x0.3125: NOBB
;; ============================
            if keyword_set(nested) then begin
                ;gcfi_nobb   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                ;    'planelog2sav/output_'+'nobb' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'
                gcfi_nobb   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                    'planelog2sav/output_'+'nobb_tmp' + '/mrg5m_wecan_c130_'+dates[n]+'.sav'
            endif
            if keyword_set(fbf) then begin
                gcfi_nobb   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                    'planelog2sav/'+'output_WECAN_C2test' + '/mrg30m_wecan_c130_'+dates[n]+'.sav'
            endif
            restore, gcfi_nobb
            tmp_utc_gc_nobb = gc.utc 
            tmp_co_gc_nobb   = gc.co*1e9
            tmp_O3_gc_nobb   = gc.o3*1e9
            tmp_pan_gc_nobb  = gc.pan*1e9
            tmp_hcho_gc_nobb = gc.ch2o*1e9
            tmp_acet_gc_nobb = gc.acet*1e9/3; 3 acetone in ppb
            tmp_benz_gc_nobb = gc.benz*1e9/6; 6C benzene in ppb
    ;        ch3oh_gc= [ch3oh_gc, gc.MOH*1e9]
            tmp_ald2_gc_nobb = gc.ald2*1e9/2; 2C ch3cho in ppb

            tmp_no_gc_nobb   = gc.no*1e9
            tmp_no2_gc_nobb  = gc.no2*1e9
            tmp_so2_gc_nobb  = gc.so2*1e9

            tmp_na = avo * gc.pres * 100./(8.31 * gc.temp) * 1e-6;; air density molec/cm3
            tmp_oh_gc_nobb  = gc.oh*tmp_na ;; v/v --> molec/cm3 

            tmp_date_gc_nobb = gc.date
            ;tmp_utc_gc_nobb  = gc.utc
            tmp_doy_gc_nobb  = gc.doy
            tmp_lat_gc_nobb  = gc.lat
            tmp_lon_gc_nobb  = gc.lon
            tmp_alt_gc_nobb  = gc.alt
            tmp_prs_gc_nobb  = gc.pres
     ;; add vocs
            tmp_c3h8_gc_nobb   = gc.c3h8*1e9/3; 3C c3h8 in ppb
            tmp_c3h6_gc_nobb   = gc.prpe*1e9/3; 3C PRPE in ppb
            tmp_c2h5oh_gc_nobb = gc.eoh*1e9/2; 2C eoh in ppb
            tmp_c5h8_gc_nobb   = gc.isop*1e9/5; 5C isoprene in ppb
            tmp_c7h8_gc_nobb   = gc.tolu*1e9/7; 7C toluene in ppb
            tmp_dms_gc_nobb    = gc.dms*1e9
            tmp_mek_gc_nobb    = gc.mek*1e9/4 ;4C  MEK in ppb
            tmp_c8h10_gc_nobb  = gc.xyle*1e9/8 ; 8C xylenes in ppb
            
            tmp_acta_gc_nobb   = gc.acta*1e9
            tmp_macr_mvk_gc_nobb = (gc.MVK+gc.MACR)*1e9  
            tmp_hcooh_gc_nobb  = gc.hcooh*1e9

            tmp_mtpa_gc_nobb  = gc.mtpa*1e9
            tmp_limo_gc_nobb  = gc.limo*1e9
            tmp_mtpo_gc_nobb  = gc.mtpo*1e9
        ;;add lumped species
            tmp_alk4_gc_nobb   = gc.alk4*1e9/4.3 ;;4.3C
            tmp_prpe_gc_nobb   = gc.prpe*1e9/3 ;;3C
            tmp_rcho_gc_nobb   = gc.rcho*1e9

            ;tmp_c2h4_gc_nobb   = gc.c2h4*1e9 
            tmp_c2h6_gc_nobb   = gc.c2h6*1e9/2 ;;2C

            co_gc_nobb = [co_gc_nobb,tmp_co_gc_nobb]
            O3_gc_nobb   = [o3_gc_nobb,tmp_o3_gc_nobb]
            pan_gc_nobb  = [pan_gc_nobb,tmp_pan_gc_nobb]
            hcho_gc_nobb = [hcho_gc_nobb,tmp_hcho_gc_nobb]
            acet_gc_nobb = [acet_gc_nobb,tmp_acet_gc_nobb]
            benz_gc_nobb = [benz_gc_nobb,tmp_benz_gc_nobb]

            ald2_gc_nobb = [ald2_gc_nobb,tmp_ald2_gc_nobb]


            no_gc_nobb   = [no_gc_nobb,tmp_no_gc_nobb]
            no2_gc_nobb  = [no2_gc_nobb,tmp_no2_gc_nobb]
            so2_gc_nobb  = [so2_gc_nobb,tmp_so2_gc_nobb]

            oh_gc_nobb  = [oh_gc_nobb,tmp_oh_gc_nobb] ;; v/v --> molec/cm3 

            date_gc_nobb = [date_gc_nobb,tmp_date_gc_nobb]
            ;utc_gc_nobb  = [utc_gc_nobb,tmp_utc_gc_nobb]
            doy_gc_nobb  = [doy_gc_nobb,tmp_doy_gc_nobb]
            lat_gc_nobb  = [lat_gc_nobb,tmp_lat_gc_nobb]
            lon_gc_nobb  = [lon_gc_nobb,tmp_lon_gc_nobb]
            alt_gc_nobb  = [alt_gc_nobb,tmp_alt_gc_nobb]
            prs_gc_nobb  = [prs_gc_nobb,tmp_prs_gc_nobb]
     ;; add vocs
            c3h8_gc_nobb   = [c3h8_gc_nobb,tmp_c3h8_gc_nobb]
            c3h6_gc_nobb   = [c3h6_gc_nobb,tmp_c3h6_gc_nobb]
            c2h5oh_gc_nobb = [c2h5oh_gc_nobb,tmp_c2h5oh_gc_nobb]
            c5h8_gc_nobb   = [c5h8_gc_nobb,tmp_c5h8_gc_nobb]
            c7h8_gc_nobb   = [c7h8_gc_nobb,tmp_c7h8_gc_nobb]
            dms_gc_nobb    = [dms_gc_nobb,tmp_dms_gc_nobb]
            mek_gc_nobb    = [mek_gc_nobb,tmp_mek_gc_nobb]
            c8h10_gc_nobb  = [c8h10_gc_nobb,tmp_c8h10_gc_nobb]

            acta_gc_nobb   = [acta_gc_nobb,tmp_acta_gc_nobb]
            macr_mvk_gc_nobb = [macr_mvk_gc_nobb,tmp_macr_mvk_gc_nobb]  
            hcooh_gc_nobb  = [hcooh_gc_nobb,tmp_hcooh_gc_nobb]

            mtpa_gc_nobb  = [mtpa_gc_nobb,tmp_mtpa_gc_nobb]
            limo_gc_nobb  = [limo_gc_nobb,tmp_limo_gc_nobb]
            mtpo_gc_nobb  = [mtpo_gc_nobb,tmp_mtpo_gc_nobb]
        ;;add lumped species
            alk4_gc_nobb  = [alk4_gc_nobb,tmp_alk4_gc_nobb]
            prpe_gc_nobb  = [prpe_gc_nobb,tmp_prpe_gc_nobb]
            rcho_gc_nobb  = [rcho_gc_nobb,tmp_rcho_gc_nobb]

            ;c2h4_gc_nobb  = [c2h4_gc_nobb,tmp_c2h4_gc_nobb]
            c2h6_gc_nobb  = [c2h6_gc_nobb,tmp_c2h6_gc_nobb]
    
            undefine,gc            
        endfor ;; loop for dates, read files
        
;; ==========================================================
;; Remove placeholder: reduced version
;; =========================================================
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
;; GEOS-Chem 0.25x0.3125: GFASs
;; ============================
    ;; GFED4
        co_gc_gfed4    = co_gc_gfed4[1:*]
        o3_gc_gfed4    = o3_gc_gfed4[1:*]
        pan_gc_gfed4   = pan_gc_gfed4[1:*]
        hcho_gc_gfed4  = hcho_gc_gfed4[1:*]
        acet_gc_gfed4  = acet_gc_gfed4[1:*]
        benz_gc_gfed4  = benz_gc_gfed4[1:*]
    ;    ch3oh_gc= ch3oh_gc[1:*]
        ald2_gc_gfed4  = ald2_gc_gfed4[1:*]

        no_gc_gfed4    = no_gc_gfed4[1:*]
        no2_gc_gfed4   = no2_gc_gfed4[1:*]
        so2_gc_gfed4   = so2_gc_gfed4[1:*]
        oh_gc_gfed4    = oh_gc_gfed4[1:*] 

        date_gc_gfed4  = date_gc_gfed4[1:*]
        ;utc_gc_gfed4   = utc_gc_gfed4[1:*]
        doy_gc_gfed4   = doy_gc_gfed4[1:*]
        lat_gc_gfed4   = lat_gc_gfed4[1:*]
        lon_gc_gfed4   = lon_gc_gfed4[1:*]
        alt_gc_gfed4   = alt_gc_gfed4[1:*]
        prs_gc_gfed4   = prs_gc_gfed4[1:*]

     ;; add vocs
        c3h8_gc_gfed4   = c3h8_gc_gfed4[1:*]
        c3h6_gc_gfed4   = c3h6_gc_gfed4[1:*]
        c2h5oh_gc_gfed4 = c2h5oh_gc_gfed4[1:*]
        c5h8_gc_gfed4   = c5h8_gc_gfed4[1:*]
        c7h8_gc_gfed4   = c7h8_gc_gfed4[1:*]
        dms_gc_gfed4    = dms_gc_gfed4[1:*]
        mek_gc_gfed4    = mek_gc_gfed4[1:*]
        c8h10_gc_gfed4  = c8h10_gc_gfed4[1:*]

        acta_gc_gfed4   = acta_gc_gfed4[1:*]
        macr_mvk_gc_gfed4 = macr_mvk_gc_gfed4[1:*]
        hcooh_gc_gfed4  = hcooh_gc_gfed4[1:*]

        mtpa_gc_gfed4  = mtpa_gc_gfed4[1:*]
        limo_gc_gfed4  = limo_gc_gfed4[1:*]
        mtpo_gc_gfed4  = mtpo_gc_gfed4[1:*]
        
        ;;add lumped species
        alk4_gc_gfed4 = alk4_gc_gfed4[1:*]
        prpe_gc_gfed4 = prpe_gc_gfed4[1:*]
        rcho_gc_gfed4 = rcho_gc_gfed4[1:*]

        c2h6_gc_gfed4 = c2h6_gc_gfed4[1:*]
    ;; FINN
        co_gc_finn    = co_gc_finn[1:*]
        o3_gc_finn    = o3_gc_finn[1:*]
        pan_gc_finn   = pan_gc_finn[1:*]
        hcho_gc_finn  = hcho_gc_finn[1:*]
        acet_gc_finn  = acet_gc_finn[1:*]
        benz_gc_finn  = benz_gc_finn[1:*]
    ;    ch3oh_gc= ch3oh_gc[1:*]
        ald2_gc_finn  = ald2_gc_finn[1:*]

        no_gc_finn    = no_gc_finn[1:*]
        no2_gc_finn   = no2_gc_finn[1:*]
        so2_gc_finn   = so2_gc_finn[1:*]
        oh_gc_finn    = oh_gc_finn[1:*] 

        date_gc_finn  = date_gc_finn[1:*]
        ;utc_gc_finn   = utc_gc_finn[1:*]
        doy_gc_finn   = doy_gc_finn[1:*]
        lat_gc_finn   = lat_gc_finn[1:*]
        lon_gc_finn   = lon_gc_finn[1:*]
        alt_gc_finn   = alt_gc_finn[1:*]
        prs_gc_finn   = prs_gc_finn[1:*]

     ;; add vocs
        c3h8_gc_finn   = c3h8_gc_finn[1:*]
        c3h6_gc_finn   = c3h6_gc_finn[1:*]
        c2h5oh_gc_finn = c2h5oh_gc_finn[1:*]
        c5h8_gc_finn   = c5h8_gc_finn[1:*]
        c7h8_gc_finn   = c7h8_gc_finn[1:*]
        dms_gc_finn    = dms_gc_finn[1:*]
        mek_gc_finn    = mek_gc_finn[1:*]
        c8h10_gc_finn  = c8h10_gc_finn[1:*]

        acta_gc_finn   = acta_gc_finn[1:*]
        macr_mvk_gc_finn = macr_mvk_gc_finn[1:*]
        hcooh_gc_finn  = hcooh_gc_finn[1:*]

        mtpa_gc_finn  = mtpa_gc_finn[1:*]
        limo_gc_finn  = limo_gc_finn[1:*]
        mtpo_gc_finn  = mtpo_gc_finn[1:*]    
    
         ;;add lumped species
        alk4_gc_finn = alk4_gc_finn[1:*]
        prpe_gc_finn = prpe_gc_finn[1:*]
        rcho_gc_finn = rcho_gc_finn[1:*]

        c2h6_gc_finn = c2h6_gc_finn[1:*]
    ;; GFAS
        co_gc_gfas   = co_gc_gfas[1:*]
        o3_gc_gfas   = o3_gc_gfas[1:*]
        pan_gc_gfas  = pan_gc_gfas[1:*]
        hcho_gc_gfas = hcho_gc_gfas[1:*]
        acet_gc_gfas = acet_gc_gfas[1:*]
        benz_gc_gfas = benz_gc_gfas[1:*]
    ;    ch3oh_gc= ch3oh_gc[1:*]
        ald2_gc_gfas = ald2_gc_gfas[1:*]

        no_gc_gfas   = no_gc_gfas[1:*]
        no2_gc_gfas  = no2_gc_gfas[1:*]
        so2_gc_gfas  = so2_gc_gfas[1:*]
        oh_gc_gfas   = oh_gc_gfas[1:*] 

        date_gc_gfas = date_gc_gfas[1:*]
        ;utc_gc_gfas  = utc_gc_gfas[1:*]
        doy_gc_gfas  = doy_gc_gfas[1:*]
        lat_gc_gfas  = lat_gc_gfas[1:*]
        lon_gc_gfas  = lon_gc_gfas[1:*]
        alt_gc_gfas  = alt_gc_gfas[1:*]
        prs_gc_gfas  = prs_gc_gfas[1:*]
        ;;temperature
        temperature_gc_gfas = temperature_gc_gfas[1:*]
     ;; add vocs
        c3h8_gc_gfas   = c3h8_gc_gfas[1:*]
        c3h6_gc_gfas   = c3h6_gc_gfas[1:*]
        c2h5oh_gc_gfas = c2h5oh_gc_gfas[1:*]
        c5h8_gc_gfas   = c5h8_gc_gfas[1:*]
        c7h8_gc_gfas   = c7h8_gc_gfas[1:*]
        dms_gc_gfas    = dms_gc_gfas[1:*]
        mek_gc_gfas    = mek_gc_gfas[1:*]
        c8h10_gc_gfas  = c8h10_gc_gfas[1:*]

        acta_gc_gfas   = acta_gc_gfas[1:*]
        macr_mvk_gc_gfas = macr_mvk_gc_gfas[1:*]
        hcooh_gc_gfas  = hcooh_gc_gfas[1:*]

        mtpa_gc_gfas  = mtpa_gc_gfas[1:*]
        limo_gc_gfas  = limo_gc_gfas[1:*]
        mtpo_gc_gfas  = mtpo_gc_gfas[1:*]
        
        ;;add lumped species
        alk4_gc_gfas = alk4_gc_gfas[1:*]
        prpe_gc_gfas = prpe_gc_gfas[1:*]
        rcho_gc_gfas = rcho_gc_gfas[1:*]

        c2h6_gc_gfas = c2h6_gc_gfas[1:*]


    ;; QFED
        co_gc_qfed    = co_gc_qfed[1:*]
        o3_gc_qfed    = o3_gc_qfed[1:*]
        pan_gc_qfed   = pan_gc_qfed[1:*]
        hcho_gc_qfed  = hcho_gc_qfed[1:*]
        acet_gc_qfed  = acet_gc_qfed[1:*]
        benz_gc_qfed  = benz_gc_qfed[1:*]
    ;    ch3oh_gc= ch3oh_gc[1:*]
        ald2_gc_qfed  = ald2_gc_qfed[1:*]

        no_gc_qfed    = no_gc_qfed[1:*]
        no2_gc_qfed   = no2_gc_qfed[1:*]
        so2_gc_qfed   = so2_gc_qfed[1:*]
        oh_gc_qfed    = oh_gc_qfed[1:*] 

        date_gc_qfed  = date_gc_qfed[1:*]
        ;utc_gc_qfed   = utc_gc_qfed[1:*]
        doy_gc_qfed   = doy_gc_qfed[1:*]
        lat_gc_qfed   = lat_gc_qfed[1:*]
        lon_gc_qfed   = lon_gc_qfed[1:*]
        alt_gc_qfed   = alt_gc_qfed[1:*]
        prs_gc_qfed   = prs_gc_qfed[1:*]

     ;; add vocs
        c3h8_gc_qfed   = c3h8_gc_qfed[1:*]
        c3h6_gc_qfed   = c3h6_gc_qfed[1:*]
        c2h5oh_gc_qfed = c2h5oh_gc_qfed[1:*]
        c5h8_gc_qfed   = c5h8_gc_qfed[1:*]
        c7h8_gc_qfed   = c7h8_gc_qfed[1:*]
        dms_gc_qfed    = dms_gc_qfed[1:*]
        mek_gc_qfed    = mek_gc_qfed[1:*]
        c8h10_gc_qfed  = c8h10_gc_qfed[1:*]

        acta_gc_qfed   = acta_gc_qfed[1:*]
        macr_mvk_gc_qfed = macr_mvk_gc_qfed[1:*]
        hcooh_gc_qfed  = hcooh_gc_qfed[1:*]

        mtpa_gc_qfed  = mtpa_gc_qfed[1:*]
        limo_gc_qfed  = limo_gc_qfed[1:*]
        mtpo_gc_qfed  = mtpo_gc_qfed[1:*]

        ;;add lumped species
        alk4_gc_qfed = alk4_gc_qfed[1:*]
        prpe_gc_qfed = prpe_gc_qfed[1:*]
        rcho_gc_qfed = rcho_gc_qfed[1:*]

        c2h6_gc_qfed = c2h6_gc_qfed[1:*]


    ;; THREEGFAS
        co_gc_threegfas    = co_gc_threegfas[1:*]
        o3_gc_threegfas    = o3_gc_threegfas[1:*]
        pan_gc_threegfas   = pan_gc_threegfas[1:*]
        hcho_gc_threegfas  = hcho_gc_threegfas[1:*]
        acet_gc_threegfas  = acet_gc_threegfas[1:*]
        benz_gc_threegfas  = benz_gc_threegfas[1:*]
    ;    ch3oh_gc= ch3oh_gc[1:*]
        ald2_gc_threegfas  = ald2_gc_threegfas[1:*]

        no_gc_threegfas    = no_gc_threegfas[1:*]
        no2_gc_threegfas   = no2_gc_threegfas[1:*]
        so2_gc_threegfas   = so2_gc_threegfas[1:*]
        oh_gc_threegfas    = oh_gc_threegfas[1:*] 

        date_gc_threegfas  = date_gc_threegfas[1:*]
        ;utc_gc_threegfas   = utc_gc_threegfas[1:*]
        doy_gc_threegfas   = doy_gc_threegfas[1:*]
        lat_gc_threegfas   = lat_gc_threegfas[1:*]
        lon_gc_threegfas   = lon_gc_threegfas[1:*]
        alt_gc_threegfas   = alt_gc_threegfas[1:*]
        prs_gc_threegfas   = prs_gc_threegfas[1:*]

     ;; add vocs
        c3h8_gc_threegfas   = c3h8_gc_threegfas[1:*]
        c3h6_gc_threegfas   = c3h6_gc_threegfas[1:*]
        c2h5oh_gc_threegfas = c2h5oh_gc_threegfas[1:*]
        c5h8_gc_threegfas   = c5h8_gc_threegfas[1:*]
        c7h8_gc_threegfas   = c7h8_gc_threegfas[1:*]
        dms_gc_threegfas    = dms_gc_threegfas[1:*]
        mek_gc_threegfas    = mek_gc_threegfas[1:*]
        c8h10_gc_threegfas  = c8h10_gc_threegfas[1:*]

        acta_gc_threegfas   = acta_gc_threegfas[1:*]
        macr_mvk_gc_threegfas = macr_mvk_gc_threegfas[1:*]
        hcooh_gc_threegfas  = hcooh_gc_threegfas[1:*]

        mtpa_gc_threegfas  = mtpa_gc_threegfas[1:*]
        limo_gc_threegfas  = limo_gc_threegfas[1:*]
        mtpo_gc_threegfas  = mtpo_gc_threegfas[1:*]

        ;;add lumped species
        alk4_gc_threegfas = alk4_gc_threegfas[1:*]
        prpe_gc_threegfas = prpe_gc_threegfas[1:*]
        rcho_gc_threegfas = rcho_gc_threegfas[1:*]

        c2h6_gc_threegfas = c2h6_gc_threegfas[1:*]

    ;; NOBB
        co_gc_nobb    = co_gc_nobb[1:*]
        o3_gc_nobb    = o3_gc_nobb[1:*]
        pan_gc_nobb   = pan_gc_nobb[1:*]
        hcho_gc_nobb  = hcho_gc_nobb[1:*]
        acet_gc_nobb  = acet_gc_nobb[1:*]
        benz_gc_nobb  = benz_gc_nobb[1:*]
    ;    ch3oh_gc= ch3oh_gc[1:*]
        ald2_gc_nobb  = ald2_gc_nobb[1:*]

        no_gc_nobb    = no_gc_nobb[1:*]
        no2_gc_nobb   = no2_gc_nobb[1:*]
        so2_gc_nobb   = so2_gc_nobb[1:*]
        oh_gc_nobb    = oh_gc_nobb[1:*] 

        date_gc_nobb  = date_gc_nobb[1:*]
        ;utc_gc_nobb   = utc_gc_nobb[1:*]
        doy_gc_nobb   = doy_gc_nobb[1:*]
        lat_gc_nobb   = lat_gc_nobb[1:*]
        lon_gc_nobb   = lon_gc_nobb[1:*]
        alt_gc_nobb   = alt_gc_nobb[1:*]
        prs_gc_nobb   = prs_gc_nobb[1:*]

     ;; add vocs
        c3h8_gc_nobb   = c3h8_gc_nobb[1:*]
        c3h6_gc_nobb   = c3h6_gc_nobb[1:*]
        c2h5oh_gc_nobb = c2h5oh_gc_nobb[1:*]
        c5h8_gc_nobb   = c5h8_gc_nobb[1:*]
        c7h8_gc_nobb   = c7h8_gc_nobb[1:*]
        dms_gc_nobb    = dms_gc_nobb[1:*]
        mek_gc_nobb    = mek_gc_nobb[1:*]
        c8h10_gc_nobb  = c8h10_gc_nobb[1:*]

        acta_gc_nobb   = acta_gc_nobb[1:*]
        macr_mvk_gc_nobb = macr_mvk_gc_nobb[1:*]
        hcooh_gc_nobb  = hcooh_gc_nobb[1:*]

        mtpa_gc_nobb  = mtpa_gc_nobb[1:*]
        limo_gc_nobb  = limo_gc_nobb[1:*]
        mtpo_gc_nobb  = mtpo_gc_nobb[1:*]    
    
         ;;add lumped species
        alk4_gc_nobb = alk4_gc_nobb[1:*]
        prpe_gc_nobb = prpe_gc_nobb[1:*]
        rcho_gc_nobb = rcho_gc_nobb[1:*]

        c2h6_gc_nobb = c2h6_gc_nobb[1:*]
;; in 5min merged, the maximum is 4.6ppm

        if n_elements(co_obs) ne n_elements(co_gc_gfed4) or $
            ;n_elements(co_obs) ne n_elements(co_gc_finn) or $
            n_elements(co_obs) ne n_elements(co_gc_gfas) or $
            n_elements(co_obs) ne n_elements(co_gc_qfed) or $
            n_elements(co_obs) ne n_elements(co_gc_threegfas) or $
            n_elements(co_obs) ne n_elements(co_gc_nobb) then stop 

;; test
    ;print, median(hcooh_gc_threegfas/hcooh_gc_gfas)
    ;print, median(acta_gc_threegfas/acta_gc_gfas)
    ;print, median(c2h5oh_gc_threegfas/c2h5oh_gc_gfas)
    ;print, median(mek_gc_threegfas/mek_gc_gfas)
    ;print, median(mtpa_gc_threegfas/mtpa_gc_gfas)
    ;print, median(macr_gc_threegfas/macr_gc_gfas)


;; test for percentile
;data_test = ch3cn_obs
;data_test = hcn_obs
;data_test = co_obs
;ind = where(data_test gt 0, ct)
;if ct gt 0 then data_test = data_test[ind]
;print, cgPercentiles(data_test, Percentiles=[0.1,0.25, 0.5, 0.75,0.9]) ; output(131ppt, 155ppt, 257ppt, 426ppt, 717ppt)
;continue

;;setting up koh of different VOCs
        temperature_mod = temperature_gc_gfas
        prs_mod = prs_gc_gfas

;; number density of air: molec/m-3
;; used to get number density of trace gas 
        scalefactor_obs = avo*(prs_obs*100)/(8.31*temperature_obs)
        scalefactor_mod = avo*(prs_mod*100)/(8.31*temperature_mod)

;;Step1: using the GC equation
        ;; CO
        tmp_koh_co_obs = GCRateConst_co(temperature_obs,prs_obs)
        tmp_koh_co_mod = GCRateConst_co(temperature_mod,prs_mod)

        ;; Propane
        tmp_koh_c3h8_1_obs = GCRateConst_c3h8(7.60E-12, 0.0E+00, -585.0, 5.87E0, 0.64E0, -816.0, temperature_obs)
        tmp_koh_c3h8_1_mod = GCRateConst_c3h8(7.60E-12, 0.0E+00, -585.0, 5.87E0, 0.64E0, -816.0, temperature_mod) 
        tmp_koh_c3h8_2_obs = GCRateConst_c3h8(7.60E-12, 0.0E+00, -585.0, 1.7E-1, -0.64E0, 816.0, temperature_obs) 
        tmp_koh_c3h8_2_mod = GCRateConst_c3h8(7.60E-12, 0.0E+00, -585.0, 1.7E-1, -0.64E0, 816.0, temperature_mod) 
        tmp_koh_c3h8_obs = (tmp_koh_c3h8_1_obs+tmp_koh_c3h8_2_obs)/2
        tmp_koh_c3h8_mod = (tmp_koh_c3h8_1_mod+tmp_koh_c3h8_2_mod)/2
        
        ;; C2H6
        tmp_koh_c2h6_obs = GCRateConst(7.66E-12, -1020, temperature_obs)
        tmp_koh_c2h6_mod = GCRateConst(7.66E-12, -1020, temperature_mod)
        
        ;; ALK4
        tmp_koh_alk4_obs = GCRateConst(9.10E-12, -405, temperature_obs)
        tmp_koh_alk4_mod = GCRateConst(9.10E-12, -405, temperature_mod)
        
        ;; PRPE
        tmp_koh_prpe_obs = GCRateConst(2.00E-12, -840, temperature_obs)
        tmp_koh_prpe_mod = GCRateConst(2.00E-12, -840, temperature_mod)
        ;; acids
        tmp_koh_hcooh_obs = GCRateConst(4.00E-13, 0, temperature_obs)
        tmp_koh_hcooh_mod = GCRateConst(4.00E-13, 0, temperature_mod)
        tmp_koh_acta_obs = GCRateConst(3.15E-14, 920, temperature_obs)
        tmp_koh_acta_mod = GCRateConst(3.15E-14, 920, temperature_mod)
        
        ;; aldehydes
        tmp_koh_ch2o_obs = GCRateConst(5.50E-12,125, temperature_obs)
        tmp_koh_ch2o_mod = GCRateConst(5.50E-12,125, temperature_mod)
        tmp_koh_ald2_obs = GCRateConst(4.63E-12,350, temperature_obs)
        tmp_koh_ald2_mod = GCRateConst(4.63E-12,350, temperature_mod)
        
        tmp_koh_rcho_obs = GCRateConst(6.0E-12, 410.0, temperature_obs)
        tmp_koh_rcho_mod = GCRateConst(6.0E-12, 410.0, temperature_mod)
        
        ;; Ketones
        tmp_koh_acet_obs = GCRateConst_acet(temperature_obs)
        tmp_koh_acet_mod = GCRateConst_acet(temperature_mod)
        tmp_koh_mek_obs = GCRateConst(1.30E-12,-25, temperature_obs)
        tmp_koh_mek_mod = GCRateConst(1.30E-12,-25, temperature_mod)
        
        ;; Aromatics
        tmp_koh_benz_obs = GCRateConst(2.33E-12,-193, temperature_obs)
        tmp_koh_benz_mod = GCRateConst(2.33E-12,-193, temperature_mod)
        tmp_koh_tolu_obs = GCRateConst(1.81E-12,338, temperature_obs)
        tmp_koh_tolu_mod = GCRateConst(1.81E-12,338, temperature_mod)
        tmp_koh_xyle_obs = GCRateConst(2.31E-11,0, temperature_obs)
        tmp_koh_xyle_mod = GCRateConst(2.31E-11,0, temperature_mod)

        
        ;; alcohols
        ;tmp_koh_moh_obs = GCRateConst(2.90E-12,-345, temperature_obs)
        ;tmp_koh_moh_mod = GCRateConst(2.90E-12,-345, temperature_mod)
        ;tmp_koh_eoh_obs = GCRateConst(3.35E-12,0, temperature_obs)
        ;tmp_koh_eoh_mod = GCRateConst(3.35E-12,0, temperature_mod)
        ;; Isoprene
        ;; MTPA
        ;tmp_koh_mtpa_obs = GCRateConst(1.21E-11, 440, temperature_obs)
        ;tmp_koh_mtpa_mod = GCRateConst(1.21E-11, 440, temperature_mod)
        
        ;;Others, may not important
        
;; ==============================
;; set up detection limit for PTR
;; not for missing values
;; ==============================   
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
        
    ;; get hhmm of each day
        utc_obs = 24.*(doy_obs mod 1.)
        utc_gc_gfed4  = 24.*(doy_gc_gfed4  mod 1.)
        utc_gc_finn  = 24.*(doy_gc_finn  mod 1.)
        utc_gc_gfas  = 24.*(doy_gc_gfas  mod 1.)
        utc_gc_qfed  = 24.*(doy_gc_qfed  mod 1.)
        utc_gc_threegfas  = 24.*(doy_gc_threegfas  mod 1.)
        utc_gc_nobb  = 24.*(doy_gc_nobb  mod 1.)

        
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

        if keyword_set(all) then s_end = 17
        if keyword_set(obs_only) then s_end = 8
        if keyword_set(each) then s_end = 0

        for s = 0, s_end do begin
    ;;all
            if keyword_set(obs_only) then begin
                case s of                     
                    ;;hcho
                    0:begin
                        ;; Emissions
                        ovoc_ptr = ch2o_obs_ptr
                        ovoc_other = ch2o_obs_toga

                        gvoc_gfas = hcho_gc_gfas
                        gvoc_gfed4 = hcho_gc_gfed4
                        gvoc_qfed  = hcho_gc_qfed
                        gvoc_finn  = hcho_gc_finn
                        gvoc_threegfas = hcho_gc_threegfas
                        gvoc_nobb = hcho_gc_nobb

                        title = 'Formaldehyde'
                        xrange = [0,6]
                        if keyword_set(filter2_nobb) and keyword_set(filter3_nobb) then xrange = [0, 2.4]
                        koh_voc_obs = tmp_koh_co_obs*scalefactor_obs*1e-15  ; dummy
                        koh_voc_mod = tmp_koh_co_mod*scalefactor_mod*1e-15  ; dummy
                    end
                    ;;ald2
                    1:begin
                        ;; Emissions
                        ovoc_ptr = ald2_obs_ptr
                        ovoc_other = ald2_obs_toga
                        gvoc_gfas = ald2_gc_gfas
                        gvoc_gfed4 = ald2_gc_gfed4
                        gvoc_qfed  = ald2_gc_qfed
                        gvoc_finn  = ald2_gc_finn
                        gvoc_threegfas = ald2_gc_threegfas
                        gvoc_nobb = ald2_gc_nobb

                        title = 'Acetaldehyde'
                        xrange = [0,3]
                        
                        if keyword_set(filter2_nobb) and keyword_set(filter3_nobb) then xrange = [0, 2.4]
                        koh_voc_obs = tmp_koh_co_obs*scalefactor_obs*1e-15  ; dummy
                        koh_voc_mod = tmp_koh_co_mod*scalefactor_mod*1e-15  ; dummy
                    end
                    ;;acet
                    2:begin
                        ;; Emissions
                        ovoc_ptr = acet_obs_ptr
                        ovoc_other = acet_obs_toga
                        gvoc_gfas = acet_gc_gfas
                        gvoc_gfed4 = acet_gc_gfed4
                        gvoc_qfed  = acet_gc_qfed
                        gvoc_finn  = acet_gc_finn
                        gvoc_threegfas = acet_gc_threegfas
                        gvoc_nobb = acet_gc_nobb
                        title = 'Acetone'
                        xrange = [0,5]
                        
                        if keyword_set(filter2_nobb) and keyword_set(filter3_nobb) then xrange = [0, 2]
                        koh_voc_obs = tmp_koh_co_obs*scalefactor_obs*1e-15  ; dummy
                        koh_voc_mod = tmp_koh_co_mod*scalefactor_mod*1e-15  ; dummy
                    end
                    ;;Formic acid
                    3:begin
                        ovoc_ptr = hcooh_obs_ptr
                        ovoc_other = hcooh_obs_cims
                        gvoc_gfas = hcooh_gc_gfas
                        gvoc_gfed4 = hcooh_gc_gfed4
                        gvoc_qfed  = hcooh_gc_qfed
                        gvoc_finn  = hcooh_gc_finn
                        gvoc_threegfas = hcooh_gc_threegfas
                        gvoc_nobb = hcooh_gc_nobb

                        title = 'HCOOH'
                        xrange = [0,15]
                        
                        if keyword_set(filter2_nobb) and keyword_set(filter3_nobb) then xrange = [0, 8]
                        koh_voc_obs = tmp_koh_co_obs*scalefactor_obs*1e-15  ; dummy
                        koh_voc_mod = tmp_koh_co_mod*scalefactor_mod*1e-15  ; dummy
                    end
                    ;;mek
                    4:begin
                        ;; Emissions
                        ovoc_ptr = mek_obs_ptr
                        ovoc_other = mek_obs_toga
                        gvoc_gfas = mek_gc_gfas
                        gvoc_gfed4 = mek_gc_gfed4
                        gvoc_qfed  = mek_gc_qfed
                        gvoc_finn  = mek_gc_finn
                        gvoc_threegfas = mek_gc_threegfas
                        gvoc_nobb = mek_gc_nobb

                        title = 'MEK'
                        xrange = [0,0.4]
                        koh_voc_obs = tmp_koh_co_obs*scalefactor_obs*1e-15  ; dummy
                        koh_voc_mod = tmp_koh_co_mod*scalefactor_mod*1e-15  ; dummy
                        
                    end

                    ;;MOH
                    5:begin
                        ovoc_ptr  = moh_obs_ptr
                        ovoc_other = moh_obs_toga
                        gvoc_gfas = mek_gc_gfas
                        gvoc_gfed4 = mek_gc_gfed4
                        gvoc_qfed  = mek_gc_qfed
                        gvoc_finn  = mek_gc_finn
                        gvoc_threegfas = mek_gc_threegfas
                        gvoc_nobb = mek_gc_nobb

                        title = 'MEK'
                        xrange = [0,0.4]

                        title = 'MOH'
                        xrange = [0,12]
                        
                        if keyword_set(filter2_nobb) and keyword_set(filter3_nobb) then xrange = [0, 6]
                        koh_voc_obs = tmp_koh_co_obs*scalefactor_obs*1e-15  ; dummy
                        koh_voc_mod = tmp_koh_co_mod*scalefactor_mod*1e-15  ; dummy
                    end
                    ;;benz
                    6:begin
                        ;; Emissions
                        ovoc_ptr = benz_obs_ptr
                        ovoc_other = benz_obs_toga
                        gvoc_gfas = benz_gc_gfas
                        gvoc_gfed4 = benz_gc_gfed4
                        gvoc_qfed  = benz_gc_qfed
                        gvoc_finn  = benz_gc_finn
                        gvoc_threegfas =benz_gc_threegfas
                        gvoc_nobb = benz_gc_nobb

                        title = 'Benzene'
                        xrange = [0,0.6]
                        
                        if keyword_set(filter2_nobb) and keyword_set(filter3_nobb) then xrange = [0, 0.2]
                        koh_voc_obs = tmp_koh_co_obs*scalefactor_obs*1e-15  ; dummy
                        koh_voc_mod = tmp_koh_co_mod*scalefactor_mod*1e-15  ; dummy
                    end
                    ;;tolu
                    7:begin
                        ;; Emissions
                        ovoc_ptr = tolu_obs_ptr
                        ovoc_other = tolu_obs_toga

                        gvoc_gfas = c7h8_gc_gfas
                        gvoc_gfed4 = c7h8_gc_gfed4
                        gvoc_qfed  = c7h8_gc_qfed
                        gvoc_finn  = c7h8_gc_finn
                        gvoc_threegfas =c7h8_gc_threegfas
                        gvoc_nobb = c7h8_gc_nobb
                        title = 'Toluene'
                        xrange = [0,0.25]
                        
                        if keyword_set(filter2_nobb) and keyword_set(filter3_nobb) then xrange = [0, 0.1]
                        koh_voc_obs = tmp_koh_co_obs*scalefactor_obs*1e-15  ; dummy
                        koh_voc_mod = tmp_koh_co_mod*scalefactor_mod*1e-15  ; dummy
                    end  
                    ;;xyle
                    8:begin
                        ;; Emissions
                        ;ovoc_ptr = c8h10_obs_ptr
                        ;ovoc_toga = c8h10_obs_toga + etbenzene_obs_toga 

                        ovoc_ptr = c8h10_obs_ptr*0.65
                        ovoc_other = c8h10_obs_toga*0.65
                        
                        gvoc_gfas = c8h10_gc_gfas
                        gvoc_gfed4 = c8h10_gc_gfed4
                        gvoc_qfed  = c8h10_gc_qfed
                        gvoc_finn  = c8h10_gc_finn
                        gvoc_threegfas =c8h10_gc_threegfas
                        gvoc_nobb = c8h10_gc_nobb
                        title = 'Xylene'
                        xrange = [0, 0.05]
                        
                        if keyword_set(filter2_nobb) and keyword_set(filter3_nobb) then xrange = [0, 0.02]
                        koh_voc_obs = tmp_koh_co_obs*scalefactor_obs*1e-15  ; dummy
                        koh_voc_mod = tmp_koh_co_mod*scalefactor_mod*1e-15  ; dummy
                    end
                    ;;Acetic acid
                    ;9:begin
                    ;    ovoc_ptr = acta_obs_ptr
                    ;    ovoc_toga = acta_obs_cims
                    ;    title = 'Acetic acid'
                    ;    xrange = [0,10]
                    ;end
                    ;;RCHO
                    ;10:begin
                    ;    ovoc_ptr = propanal_obs_ptr + butanal_obs_ptr
                    ;    ovoc_toga = propanal_obs_toga + butanal_obs_toga
                    ;    title = 'RCHO'
                    ;    xrange = [0,1.5]
                    ;end
                endcase 
            endif
    ;;set for individual species, e.g., CO, benzene
            if keyword_set(each) then begin
                if keyword_set(ind_co) then begin
                    ovoc_ptr = co_obs
                    ovoc_other = co_obs

                    gvoc_gfed4 = co_gc_gfed4
                    gvoc_finn = co_gc_finn
                    gvoc_gfas = co_gc_gfas
                    gvoc_threegfas = co_gc_threegfas
                    gvoc_qfed = co_gc_qfed
                    gvoc_nobb = co_gc_nobb
                    title = 'CO'
                    xrange = [0,600]
                    if keyword_set(boise) then xrange = [0,300]
                    if keyword_set(CentralValley) and keyword_set(all) then xrange = [0,800]

                    if keyword_set(young) then xrange = [0,1200]
                    if keyword_set(intermed) then xrange = [0,700]
                    if keyword_set(aged) then xrange = [0,700]
                    if keyword_set(nonsmoke) then xrange = [40,120]
                    if keyword_set(normalized) then xrange = [0,2]      
                    
                    if keyword_set(filter1) then xrange = [0,200]
                    if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter4_nobb) or keyword_set(filter5_nobb) and dd eq 0 then xrange = [0,300]
                    if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter4_nobb) or keyword_set(filter5_nobb) and dd ne 0 then xrange = [0,1000]
                    
                    
                    koh_voc_obs = tmp_koh_co_obs*scalefactor_obs*1e-15
                    koh_voc_mod = tmp_koh_co_mod*scalefactor_mod*1e-15
                        
                    if keyword_set(OHrplot) then xrange = [0,10]
                
                endif
                if keyword_set(ind_benz) then begin
                    ovoc_ptr = benz_obs_ptr
                    ovoc_other = benz_obs_toga
                    
                    gvoc_gfed4 = benz_gc_gfed4
                    gvoc_finn = benz_gc_finn
                    gvoc_gfas = benz_gc_gfas
                    gvoc_qfed = benz_gc_qfed
                    gvoc_nobb = benz_gc_nobb
                    title = 'Benzene'
                    xrange = [0,0.6]
                    if keyword_set(boise) then xrange = [0,0.2]
                    if keyword_set(CentralValley) then xrange = [0,1]

                    if keyword_set(young) then xrange = [0,2]
                    if keyword_set(intermed) then xrange = [0,0.8]
                    if keyword_set(aged) then xrange = [0,0.5]
                    ;if keyword_set(LT4) then xrange = [0,0.6]
                    if keyword_set(nonsmoke) then xrange = [0,0.1]
                    if keyword_set(normalized) then begin
                        xrange = [0,2]
                        if keyword_set(aged) then xrange = [0,3]
                    endif
                endif        
                if keyword_set(ind_o3) then begin
                    ovoc_ptr = o3_obs
                    ovoc_other = o3_obs
                    
                    gvoc_gfed4 = o3_gc_gfed4
                    gvoc_finn = o3_gc_finn
                    gvoc_gfas = o3_gc_gfas
                    gvoc_qfed = o3_gc_qfed
                    gvoc_nobb = o3_gc_nobb
                    title = 'O3'
                    xrange = [0,100]
                    if keyword_set(boise) then xrange = [0,100]
                    if keyword_set(CentralValley) then xrange = [0,100]

                    if keyword_set(young) then xrange = [0,100]
                    if keyword_set(intermed) then xrange = [0,100]
                    if keyword_set(aged) then xrange = [0,100]
                    ;if keyword_set(LT4) then xrange = [0,0.6]
                    if keyword_set(nonsmoke) then xrange = [0,100]
                    if keyword_set(normalized) then begin
                        xrange = [0,100]
                        if keyword_set(aged) then xrange = [0,100]
                    endif
                endif  
            endif
            
    ;;all
            if keyword_set(all) then begin
                case s of                     
                    ;;c2h6(AWAS, discrete data)
                    ;;hcho
                    0:begin
                        ;; Emissions
                        ovoc_ptr = ch2o_obs_ptr 
                        ovoc_other = ch2o_obs_toga 
                        
                        gvoc_gfas = hcho_gc_gfas
                        gvoc_gfed4 = hcho_gc_gfed4
                        gvoc_qfed  = hcho_gc_qfed
                        gvoc_finn  = hcho_gc_finn
                        gvoc_threegfas = hcho_gc_threegfas
                        gvoc_nobb = hcho_gc_nobb
                        
                        ;; OH reactivity
                        ;; Units
                        ;; scalefactor is used to convert ppm to molec/m3, 1e-6 from m3, 1e-9 from ppb
                        koh_voc_obs = tmp_koh_ch2o_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_ch2o_mod*scalefactor_mod*1e-15
                        
                        
                        title = 'Formaldehyde'
                        
                        xrange = [0,6]
                        
                        if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) then xrange = [0,2]    
                        
                        if keyword_set(OHrplot) then xrange = [0,1.5]
    
                    end
                    ;;ald2
                    1:begin
                        ;; Emissions
                        ovoc_ptr = ald2_obs_ptr
                        ovoc_other = ald2_obs_toga
                        
                        gvoc_gfas = ald2_gc_gfas
                        gvoc_gfed4= ald2_gc_gfed4
                        gvoc_qfed = ald2_gc_qfed
                        gvoc_finn = ald2_gc_finn
                        gvoc_threegfas = ald2_gc_threegfas
                        gvoc_nobb = ald2_gc_nobb
                        
                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_ald2_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_ald2_mod*scalefactor_mod*1e-15
                        
                        title = 'Acetaldehyde'
                        xrange = [0,3]

                        if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) then xrange = [0,0.6]    
                        if keyword_set(OHrplot) then xrange = [0,0.8]
                    end
                    ;;acet
                    2:begin
                        ;; Emissions
                        ovoc_ptr = acet_obs_ptr 
                        ovoc_other = acet_obs_toga
                        
                        gvoc_gfas = acet_gc_gfas 
                        gvoc_gfed4 = acet_gc_gfed4
                        gvoc_qfed = acet_gc_qfed
                        gvoc_finn = acet_gc_finn
                        gvoc_threegfas =  acet_gc_threegfas

                        gvoc_nobb = acet_gc_nobb
                        
                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_acet_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_acet_mod*scalefactor_mod*1e-15
                        title = 'Acetone'
                        xrange = [0,4]


                        if keyword_set(filter2_nobb) or keyword_set(filter3_nobb)  or keyword_set(filter5_nobb) then xrange = [0,2]  
                        IF keyword_set(OHrplot) then xrange = [0,0.025]
                    end
                    
                    ;;c3h8
                    3:begin
                        ;; Emissions
                        ovoc_ptr = propane_obs_toga
                        ovoc_other = propane_obs_toga
                        
                        gvoc_gfas = c3h8_gc_gfas
                        gvoc_gfed4 = c3h8_gc_gfed4
                        gvoc_qfed = c3h8_gc_qfed
                        gvoc_finn = c3h8_gc_finn
                        gvoc_threegfas = c3h8_gc_threegfas

                        gvoc_nobb = c3h8_gc_nobb
                    
                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_c3h8_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_c3h8_mod*scalefactor_mod*1e-15
                        
                        title = 'Propane'
                        xrange = [0,1]

                        yrange = [0,7]
                        
                        
                        if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) then xrange = [0,0.5]    
                        
                        if keyword_set(OHrplot) then xrange = [0,0.012]
                    end

                    ;;mek
                    4:begin
                        ;; Emissions
                        ovoc_ptr = mek_obs_ptr
                        ovoc_other = mek_obs_toga
                        
                        gvoc_gfas = mek_gc_gfas
                        gvoc_gfed4 = mek_gc_gfed4
                        gvoc_qfed  = mek_gc_qfed
                        gvoc_finn  = mek_gc_finn
                        gvoc_threegfas = mek_gc_threegfas
                        gvoc_nobb = mek_gc_nobb
                        
                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_mek_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_mek_mod*scalefactor_mod*1e-15
                        
                        title = 'MEK'
                        
                        xrange = [0,0.4]

                        if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) then xrange = [0,0.2]   
                        
                        if keyword_set(OHrplot) then xrange = [0,0.015]

                    end
                    
                    ;;eoh
                    5:begin
                        ovoc_ptr = benz_obs_ptr
                        ovoc_other = benz_obs_toga
                        
                        gvoc_gfas = benz_gc_gfas 
                        gvoc_gfed4 = benz_gc_gfed4
                        gvoc_qfed = benz_gc_qfed
                        gvoc_finn = benz_gc_finn
                        gvoc_threegfas = benz_gc_threegfas
                        gvoc_nobb = benz_gc_nobb
                        
                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_benz_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_benz_mod*scalefactor_mod*1e-15
                        
                        title = 'Benzene'

                        xrange = [0,0.6]
                        
                        if keyword_set(OHrplot) then xrange = [0,15]

                    end
                    ;;benz
                    6:begin
                        ;; Emissions
                        ;; test
                        ovoc_ptr = benz_obs_ptr
                        ovoc_other = benz_obs_toga
                        
                        gvoc_gfas = benz_gc_gfas 
                        gvoc_gfed4 = benz_gc_gfed4
                        gvoc_qfed = benz_gc_qfed
                        gvoc_finn = benz_gc_finn
                        gvoc_threegfas = benz_gc_threegfas
                        gvoc_nobb = benz_gc_nobb
                        
                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_benz_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_benz_mod*scalefactor_mod*1e-15
                        
                        title = 'Benzene'
                        
                        xrange = [0,0.6]

                        if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) then xrange = [0,0.08]    
                        if keyword_set(OHrplot) then xrange = [0,0.015]

                    end
                    ;;tolu
                    7:begin
                        ;; Emissions
                        ovoc_ptr = tolu_obs_ptr
                        ovoc_other = tolu_obs_toga
                        
                        gvoc_gfas = c7h8_gc_gfas 
                        gvoc_gfed4 = c7h8_gc_gfed4
                        gvoc_qfed = c7h8_gc_qfed
                        gvoc_finn = c7h8_gc_finn
                        gvoc_threegfas = c7h8_gc_threegfas
                        gvoc_nobb = c7h8_gc_nobb
                        
                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_tolu_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_tolu_mod*scalefactor_mod*1e-15
                        
                        title = 'Toluene'
                        xrange = [0,0.25]

                        if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) then xrange = [0,0.08]                            
                        if keyword_set(OHrplFot) then xrange = [0,0.03]


                    end  
                    ;;xyle
                    8:begin
                        ;; Emissions
                        ovoc_ptr = Xylenes_obs_ptr ;; PTR would be different 
                        ovoc_other = Xylenes_obs_toga ;; PTR would be different 
                        
                        gvoc_gfas = c8h10_gc_gfas
                        gvoc_gfed4 = c8h10_gc_gfed4
                        gvoc_qfed  = c8h10_gc_qfed
                        gvoc_finn  = c8h10_gc_finn
                        gvoc_threegfas = c8h10_gc_threegfas
                        gvoc_nobb = c8h10_gc_nobb
                        
                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_xyle_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_xyle_mod*scalefactor_mod*1e-15
                        
                        title = 'Xylene'
                        xrange = [0, 0.025] 
                        
                        if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) then xrange = [0,0.002]  

                        if keyword_set(OHrplot) then xrange = [0, 0.03]
                    end

                    ;;CO
                    9:begin
                        ovoc_ptr = co_obs
                        ovoc_other = co_obs
                        
                        gvoc_gfed4 = co_gc_gfed4
                        gvoc_finn = co_gc_finn
                        gvoc_gfas = co_gc_gfas
                        gvoc_qfed = co_gc_qfed
                        gvoc_threegfas = co_gc_threegfas

                        gvoc_nobb = co_gc_nobb
                        
                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_co_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_co_mod*scalefactor_mod*1e-15
            
                        title = 'CO'
                        xrange = [0,600]
                        if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) then xrange = [0,300]  

                        if keyword_set(OHrplot) then xrange = [0,15]
                        
                    end
                    
                    ;;ethane
                    10:begin
                        ovoc_ptr = c2h6_obs_awas
                        ovoc_other = c2h6_obs_awas
                        
                        gvoc_gfed4 = c2h6_gc_gfed4
                        gvoc_finn = c2h6_gc_finn
                        gvoc_gfas = c2h6_gc_gfas
                        gvoc_qfed = c2h6_gc_qfed
                        gvoc_threegfas = c2h6_gc_threegfas

                        gvoc_nobb = c2h6_gc_nobb
                        
                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_c2h6_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_c2h6_mod*scalefactor_mod*1e-15
            
                        title = 'C2H6'
                        xrange = [0, 8]
                        if keyword_set(OHrplot) then xrange = [0,0.1]
                    end

                    ;;ALK4
                    11:begin
                        ovoc_ptr = alk4_obs_awas
                        ovoc_other = alk4_obs_awas
                        
                        gvoc_gfed4 = alk4_gc_gfed4
                        gvoc_finn = alk4_gc_finn
                        gvoc_gfas = alk4_gc_gfas

                        gvoc_qfed = alk4_gc_qfed
                        gvoc_threegfas = alk4_gc_threegfas

                        gvoc_nobb = alk4_gc_nobb
                        
                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_alk4_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_alk4_mod*scalefactor_mod*1e-15
            
                        title = 'ALK4'
                        if keyword_set(OHrplot) then xrange = [0,1.5]
                    end
                    
                    ;;RCHO
                    12:begin
                        ovoc_ptr = propanal_obs_toga + butanal_obs_toga
                        ovoc_other = propanal_obs_toga + butanal_obs_toga
                        
                        gvoc_gfed4 = rcho_gc_gfed4
                        gvoc_finn = rcho_gc_finn
                        gvoc_gfas = rcho_gc_gfas 
                        gvoc_qfed = rcho_gc_qfed
                        gvoc_threegfas = rcho_gc_threegfas
                        gvoc_nobb = rcho_gc_nobb
                        
                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_rcho_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_rcho_mod*scalefactor_mod*1e-15
                        
                        title = 'RCHO'
                        ;xrange = [0,2]
                        ;xrange = [0,0.4]
                        xrange = [0,0.2]

                        if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) then xrange = [0,0.10]  

                        if keyword_set(OHrplot) then xrange = [0,0.15]
                    end
            
                    ;;PRPE
                    13:begin
                        ovoc_ptr = prpe_obs_awas
                        ovoc_other = prpe_obs_awas
                        
                        gvoc_gfed4 = prpe_gc_gfed4
                        gvoc_finn = prpe_gc_finn
                        gvoc_gfas = prpe_gc_gfas
                        gvoc_qfed = prpe_gc_qfed
                        gvoc_threegfas = prpe_gc_threegfas
                        gvoc_nobb = prpe_gc_nobb
                        
                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_prpe_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_prpe_mod*scalefactor_mod*1e-15
            
                        title = 'PRPE'
                        
                        xrange = [0,5]
                        if keyword_set(OHrplot) then xrange = [0,1.5]
                    end
                    ;;PRPE
                    ;; could be too low since only one trace gas is implemented 
                    14:begin
                        ovoc_ptr = prpe_obs_toga
                        ovoc_other = prpe_obs_toga
                        
                        gvoc_gfed4 = prpe_gc_gfed4
                        gvoc_finn = prpe_gc_finn
                        gvoc_gfas = prpe_gc_gfas
                        gvoc_qfed = prpe_gc_qfed
                        gvoc_threegfas = prpe_gc_threegfas
                        gvoc_nobb = prpe_gc_nobb
                        
                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_prpe_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_prpe_mod*scalefactor_mod*1e-15
            
                        title = 'PRPE'
                        xrange = [0,800]
                        if keyword_set(OHrplot) then xrange = [0,0.4]
                    end
                    ;;Formic acid
                    15:begin
                        ovoc_ptr = hcooh_obs_cims
                        ovoc_other = hcooh_obs_ptr
                        
                        gvoc_gfed4 = hcooh_gc_gfed4
                        gvoc_finn = hcooh_gc_finn
                        gvoc_gfas = hcooh_gc_gfas
                        gvoc_qfed = hcooh_gc_qfed
                        gvoc_threegfas = hcooh_gc_threegfas
                        gvoc_nobb = hcooh_gc_nobb
                        
                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_hcooh_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_hcooh_mod*scalefactor_mod*1e-15
            
                        title = 'HCOOH'
                        xrange = [0,15]
                        if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) then xrange = [0,1]  

                        
                        if keyword_set(OHrplot) then xrange = [0,0.15]
                        
                    end
                    ;;Acetic acid
                    16:begin
                        ovoc_ptr = acta_obs_ptr
                        ovoc_other = acta_obs_cims
                        
                        gvoc_gfed4 = acta_gc_gfed4
                        gvoc_finn = acta_gc_finn
                        gvoc_gfas = acta_gc_gfas
                        gvoc_qfed = acta_gc_qfed
                        gvoc_threegfas = acta_gc_threegfas
                        gvoc_nobb = acta_gc_nobb
                        
                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_acta_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_acta_mod*scalefactor_mod*1e-15
            
                        title = 'Acetic'
                        xrange = [0,9]
                        if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) then xrange = [0,0.5]  

                        if keyword_set(OHrplot) then xrange = [0,0.15]

                    end
                    
                    ;;RCHO
                    17:begin
                        ovoc_ptr = propanal_obs_toga + butanal_obs_toga
                        ovoc_other = propanal_obs_toga + butanal_obs_toga
                        
                        gvoc_gfed4 = rcho_gc_gfed4
                        gvoc_finn = rcho_gc_finn
                        gvoc_gfas = rcho_gc_gfas 
                        gvoc_qfed = rcho_gc_qfed
                        gvoc_threegfas = rcho_gc_threegfas
                        gvoc_nobb = rcho_gc_nobb
                        
                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_rcho_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_rcho_mod*scalefactor_mod*1e-15
                        
                        title = 'RCHO'
                        ;xrange = [0,2]
                        ;xrange = [0,0.4]
                        xrange = [0,0.2]

                        if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) then xrange = [0,0.10]  

                        if keyword_set(OHrplot) then xrange = [0,0.15]

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
            zvar_gc_gfas   = alt_gc_gfas
            zvar_gc_gfed4   = alt_gc_gfed4
            zvar_gc_qfed   = alt_gc_qfed
            zvar_gc_finn   = alt_gc_finn
            zvar_gc_threegfas   = alt_gc_threegfas
            zvar_gc_nobb  = alt_gc_nobb
            
    ;test,lxu
    if keyword_set(prs) then begin
            zvar_obs  = prs_obs
            zvar_gc_gfas   = prs_gc_gfas
            zvar_gc_gfed4   = prs_gc_gfed4
            zvar_gc_qfed   = prs_gc_qfed
            zvar_gc_finn   = prs_gc_finn  
            zvar_gc_threegfas   = prs_gc_threegfas
            zvar_gc_nobb   = prs_gc_nobb
    endif

    ;for time series, reset the time array every loop (diff)
            xtime_obs  = doy_obs
            xtime_gc   = doy_gc_gfas

    ;for voc .vs. co
            oco_tmp=co_obs
            gco_gfas=co_gc_gfas
            gco_gfed4=co_gc_gfed4
            gco_qfed=co_gc_qfed     
            gco_finn=co_gc_finn
            gco_threegfas=co_gc_threegfas
            gco_nobb=co_gc_nobb

    ;;add for OVOC VS. ISOP
            oisop_ptr = isop_obs_ptr
            oisop_toga = isop_obs_toga
            gisop_gfas = c5h8_gc_gfas
            gisop_gfed4 = c5h8_gc_gfed4
            gisop_qfed = c5h8_gc_qfed
            gisop_finn = c5h8_gc_finn
            gisop_threegfas = c5h8_gc_threegfas
            gisop_nobb = c5h8_gc_nobb

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
            lat_gc = lat_gc_gfas
            lon_gc = lon_gc_gfas
            temperature_gc = temperature_gc_gfas
    ;;setting tmp value
            tmp_isop_ptr = oisop_ptr
            tmp_isop_toga = oisop_toga
            tmp_lat_obs = lat_obs
            tmp_lon_obs = lon_obs
            tmp_lat_gc  = lat_gc
            tmp_lon_gc  = lon_gc
            tmp_temperature_gc = temperature_gc
            ;cloud
            tmp_rhum = rhum
            ;terperature
            tmp_temperature_obs=temperature_obs

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
    ;;all
            ;i = WHERE(ovoc_ptr LE 0, count)
            ;if count gt 0 then ovoc_ptr[i] = !VALUES.F_NAN    
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
                ;ovoc_ptr = ovoc_ptr[remove_ind]
                ;ovoc_toga = ovoc_toga[remove_ind] 
                ovoc_ptr = ovoc_ptr[remove_ind]
                ovoc_other = ovoc_other[remove_ind]

                gvoc_gfas = gvoc_gfas[remove_ind]
                gvoc_qfed = gvoc_qfed[remove_ind]
                gvoc_gfed4 = gvoc_gfed4[remove_ind]
                gvoc_finn = gvoc_finn[remove_ind]
                gvoc_threegfas = gvoc_threegfas[remove_ind]
                gvoc_nobb = gvoc_nobb[remove_ind]
                
                ;for vertical profile
                zvar_obs=zvar_obs[remove_ind]
                zvar_gc_gfed4=zvar_gc_gfed4[remove_ind]
                zvar_gc_finn = zvar_gc_finn[remove_ind]
                zvar_gc_gfas=zvar_gc_gfas[remove_ind]
                zvar_gc_qfed=zvar_gc_qfed[remove_ind]
                zvar_gc_threegfas=zvar_gc_threegfas[remove_ind]
                zvar_gc_nobb = zvar_gc_nobb[remove_ind]

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
                gco_gfed4=gco_gfed4[remove_ind]
                gco_finn=gco_finn[remove_ind]
                gco_gfas=gco_gfas[remove_ind]
                gco_qfed=gco_qfed[remove_ind]
                gco_threegfas=gco_threegfas[remove_ind]
                gco_nobb=gco_nobb[remove_ind]

                ;for voc vs. isop
                tmp_isop_ptr = tmp_isop_ptr[remove_ind]
                tmp_isop_toga = tmp_isop_toga[remove_ind]
                gisop_gfed4 = gisop_gfed4[remove_ind]
                gisop_finn = gisop_finn[remove_ind]
                gisop_gfas = gisop_gfas[remove_ind]
                gisop_qfed = gisop_qfed[remove_ind]
                gisop_threegfas = gisop_threegfas[remove_ind]
                gisop_nobb = gisop_nobb[remove_ind]

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
                ;;temperature
                tmp_temperature_obs=tmp_temperature_obs[remove_ind]
                tmp_temperature_gc=tmp_temperature_gc[remove_ind]
            endif
            ;;set NA 
            if remove_ind[0] eq -1 then gvoc_gfed4 = !values.f_nan
            if remove_ind[0] eq -1 then gvoc_finn = !values.f_nan
            if remove_ind[0] eq -1 then gvoc_gfas = !values.f_nan
            if remove_ind[0] eq -1 then gvoc_qfed = !values.f_nan
            if remove_ind[0] eq -1 then gvoc_threegfas = !values.f_nan

            if remove_ind[0] eq -1 then gvoc_nobb = !values.f_nan


            if remove_ind[0] eq -1 then ovoc_ptr = !values.f_nan
            if remove_ind[0] eq -1 then ovoc_other = !values.f_nan

;///////////////////////////////////////////////////////////////////////////////////////////
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;++                                      STEP 5                                           ++
;++                                  DO THE PLOTTING                                      ++
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;///////////////////////////////////////////////////////////////////////////////////////////W
    ;; Vertical profiles
    ;; pressure
    ;; just use one altitude
            zvar_gc    = zvar_gc_gfas                              
            if keyword_set(prs) then begin
                    oVert_ptr = bin_vert(Data=ovoc_ptr, ZData = zvar_obs, $
                                        ZEdges=PEdges,/Press)
                                        
                    oVert_other = bin_vert(Data=ovoc_other, ZData = zvar_obs, $
                                        ZEdges=PEdges,/Press)
                                        
                    gVert_gfed4  = bin_vert(Data=gvoc_gfed4, ZData=zvar_gc, $
                                        ZEdges=PEdges,/Press)
                    gVert_finn  = bin_vert(Data=gvoc_finn, ZData=zvar_gc, $
                                        ZEdges=PEdges,/Press)

                    gVert_gfas  = bin_vert(Data=gvoc_gfas, ZData=zvar_gc, $
                                        ZEdges=PEdges,/Press)

                                        
                    gVert_qfed  = bin_vert(Data=gvoc_qfed, ZData=zvar_gc, $
                                        ZEdges=PEdges,/Press)
                    gVert_threegfas = bin_vert(Data=gvoc_threegfas, ZData=zvar_gc, $
                                        ZEdges=PEdges,/Press)
                    gVert_nobb  = bin_vert(Data=gvoc_nobb, ZData=zvar_gc, $
                                        ZEdges=PEdges,/Press)
            endif
            if keyword_set(average) then begin
                obs_ptr = oVert_ptr.DataMean
                obs_other = oVert_other.DataMean
                
                gfed4 = gVert_gfed4.DataMean
                finn  = gVert_finn.DataMean
                gfas  = gVert_gfas.DataMean
                

                qfed  = gVert_qfed.DataMean
                threegfas  = gVert_threegfas.DataMean

                nobb  = gVert_nobb.DataMean
            endif
            if keyword_set(med) then begin
                obs_ptr = oVert_ptr.DataMed
                obs_other = oVert_other.DataMed

                gfed4 = gVert_gfed4.DataMed
                finn  = gVert_finn.DataMed
                gfas  = gVert_gfas.DataMed
                
  
                qfed  = gVert_qfed.DataMed
                nobb  = gVert_nobb.DataMed
                threegfas  = gVert_threegfas.DataMed
            endif
;if s eq 10 then print, obs
            z_obs_ptr = oVert_ptr.zMean
            z_obs_other = oVert_other.zMean
            
            z_gfed4 = gVert_gfed4.zMean
            z_finn = gVert_finn.zMean
            z_gfas = gVert_gfas.zMean
            z_qfed = gVert_qfed.zMean
            z_threegfas = gVert_threegfas.zMean

            z_nobb = gVert_nobb.zMean

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
            
            ;; gfed4
            gfed4_q10 = gVert_gfed4.DataQ10
            gfed4_q90 = gVert_gfed4.DataQ90
            gfed4_q25 = gVert_gfed4.DataQ25
            gfed4_q75 = gVert_gfed4.DataQ75
            gfed4_numpts = gVert_gfed4.NumPts

            ;; finn
            finn_q10 = gVert_finn.DataQ10
            finn_q90 = gVert_finn.DataQ90
            finn_q25 = gVert_finn.DataQ25
            finn_q75 = gVert_finn.DataQ75
            finn_numpts = gVert_finn.NumPts
            
            ;; gfas
            gfas_q10 = gVert_gfas.DataQ10
            gfas_q90 = gVert_gfas.DataQ90
            gfas_q25 = gVert_gfas.DataQ25
            gfas_q75 = gVert_gfas.DataQ75
            gfas_numpts = gVert_gfas.NumPts
            

            
            ;; qfed
            qfed_q10 = gVert_qfed.DataQ10
            qfed_q90 = gVert_qfed.DataQ90
            qfed_q25 = gVert_qfed.DataQ25
            qfed_q75 = gVert_qfed.DataQ75
            qfed_numpts = gVert_qfed.NumPts
            
            ;; threegfas
            threegfas_q10 = gVert_threegfas.DataQ10
            threegfas_q90 = gVert_threegfas.DataQ90
            threegfas_q25 = gVert_threegfas.DataQ25
            threegfas_q75 = gVert_threegfas.DataQ75
            threegfas_numpts = gVert_threegfas.NumPts

            ;; nobb
            nobb_q10 = gVert_nobb.DataQ10
            nobb_q90 = gVert_nobb.DataQ90
            nobb_q25 = gVert_nobb.DataQ25
            nobb_q75 = gVert_nobb.DataQ75   
            nobb_numpts = gVert_nobb.NumPts

            ;; delete data point less than 20 for full and 10 for nobb
            ;; PTR
            if plume eq 0 then begin
                numpt_thresh = 10
            endif else begin
                if keyword_set(each) then numpt_thresh = 6 ; slightly modify the trehshold
                if keyword_set(all) then numpt_thresh = 5
            endelse
            
            ind_num = where(obs_ptr_numpts le numpt_thresh, ct)
            if ct gt 0 then obs_ptr[ind_num] = !VALUES.F_NAN
            
            ;; other
            ind_num = where(obs_other_numpts le numpt_thresh, ct)
            if ct gt 0 then obs_other[ind_num] = !VALUES.F_NAN
            
            ind =  where(finite(obs_ptr),ct)
            if ct gt 0 then begin
                obs_ptr = obs_ptr[ind]
                obs_other = obs_other[ind]

                gfed4 = gfed4[ind]
                finn  = finn[ind]
                gfas  = gfas[ind]
                qfed  = qfed[ind]
                threegfas  = threegfas[ind]
                nobb  = nobb[ind]

                z_obs_ptr = z_obs_ptr[ind]
                z_obs_other = z_obs_other[ind]
                
                z_gfed4 = z_gfed4[ind]
                z_finn = z_finn[ind]
                z_gfas = z_gfas[ind]
                
                
                z_qfed = z_qfed[ind]
                z_threegfas = z_threegfas[ind]

                z_nobb = z_nobb[ind]
                
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

                gfed4_q10 = gfed4_q10[ind]
                gfed4_q90 = gfed4_q90[ind]
                gfed4_q25 = gfed4_q25[ind]
                gfed4_q75 = gfed4_q75[ind]
                gfed4_numpts = gfed4_numpts[ind]
                
                finn_q10 = finn_q10[ind]
                finn_q90 = finn_q90[ind]
                finn_q25 = finn_q25[ind]
                finn_q75 = finn_q75[ind]
                finn_numpts = finn_numpts[ind]


                gfas_q10 = gfas_q10[ind]
                gfas_q90 = gfas_q90[ind]
                gfas_q25 = gfas_q25[ind]
                gfas_q75 = gfas_q75[ind]
                gfas_numpts = gfas_numpts[ind]


                qfed_q10 = qfed_q10[ind]
                qfed_q90 = qfed_q90[ind]
                qfed_q25 = qfed_q25[ind]
                qfed_q75 = qfed_q75[ind]
                qfed_numpts = qfed_numpts[ind]

                threegfas_q10 = threegfas_q10[ind]
                threegfas_q90 = threegfas_q90[ind]
                threegfas_q25 = threegfas_q25[ind]
                threegfas_q75 = threegfas_q75[ind]
                threegfas_numpts = threegfas_numpts[ind]
                
                nobb_q10 = nobb_q10[ind]
                nobb_q90 = nobb_q90[ind]
                nobb_q25 = nobb_q25[ind]
                nobb_q75 = nobb_q75[ind]
                nobb_numpts = nobb_numpts[ind]
            endif
            
            
            if plume eq 0 then yrange = [1000,400]
            ;if plume eq 1 then yrange = [850,400]
            if plume eq 1 then yrange = [1000,400]

            if keyword_set(all) or keyword_set(obs_only) then begin
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

            
            
            if (keyword_set(obs_major) or plume eq 1) and s eq 15  then obs_ptr = obs_other ;; set for HCOOH
            ;if keyword_set(OHrplot) then begin
                if keyword_set(all) or keyword_set(obs_only) or keyword_set(each) then begin
                    if s eq 0 or s eq 3 or s eq 6 or s eq 9 or s eq 12 or s eq 15 then begin                
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
                ;endif
                endif 
                
            ;if keyword_set(OHrplot) then begin
                if keyword_set(chemistry) then begin
                    if s eq 0 or s eq 2 then begin                
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
                ;endif
                endif 
                
                
                ;; other measurements
                if keyword_set(codeployed) and not keyword_set(obs_only) and not keyword_set(obs_major) then begin
                    if (s ne 15) then oplot, obs_other, z_obs_other, col=3,thick=thick1,line=0
                    if (s eq 15) then oplot, obs_other, z_obs_other, col=9,thick=thick1,line=0
                endif
                
                if keyword_set(codeployed) and keyword_set(obs_only) then begin
                    if (s ne 3) then oplot, obs_other, z_obs_other, col=3,thick=thick1,line=0
                    if (s eq 3) then oplot, obs_other, z_obs_other, col=9,thick=thick1,line=0
                endif
                
                if not keyword_set(obs_only) then begin
                    ;; only set sepcies implemented in the model
                    if ((s ne 15) and (s ne 16) and (s ne 17)) then oplot,gfed4,z_gfed4,col=2,thick=thick1,line=0


                    ;if s ne 4 then begin
                        oplot,gfas,z_gfas,col=4,thick=thick1,line=0

                    if (s eq 0 or s eq 1 or s eq 2 or s eq 3 or s eq 4) and not keyword_set(chemistry) then begin
                        oplot,qfed,z_qfed,col=5,thick=thick1,line=0

                    endif

                    if not keyword_set(filter2_nobb) then begin
                        if s ne 4 then begin
                            oplot,threegfas,z_threegfas,col=12,thick=thick1,line=1

                        endif
                    endif

                    if not (keyword_set(filter2_nobb) and keyword_set(filter2_nobb)) then $
                    oplot,nobb,z_nobb,col=7,thick=thick1,line=1

                endif
            
            if keyword_set(all) or keyword_set(obs_only) then begin
                charsize =1.2
                charthick =2.8
            endif 
            if keyword_set(each) then begin
                charsize =2
                charthick =6
            endif 
            ;; change the postition later.
            ;; xyouts for NumPoints for each layers            
            for num = 0, n_elements(obs_ptr_numpts)-1 do begin
                if z_obs_ptr[num] ge yrange[0] then continue            
                xyouts, 1.015*max(xrange),z_obs_ptr[num], STRTRIM(string(obs_ptr_numpts[num]), 1), /data, col=1, charsize=charsize, $
                    charthick=charthick, alignment=0.      
            endfor
     

            if dd eq 0 then begin
                spec_obs_total   = [0, obs_ptr] 
                spec_gfas_total  = [1, gfas]
                spec_pres        = [2, z_obs_ptr]
                print, 'currently saving out data...'
                if plume eq 0 then write_csv, './tmp/WE-CAN_full_' + title + '.csv' ,  spec_pres, spec_obs_total, spec_gfas_total
                if plume eq 1 then write_csv, './tmp/WE-CAN_nobb_' + title + '.csv' ,  spec_pres, spec_obs_total, spec_gfas_total
            endif
            
;;; task holder
        endfor ;; for s
    endfor ;; dd loop
;;; print out data in one csv file
    close_device
    ;DEVICE, /CLOSE
    print,'done!'
end
