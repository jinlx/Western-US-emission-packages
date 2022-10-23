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
 ;; We may need to change file directory after finished the test.

;; v15: add O3, NO2, PAN in individual VOCs, may need to give them a better position later
;; v15_v2: add OVOCS and NMHCs, methonal, formic acid and acetic acids are waiting to be added
;; fix no, no2
;; detailed pedges from sree
;; v16, try ch3cn

;;filter1:filter STE, urban plumes
;;filter2:use 25 percentile of acetonitrile
;;fitler3: use acn/co ratio (2.01ppb/ppm) lt the threshold: nonbb (doesn't work well)

; comp_mod_obs_vpro_firexaq_clean,plume=0,/all,/nested,/save,/prs,/med,/errorbar,/wus,/whole,/obs_major,/test
; comp_mod_obs_vpro_firexaq_clean,plume=0,/each,/ind_co,/nested,/save,/prs,/med,/errorbar,/wus,/whole,/obs_major,/test

; comp_mod_obs_vpro_firexaq_clean,plume=1,/filter2_nobb,/filter3_nobb,/all,/nested,/save,/prs,/med,/errorbar,/wus,/whole,/obs_major,/test
@bin_vert
@GCRateConst_co
@GCRateConst_acet
@GCRateConst_c3h8
@GCRateConst
pro comp_mod_obs_vpro_firexaq_clean,$
        plume=plume,$
        filter2_nobb=filter2_nobb,filter2_bb=filter2_bb,$
        filter3_nobb=filter3_nobb,filter3_bb=filter3_bb,$
        all=all,obs_only = obs_only, each=each,measurements=measurements, $
        ind_co=ind_co,ind_benz=ind_benz,ind_o3=ind_o3,ind_pan=ind_pan,ind_no2=ind_no2, $;; it should be used with 'each' keyword_set
        test=test,save=save, $
        nested=nested,fbf=fbf,$
        prs=prs, $
        med=med,average=average,$
        errorbar=errorbar, $
        wus=wus,sus=sus, $
        whole=whole, $ ; meaningless but just be consistent with WE-CAN script
        codeployed=codeployed, obs_major=obs_major;, $
        
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
    fi = 'emipass_ERs'
    if keyword_set(primary) then fi = fi + '_primary'
    if keyword_set(secondary) then fi = fi + '_secondary'
    if keyword_set(filter1) then fi = fi + '_rm'
    if keyword_set(filter2) then fi = fi + '_keep'
    if keyword_set(aged) then fi = fi + '_aged'
    if keyword_set(fresh) then fi = fi + '_fresh'
    if keyword_set(big) then fi = fi + '_big'
    if keyword_set(small) then fi = fi + '_small'
    if keyword_set(cloud1) then fi = fi +'_rmdp'
    if keyword_set(cloud2) then fi = fi +'_rmflight'
    if keyword_set(n_criteria) then fi = fi +'_n_criteria'
    if keyword_set(nested) then fi = fi + '_nested'
    if keyword_set(fbf) then fi = fi + '_fbf'
  if keyword_set(test) then fi = './test_vp'  
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
    ;dates_all=  ['20190722', '20190724', '20190725', '20190729', '20190730',$
    ;                '20190802', '20190803', '20190806', '20190807', '20190808', $
    ;                '20190812', '20190813', '20190815', '20190816', '20190819', $
    ;                '20190821', '20190823', '20190826', '20190830', '20190831', $ 
    ;                '20190903', '20190905']
    dates_all=  ['20190722', '20190724', '20190725', '20190729', '20190730',$
                    '20190802', '20190803', '20190806', '20190807', '20190808', $
                    '20190812', '20190813', '20190815', '20190816', '20190819', $
                    '20190821', '20190823', '20190826', '20190830', '20190831', $ 
                    '20190903']
    
    avo= 6.022e23  ;; molec/mol
    
    for dd = 0, n_elements(dates_all) do begin
        ;;only consider whole flight
        if dd ne 0 then break
    
;; ==============
;; FIREXAQ DATA
;; ==============
        utc_obs = [0]
        doy_obs = [0]
        lat_obs  = [0]  
        lon_obs  = [0]
        prs_obs = [0]
        alt_obs = [0]
        rh_obs  = [0]
        
        temp_obs = [0]
        
;; using mce odr as the transect time
        transect_mce_odr_obs = [0]
        transect_smoke_age_obs=  [0]
;; CL measurement
        o3_obs_cl = [0]
        no_obs_cl = [0]
        no2_obs_cl = [0]
        noy_obs_cl = [0]
;; Others
        o3_obs_roze = [0]
        no_obs_lif  = [0]
        no2_obs_aces = [0]
        no2_obs_canoe = [0]
        hno2_obs_cims = [0]
        hno2_obs_aces = [0]
        hono_obs_saga = [0]
        hno3_obs_cit =  [0]
;; CIMS
        pan_obs_cims = [0]
        ppn_obs_cims = [0]
        apan_obs_cims= [0]
        pbn_obs_cims = [0]
        n2o5_obs_cims= [0]
;; others
        n2o_obs_lgr = [0]
;; DACOM
        co_obs_dacom = [0]
        co_obs_lgr   = [0]
        ch4_obs_dacom = [0]
        co2_obs_700   = [0]
        c2h6_obs_cams = [0]
        ch2o_obs_cams = [0]
        ch2o_obs_isaf = [0]
        ch2o_obs_uioptr  = [0]
;; WAS
        ocs_obs_was  =  [0]
        dms_obs_was   = [0]
        ;; RNO2, CH3X
        c2h6_obs_was = [0]
        c2h4_obs_was = [0]
        c2h2_obs_was = [0]
        c3h6_obs_was = [0]
        c3h8_obs_was = [0]
        propadiene_obs_was = [0]
        propyne_obs_was = [0]
        iButane_obs_was = [0]
        nButane_obs_was = [0]
        x1Butene_obs_was= [0]
        iButene_obs_was = [0]
        t2Butene_obs_was= [0]
        c2Butene_obs_was= [0]
        x13Butadiene_obs_was = [0]
        x12Butadiene_obs_was = [0]
        x1Buten3yne_obs_was  = [0]
        x13Butadyine_obs_was = [0]
        x1Butyne_obs_was = [0]
        x2Butyne_obs_was = [0]
        iPentane_obs_was = [0]
        nPentane_obs_was = [0]
        Isoprene_obs_was = [0]
        x1Pentene_obs_was = [0]
        t2Pentene_obs_was = [0]
        c2Pentene_obs_was = [0]
        threeMe1Butene_obs_was = [0]
        twoMe1Butene_obs_was   = [0]
        twoMe2Butene_obs_was  = [0]
        x13Pentadienes_obs_was = [0]
        x3Me1PenteneAnd4Me1Pentene_obs_was = [0]
        x1Hexene_obs_was = [0]
        x1Heptene_obs_was = [0]
        x1Octene_obs_was = [0]
        x1Nonene_obs_was = [0]
        x1Decene_obs_was = [0]
        nHexane_obs_was  = [0]
        nHeptane_obs_was = [0]
        nOctane_obs_was  = [0]
        nNonane_obs_was  = [0]
        nDecane_obs_was  = [0]
        nUndecane_obs_was = [0]
        x22Dimebutane_obs_was = [0]
        x23Dimebutane_obs_was = [0]
        x2MePentane_obs_was   = [0]
        x3MePentane_obs_was   = [0]
        x2MeHexane_obs_was    = [0]
        x3MeHexane_obs_was    = [0]
        x23DimePentane_obs_was= [0]
        x224TrimePentane_obs_was = [0]
        x234TrimePentane_obs_was = [0]
        CycPentane_obs_was       = [0]
        MeCycPentane_obs_was     = [0]
        CycHexane_obs_was        = [0]
        MeCycHexane_obs_was      = [0]
        CycPentene_obs_was       = [0]
        Benzene_obs_was = [0]
        Toluene_obs_was = [0]
        EthBenzene_obs_was = [0]
        mpXylene_obs_was = [0]
        oXylene_obs_was = [0]
        Styrene_obs_was = [0]
        EthynylBenzene_obs_was = [0]
        iPropBenzene_obs_was = [0]
        nPropBenzene_obs_was = [0]
        x3EthToluene_obs_was = [0]
        x4EthToluene_obs_was = [0]
        x2EthToluene_obs_was = [0]
        x135rimeBenzene_obs_was = [0]
        x124rimeBenzene_obs_was = [0]
        ClBenzene_obs_was = [0]
        aPinene_obs_was = [0]
        bPinene_obs_was = [0]
        Tricyclene_obs_was = [0]
        Camphene_obs_was = [0]
        Myrcene_obs_was = [0]
        Limonene_obs_was = [0]
        Furan_obs_was = [0]
        x2MeFuran_obs_was = [0]
        x3MeFuran_obs_was = [0]
        BenzFuran_obs_was = [0]
        iButanal_obs_was = [0]
        Butanal_obs_was = [0]
        AcetonePropanal_obs_was = [0]
        Acetone_obs_was = [0]
        Propanal_obs_was = [0]
        MEK_obs_was = [0]
        MAC_obs_was = [0]
        MVK_obs_was = [0]
        Acrolein_obs_was = [0]
        iPropanol_obs_was = [0]
        Nitromethane_obs_was = [0]
        Acrylonitrile_obs_was = [0]
        PropNitrile_obs_was = [0]
        MeAcetate_obs_was = [0]
        
        alk4_obs_was = [0]
        prpe_obs_was = [0]
        
;; TOGA measurement
        ;; CH3X
        DMS_obs_toga = [0]
        Propane_obs_toga = [0]
        iButane_obs_toga = [0]
        nButane_obs_toga = [0]
        iPentane_obs_toga = [0]
        nPentane_obs_toga = [0]
        x2MePentane_obs_toga = [0]
        x3MePentane_obs_toga = [0]
        nHexane_obs_toga = [0]
        x224TrimePentane_obs_toga = [0]
        nHeptane_obs_toga = [0]
        nOctane_obs_toga = [0]
        Propene_obs_toga = [0]
        iButene1Butene_obs_toga = [0]
        Isoprene_obs_toga = [0]
        Tricyclene_obs_toga = [0]
        aPinene_obs_toga = [0]
        Camphene_obs_toga = [0]
        bPineneMyrcene_obs_toga = [0]
        LimoneneD3Carene_obs_toga = [0]
        Benzene_obs_toga = [0]
        Toluene_obs_toga = [0]
        Xylenes_obs_toga = [0]
        C8Aromatics_obs_toga = [0]
        EthBenzene_obs_toga = [0]
        mpXylene_obs_toga = [0]
        oXylene_obs_toga = [0]
        Styrene_obs_toga = [0]
        EthynylBenzene_obs_toga = [0]
        CH2O_obs_toga = [0]
        CH3CHO_obs_toga = [0]
        Propanal_obs_toga = [0]
        Butanal_obs_toga = [0]
        iButanal_obs_toga = [0] 
        Acrolein_obs_toga = [0]
        x2Butenals_obs_toga = [0]
        Acetone_obs_toga = [0]
        MEK_obs_toga = [0]
        CH3OH_obs_toga = [0]
        C2H5OH_obs_toga = [0]
        iPropanol_obs_toga = [0] 
        MBO_obs_toga = [0]
        MAC_obs_toga = [0] 
        MVK_obs_toga = [0]
        MeFormate_obs_toga = [0] 
        MeAcetate_obs_toga = [0]
        Furan_obs_toga = [0]
        x2MeFuran_obs_toga = [0]
        x3MeFuran_obs_toga = [0]
        Furfural_obs_toga = [0]
        HCN_obs_toga = [0]
        CH3CN_obs_toga = [0]
        PropNitrile_obs_toga = [0]
        Acrylonitrile_obs_toga = [0]
        MeAcrylonitrile_obs_toga = [0]
        Pyrrole_obs_toga = [0]
        Nitromethane_obs_toga = [0]
        MeONO2_obs_toga = [0]
        EthONO2_obs_toga = [0]
        iPropONO2_obs_toga = [0]
        x2ButONO2iButONO2_obs_toga = [0]
;; IWAS
        ;; CH3X
        Ethane_obs_iwas  = [0]
        Propane_obs_iwas = [0]
        nButane_obs_iwas = [0]
        iButane_obs_iwas = [0]
        nPentane_obs_iwas= [0]
        iPentane_obs_iwas= [0]
        nHexane_obs_iwas = [0]
        x2MePentane_obs_iwas = [0]
        x3MePentane_obs_iwas = [0]
        x22DiMeButane_obs_iwas = [0]
        x24DiMePentane_obs_iwas= [0]
        nOctane_obs_iwas = [0]
        x224TriMePentane_obs_iwas = [0]
        nNonane_obs_iwas = [0]
        nDecane_obs_iwas = [0]
        MeCycPentane_obs_iwas = [0]
        CycHexane_obs_iwas = [0]
        MeCycHexane_obs_iwas = [0]
        Ethyne_obs_iwas = [0]
        Ethene_obs_iwas = [0]
        Propene_obs_iwas = [0]
        x1Butene_obs_iwas = [0]
        c2Butene_obs_iwas = [0]
        t2butene_obs_iwas = [0]
        iButene_obs_iwas = [0]
        x1Pentene_obs_iwas = [0]
        c2Pentene_obs_iwas = [0]
        t2Pentene_obs_iwas = [0]
        x2Me1Butene_obs_iwas = [0]
        x3Me1Butene_obs_iwas = [0]
        t13Pentadiene_obs_iwas = [0]
        Isoprene_obs_iwas = [0]
        aPinene_obs_iwas = [0]
        Benzene_obs_iwas = [0]
        Toluene_obs_iwas = [0]
        EthBenzene_obs_iwas = [0]
        oXylene_obs_iwas = [0]
        mpXylene_obs_iwas = [0]
        Acetone_obs_iwas = [0]
        MEK_obs_iwas = [0]
        MeFormate_obs_iwas = [0]
        Furan_obs_iwas = [0]
        CH3CN_obs_iwas = [0]
        Acrylonitrile_obs_iwas = [0]
        
        
        alk4_obs_iwas = [0]
        prpe_obs_iwas = [0]
;; PTR
        HCN_obs_noaaptr = [0]
        CH2O_obs_noaaptr = [0]
        CH3OH_obs_noaaptr = [0]
        CH3CN_obs_noaaptr = [0]
        HNCO_obs_noaaptr = [0]
        CH3CHO_obs_noaaptr = [0]
        C2H5OH_obs_noaaptr = [0]
        HCOOH_obs_noaaptr = [0]
        Acrylonitrile_obs_noaaptr = [0]
        Acrolein_obs_noaaptr = [0]
        AcetonePropanal_obs_noaaptr = [0]
        Acetone_obs_noaaptr = [0]
        Propanal_obs_noaaptr = [0]
        GlycolaldehydeCH3COOH_obs_noaaptr = [0]
        CH3NO2_obs_noaaptr = [0]
        DMS_obs_noaaptr = [0]
        C4H5N_obs_noaaptr = [0]
        Furan_obs_noaaptr = [0]
        MVKMAC_obs_noaaptr = [0]
        C4Carbonyls_obs_noaaptr = [0]
        MEK_obs_noaaptr = [0]
        Butanal_obs_noaaptr = [0]
        C3H6O2_obs_noaaptr = [0]
        Benzene_obs_noaaptr = [0]
        x2MeFuranx3MeFuran_obs_noaaptr = [0]
        x2Furanone_obs_noaaptr = [0]
        x23Butanedione_obs_noaaptr = [0]
        Toluene_obs_noaaptr = [0]
        Phenol_obs_noaaptr = [0]
        Furfural_obs_noaaptr = [0]
        DimeFurans_obs_noaaptr = [0]
        MaleicAnhyd_obs_noaaptr = [0]
        BenzNitrile_obs_noaaptr = [0] 
        Styrene_obs_noaaptr = [0] 
        Benzaldehyde_obs_noaaptr = [0]
        Xylenes_obs_noaaptr = [0]
        C8Aromatics_obs_noaaptr = [0]
        C7H8O_obs_noaaptr = [0]
        Catecholx5MeFurfural_obs_noaaptr = [0]
        BenzFuran_obs_noaaptr = [0]
        C9Aromatics_obs_noaaptr = [0]
        C6H4O3_obs_noaaptr = [0]
        Guaiacol_obs_noaaptr = [0]
        Naphthalene_obs_noaaptr = [0]
        Monoterpenes_obs_noaaptr = [0]
        Creosols_obs_noaaptr = [0]
        Syringol_obs_noaaptr = [0]

;; CIMS
        ;; Some unknown VOCs
        HCN_obs_cims = [0]
        HCOOH_obs_cims = [0]
        HNCO_obs_cims = [0]

    ;; empty column for time index
        ind_enter_time = [0]
        ind_exit_time = [0]

;; ==============
;; GEOS-Chem DATA
;; ==============

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
        temp_gc_gfas = [0]

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
        temp_gc_threegfas = [0]

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

        mtpa_gc_threegfas   = [0]
        limo_gc_threegfas   = [0]
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
        temp_gc_nobb = [0]

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

        mtpa_gc_nobb   = [0]
        limo_gc_nobb   = [0]
        mtpo_gc_nobb   = [0]
        
        c2h4_gc_nobb   = [0]
        c2h6_gc_nobb   = [0]
    
    ;;add lumped species
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
    ;; GC file for all data
            str_date=dates[n]
            print,'processing '+str_date

;; ======================
;; FIREXAQ observation
;; ======================
;; -------------------------------------------------------------------------------------
            ;dc8fi = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
            ;    'mrg2sav/R4_merges/wecan-mrg60-c130_merge_'+dates[n]+'_R4.sav'

            ;dc8fi = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
            ;    'mrg2sav/R4_merges_avg/wecan-mrg1m-c130_merge_'+dates[n]+'_R4.sav'
            
            dc8fi = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                'mrg2sav/FIREX-AQ/firexaq_avg/firexaq-mrg5m-dc8_merge_'+dates[n]+'_R1.sav'
           
            restore, dc8fi
;; ==============
;; FIREXAQ DATA
;; ==============
            tmp_utc_obs = dc8.Time_Stop
            tmp_doy_obs = dc8.Day_Of_Year_YANG
            tmp_lat_obs  = dc8.Latitude_YANG
            tmp_lon_obs  = dc8.Longitude_YANG
            tmp_prs_obs = dc8.Static_Pressure_YANG
            tmp_alt_obs = dc8.MSL_GPS_Altitude_YANG
            tmp_rh_obs  = dc8.Relative_Humidity_YANG
            tmp_temp_obs = dc8.Static_Air_Temp_YANG +  273.15;; degree C to Kalvin
            
            tmp_transect_mce_odr_obs = dc8.transect_number_SCHWARZ
            tmp_transect_smoke_age_obs = dc8.transect_smoke_age_SCHWARZ

;; CL measurement
            tmp_o3_obs_cl = dc8.O3_CL_RYERSON
            tmp_no_obs_cl = dc8.NO_CL_RYERSON
            tmp_no2_obs_cl = dc8.NO2_CL_RYERSON
            tmp_noy_obs_cl = dc8.NOy_CL_RYERSON
;; Others
            tmp_o3_obs_roze = dc8.O3_ROZE_HANISCO
            tmp_no_obs_lif  = dc8.NO_LIF_ROLLINS/1e3
            tmp_no2_obs_aces = dc8.NO2_ACES_WOMACK
            tmp_no2_obs_canoe = dc8.NO2_CANOE_STCLAIR/1e3
            tmp_hno2_obs_cims = dc8.HNO2_NOAACIMS_VERES/1e3
            tmp_hno2_obs_aces = dc8.HNO2_ACES_WOMACK
            tmp_hono_obs_saga = dc8.HONO_SAGA_DIBB/1e3
;; remaining to be fixed, illegal string
            ;tmp_hno3_obs_cit =  dc8.HNO3+submicron-NO3_SAGA_DIBB/1e3
;; CIMS
            tmp_pan_obs_cims = dc8.PAN_GTCIMS_HUEY/1e3
            tmp_ppn_obs_cims = dc8.PPN_GTCIMS_HUEY/1e3
            tmp_apan_obs_cims= dc8.APAN_GTCIMS_HUEY/1e3
            tmp_pbn_obs_cims = dc8.PBN_GTCIMS_HUEY/1e3
            tmp_n2o5_obs_cims= dc8.N2O5_NOAACIMS_VERES/1e3
;; others
            tmp_n2o_obs_lgr = dc8.N2O_LGR_ppb_PEISCHL
;; DACOM
            tmp_co_obs_dacom = dc8.CO_DACOM_DISKIN
            tmp_co_obs_lgr   = dc8.CO_LGR_ppb_PEISCHL
            tmp_ch4_obs_dacom = dc8.CH4_DACOM_DISKIN
            tmp_co2_obs_700   = dc8.CO2_7000_ppm_DISKIN
            tmp_c2h6_obs_cams = dc8.C2H6_CAMS_pptv_FRIED/1e3
            tmp_ch2o_obs_cams = dc8.CH2O_CAMS_pptv_FRIED/1e3
            tmp_ch2o_obs_isaf = dc8.CH2O_ISAF_HANISCO/1e3
            tmp_ch2o_obs_uioptr  = dc8.CH2O_UIOPTR_ppbV_WISTHALER
;; WAS
            tmp_ocs_obs_was  =  dc8.OCS_WAS_BLAKE/1e3
            tmp_dms_obs_was   = dc8.DMS_WAS_BLAKE/1e3
            ;; RNO2, CH3X
            tmp_c2h6_obs_was = dc8.Ethane_WAS_BLAKE/1e3
            tmp_c2h4_obs_was = dc8.Ethene_WAS_BLAKE/1e3
            tmp_c2h2_obs_was = dc8.Ethyne_WAS_BLAKE/1e3
            tmp_c3h6_obs_was = dc8.Propene_WAS_BLAKE/1e3
            tmp_c3h8_obs_was = dc8.Propane_WAS_BLAKE/1e3
            tmp_propadiene_obs_was = dc8.Propadiene_WAS_BLAKE/1e3
            tmp_propyne_obs_was = dc8.Propyne_WAS_BLAKE/1e3
            tmp_iButane_obs_was = dc8.iButane_WAS_BLAKE/1e3
            tmp_nButane_obs_was = dc8.nButane_WAS_BLAKE/1e3
            tmp_x1Butene_obs_was= dc8.x1Butene_WAS_BLAKE/1e3
            tmp_iButene_obs_was = dc8.iButene_WAS_BLAKE/1e3
            tmp_t2Butene_obs_was= dc8.t2Butene_WAS_BLAKE/1e3
            tmp_c2Butene_obs_was= dc8.c2Butene_WAS_BLAKE/1e3
            tmp_x13Butadiene_obs_was = dc8.x13Butadiene_WAS_BLAKE/1e3
            tmp_x12Butadiene_obs_was = dc8.x12Butadiene_WAS_BLAKE/1e3
            tmp_x1Buten3yne_obs_was  = dc8.x1Buten3yne_WAS_BLAKE/1e3
            tmp_x13Butadyine_obs_was = dc8.x13Butadyine_WAS_BLAKE/1e3
            tmp_x1Butyne_obs_was = dc8.x1Butyne_WAS_BLAKE/1e3
            tmp_x2Butyne_obs_was = dc8.x2Butyne_WAS_BLAKE/1e3
            tmp_iPentane_obs_was = dc8.iPentane_WAS_BLAKE/1e3
            tmp_nPentane_obs_was = dc8.nPentane_WAS_BLAKE/1e3
            tmp_Isoprene_obs_was = dc8.Isoprene_WAS_BLAKE/1e3
            tmp_x1Pentene_obs_was = dc8.x1Pentene_WAS_BLAKE/1e3
            tmp_t2Pentene_obs_was = dc8.t2Pentene_WAS_BLAKE/1e3
            tmp_c2Pentene_obs_was = dc8.c2Pentene_WAS_BLAKE/1e3
;; test, may need to change it later
            ;tmp_threeMe1Butene_obs_was = dc8.3Me1Butene_WAS_BLAKE/1e3
            ;tmp_twoMe1Butene_obs_was   = dc8.2Me1Butene_WAS_BLAKE/1e3
            ;tmp_twoMe2Butene_obs_was  = dc8.2Me2Butene_WAS_BLAKE/1e3
            tmp_x13Pentadienes_obs_was = dc8.x13Pentadienes_WAS_BLAKE/1e3
            tmp_x3Me1PenteneAnd4Me1Pentene_obs_was = dc8.x3Me1PenteneAnd4Me1Pentene_WAS_BLAKE/1e3
            tmp_x1Hexene_obs_was = dc8.x1Hexene_WAS_BLAKE/1e3
            tmp_x1Heptene_obs_was = dc8.x1Heptene_WAS_BLAKE/1e3
            tmp_x1Octene_obs_was = dc8.x1Octene_WAS_BLAKE/1e3
            tmp_x1Nonene_obs_was = dc8.x1Nonene_WAS_BLAKE/1e3
            tmp_x1Decene_obs_was = dc8.x1Decene_WAS_BLAKE/1e3
            tmp_nHexane_obs_was  = dc8.nHexane_WAS_BLAKE/1e3
            tmp_nHeptane_obs_was = dc8.nHeptane_WAS_BLAKE/1e3
            tmp_nOctane_obs_was  = dc8.nOctane_WAS_BLAKE/1e3
            tmp_nNonane_obs_was  = dc8.nNonane_WAS_BLAKE/1e3
            tmp_nDecane_obs_was  = dc8.nDecane_WAS_BLAKE/1e3
            tmp_nUndecane_obs_was = dc8.nUndecane_WAS_BLAKE/1e3
            tmp_x22Dimebutane_obs_was = dc8.x22Dimebutane_WAS_BLAKE/1e3
            tmp_x23Dimebutane_obs_was = dc8.x23Dimebutane_WAS_BLAKE/1e3
            tmp_x2MePentane_obs_was   = dc8.x2MePentane_WAS_BLAKE/1e3
            tmp_x3MePentane_obs_was   = dc8.x3MePentane_WAS_BLAKE/1e3
            tmp_x2MeHexane_obs_was    = dc8.x2MeHexane_WAS_BLAKE/1e3
            tmp_x3MeHexane_obs_was    = dc8.x3MeHexane_WAS_BLAKE/1e3
            tmp_x23DimePentane_obs_was= dc8.x23DimePentane_BLAKE/1e3
            tmp_x224TrimePentane_obs_was = dc8.x224TrimePentane_WAS_BLAKE/1e3
            tmp_x234TrimePentane_obs_was = dc8.x234TrimePentane_WAS_BLAKE/1e3
            tmp_CycPentane_obs_was       = dc8.CycPentane_WAS_BLAKE/1e3
            tmp_MeCycPentane_obs_was     = dc8.MeCycPentane_WAS_BLAKE/1e3
            tmp_CycHexane_obs_was        = dc8.CycHexane_WAS_BLAKE/1e3
            tmp_MeCycHexane_obs_was      = dc8.MeCycHexane_WAS_BLAKE/1e3
            tmp_CycPentene_obs_was       = dc8.CycPentene_WAS_BLAKE/1e3
            tmp_Benzene_obs_was = dc8.Benzene_WAS_BLAKE/1e3
            tmp_Toluene_obs_was = dc8.Toluene_WAS_BLAKE/1e3
            tmp_EthBenzene_obs_was = dc8.EthBenzene_WAS_BLAKE/1e3
            tmp_mpXylene_obs_was = dc8.mpXylene_WAS_BLAKE/1e3
            tmp_oXylene_obs_was = dc8.oXylene_WAS_BLAKE/1e3
            tmp_Styrene_obs_was = dc8.Styrene_WAS_BLAKE/1e3
            tmp_EthynylBenzene_obs_was = dc8.EthynylBenzene_WAS_BLAKE/1e3
            tmp_iPropBenzene_obs_was = dc8.iPropBenzene_WAS_BLAKE/1e3
            tmp_nPropBenzene_obs_was = dc8.nPropBenzene_WAS_BLAKE/1e3
            tmp_x3EthToluene_obs_was = dc8.x3EthToluene_WAS_BLAKE/1e3
            tmp_x4EthToluene_obs_was = dc8.x4EthToluene_WAS_BLAKE/1e3
            tmp_x2EthToluene_obs_was = dc8.x2EthToluene_WAS_BLAKE/1e3
            tmp_x135rimeBenzene_obs_was = dc8.x135rimeBenzene_WAS_BLAKE/1e3
            tmp_x124rimeBenzene_obs_was = dc8.x124rimeBenzene_WAS_BLAKE/1e3
            tmp_ClBenzene_obs_was = dc8.ClBenzene_WAS_BLAKE/1e3
            tmp_aPinene_obs_was = dc8.aPinene_WAS_BLAKE/1e3
            tmp_bPinene_obs_was = dc8.bPinene_WAS_BLAKE/1e3
            tmp_Tricyclene_obs_was = dc8.Tricyclene_WAS_BLAKE/1e3
            tmp_Camphene_obs_was = dc8.Camphene_WAS_BLAKE/1e3
            tmp_Myrcene_obs_was = dc8.Myrcene_WAS_BLAKE/1e3
            tmp_Limonene_obs_was = dc8.Limonene_WAS_BLAKE/1e3
            tmp_Furan_obs_was = dc8.Furan_WAS_BLAKE/1e3
            tmp_x2MeFuran_obs_was = dc8.x2MeFuran_WAS_BLAKE/1e3
            tmp_x3MeFuran_obs_was = dc8.x3MeFuran_WAS_BLAKE/1e3
            tmp_BenzFuran_obs_was = dc8.BenzFuran_WAS_BLAKE/1e3
            tmp_iButanal_obs_was = dc8.iButanal_WAS_BLAKE/1e3
            tmp_Butanal_obs_was = dc8.Butanal_WAS_BLAKE/1e3
            tmp_AcetonePropanal_obs_was = dc8.AcetonePropanal_WAS_BLAKE/1e3
            tmp_Acetone_obs_was = dc8.AcetonePropanal_WAS_BLAKE/1e3  *0.78
            tmp_Propanal_obs_was = dc8.AcetonePropanal_WAS_BLAKE/1e3 *0.22
            tmp_MEK_obs_was = dc8.MEK_WAS_BLAKE/1e3
            tmp_MAC_obs_was = dc8.MAC_WAS_BLAKE/1e3
            tmp_MVK_obs_was = dc8.MVK_WAS_BLAKE/1e3
            tmp_Acrolein_obs_was = dc8.Acrolein_WAS_BLAKE/1e3
            tmp_iPropanol_obs_was = dc8.iPropanol_WAS_BLAKE/1e3
            tmp_Nitromethane_obs_was = dc8.Nitromethane_WAS_BLAKE/1e3
            tmp_Acrylonitrile_obs_was = dc8.Acrylonitrile_WAS_BLAKE/1e3
            tmp_PropNitrile_obs_was = dc8.PropNitrile_WAS_BLAKE/1e3
            tmp_MeAcetate_obs_was = dc8.MeAcetate_WAS_BLAKE/1e3

;; TOGA measurement
            ;; CH3X
            tmp_DMS_obs_toga = dc8.DMS_TOGA_APEL/1E3
            tmp_Propane_obs_toga = dc8.Propane_TOGA_APEL/1E3

            tmp_iButane_obs_toga = dc8.iButane_TOGA_APEL/1E3
            tmp_nButane_obs_toga = dc8.nButane_TOGA_APEL/1E3
            tmp_iPentane_obs_toga = dc8.iPentane_TOGA_APEL/1E3
            tmp_nPentane_obs_toga = dc8.nPentane_TOGA_APEL/1E3
            tmp_x2MePentane_obs_toga = dc8.x2MePentane_TOGA_APEL/1E3
            tmp_x3MePentane_obs_toga = dc8.x3MePentane_TOGA_APEL/1E3
            tmp_nHexane_obs_toga = dc8.nHexane_TOGA_APEL/1E3
            tmp_x224TrimePentane_obs_toga = dc8.x224TrimePentane_TOGA_APEL/1E3
            tmp_nHeptane_obs_toga = dc8.nHeptane_TOGA_APEL/1E3
            tmp_nOctane_obs_toga = dc8.nOctane_TOGA_APEL/1E3
            tmp_Propene_obs_toga = dc8.Propene_TOGA_APEL/1E3
            tmp_iButene1Butene_obs_toga = dc8.iButene1Butene_TOGA_APEL/1E3
            tmp_Isoprene_obs_toga = dc8.Isoprene_TOGA_APEL/1E3
            tmp_Tricyclene_obs_toga = dc8.Tricyclene_TOGA_APEL/1E3
            tmp_aPinene_obs_toga = dc8.aPinene_TOGA_APEL/1E3
            tmp_Camphene_obs_toga = dc8.Camphene_TOGA_APEL/1E3
            tmp_bPineneMyrcene_obs_toga = dc8.bPineneMyrcene_TOGA_APEL/1E3
            tmp_LimoneneD3Carene_obs_toga = dc8.LimoneneD3Carene_TOGA_APEL/1E3
            tmp_Benzene_obs_toga = dc8.Benzene_TOGA_APEL/1E3
            tmp_Toluene_obs_toga = dc8.Toluene_TOGA_APEL/1E3
            tmp_Xylenes_obs_toga = dc8.mpXylene_TOGA_APEL/1E3 + dc8.oXylene_TOGA_APEL/1E3
            tmp_C8Aromatics_obs_toga = dc8.mpXylene_TOGA_APEL/1E3 + dc8.oXylene_TOGA_APEL/1E3 + dc8.EthBenzene_TOGA_APEL/1E3 
            tmp_EthBenzene_obs_toga = dc8.EthBenzene_TOGA_APEL/1E3
            tmp_mpXylene_obs_toga = dc8.mpXylene_TOGA_APEL/1E3
            tmp_oXylene_obs_toga = dc8.oXylene_TOGA_APEL/1E3
            tmp_Styrene_obs_toga = dc8.Styrene_TOGA_APEL/1E3
            tmp_EthynylBenzene_obs_toga = dc8.EthynylBenzene_TOGA_APEL/1E3
            tmp_CH2O_obs_toga = dc8.CH2O_TOGA_APEL/1E3
            tmp_CH3CHO_obs_toga = dc8.CH3CHO_TOGA_APEL/1E3
            tmp_Propanal_obs_toga = dc8.Propanal_TOGA_APEL/1E3
            tmp_Butanal_obs_toga = dc8.Butanal_TOGA_APEL/1E3
            tmp_iButanal_obs_toga = dc8. iButanal_TOGA_APEL/1E3
            tmp_Acrolein_obs_toga = dc8.Acrolein_TOGA_APEL/1E3
            tmp_x2Butenals_obs_toga = dc8.x2Butenals_TOGA_APEL/1E3
            tmp_Acetone_obs_toga = dc8.Acetone_TOGA_APEL/1E3
            tmp_MEK_obs_toga = dc8.MEK_TOGA_APEL/1E3
            tmp_CH3OH_obs_toga = dc8.CH3OH_TOGA_APEL/1E3
            tmp_C2H5OH_obs_toga = dc8.C2H5OH_TOGA_APEL/1E3
            tmp_iPropanol_obs_toga = dc8.iPropanol_TOGA_APEL/1E3
            tmp_MBO_obs_toga = dc8.MBO_TOGA_APEL/1E3
            tmp_MAC_obs_toga = dc8.MAC_TOGA_APEL/1E3
            tmp_MVK_obs_toga = dc8.MVK_TOGA_APEL/1E3
            tmp_MeFormate_obs_toga = dc8.MeFormate_TOGA_APEL/1E3
            tmp_MeAcetate_obs_toga = dc8.MeAcetate_TOGA_APEL/1E3
            tmp_Furan_obs_toga = dc8.Furan_TOGA_APEL/1E3
            tmp_x2MeFuran_obs_toga = dc8.x2MeFuran_TOGA_APEL/1E3
            tmp_x3MeFuran_obs_toga = dc8.x3MeFuran_TOGA_APEL/1E3
            tmp_Furfural_obs_toga = dc8.Furfural_TOGA_APEL/1E3
            tmp_HCN_obs_toga = dc8.HCN_TOGA_APEL/1E3
            tmp_CH3CN_obs_toga = dc8.CH3CN_TOGA_APEL/1E3
            tmp_PropNitrile_obs_toga = dc8.PropNitrile_TOGA_APEL/1E3
            tmp_Acrylonitrile_obs_toga = dc8.Acrylonitrile_TOGA_APEL/1E3
            tmp_MeAcrylonitrile_obs_toga = dc8.MeAcrylonitrile_TOGA_APEL/1E3
            tmp_Pyrrole_obs_toga = dc8.Pyrrole_TOGA_APEL/1E3
            tmp_Nitromethane_obs_toga = dc8.Nitromethane_TOGA_APEL/1E3
            tmp_MeONO2_obs_toga = dc8.MeONO2_TOGA_APEL/1E3
            tmp_EthONO2_obs_toga = dc8.EthONO2_TOGA_APEL/1E3
            tmp_iPropONO2_obs_toga = dc8.iPropONO2_TOGA_APEL/1E3
            tmp_x2ButONO2iButONO2_obs_toga = dc8.x2ButONO2iButONO2_TOGA_APEL/1E3

;; IWAS
            ;; CH3X
            tmp_Ethane_obs_iwas  = dc8.Ethane_NOAAiWAS_GILMAN
            tmp_Propane_obs_iwas = dc8.Propane_NOAAiWAS_GILMAN
            tmp_nButane_obs_iwas = dc8.nButane_NOAAiWAS_GILMAN
            tmp_iButane_obs_iwas = dc8.iButane_NOAAiWAS_GILMAN
            tmp_nPentane_obs_iwas= dc8.nPentane_NOAAiWAS_GILMAN
            tmp_iPentane_obs_iwas= dc8.iPentane_NOAAiWAS_GILMAN
            tmp_nHexane_obs_iwas = dc8.nHexane_NOAAiWAS_GILMAN
            tmp_x2MePentane_obs_iwas = dc8.x2MePentane_NOAAiWAS_GILMAN
            tmp_x3MePentane_obs_iwas = dc8.x3MePentane_NOAAiWAS_GILMAN
            tmp_x22DiMeButane_obs_iwas = dc8.x22DiMeButane_NOAAiWAS_GILMAN
            tmp_x24DiMePentane_obs_iwas= dc8.x24DiMePentane_NOAAiWAS_GILMAN
            tmp_nOctane_obs_iwas = dc8.nOctane_NOAAiWAS_GILMAN
            tmp_x224TriMePentane_obs_iwas = dc8.x224TriMePentane_NOAAiWAS_GILMAN
            tmp_nNonane_obs_iwas = dc8.nNonane_NOAAiWAS_GILMAN
            tmp_nDecane_obs_iwas = dc8.nDecane_NOAAiWAS_GILMAN
            tmp_MeCycPentane_obs_iwas = dc8.MeCycPentane_NOAAiWAS_GILMAN
            tmp_CycHexane_obs_iwas = dc8.CycHexane_NOAAiWAS_GILMAN
            tmp_MeCycHexane_obs_iwas = dc8.MeCycHexane_NOAAiWAS_GILMAN
            tmp_Ethyne_obs_iwas = dc8.Ethyne_NOAAiWAS_GILMAN
            tmp_Ethene_obs_iwas = dc8.Ethene_NOAAiWAS_GILMAN
            tmp_Propene_obs_iwas = dc8.Propene_NOAAiWAS_GILMAN
            tmp_x1Butene_obs_iwas = dc8.x1Butene_NOAAiWAS_GILMAN
            tmp_c2Butene_obs_iwas = dc8.c2Butene_NOAAiWAS_GILMAN
            tmp_t2butene_obs_iwas = dc8.t2butene_NOAAiWAS_GILMAN
            tmp_iButene_obs_iwas = dc8.iButene_NOAAiWAS_GILMAN
            tmp_x1Pentene_obs_iwas = dc8.x1Pentene_NOAAiWAS_GILMAN
            tmp_c2Pentene_obs_iwas = dc8.c2Pentene_NOAAiWAS_GILMAN
            tmp_t2Pentene_obs_iwas = dc8.t2Pentene_NOAAiWAS_GILMAN
            tmp_x2Me1Butene_obs_iwas = dc8.x2Me1Butene_NOAAiWAS_GILMAN
            tmp_x3Me1Butene_obs_iwas = dc8.x3Me1Butene_NOAAiWAS_GILMAN
            tmp_t13Pentadiene_obs_iwas = dc8.t13Pentadiene_NOAAiWAS_GILMAN
            tmp_Isoprene_obs_iwas = dc8.Isoprene_NOAAiWAS_GILMAN
            tmp_aPinene_obs_iwas = dc8.aPinene_NOAAiWAS_GILMAN
            tmp_Benzene_obs_iwas = dc8.Benzene_NOAAiWAS_GILMAN
            tmp_Toluene_obs_iwas = dc8.Toluene_NOAAiWAS_GILMAN
            tmp_EthBenzene_obs_iwas = dc8.EthBenzene_NOAAiWAS_GILMAN
            tmp_oXylene_obs_iwas = dc8.oXylene_NOAAiWAS_GILMAN
            tmp_mpXylene_obs_iwas = dc8.mpXylene_NOAAiWAS_GILMAN
            tmp_Acetone_obs_iwas = dc8.Acetone_NOAAiWAS_GILMAN
            tmp_MEK_obs_iwas = dc8.MEK_NOAAiWAS_GILMAN
            tmp_MeFormate_obs_iwas = dc8.MeFormate_NOAAiWAS_GILMAN
            tmp_Furan_obs_iwas = dc8.Furan_NOAAiWAS_GILMAN
            tmp_CH3CN_obs_iwas = dc8.CH3CN_NOAAiWAS_GILMAN
            tmp_Acrylonitrile_obs_iwas = dc8.Acrylonitrile_NOAAiWAS_GILMAN
;; PTR
            tmp_HCN_obs_noaaptr = dc8.HCN_PTR_WARNEKE_WISTHALER
            tmp_CH2O_obs_noaaptr = dc8.CH2O_PTR_WARNEKE_WISTHALER
            tmp_CH3OH_obs_noaaptr = dc8.CH3OH_PTR_WARNEKE_WISTHALER
            tmp_CH3CN_obs_noaaptr = dc8.CH3CN_PTR_WARNEKE_WISTHALER
            tmp_HNCO_obs_noaaptr = dc8.HNCO_PTR_WARNEKE_WISTHALER
            tmp_CH3CHO_obs_noaaptr = dc8.CH3CHO_PTR_WARNEKE_WISTHALER
            tmp_C2H5OH_obs_noaaptr = dc8.C2H5OH_PTR_WARNEKE_WISTHALER
            tmp_HCOOH_obs_noaaptr = dc8.HCOOH_PTR_WARNEKE_WISTHALER
            tmp_Acrylonitrile_obs_noaaptr = dc8.Acrylonitrile_PTR_WARNEKE_WISTHALER
            tmp_Acrolein_obs_noaaptr = dc8.Acrolein_PTR_WARNEKE_WISTHALER
            tmp_AcetonePropanal_obs_noaaptr = dc8.AcetonePropanal_PTR_WARNEKE_WISTHALER
            tmp_Acetone_obs_noaaptr  = dc8.AcetonePropanal_PTR_WARNEKE_WISTHALER*0.78
            tmp_Propanal_obs_noaaptr = dc8.AcetonePropanal_PTR_WARNEKE_WISTHALER*0.22
            tmp_GlycolaldehydeCH3COOH_obs_noaaptr = dc8.Glycolaldehyde_PTR_WARNEKE_WISTHALER
            tmp_CH3NO2_obs_noaaptr = dc8.CH3NO2_PTR_WARNEKE_WISTHALER
            tmp_DMS_obs_noaaptr = dc8.DMS_PTR_WARNEKE_WISTHALER
            tmp_C4H5N_obs_noaaptr = dc8.C4H5N_PTR_WARNEKE_WISTHALER
            tmp_Furan_obs_noaaptr = dc8.Furan_PTR_WARNEKE_WISTHALER
            tmp_MVKMAC_obs_noaaptr = dc8.MVKMAC_PTR_WARNEKE_WISTHALER
            tmp_C4Carbonyls_obs_noaaptr = dc8.C4Carbonyls_PTR_WARNEKE_WISTHALER
            tmp_MEK_obs_noaaptr = dc8.C4Carbonyls_PTR_WARNEKE_WISTHALER * 0.8
            tmp_Butanal_obs_noaaptr = dc8.C4Carbonyls_PTR_WARNEKE_WISTHALER * 0.2
            tmp_C3H6O2_obs_noaaptr = dc8.C3H6O2_PTR_WARNEKE_WISTHALER
            tmp_Benzene_obs_noaaptr = dc8.Benzene_PTR_WARNEKE_WISTHALER
            tmp_x2MeFuranx3MeFuran_obs_noaaptr = dc8.x2MeFuranx3MeFuran_PTR_WARNEKE_WISTHALER
            tmp_x2Furanone_obs_noaaptr = dc8.x2Furanone_PTR_WARNEKE_WISTHALER
            tmp_x23Butanedione_obs_noaaptr = dc8.x23Butanedione_PTR_WARNEKE_WISTHALER
            tmp_Toluene_obs_noaaptr = dc8.Toluene_PTR_WARNEKE_WISTHALER
            tmp_Phenol_obs_noaaptr = dc8.Phenol_PTR_WARNEKE_WISTHALER
            tmp_Furfural_obs_noaaptr = dc8.Furfural_PTR_WARNEKE_WISTHALER
            tmp_DimeFurans_obs_noaaptr = dc8.DimeFurans_PTR_WARNEKE_WISTHALER
            tmp_MaleicAnhyd_obs_noaaptr = dc8.MaleicAnhyd_PTR_WARNEKE_WISTHALER
            tmp_BenzNitrile_obs_noaaptr = dc8.BenzNitrile_PTR_WARNEKE_WISTHALER 
            tmp_Styrene_obs_noaaptr = dc8.Styrene_PTR_WARNEKE_WISTHALER 
            tmp_Benzaldehyde_obs_noaaptr = dc8.Benzaldehyde_PTR_WARNEKE_WISTHALER
            tmp_Xylenes_obs_noaaptr = dc8.C8Aromatics_PTR_WARNEKE_WISTHALER*0.65
            tmp_C8Aromatics_obs_noaaptr = dc8.C8Aromatics_PTR_WARNEKE_WISTHALER
            tmp_C7H8O_obs_noaaptr = dc8.C7H8O_PTR_WARNEKE_WISTHALER
            tmp_Catecholx5MeFurfural_obs_noaaptr = dc8.Catecholx5MeFurfural_PTR_WARNEKE_WISTHALER
            tmp_BenzFuran_obs_noaaptr = dc8.BenzFuran_PTR_WARNEKE_WISTHALER
            tmp_C9Aromatics_obs_noaaptr = dc8.C9Aromatics_PTR_WARNEKE_WISTHALER
            tmp_C6H4O3_obs_noaaptr = dc8.C6H4O3_PTR_WARNEKE_WISTHALER
            tmp_Guaiacol_obs_noaaptr = dc8.Guaiacol_PTR_WARNEKE_WISTHALER
            tmp_Naphthalene_obs_noaaptr = dc8.Naphthalene_PTR_WARNEKE_WISTHALER
            tmp_Monoterpenes_obs_noaaptr = dc8.Monoterpenes_PTR_WARNEKE_WISTHALER
            tmp_Creosols_obs_noaaptr = dc8.Cresols_PTR_WARNEKE_WISTHALER
            tmp_Syringol_obs_noaaptr = dc8.Syringol_PTR_WARNEKE_WISTHALER

;; CIMS
            ;; Some unknown VOCs
            tmp_HCN_obs_cims = dc8.HCN_NOAACIMS_VERES/1e3
            tmp_HCOOH_obs_cims = dc8.HCOOH_NOAACIMS_VERES/1e3
            tmp_HNCO_obs_cims = dc8.HNCO_NOAACIMS_VERES/1e3

    ;; Basic settings
            utc_obs = [utc_obs, tmp_utc_obs]
            doy_obs = [doy_obs, tmp_doy_obs]
            lat_obs  = [lat_obs, tmp_lat_obs]  
            lon_obs  = [lon_obs, tmp_lon_obs]
            prs_obs = [prs_obs, tmp_prs_obs]
            alt_obs = [alt_obs, tmp_alt_obs]
            rh_obs  = [rh_obs, tmp_rh_obs]
            temp_obs = [temp_obs, tmp_temp_obs]
            transect_mce_odr_obs = [transect_mce_odr_obs, tmp_transect_mce_odr_obs]
            transect_smoke_age_obs = [transect_smoke_age_obs, tmp_transect_smoke_age_obs]

    ;; CL measurement
            o3_obs_cl = [o3_obs_cl, tmp_o3_obs_cl]
            no_obs_cl = [no_obs_cl, tmp_no_obs_cl]
            no2_obs_cl = [no2_obs_cl, tmp_no2_obs_cl]
            noy_obs_cl = [noy_obs_cl, tmp_noy_obs_cl]
    ;; Others
            o3_obs_roze = [o3_obs_roze, tmp_o3_obs_roze]
            no_obs_lif  = [no_obs_lif, tmp_no_obs_lif]
            no2_obs_aces = [no2_obs_aces, tmp_no2_obs_aces]
            no2_obs_canoe = [no2_obs_canoe, tmp_no2_obs_canoe]
            hno2_obs_cims = [hno2_obs_cims, tmp_hno2_obs_cims]
            hno2_obs_aces = [hno2_obs_aces, tmp_hno2_obs_aces]
            hono_obs_saga = [hono_obs_saga, tmp_hono_obs_saga]
;; illegal string
            ;hno3_obs_cit =  [hno3_obs_cit, tmp_hno3_obs_cit]
    ;; CIMS
            pan_obs_cims = [pan_obs_cims, tmp_pan_obs_cims]
            ppn_obs_cims = [ppn_obs_cims, tmp_ppn_obs_cims]
            apan_obs_cims= [apan_obs_cims, tmp_apan_obs_cims]
            pbn_obs_cims = [pbn_obs_cims, tmp_pbn_obs_cims]
            n2o5_obs_cims= [n2o5_obs_cims, tmp_n2o5_obs_cims]
    ;; others
            n2o_obs_lgr = [n2o_obs_lgr, tmp_n2o_obs_lgr]
    ;; DACOM
            co_obs_dacom = [co_obs_dacom, tmp_co_obs_dacom]
            co_obs_lgr   = [co_obs_lgr, tmp_co_obs_lgr]
            ch4_obs_dacom = [ch4_obs_dacom, tmp_ch4_obs_dacom]
            co2_obs_700   = [co2_obs_700, tmp_co2_obs_700]
            c2h6_obs_cams = [c2h6_obs_cams, tmp_c2h6_obs_cams]
            ch2o_obs_cams = [ch2o_obs_cams, tmp_ch2o_obs_cams]
            ch2o_obs_isaf = [ch2o_obs_isaf, tmp_ch2o_obs_isaf]
            ch2o_obs_uioptr  = [ch2o_obs_uioptr, tmp_ch2o_obs_uioptr]
    ;; WAS
            ocs_obs_was  =  [ocs_obs_was, tmp_ocs_obs_was]
            dms_obs_was   = [dms_obs_was, tmp_dms_obs_was]
            ;; RNO2, CH3X
            c2h6_obs_was = [c2h6_obs_was, tmp_c2h6_obs_was]
            c2h4_obs_was = [c2h4_obs_was, tmp_c2h4_obs_was]
            c2h2_obs_was = [c2h2_obs_was, tmp_c2h2_obs_was]
            c3h6_obs_was = [c3h6_obs_was, tmp_c3h6_obs_was]
            c3h8_obs_was = [c3h8_obs_was, tmp_c3h8_obs_was]
            propadiene_obs_was = [propadiene_obs_was, tmp_propadiene_obs_was]
            propyne_obs_was = [propyne_obs_was, tmp_propyne_obs_was]
            iButane_obs_was = [iButane_obs_was, tmp_iButane_obs_was]
            nButane_obs_was = [nButane_obs_was, tmp_nButane_obs_was]
            x1Butene_obs_was= [x1Butene_obs_was, tmp_x1Butene_obs_was]
            iButene_obs_was = [iButene_obs_was, tmp_iButene_obs_was]
            t2Butene_obs_was= [t2Butene_obs_was, tmp_t2Butene_obs_was]
            c2Butene_obs_was= [c2Butene_obs_was, tmp_c2Butene_obs_was]
            x13Butadiene_obs_was = [x13Butadiene_obs_was, tmp_x13Butadiene_obs_was]
            x12Butadiene_obs_was = [x12Butadiene_obs_was, tmp_x12Butadiene_obs_was]
            x1Buten3yne_obs_was  = [x1Buten3yne_obs_was, tmp_x1Buten3yne_obs_was]
            x13Butadyine_obs_was = [x13Butadyine_obs_was, tmp_x13Butadyine_obs_was]
            x1Butyne_obs_was = [x1Butyne_obs_was, tmp_x1Butyne_obs_was]
            x2Butyne_obs_was = [x2Butyne_obs_was, tmp_x2Butyne_obs_was]
            iPentane_obs_was = [iPentane_obs_was, tmp_iPentane_obs_was]
            nPentane_obs_was = [nPentane_obs_was, tmp_nPentane_obs_was]
            Isoprene_obs_was = [Isoprene_obs_was, tmp_Isoprene_obs_was]
            x1Pentene_obs_was = [x1Pentene_obs_was, tmp_x1Pentene_obs_was]
            t2Pentene_obs_was = [t2Pentene_obs_was, tmp_t2Pentene_obs_was]
            c2Pentene_obs_was = [c2Pentene_obs_was, tmp_c2Pentene_obs_was]
;; illegal string
            ;threeMe1Butene_obs_was = [threeMe1Butene_obs_was, tmp_threeMe1Butene_obs_was]
            ;twoMe1Butene_obs_was   = [twoMe1Butene_obs_was, tmp_twoMe1Butene_obs_was]
            ;twoMe2Butene_obs_was  = [twoMe2Butene_obs_was, tmp_twoMe2Butene_obs_was]
            x13Pentadienes_obs_was = [x13Pentadienes_obs_was, tmp_x13Pentadienes_obs_was]
            x3Me1PenteneAnd4Me1Pentene_obs_was = [x3Me1PenteneAnd4Me1Pentene_obs_was, tmp_x3Me1PenteneAnd4Me1Pentene_obs_was]
            x1Hexene_obs_was = [x1Hexene_obs_was, tmp_x1Hexene_obs_was]
            x1Heptene_obs_was = [x1Heptene_obs_was, tmp_x1Heptene_obs_was]
            x1Octene_obs_was = [x1Octene_obs_was, tmp_x1Octene_obs_was]
            x1Nonene_obs_was = [x1Nonene_obs_was, tmp_x1Nonene_obs_was]
            x1Decene_obs_was = [x1Decene_obs_was, tmp_x1Decene_obs_was]
            nHexane_obs_was  = [nHexane_obs_was, tmp_nHexane_obs_was]
            nHeptane_obs_was = [nHeptane_obs_was, tmp_nHeptane_obs_was]
            nOctane_obs_was  = [nOctane_obs_was, tmp_nOctane_obs_was]
            nNonane_obs_was  = [nNonane_obs_was, tmp_nNonane_obs_was]
            nDecane_obs_was  = [nDecane_obs_was, tmp_nDecane_obs_was]
            nUndecane_obs_was = [nUndecane_obs_was, tmp_nUndecane_obs_was]
            x22Dimebutane_obs_was = [x22Dimebutane_obs_was, tmp_x22Dimebutane_obs_was]
            x23Dimebutane_obs_was = [x23Dimebutane_obs_was, tmp_x23Dimebutane_obs_was]
            x2MePentane_obs_was   = [x2MePentane_obs_was, tmp_x2MePentane_obs_was]
            x3MePentane_obs_was   = [x3MePentane_obs_was, tmp_x3MePentane_obs_was]
            x2MeHexane_obs_was    = [x2MeHexane_obs_was, tmp_x2MeHexane_obs_was]
            x3MeHexane_obs_was    = [x3MeHexane_obs_was, tmp_x3MeHexane_obs_was]
            x23DimePentane_obs_was= [x23DimePentane_obs_was, tmp_x23DimePentane_obs_was]
            x224TrimePentane_obs_was = [x224TrimePentane_obs_was, tmp_x224TrimePentane_obs_was]
            x234TrimePentane_obs_was = [x234TrimePentane_obs_was, tmp_x234TrimePentane_obs_was]
            CycPentane_obs_was       = [CycPentane_obs_was, tmp_CycPentane_obs_was]
            MeCycPentane_obs_was     = [MeCycPentane_obs_was, tmp_MeCycPentane_obs_was]
            CycHexane_obs_was        = [CycHexane_obs_was, tmp_CycHexane_obs_was]
            MeCycHexane_obs_was      = [MeCycHexane_obs_was, tmp_MeCycHexane_obs_was]
            CycPentene_obs_was       = [CycPentene_obs_was, tmp_CycPentene_obs_was]
            Benzene_obs_was = [Benzene_obs_was, tmp_Benzene_obs_was]
            Toluene_obs_was = [Toluene_obs_was, tmp_Toluene_obs_was]
            EthBenzene_obs_was = [EthBenzene_obs_was, tmp_EthBenzene_obs_was]
            mpXylene_obs_was = [mpXylene_obs_was, tmp_mpXylene_obs_was]
            oXylene_obs_was = [oXylene_obs_was,tmp_oXylene_obs_was]
            Styrene_obs_was = [Styrene_obs_was, tmp_Styrene_obs_was]
            EthynylBenzene_obs_was = [EthynylBenzene_obs_was, tmp_EthynylBenzene_obs_was]
            iPropBenzene_obs_was = [iPropBenzene_obs_was, tmp_iPropBenzene_obs_was]
            nPropBenzene_obs_was = [nPropBenzene_obs_was, tmp_nPropBenzene_obs_was]
            x3EthToluene_obs_was = [x3EthToluene_obs_was, tmp_x3EthToluene_obs_was]
            x4EthToluene_obs_was = [x4EthToluene_obs_was, tmp_x4EthToluene_obs_was]
            x2EthToluene_obs_was = [x2EthToluene_obs_was, tmp_x2EthToluene_obs_was]
            x135rimeBenzene_obs_was = [x135rimeBenzene_obs_was, tmp_x135rimeBenzene_obs_was]
            x124rimeBenzene_obs_was = [x124rimeBenzene_obs_was, tmp_x124rimeBenzene_obs_was]
            ClBenzene_obs_was = [ClBenzene_obs_was, tmp_ClBenzene_obs_was]
            aPinene_obs_was = [aPinene_obs_was, tmp_aPinene_obs_was]
            bPinene_obs_was = [bPinene_obs_was, tmp_bPinene_obs_was]
            Tricyclene_obs_was = [Tricyclene_obs_was, tmp_Tricyclene_obs_was]
            Camphene_obs_was = [Camphene_obs_was, tmp_Camphene_obs_was]
            Myrcene_obs_was = [Myrcene_obs_was, tmp_Myrcene_obs_was]
            Limonene_obs_was = [Limonene_obs_was, tmp_Limonene_obs_was]
            Furan_obs_was = [Furan_obs_was, tmp_Furan_obs_was]
            x2MeFuran_obs_was = [x2MeFuran_obs_was, tmp_x2MeFuran_obs_was]
            x3MeFuran_obs_was = [x3MeFuran_obs_was, tmp_x3MeFuran_obs_was]
            BenzFuran_obs_was = [BenzFuran_obs_was, tmp_BenzFuran_obs_was]
            iButanal_obs_was = [iButanal_obs_was, tmp_iButanal_obs_was]
            Butanal_obs_was = [Butanal_obs_was, tmp_Butanal_obs_was]
            AcetonePropanal_obs_was = [AcetonePropanal_obs_was, tmp_AcetonePropanal_obs_was]
            Acetone_obs_was = [Acetone_obs_was, tmp_Acetone_obs_was]
            Propanal_obs_was = [Propanal_obs_was, tmp_Propanal_obs_was]
            MEK_obs_was = [MEK_obs_was, tmp_MEK_obs_was]
            MAC_obs_was = [MAC_obs_was, tmp_MAC_obs_was]
            MVK_obs_was = [MVK_obs_was, tmp_MVK_obs_was]
            Acrolein_obs_was = [Acrolein_obs_was, tmp_Acrolein_obs_was]
            iPropanol_obs_was = [iPropanol_obs_was, tmp_iPropanol_obs_was]
            Nitromethane_obs_was = [Nitromethane_obs_was, tmp_Nitromethane_obs_was]
            Acrylonitrile_obs_was = [Acrylonitrile_obs_was, tmp_Acrylonitrile_obs_was]
            PropNitrile_obs_was = [PropNitrile_obs_was, tmp_PropNitrile_obs_was]
            MeAcetate_obs_was = [MeAcetate_obs_was, tmp_MeAcetate_obs_was]
    ;; TOGA measurement
            ;; CH3X
            DMS_obs_toga = [DMS_obs_toga, tmp_DMS_obs_toga]
            Propane_obs_toga = [Propane_obs_toga, tmp_Propane_obs_toga]
            iButane_obs_toga = [iButane_obs_toga, tmp_iButane_obs_toga]
            nButane_obs_toga = [nButane_obs_toga, tmp_nButane_obs_toga]
            iPentane_obs_toga = [iPentane_obs_toga, tmp_iPentane_obs_toga]
            nPentane_obs_toga = [nPentane_obs_toga, tmp_nPentane_obs_toga]
            x2MePentane_obs_toga = [x2MePentane_obs_toga, tmp_x2MePentane_obs_toga]
            x3MePentane_obs_toga = [x3MePentane_obs_toga, tmp_x3MePentane_obs_toga]
            nHexane_obs_toga = [nHexane_obs_toga, tmp_nHexane_obs_toga]
            x224TrimePentane_obs_toga = [x224TrimePentane_obs_toga, tmp_x224TrimePentane_obs_toga]
            nHeptane_obs_toga = [nHeptane_obs_toga, tmp_nHeptane_obs_toga]
            nOctane_obs_toga = [nOctane_obs_toga, tmp_nOctane_obs_toga]
            Propene_obs_toga = [Propene_obs_toga, tmp_Propene_obs_toga]
            iButene1Butene_obs_toga = [iButene1Butene_obs_toga, tmp_iButene1Butene_obs_toga]
            Isoprene_obs_toga = [Isoprene_obs_toga, tmp_Isoprene_obs_toga]
            Tricyclene_obs_toga = [Tricyclene_obs_toga, tmp_Tricyclene_obs_toga]
            aPinene_obs_toga = [aPinene_obs_toga, tmp_aPinene_obs_toga]
            Camphene_obs_toga = [Camphene_obs_toga, tmp_Camphene_obs_toga]
            bPineneMyrcene_obs_toga = [bPineneMyrcene_obs_toga, tmp_bPineneMyrcene_obs_toga]
            LimoneneD3Carene_obs_toga = [LimoneneD3Carene_obs_toga, tmp_LimoneneD3Carene_obs_toga]
            Benzene_obs_toga = [Benzene_obs_toga, tmp_Benzene_obs_toga]
            Toluene_obs_toga = [Toluene_obs_toga, tmp_Toluene_obs_toga]
            Xylenes_obs_toga = [Xylenes_obs_toga, tmp_Xylenes_obs_toga] 
            C8Aromatics_obs_toga = [C8Aromatics_obs_toga, tmp_C8Aromatics_obs_toga]
            EthBenzene_obs_toga  = [EthBenzene_obs_toga, tmp_EthBenzene_obs_toga]
            mpXylene_obs_toga = [mpXylene_obs_toga, tmp_mpXylene_obs_toga]
            oXylene_obs_toga = [oXylene_obs_toga, tmp_oXylene_obs_toga]
            Styrene_obs_toga = [Styrene_obs_toga, tmp_Styrene_obs_toga]
            EthynylBenzene_obs_toga = [EthynylBenzene_obs_toga, tmp_EthynylBenzene_obs_toga]
            CH2O_obs_toga = [CH2O_obs_toga, tmp_CH2O_obs_toga]
            CH3CHO_obs_toga = [CH3CHO_obs_toga, tmp_CH3CHO_obs_toga]
            Propanal_obs_toga = [Propanal_obs_toga, tmp_Propanal_obs_toga]
            Butanal_obs_toga = [Butanal_obs_toga, tmp_Butanal_obs_toga]
            iButanal_obs_toga = [iButanal_obs_toga, tmp_iButanal_obs_toga] 
            Acrolein_obs_toga = [Acrolein_obs_toga, tmp_Acrolein_obs_toga]
            x2Butenals_obs_toga = [x2Butenals_obs_toga, tmp_x2Butenals_obs_toga]
            Acetone_obs_toga = [Acetone_obs_toga, tmp_Acetone_obs_toga]
            MEK_obs_toga = [MEK_obs_toga, tmp_MEK_obs_toga]
            CH3OH_obs_toga = [CH3OH_obs_toga, tmp_CH3OH_obs_toga]
            C2H5OH_obs_toga = [C2H5OH_obs_toga, tmp_C2H5OH_obs_toga]
            iPropanol_obs_toga = [iPropanol_obs_toga, tmp_iPropanol_obs_toga] 
            MBO_obs_toga = [MBO_obs_toga, tmp_MBO_obs_toga]
            MAC_obs_toga = [MAC_obs_toga, tmp_MAC_obs_toga] 
            MVK_obs_toga = [MVK_obs_toga, tmp_MVK_obs_toga]
            MeFormate_obs_toga = [MeFormate_obs_toga, tmp_MeFormate_obs_toga] 
            MeAcetate_obs_toga = [MeAcetate_obs_toga, tmp_MeAcetate_obs_toga]
            Furan_obs_toga = [Furan_obs_toga, tmp_Furan_obs_toga]
            x2MeFuran_obs_toga = [x2MeFuran_obs_toga, tmp_x2MeFuran_obs_toga]
            x3MeFuran_obs_toga = [x3MeFuran_obs_toga, tmp_x3MeFuran_obs_toga]
            Furfural_obs_toga = [Furfural_obs_toga, tmp_Furfural_obs_toga]
            HCN_obs_toga = [HCN_obs_toga, tmp_HCN_obs_toga]
            CH3CN_obs_toga = [CH3CN_obs_toga, tmp_CH3CN_obs_toga]
            PropNitrile_obs_toga = [PropNitrile_obs_toga, tmp_PropNitrile_obs_toga]
            Acrylonitrile_obs_toga = [Acrylonitrile_obs_toga, tmp_Acrylonitrile_obs_toga]
            MeAcrylonitrile_obs_toga = [MeAcrylonitrile_obs_toga, tmp_MeAcrylonitrile_obs_toga]
            Pyrrole_obs_toga = [Pyrrole_obs_toga, tmp_Pyrrole_obs_toga]
            Nitromethane_obs_toga = [Nitromethane_obs_toga, tmp_Nitromethane_obs_toga]
            MeONO2_obs_toga = [MeONO2_obs_toga, tmp_MeONO2_obs_toga]
            EthONO2_obs_toga = [EthONO2_obs_toga, tmp_EthONO2_obs_toga]
            iPropONO2_obs_toga = [iPropONO2_obs_toga, tmp_iPropONO2_obs_toga]
            x2ButONO2iButONO2_obs_toga = [x2ButONO2iButONO2_obs_toga, tmp_x2ButONO2iButONO2_obs_toga]
    ;; IWAS
            ;; CH3X
            Ethane_obs_iwas  = [Ethane_obs_iwas, tmp_Ethane_obs_iwas]
            Propane_obs_iwas = [Propane_obs_iwas, tmp_Propane_obs_iwas]
            nButane_obs_iwas = [nButane_obs_iwas, tmp_nButane_obs_iwas]
            iButane_obs_iwas = [iButane_obs_iwas, tmp_iButane_obs_iwas]
            nPentane_obs_iwas= [nPentane_obs_iwas, tmp_nPentane_obs_iwas]
            iPentane_obs_iwas= [iPentane_obs_iwas, tmp_iPentane_obs_iwas]
            nHexane_obs_iwas = [nHexane_obs_iwas, tmp_nHexane_obs_iwas]
            x2MePentane_obs_iwas = [x2MePentane_obs_iwas, tmp_x2MePentane_obs_iwas]
            x3MePentane_obs_iwas = [x3MePentane_obs_iwas, tmp_x3MePentane_obs_iwas]
            x22DiMeButane_obs_iwas = [x22DiMeButane_obs_iwas, tmp_x22DiMeButane_obs_iwas]
            x24DiMePentane_obs_iwas= [x24DiMePentane_obs_iwas, tmp_x24DiMePentane_obs_iwas]
            nOctane_obs_iwas = [nOctane_obs_iwas, tmp_nOctane_obs_iwas]
            x224TriMePentane_obs_iwas = [x224TriMePentane_obs_iwas, tmp_x224TriMePentane_obs_iwas]
            nNonane_obs_iwas = [nNonane_obs_iwas, tmp_nNonane_obs_iwas]
            nDecane_obs_iwas = [nDecane_obs_iwas, tmp_nDecane_obs_iwas]
            MeCycPentane_obs_iwas = [MeCycPentane_obs_iwas, tmp_MeCycPentane_obs_iwas]
            CycHexane_obs_iwas = [CycHexane_obs_iwas, tmp_CycHexane_obs_iwas]
            MeCycHexane_obs_iwas = [MeCycHexane_obs_iwas, tmp_MeCycHexane_obs_iwas]
            Ethyne_obs_iwas = [Ethyne_obs_iwas, tmp_Ethyne_obs_iwas]
            Ethene_obs_iwas = [Ethene_obs_iwas, tmp_Ethene_obs_iwas]
            Propene_obs_iwas = [Propene_obs_iwas, tmp_Propene_obs_iwas]
            x1Butene_obs_iwas = [x1Butene_obs_iwas, tmp_x1Butene_obs_iwas]
            c2Butene_obs_iwas = [c2Butene_obs_iwas, tmp_c2Butene_obs_iwas]
            t2butene_obs_iwas = [t2butene_obs_iwas, tmp_t2butene_obs_iwas]
            iButene_obs_iwas = [iButene_obs_iwas, tmp_iButene_obs_iwas]
            x1Pentene_obs_iwas = [x1Pentene_obs_iwas, tmp_x1Pentene_obs_iwas]
            c2Pentene_obs_iwas = [c2Pentene_obs_iwas, tmp_c2Pentene_obs_iwas]
            t2Pentene_obs_iwas = [t2Pentene_obs_iwas, tmp_t2Pentene_obs_iwas]
            x2Me1Butene_obs_iwas = [x2Me1Butene_obs_iwas, tmp_x2Me1Butene_obs_iwas]
            x3Me1Butene_obs_iwas = [x2Me1Butene_obs_iwas, tmp_x2Me1Butene_obs_iwas]
            t13Pentadiene_obs_iwas = [t13Pentadiene_obs_iwas, tmp_t13Pentadiene_obs_iwas]
            Isoprene_obs_iwas = [Isoprene_obs_iwas, tmp_Isoprene_obs_iwas]
            aPinene_obs_iwas = [aPinene_obs_iwas, tmp_aPinene_obs_iwas]
            Benzene_obs_iwas = [Benzene_obs_iwas, tmp_Benzene_obs_iwas]
            Toluene_obs_iwas = [Toluene_obs_iwas, tmp_Toluene_obs_iwas]
            EthBenzene_obs_iwas = [EthBenzene_obs_iwas, tmp_EthBenzene_obs_iwas]
            oXylene_obs_iwas = [oXylene_obs_iwas, tmp_oXylene_obs_iwas]
            mpXylene_obs_iwas = [mpXylene_obs_iwas, tmp_mpXylene_obs_iwas]
            Acetone_obs_iwas = [Acetone_obs_iwas, tmp_Acetone_obs_iwas]
            MEK_obs_iwas = [MEK_obs_iwas, tmp_MEK_obs_iwas]
            MeFormate_obs_iwas = [MeFormate_obs_iwas, tmp_MeFormate_obs_iwas]
            Furan_obs_iwas = [Furan_obs_iwas, tmp_Furan_obs_iwas]
            CH3CN_obs_iwas = [CH3CN_obs_iwas, tmp_CH3CN_obs_iwas]
            Acrylonitrile_obs_iwas = [Acrylonitrile_obs_iwas, tmp_Acrylonitrile_obs_iwas]
    ;; NOAAPTR
            HCN_obs_noaaptr = [HCN_obs_noaaptr, tmp_HCN_obs_noaaptr]
            CH2O_obs_noaaptr = [CH2O_obs_noaaptr, tmp_CH2O_obs_noaaptr]
            CH3OH_obs_noaaptr = [CH3OH_obs_noaaptr, tmp_CH3OH_obs_noaaptr]
            CH3CN_obs_noaaptr = [CH3CN_obs_noaaptr, tmp_CH3CN_obs_noaaptr]
            HNCO_obs_noaaptr = [HNCO_obs_noaaptr, tmp_HNCO_obs_noaaptr]
            CH3CHO_obs_noaaptr = [CH3CHO_obs_noaaptr, tmp_CH3CHO_obs_noaaptr]
            C2H5OH_obs_noaaptr = [C2H5OH_obs_noaaptr, tmp_C2H5OH_obs_noaaptr]
            HCOOH_obs_noaaptr = [HCOOH_obs_noaaptr, tmp_HCOOH_obs_noaaptr]
            Acrylonitrile_obs_noaaptr = [Acrylonitrile_obs_noaaptr, tmp_Acrylonitrile_obs_noaaptr]
            Acrolein_obs_noaaptr = [Acrolein_obs_noaaptr, tmp_Acrolein_obs_noaaptr]
            AcetonePropanal_obs_noaaptr = [AcetonePropanal_obs_noaaptr, tmp_AcetonePropanal_obs_noaaptr]
            Acetone_obs_noaaptr = [Acetone_obs_noaaptr, tmp_Acetone_obs_noaaptr]
            Propanal_obs_noaaptr = [Propanal_obs_noaaptr, tmp_Propanal_obs_noaaptr]
            GlycolaldehydeCH3COOH_obs_noaaptr = [GlycolaldehydeCH3COOH_obs_noaaptr, tmp_GlycolaldehydeCH3COOH_obs_noaaptr]
            CH3NO2_obs_noaaptr = [CH3NO2_obs_noaaptr, tmp_CH3NO2_obs_noaaptr]
            DMS_obs_noaaptr = [DMS_obs_noaaptr, tmp_DMS_obs_noaaptr]
            C4H5N_obs_noaaptr = [C4H5N_obs_noaaptr, tmp_C4H5N_obs_noaaptr]
            Furan_obs_noaaptr = [Furan_obs_noaaptr, tmp_Furan_obs_noaaptr]
            MVKMAC_obs_noaaptr = [MVKMAC_obs_noaaptr, tmp_MVKMAC_obs_noaaptr]
            C4Carbonyls_obs_noaaptr = [C4Carbonyls_obs_noaaptr, tmp_C4Carbonyls_obs_noaaptr]
            MEK_obs_noaaptr = [MEK_obs_noaaptr, tmp_MEK_obs_noaaptr]
            Butanal_obs_noaaptr = [Butanal_obs_noaaptr, tmp_Butanal_obs_noaaptr]
            C3H6O2_obs_noaaptr = [C3H6O2_obs_noaaptr, tmp_C3H6O2_obs_noaaptr]
            Benzene_obs_noaaptr = [Benzene_obs_noaaptr, tmp_Benzene_obs_noaaptr]
            x2MeFuranx3MeFuran_obs_noaaptr = [x2MeFuranx3MeFuran_obs_noaaptr, tmp_x2MeFuranx3MeFuran_obs_noaaptr]
            x2Furanone_obs_noaaptr = [x2Furanone_obs_noaaptr, tmp_x2Furanone_obs_noaaptr]
            x23Butanedione_obs_noaaptr = [x23Butanedione_obs_noaaptr, tmp_x23Butanedione_obs_noaaptr]
            Toluene_obs_noaaptr = [Toluene_obs_noaaptr, tmp_Toluene_obs_noaaptr]
            Phenol_obs_noaaptr = [Phenol_obs_noaaptr, tmp_Phenol_obs_noaaptr]
            Furfural_obs_noaaptr = [Furfural_obs_noaaptr, tmp_Furfural_obs_noaaptr]
            DimeFurans_obs_noaaptr = [DimeFurans_obs_noaaptr, tmp_DimeFurans_obs_noaaptr]
            MaleicAnhyd_obs_noaaptr = [MaleicAnhyd_obs_noaaptr, tmp_MaleicAnhyd_obs_noaaptr]
            BenzNitrile_obs_noaaptr = [BenzNitrile_obs_noaaptr, tmp_BenzNitrile_obs_noaaptr] 
            Styrene_obs_noaaptr = [Styrene_obs_noaaptr, tmp_Styrene_obs_noaaptr] 
            Benzaldehyde_obs_noaaptr = [Benzaldehyde_obs_noaaptr, tmp_Benzaldehyde_obs_noaaptr]
            Xylenes_obs_noaaptr = [Xylenes_obs_noaaptr, tmp_Xylenes_obs_noaaptr]
            C8Aromatics_obs_noaaptr = [C8Aromatics_obs_noaaptr, tmp_C8Aromatics_obs_noaaptr]
            C7H8O_obs_noaaptr = [C7H8O_obs_noaaptr, tmp_C7H8O_obs_noaaptr]
            Catecholx5MeFurfural_obs_noaaptr = [Catecholx5MeFurfural_obs_noaaptr, tmp_Catecholx5MeFurfural_obs_noaaptr]
            BenzFuran_obs_noaaptr = [BenzFuran_obs_noaaptr, tmp_BenzFuran_obs_noaaptr]
            C9Aromatics_obs_noaaptr = [C9Aromatics_obs_noaaptr, tmp_C9Aromatics_obs_noaaptr]
            C6H4O3_obs_noaaptr = [C6H4O3_obs_noaaptr, tmp_C6H4O3_obs_noaaptr]
            Guaiacol_obs_noaaptr = [Guaiacol_obs_noaaptr, tmp_Guaiacol_obs_noaaptr]
            Naphthalene_obs_noaaptr = [Naphthalene_obs_noaaptr, tmp_Naphthalene_obs_noaaptr]
            Monoterpenes_obs_noaaptr = [Monoterpenes_obs_noaaptr, tmp_Monoterpenes_obs_noaaptr]
            Creosols_obs_noaaptr = [Creosols_obs_noaaptr, tmp_Creosols_obs_noaaptr]
            Syringol_obs_noaaptr = [Syringol_obs_noaaptr, tmp_Syringol_obs_noaaptr]

    ;; CIMS
            ;; Some unknown VOCs
            HCN_obs_cims = [HCN_obs_cims, tmp_HCN_obs_cims]
            HCOOH_obs_cims = [HCOOH_obs_cims, tmp_HCOOH_obs_cims]
            HNCO_obs_cims = [HNCO_obs_cims, tmp_HNCO_obs_cims]

            undefine,dc8    
        
;;--------------------------------------------------------------------------------------------------------------
;; ============================
;; GEOS-Chem 0.25x0.3125: GFAS
;; ============================
            if keyword_set(nested) then begin
                ;gcfi_gfas   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                ;    'planelog2sav/output_'+'gfas' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'
                ;gcfi_gfas   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                ;    'planelog2sav/output_'+'gfas_tmp' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'         
                gcfi_gfas   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                    'planelog2sav/firexaq_test' + '/mrg5m_firexaq_dc8_'+dates[n]+'.sav'   
            endif
            if keyword_set(fbf) then begin
                ;gcfi_gfas   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                ;    'planelog2sav/output_'+'gfas_4x5' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'
                gcfi_gfas   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                    'planelog2sav/output_'+'gfas_4x5_opt' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'
            endif
            restore, gcfi_gfas

            
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
            tmp_utc_gc_gfas  = gc.utc
            tmp_doy_gc_gfas  = gc.doy
            tmp_lat_gc_gfas  = gc.lat
            tmp_lon_gc_gfas  = gc.lon
            tmp_alt_gc_gfas  = gc.alt
            tmp_prs_gc_gfas  = gc.pres
            tmp_temp_gc_gfas = gc.temp

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
            
            tmp_c2h6_gc_gfas = gc.c2h6*1e9
        ;;add lumped species
            tmp_alk4_gc_gfas   = gc.alk4*1e9/4.3 ; 4.3C
            tmp_prpe_gc_gfas   = gc.prpe*1e9/3   ; 3C 
            tmp_rcho_gc_gfas   = gc.rcho*1e9
            
    ;; concert hhmm(utc) into ss(utc)      
            ;print,tmp_utc_gc_gfas
            hh = floor(tmp_utc_gc_gfas/100)
            ind = where(hh gt 12, ct)
            if ct gt 0 then hh[ind] = hh[ind]
            ind = where(hh lt 12, ct)
            if ct gt 0 then hh[ind] = hh[ind] + 24
            
            mm = tmp_utc_gc_gfas- floor(tmp_utc_gc_gfas/100)*100
            tmp_utc_gc_gfas = float(hh)*60*60+float(mm)*60
            
            
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
            temp_gc_gfas = [temp_gc_gfas,tmp_temp_gc_gfas]

            
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
            
            c2h6_gc_gfas  = [c2h6_gc_gfas, tmp_c2h6_gc_gfas]
        ;;add lumped species
            alk4_gc_gfas  = [alk4_gc_gfas,tmp_alk4_gc_gfas]
            prpe_gc_gfas  = [prpe_gc_gfas,tmp_prpe_gc_gfas]
            rcho_gc_gfas  = [rcho_gc_gfas,tmp_rcho_gc_gfas]

            undefine,gc
            
;; ============================
;; GEOS-Chem 0.25x0.3125: THREEGFAS
;; ============================
            if keyword_set(nested) then begin
                ;gcfi_threegfas   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                ;    'planelog2sav/output_'+'threegfas' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'
                ;gcfi_threegfas   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                ;    'planelog2sav/output_'+'threegfas_tmp' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'         
                gcfi_threegfas   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                    'planelog2sav/output_3Xgfas_firexaq' + '/mrg5m_firexaq_dc8_'+dates[n]+'.sav'   
            endif
            if keyword_set(fbf) then begin
                ;gcfi_threegfas   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                ;    'planelog2sav/output_'+'threegfas_4x5' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'
                gcfi_threegfas   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                    'planelog2sav/output_'+'threegfas_4x5_opt' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'
            endif
            restore, gcfi_threegfas

            
            tmp_co_gc_threegfas   = gc.co*1e9
            tmp_O3_gc_threegfas   = gc.o3*1e9
            tmp_pan_gc_threegfas  = gc.pan*1e9
            tmp_hcho_gc_threegfas = gc.ch2o*1e9
            tmp_acet_gc_threegfas = gc.acet*1e9/3; 3C acetone in ppb
            tmp_benz_gc_threegfas = gc.benz*1e9/6; 6C benzene in ppb
    ;        ch3oh_gc= [ch3oh_gc, gc.MOH*1e9]
            tmp_ald2_gc_threegfas = gc.ald2*1e9/2; 2C ch3cho in ppb

            tmp_no_gc_threegfas   = gc.no*1e9
            tmp_no2_gc_threegfas  = gc.no2*1e9
            tmp_so2_gc_threegfas  = gc.so2*1e9

            tmp_na = avo * gc.pres * 100./(8.31 * gc.temp) * 1e-6;; air density molec/cm3
            tmp_oh_gc_threegfas  = gc.oh*tmp_na ;; v/v --> molec/cm3 

            tmp_date_gc_threegfas = gc.date
            tmp_utc_gc_threegfas  = gc.utc
            tmp_doy_gc_threegfas  = gc.doy
            tmp_lat_gc_threegfas  = gc.lat
            tmp_lon_gc_threegfas  = gc.lon
            tmp_alt_gc_threegfas  = gc.alt
            tmp_prs_gc_threegfas  = gc.pres
            tmp_temp_gc_threegfas = gc.temp

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
            
            tmp_c2h6_gc_threegfas = gc.c2h6*1e9
        ;;add lumped species
            tmp_alk4_gc_threegfas   = gc.alk4*1e9/4.3 ; 4.3C
            tmp_prpe_gc_threegfas   = gc.prpe*1e9/3   ; 3C 
            tmp_rcho_gc_threegfas   = gc.rcho*1e9

    ;; concert hhmm(utc) into ss(utc)      
            ;print,tmp_utc_gc_threegfas
            hh = floor(tmp_utc_gc_threegfas/100)
            ind = where(hh gt 12, ct)
            if ct gt 0 then hh[ind] = hh[ind]
            ind = where(hh lt 12, ct)
            if ct gt 0 then hh[ind] = hh[ind] + 24
            
            mm = tmp_utc_gc_threegfas- floor(tmp_utc_gc_threegfas/100)*100
            tmp_utc_gc_threegfas = float(hh)*60*60+float(mm)*60
            
            
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
            temp_gc_threegfas = [temp_gc_threegfas,tmp_temp_gc_threegfas]

            
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
            
            c2h6_gc_threegfas  = [c2h6_gc_threegfas, tmp_c2h6_gc_threegfas]
        ;;add lumped species
            alk4_gc_threegfas  = [alk4_gc_threegfas,tmp_alk4_gc_threegfas]
            prpe_gc_threegfas  = [prpe_gc_threegfas,tmp_prpe_gc_threegfas]
            rcho_gc_threegfas  = [rcho_gc_threegfas,tmp_rcho_gc_threegfas]

            undefine,gc
;; ============================
;; GEOS-Chem 0.25x0.3125: NOBB
;; ============================
            if keyword_set(nested) then begin
                ;gcfi_nobb   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                ;    'planelog2sav/output_'+'nobb' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'
                ;gcfi_nobb   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                ;    'planelog2sav/output_'+'nobb_tmp' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'         
                gcfi_nobb   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                    'planelog2sav/output_nobb_firexaq' + '/mrg5m_firexaq_dc8_'+dates[n]+'.sav'   
            endif
            if keyword_set(fbf) then begin
                ;gcfi_nobb   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                ;    'planelog2sav/output_'+'nobb_4x5' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'
                gcfi_nobb   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                    'planelog2sav/output_'+'nobb_4x5_opt' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'
            endif
            restore, gcfi_nobb

            
            tmp_co_gc_nobb   = gc.co*1e9
            tmp_O3_gc_nobb   = gc.o3*1e9
            tmp_pan_gc_nobb  = gc.pan*1e9
            tmp_hcho_gc_nobb = gc.ch2o*1e9
            tmp_acet_gc_nobb = gc.acet*1e9/3; 3C acetone in ppb
            tmp_benz_gc_nobb = gc.benz*1e9/6; 6C benzene in ppb
    ;        ch3oh_gc= [ch3oh_gc, gc.MOH*1e9]
            tmp_ald2_gc_nobb = gc.ald2*1e9/2; 2C ch3cho in ppb

            tmp_no_gc_nobb   = gc.no*1e9
            tmp_no2_gc_nobb  = gc.no2*1e9
            tmp_so2_gc_nobb  = gc.so2*1e9

            tmp_na = avo * gc.pres * 100./(8.31 * gc.temp) * 1e-6;; air density molec/cm3
            tmp_oh_gc_nobb  = gc.oh*tmp_na ;; v/v --> molec/cm3 

            tmp_date_gc_nobb = gc.date
            tmp_utc_gc_nobb  = gc.utc
            tmp_doy_gc_nobb  = gc.doy
            tmp_lat_gc_nobb  = gc.lat
            tmp_lon_gc_nobb  = gc.lon
            tmp_alt_gc_nobb  = gc.alt
            tmp_prs_gc_nobb  = gc.pres
            tmp_temp_gc_nobb = gc.temp

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
            
            tmp_c2h6_gc_nobb = gc.c2h6*1e9
        ;;add lumped species
            tmp_alk4_gc_nobb   = gc.alk4*1e9/4.3 ; 4.3C
            tmp_prpe_gc_nobb   = gc.prpe*1e9/3   ; 3C 
            tmp_rcho_gc_nobb   = gc.rcho*1e9

    ;; concert hhmm(utc) into ss(utc)      
            ;print,tmp_utc_gc_nobb
            hh = floor(tmp_utc_gc_nobb/100)
            ind = where(hh gt 12, ct)
            if ct gt 0 then hh[ind] = hh[ind]
            ind = where(hh lt 12, ct)
            if ct gt 0 then hh[ind] = hh[ind] + 24
            
            mm = tmp_utc_gc_nobb- floor(tmp_utc_gc_nobb/100)*100
            tmp_utc_gc_nobb = float(hh)*60*60+float(mm)*60
            
            
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
            temp_gc_nobb = [temp_gc_nobb,tmp_temp_gc_nobb]

            
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
            
            c2h6_gc_nobb  = [c2h6_gc_nobb, tmp_c2h6_gc_nobb]
        ;;add lumped species
            alk4_gc_nobb  = [alk4_gc_nobb,tmp_alk4_gc_nobb]
            prpe_gc_nobb  = [prpe_gc_nobb,tmp_prpe_gc_nobb]
            rcho_gc_nobb  = [rcho_gc_nobb,tmp_rcho_gc_nobb]

            undefine,gc

            
        endfor
;; ==========================================================
;; Remove placeholder
;; ==========================================================
        utc_obs = utc_obs[1:*]
        doy_obs = doy_obs[1:*]
        lat_obs  = lat_obs[1:*]  
        lon_obs  = lon_obs[1:*]
        prs_obs = prs_obs[1:*]
        alt_obs = alt_obs[1:*]
        rh_obs  = rh_obs[1:*]
        temp_obs = temp_obs[1:*]
        transect_mce_odr_obs = transect_mce_odr_obs[1:*]
        transect_smoke_age_obs = transect_smoke_age_obs[1:*]
        ;; play with lon
        ind = where(lon_obs gt 180,ct)
        if ct gt 0 then lon_obs[ind] = lon_obs[ind] - 360
        
        ;; get local time in seconds of the day
        timediff = lon_obs/15
        lct_obs = utc_obs + timediff*60*60
        
        ;; use local time to get hh, mm, and ss
        local_time_hh = floor(lct_obs/60/60)
        local_time_mm = floor((lct_obs - local_time_hh*60*60)/60)
        local_time_ss = lct_obs - local_time_hh*60*60 - local_time_mm*60        

;; CL measurement
        o3_obs_cl = o3_obs_cl[1:*]
        no_obs_cl = no_obs_cl[1:*]
        no2_obs_cl = no2_obs_cl[1:*]
        noy_obs_cl = noy_obs_cl[1:*]
;; Others
        o3_obs_roze = o3_obs_roze[1:*]
        no_obs_lif  = no_obs_lif[1:*]
        no2_obs_aces = no2_obs_aces[1:*]
        no2_obs_canoe = no2_obs_canoe[1:*]
        hno2_obs_cims = hno2_obs_cims[1:*]
        hno2_obs_aces = hno2_obs_aces[1:*]
        hono_obs_saga = hono_obs_saga[1:*]
;; illegal string, lxu
        ;hno3_obs_cit =  hno3_obs_cit[1:*]
;; CIMS
        pan_obs_cims = pan_obs_cims[1:*]
        ppn_obs_cims = ppn_obs_cims[1:*]
        apan_obs_cims= apan_obs_cims[1:*]
        pbn_obs_cims = pbn_obs_cims[1:*]
        n2o5_obs_cims= n2o5_obs_cims[1:*]
;; others
        n2o_obs_lgr = n2o_obs_lgr[1:*]
;; DACOM
        co_obs_dacom = co_obs_dacom[1:*]
        co_obs_lgr   = co_obs_lgr[1:*]
        ch4_obs_dacom = ch4_obs_dacom[1:*]
        co2_obs_700   = co2_obs_700[1:*]
        c2h6_obs_cams = c2h6_obs_cams[1:*]
        ch2o_obs_cams = ch2o_obs_cams[1:*]
        ch2o_obs_isaf = ch2o_obs_isaf[1:*]
        ch2o_obs_uioptr  = ch2o_obs_uioptr[1:*]
;; WAS
        ocs_obs_was  =  ocs_obs_was[1:*]
        dms_obs_was   = dms_obs_was[1:*]
        ;; RNO2, CH3X
        c2h6_obs_was = c2h6_obs_was[1:*]
        c2h4_obs_was = c2h4_obs_was[1:*]
        c2h2_obs_was = c2h2_obs_was[1:*]
        c3h6_obs_was = c3h6_obs_was[1:*]
        c3h8_obs_was = c3h8_obs_was[1:*]
        propadiene_obs_was = propadiene_obs_was[1:*]
        propyne_obs_was = propyne_obs_was[1:*]
        iButane_obs_was = iButane_obs_was[1:*]
        nButane_obs_was = nButane_obs_was[1:*]
        x1Butene_obs_was= x1Butene_obs_was[1:*]
        iButene_obs_was = iButene_obs_was[1:*]
        t2Butene_obs_was= t2Butene_obs_was[1:*]
        c2Butene_obs_was= c2Butene_obs_was[1:*]
        x13Butadiene_obs_was = x13Butadiene_obs_was[1:*]
        x12Butadiene_obs_was = x12Butadiene_obs_was[1:*]
        x1Buten3yne_obs_was  = x1Buten3yne_obs_was[1:*]
        x13Butadyine_obs_was = x13Butadyine_obs_was[1:*]
        x1Butyne_obs_was = x1Butyne_obs_was[1:*]
        x2Butyne_obs_was = x2Butyne_obs_was[1:*]
        iPentane_obs_was = iPentane_obs_was[1:*]
        nPentane_obs_was = nPentane_obs_was[1:*]
        Isoprene_obs_was = Isoprene_obs_was[1:*]
        x1Pentene_obs_was = x1Pentene_obs_was[1:*]
        t2Pentene_obs_was = t2Pentene_obs_was[1:*]
        c2Pentene_obs_was = c2Pentene_obs_was[1:*]
;; illegal string, waiting for being fixed
        ;threeMe1Butene_obs_was = threeMe1Butene_obs_was[1:*]
        ;twoMe1Butene_obs_was   = twoMe1Butene_obs_was[1:*]
        ;twoMe2Butene_obs_was  = twoMe2Butene_obs_was[1:*]
        x13Pentadienes_obs_was = x13Pentadienes_obs_was[1:*]
        x3Me1PenteneAnd4Me1Pentene_obs_was = x3Me1PenteneAnd4Me1Pentene_obs_was[1:*]
        x1Hexene_obs_was = x1Hexene_obs_was[1:*]
        x1Heptene_obs_was = x1Heptene_obs_was[1:*]
        x1Octene_obs_was = x1Octene_obs_was[1:*]
        x1Nonene_obs_was = x1Nonene_obs_was[1:*]
        x1Decene_obs_was = x1Decene_obs_was[1:*]
        nHexane_obs_was  = nHexane_obs_was[1:*]
        nHeptane_obs_was = nHeptane_obs_was[1:*]
        nOctane_obs_was  = nOctane_obs_was[1:*]
        nNonane_obs_was  = nNonane_obs_was[1:*]
        nDecane_obs_was  = nDecane_obs_was[1:*]
        nUndecane_obs_was = nUndecane_obs_was[1:*]
        x22Dimebutane_obs_was = x22Dimebutane_obs_was[1:*]
        x23Dimebutane_obs_was = x23Dimebutane_obs_was[1:*]
        x2MePentane_obs_was   = x2MePentane_obs_was[1:*]
        x3MePentane_obs_was   = x3MePentane_obs_was[1:*]
        x2MeHexane_obs_was    = x2MeHexane_obs_was[1:*]
        x3MeHexane_obs_was    = x3MeHexane_obs_was[1:*]
        x23DimePentane_obs_was= x23DimePentane_obs_was[1:*]
        x224TrimePentane_obs_was = x224TrimePentane_obs_was[1:*]
        x234TrimePentane_obs_was = x234TrimePentane_obs_was[1:*]
        CycPentane_obs_was       = CycPentane_obs_was[1:*]
        MeCycPentane_obs_was     = MeCycPentane_obs_was[1:*]
        CycHexane_obs_was        = CycHexane_obs_was[1:*]
        MeCycHexane_obs_was      = MeCycHexane_obs_was[1:*]
        CycPentene_obs_was       = CycPentene_obs_was[1:*]
        Benzene_obs_was = Benzene_obs_was[1:*]
        Toluene_obs_was = Toluene_obs_was[1:*]
        EthBenzene_obs_was = EthBenzene_obs_was[1:*]
        mpXylene_obs_was = mpXylene_obs_was[1:*]
        oXylene_obs_was = oXylene_obs_was[1:*]
        Styrene_obs_was = Styrene_obs_was[1:*]
        EthynylBenzene_obs_was = EthynylBenzene_obs_was[1:*]
        iPropBenzene_obs_was = iPropBenzene_obs_was[1:*]
        nPropBenzene_obs_was = nPropBenzene_obs_was[1:*]
        x3EthToluene_obs_was = x3EthToluene_obs_was[1:*]
        x4EthToluene_obs_was = x4EthToluene_obs_was[1:*]
        x2EthToluene_obs_was = x2EthToluene_obs_was[1:*]
        x135rimeBenzene_obs_was = x135rimeBenzene_obs_was[1:*]
        x124rimeBenzene_obs_was = x124rimeBenzene_obs_was[1:*]
        ClBenzene_obs_was = ClBenzene_obs_was[1:*]
        aPinene_obs_was = aPinene_obs_was[1:*]
        bPinene_obs_was = bPinene_obs_was[1:*]
        Tricyclene_obs_was = Tricyclene_obs_was[1:*]
        Camphene_obs_was = Camphene_obs_was[1:*]
        Myrcene_obs_was = Myrcene_obs_was[1:*]
        Limonene_obs_was = Limonene_obs_was[1:*]
        Furan_obs_was = Furan_obs_was[1:*]
        x2MeFuran_obs_was = x2MeFuran_obs_was[1:*]
        x3MeFuran_obs_was = x3MeFuran_obs_was[1:*]
        BenzFuran_obs_was = BenzFuran_obs_was[1:*]
        iButanal_obs_was = iButanal_obs_was[1:*]
        Butanal_obs_was = Butanal_obs_was[1:*]
        AcetonePropanal_obs_was = AcetonePropanal_obs_was[1:*]
        Acetone_obs_was = Acetone_obs_was[1:*]
        Propanal_obs_was = Propanal_obs_was[1:*]

        ;; RCHO
        rcho_obs_was = iButanal_obs_was + Butanal_obs_was + Propanal_obs_was
        
        MEK_obs_was = MEK_obs_was[1:*]
        MAC_obs_was = MAC_obs_was[1:*]
        MVK_obs_was = MVK_obs_was[1:*]
        Acrolein_obs_was = Acrolein_obs_was[1:*]
        iPropanol_obs_was = iPropanol_obs_was[1:*]
        Nitromethane_obs_was = Nitromethane_obs_was[1:*]
        Acrylonitrile_obs_was = Acrylonitrile_obs_was[1:*]
        PropNitrile_obs_was = PropNitrile_obs_was[1:*]
        MeAcetate_obs_was = MeAcetate_obs_was[1:*]
        ;; ALK4 and PRPE
        ;; ALK4
        alk4_obs_was = iButane_obs_was + nButane_obs_was + $
                        ;x1Butene_obs_was + iButene_obs_was +  t2Butene_obs_was + c2Butene_obs_was + $
                        iPentane_obs_was + nPentane_obs_was + $
                        nHexane_obs_was  + x22Dimebutane_obs_was + x23Dimebutane_obs_was + x2MePentane_obs_was + x3MePentane_obs_was + $
                        nHeptane_obs_was  + x2MeHexane_obs_was + x3MeHexane_obs_was + x23DimePentane_obs_was + $
                        nOctane_obs_was  + x224TrimePentane_obs_was + x234TrimePentane_obs_was + $
                        nNonane_obs_was + $
                        nDecane_obs_was + $
                        nUndecane_obs_was 
        ;; PRPE
        prpe_obs_was = c3h6_obs_was + $
                        x1Butene_obs_was + iButene_obs_was +  t2Butene_obs_was + c2Butene_obs_was + $
                        x1Pentene_obs_was + t2Pentene_obs_was + c2Pentene_obs_was + $
                        ;threeMe1Butene_obs_was + twoMe1Butene_obs_was + twoMe2Butene_obs_was + $ ; ilegal string, waiting for being fixed
                        x3Me1PenteneAnd4Me1Pentene_obs_was + x1Hexene_obs_was + $
                        x1Heptene_obs_was + $
                        x1Octene_obs_was  + $
                        x1Nonene_obs_was + $
                        x1Decene_obs_was 

    
;; TOGA measurement
        ;; CH3X
        DMS_obs_toga = DMS_obs_toga[1:*]
        Propane_obs_toga = Propane_obs_toga[1:*]
        iButane_obs_toga = iButane_obs_toga[1:*]
        nButane_obs_toga = nButane_obs_toga[1:*]
        iPentane_obs_toga = iPentane_obs_toga[1:*]
        nPentane_obs_toga = nPentane_obs_toga[1:*]
        x2MePentane_obs_toga = x2MePentane_obs_toga[1:*]
        x3MePentane_obs_toga = x3MePentane_obs_toga[1:*]
        nHexane_obs_toga = nHexane_obs_toga[1:*]
        x224TrimePentane_obs_toga = x224TrimePentane_obs_toga[1:*]
        nHeptane_obs_toga = nHeptane_obs_toga[1:*]
        nOctane_obs_toga = nOctane_obs_toga[1:*]
        Propene_obs_toga = Propene_obs_toga[1:*]
        iButene1Butene_obs_toga = iButene1Butene_obs_toga[1:*]
        Isoprene_obs_toga = Isoprene_obs_toga[1:*]
        Tricyclene_obs_toga = Tricyclene_obs_toga[1:*]
        aPinene_obs_toga = aPinene_obs_toga[1:*]
        Camphene_obs_toga = Camphene_obs_toga[1:*]
        bPineneMyrcene_obs_toga = bPineneMyrcene_obs_toga[1:*]
        LimoneneD3Carene_obs_toga = LimoneneD3Carene_obs_toga[1:*]
        Benzene_obs_toga = Benzene_obs_toga[1:*]
        Toluene_obs_toga = Toluene_obs_toga[1:*]
        Xylenes_obs_toga = Xylenes_obs_toga[1:*]
        C8Aromatics_obs_toga= C8Aromatics_obs_toga[1:*]
        EthBenzene_obs_toga = EthBenzene_obs_toga[1:*]
        mpXylene_obs_toga = mpXylene_obs_toga[1:*]
        oXylene_obs_toga = oXylene_obs_toga[1:*]
        Styrene_obs_toga = Styrene_obs_toga[1:*]
        EthynylBenzene_obs_toga = EthynylBenzene_obs_toga[1:*]
        CH2O_obs_toga = CH2O_obs_toga[1:*]
        CH3CHO_obs_toga = CH3CHO_obs_toga[1:*]
        Propanal_obs_toga = Propanal_obs_toga[1:*]
        Butanal_obs_toga = Butanal_obs_toga[1:*]
        iButanal_obs_toga = iButanal_obs_toga[1:*] 
        Acrolein_obs_toga = Acrolein_obs_toga[1:*]
        x2Butenals_obs_toga = x2Butenals_obs_toga[1:*]
        
        
        ;; RCHO from TOGA
        rcho_obs_toga = Propanal_obs_toga + Butanal_obs_toga + iButanal_obs_toga
        
        Acetone_obs_toga = Acetone_obs_toga[1:*]
        MEK_obs_toga = MEK_obs_toga[1:*]
        CH3OH_obs_toga = CH3OH_obs_toga[1:*]
        C2H5OH_obs_toga = C2H5OH_obs_toga[1:*]
        iPropanol_obs_toga = iPropanol_obs_toga[1:*] 
        MBO_obs_toga = MBO_obs_toga[1:*]
        MAC_obs_toga = MAC_obs_toga[1:*] 
        MVK_obs_toga = MVK_obs_toga[1:*]
        MeFormate_obs_toga = MeFormate_obs_toga[1:*] 
        MeAcetate_obs_toga = MeAcetate_obs_toga[1:*]
        Furan_obs_toga = Furan_obs_toga[1:*]
        x2MeFuran_obs_toga = x2MeFuran_obs_toga[1:*]
        x3MeFuran_obs_toga = x3MeFuran_obs_toga[1:*]
        Furfural_obs_toga = Furfural_obs_toga[1:*]
        HCN_obs_toga = HCN_obs_toga[1:*]
        CH3CN_obs_toga = CH3CN_obs_toga[1:*]
        PropNitrile_obs_toga = PropNitrile_obs_toga[1:*]
        Acrylonitrile_obs_toga = Acrylonitrile_obs_toga[1:*]
        MeAcrylonitrile_obs_toga = MeAcrylonitrile_obs_toga[1:*]
        Pyrrole_obs_toga = Pyrrole_obs_toga[1:*]
        Nitromethane_obs_toga = Nitromethane_obs_toga[1:*]
        MeONO2_obs_toga = MeONO2_obs_toga[1:*]
        EthONO2_obs_toga = EthONO2_obs_toga[1:*]
        iPropONO2_obs_toga = iPropONO2_obs_toga[1:*]
        x2ButONO2iButONO2_obs_toga = x2ButONO2iButONO2_obs_toga[1:*]
;; IWAS
        ;; CH3X
        Ethane_obs_iwas  = Ethane_obs_iwas[1:*]
        Propane_obs_iwas = Propane_obs_iwas[1:*]
        nButane_obs_iwas = nButane_obs_iwas[1:*]
        iButane_obs_iwas = iButane_obs_iwas[1:*]
        nPentane_obs_iwas= nPentane_obs_iwas[1:*]
        iPentane_obs_iwas= iPentane_obs_iwas[1:*]
        nHexane_obs_iwas = nHexane_obs_iwas[1:*]
        x2MePentane_obs_iwas = x2MePentane_obs_iwas[1:*]
        x3MePentane_obs_iwas = x3MePentane_obs_iwas[1:*]
        x22DiMeButane_obs_iwas = x22DiMeButane_obs_iwas[1:*]
        x24DiMePentane_obs_iwas= x24DiMePentane_obs_iwas[1:*]
        nOctane_obs_iwas = nOctane_obs_iwas[1:*]
        x224TriMePentane_obs_iwas = x224TriMePentane_obs_iwas[1:*]
        nNonane_obs_iwas = nNonane_obs_iwas[1:*]
        nDecane_obs_iwas = nDecane_obs_iwas[1:*]
        MeCycPentane_obs_iwas = MeCycPentane_obs_iwas[1:*]
        CycHexane_obs_iwas = CycHexane_obs_iwas[1:*]
        MeCycHexane_obs_iwas = MeCycHexane_obs_iwas[1:*]
        Ethyne_obs_iwas = Ethyne_obs_iwas[1:*]
        Ethene_obs_iwas = Ethene_obs_iwas[1:*]
        Propene_obs_iwas = Propene_obs_iwas[1:*]
        x1Butene_obs_iwas = x1Butene_obs_iwas[1:*]
        c2Butene_obs_iwas = c2Butene_obs_iwas[1:*]
        t2butene_obs_iwas = t2butene_obs_iwas[1:*]
        iButene_obs_iwas = iButene_obs_iwas[1:*]
        x1Pentene_obs_iwas = x1Pentene_obs_iwas[1:*]
        c2Pentene_obs_iwas = c2Pentene_obs_iwas[1:*]
        t2Pentene_obs_iwas = t2Pentene_obs_iwas[1:*]
        x2Me1Butene_obs_iwas = x2Me1Butene_obs_iwas[1:*]
        x3Me1Butene_obs_iwas = x3Me1Butene_obs_iwas[1:*]
        t13Pentadiene_obs_iwas = t13Pentadiene_obs_iwas[1:*]
        Isoprene_obs_iwas = Isoprene_obs_iwas[1:*]
        aPinene_obs_iwas = aPinene_obs_iwas[1:*]
        Benzene_obs_iwas = Benzene_obs_iwas[1:*]
        Toluene_obs_iwas = Toluene_obs_iwas[1:*]
        EthBenzene_obs_iwas = EthBenzene_obs_iwas[1:*]
        oXylene_obs_iwas = oXylene_obs_iwas[1:*]
        mpXylene_obs_iwas = mpXylene_obs_iwas[1:*]
        Acetone_obs_iwas = Acetone_obs_iwas[1:*]
        MEK_obs_iwas = MEK_obs_iwas[1:*]
        MeFormate_obs_iwas = MeFormate_obs_iwas[1:*]
        Furan_obs_iwas = Furan_obs_iwas[1:*]
        CH3CN_obs_iwas = CH3CN_obs_iwas[1:*]
        Acrylonitrile_obs_iwas = Acrylonitrile_obs_iwas[1:*]
        
        ;; PTR and ALK4
        ;; ALK4
        alk4_obs_iwas = nButane_obs_iwas + iButane_obs_iwas + $
                        nPentane_obs_iwas + iPentane_obs_iwas + $
                        nHexane_obs_iwas + x2MePentane_obs_iwas + x3MePentane_obs_iwas + x22DiMeButane_obs_iwas + $
                        x24DiMePentane_obs_iwas + $
                        nOctane_obs_iwas + x224TriMePentane_obs_iwas + $
                        nNonane_obs_iwas + $
                        nDecane_obs_iwas 
        
        
        ;; PRPE
        prpe_obs_iwas = Propene_obs_iwas + $
                        x1Butene_obs_iwas + c2Butene_obs_iwas + t2Butene_obs_iwas +  iButene_obs_iwas + $
                        x1Pentene_obs_iwas + c2Pentene_obs_iwas + t2Pentene_obs_iwas + x2Me1Butene_obs_iwas + x3Me1Butene_obs_iwas 
        
;; NOAAPTR
        HCN_obs_noaaptr = HCN_obs_noaaptr[1:*]
        CH2O_obs_noaaptr = CH2O_obs_noaaptr[1:*]
        CH3OH_obs_noaaptr = CH3OH_obs_noaaptr[1:*]
        CH3CN_obs_noaaptr = CH3CN_obs_noaaptr[1:*]
        HNCO_obs_noaaptr = HNCO_obs_noaaptr[1:*]
        CH3CHO_obs_noaaptr = CH3CHO_obs_noaaptr[1:*]
        C2H5OH_obs_noaaptr = C2H5OH_obs_noaaptr[1:*]
        HCOOH_obs_noaaptr = HCOOH_obs_noaaptr[1:*]
        Acrylonitrile_obs_noaaptr = Acrylonitrile_obs_noaaptr[1:*]
        Acrolein_obs_noaaptr = Acrolein_obs_noaaptr[1:*]
        AcetonePropanal_obs_noaaptr = AcetonePropanal_obs_noaaptr[1:*]
        Acetone_obs_noaaptr = Acetone_obs_noaaptr[1:*]
        Propanal_obs_noaaptr = Propanal_obs_noaaptr[1:*]

        ;; RCHO
        rcho_obs_ptr = Propanal_obs_noaaptr
        
        GlycolaldehydeCH3COOH_obs_noaaptr = GlycolaldehydeCH3COOH_obs_noaaptr[1:*]
        CH3NO2_obs_noaaptr = CH3NO2_obs_noaaptr[1:*]
        DMS_obs_noaaptr = DMS_obs_noaaptr[1:*]
        C4H5N_obs_noaaptr = C4H5N_obs_noaaptr[1:*]
        Furan_obs_noaaptr = Furan_obs_noaaptr[1:*]
        MVKMAC_obs_noaaptr = MVKMAC_obs_noaaptr[1:*]
        C4Carbonyls_obs_noaaptr = C4Carbonyls_obs_noaaptr[1:*]
        MEK_obs_noaaptr = MEK_obs_noaaptr[1:*]
        Butanal_obs_noaaptr = Butanal_obs_noaaptr[1:*]
        C3H6O2_obs_noaaptr = C3H6O2_obs_noaaptr[1:*]
        Benzene_obs_noaaptr = Benzene_obs_noaaptr[1:*]
        x2MeFuranx3MeFuran_obs_noaaptr = x2MeFuranx3MeFuran_obs_noaaptr[1:*]
        x2Furanone_obs_noaaptr = x2Furanone_obs_noaaptr[1:*]
        x23Butanedione_obs_noaaptr = x23Butanedione_obs_noaaptr[1:*]
        Toluene_obs_noaaptr = Toluene_obs_noaaptr[1:*]
        Phenol_obs_noaaptr = Phenol_obs_noaaptr[1:*]
        Furfural_obs_noaaptr = Furfural_obs_noaaptr[1:*]
        DimeFurans_obs_noaaptr = DimeFurans_obs_noaaptr[1:*]
        MaleicAnhyd_obs_noaaptr = MaleicAnhyd_obs_noaaptr[1:*]
        BenzNitrile_obs_noaaptr = BenzNitrile_obs_noaaptr[1:*] 
        Styrene_obs_noaaptr = Styrene_obs_noaaptr[1:*] 
        Benzaldehyde_obs_noaaptr = Benzaldehyde_obs_noaaptr[1:*]
        Xylenes_obs_noaaptr = Xylenes_obs_noaaptr[1:*]
        C8Aromatics_obs_noaaptr = C8Aromatics_obs_noaaptr[1:*]
        C7H8O_obs_noaaptr = C7H8O_obs_noaaptr[1:*]
        Catecholx5MeFurfural_obs_noaaptr = Catecholx5MeFurfural_obs_noaaptr[1:*]
        BenzFuran_obs_noaaptr = BenzFuran_obs_noaaptr[1:*]
        C9Aromatics_obs_noaaptr = C9Aromatics_obs_noaaptr[1:*]
        C6H4O3_obs_noaaptr = C6H4O3_obs_noaaptr[1:*]
        Guaiacol_obs_noaaptr = Guaiacol_obs_noaaptr[1:*]
        Naphthalene_obs_noaaptr = Naphthalene_obs_noaaptr[1:*]
        Monoterpenes_obs_noaaptr = Monoterpenes_obs_noaaptr[1:*]
        Creosols_obs_noaaptr = Creosols_obs_noaaptr[1:*]
        Syringol_obs_noaaptr = Syringol_obs_noaaptr[1:*]

;; CIMS
        ;; Some unknown VOCs
        HCN_obs_cims = HCN_obs_cims[1:*]
        HCOOH_obs_cims = HCOOH_obs_cims[1:*]
        HNCO_obs_cims = HNCO_obs_cims[1:*]



    
;; ============================
;; GEOS-Chem 0.25x0.3125: GFAS
;; ============================
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
        temp_gc_gfas = temp_gc_gfas[1:*]

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
        
        c2h6_gc_gfas = c2h6_gc_gfas[1:*]
        ;;add lumped species
        alk4_gc_gfas = alk4_gc_gfas[1:*]
        prpe_gc_gfas = prpe_gc_gfas[1:*]
        rcho_gc_gfas = rcho_gc_gfas[1:*]


;; ============================
;; GEOS-Chem 0.25x0.3125: THREEGFAS
;; ============================
    ;; THREEGFAS
        co_gc_threegfas   = co_gc_threegfas[1:*]
        o3_gc_threegfas   = o3_gc_threegfas[1:*]
        pan_gc_threegfas  = pan_gc_threegfas[1:*]
        hcho_gc_threegfas = hcho_gc_threegfas[1:*]
        acet_gc_threegfas = acet_gc_threegfas[1:*]
        benz_gc_threegfas = benz_gc_threegfas[1:*]
    ;    ch3oh_gc= ch3oh_gc[1:*]
        ald2_gc_threegfas = ald2_gc_threegfas[1:*]

        no_gc_threegfas   = no_gc_threegfas[1:*]
        no2_gc_threegfas  = no2_gc_threegfas[1:*]
        so2_gc_threegfas  = so2_gc_threegfas[1:*]
        oh_gc_threegfas   = oh_gc_threegfas[1:*] 

        date_gc_threegfas = date_gc_threegfas[1:*]
        ;utc_gc_threegfas  = utc_gc_threegfas[1:*]
        doy_gc_threegfas  = doy_gc_threegfas[1:*]
        lat_gc_threegfas  = lat_gc_threegfas[1:*]
        lon_gc_threegfas  = lon_gc_threegfas[1:*]
        alt_gc_threegfas  = alt_gc_threegfas[1:*]
        prs_gc_threegfas  = prs_gc_threegfas[1:*]
        temp_gc_threegfas = temp_gc_threegfas[1:*]

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
        
        c2h6_gc_threegfas = c2h6_gc_threegfas[1:*]
        ;;add lumped species
        alk4_gc_threegfas = alk4_gc_threegfas[1:*]
        prpe_gc_threegfas = prpe_gc_threegfas[1:*]
        rcho_gc_threegfas = rcho_gc_threegfas[1:*]

;; ============================
;; GEOS-Chem 0.25x0.3125: NOBB
;; ============================
    ;; NOBB
        co_gc_nobb   = co_gc_nobb[1:*]
        o3_gc_nobb   = o3_gc_nobb[1:*]
        pan_gc_nobb  = pan_gc_nobb[1:*]
        hcho_gc_nobb = hcho_gc_nobb[1:*]
        acet_gc_nobb = acet_gc_nobb[1:*]
        benz_gc_nobb = benz_gc_nobb[1:*]
    ;    ch3oh_gc= ch3oh_gc[1:*]
        ald2_gc_nobb = ald2_gc_nobb[1:*]

        no_gc_nobb   = no_gc_nobb[1:*]
        no2_gc_nobb  = no2_gc_nobb[1:*]
        so2_gc_nobb  = so2_gc_nobb[1:*]
        oh_gc_nobb   = oh_gc_nobb[1:*] 

        date_gc_nobb = date_gc_nobb[1:*]
        ;utc_gc_nobb  = utc_gc_nobb[1:*]
        doy_gc_nobb  = doy_gc_nobb[1:*]
        lat_gc_nobb  = lat_gc_nobb[1:*]
        lon_gc_nobb  = lon_gc_nobb[1:*]
        alt_gc_nobb  = alt_gc_nobb[1:*]
        prs_gc_nobb  = prs_gc_nobb[1:*]
        temp_gc_nobb = temp_gc_nobb[1:*]

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
        
        c2h6_gc_nobb = c2h6_gc_nobb[1:*]
        ;;add lumped species
        alk4_gc_nobb = alk4_gc_nobb[1:*]
        prpe_gc_nobb = prpe_gc_nobb[1:*]
        rcho_gc_nobb = rcho_gc_nobb[1:*]
    endfor
;; ==============================
;; FIREX-AQ RL propane
;; ==============================
    for dd = 0, n_elements(dates_all) do begin
        ;;only consider whole flight
        if dd ne 0 then break
        Propane_obs_toga = [0]
        Propane_obs_iwas = [0]
    ; smart loop for each fligiht and total flights from lu
        if dd eq 0 then dates = dates_all else dates=dates_all[dd-1]
        dates = STRTRIM(string(dates),1)
        for n=0,n_elements(dates)-1 do begin
    ;; GC file for all data
            str_date=dates[n]
            print,'processing '+str_date
            dc8fi = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                'mrg2sav/FIREX-AQ/firexaq_avg/firexaq-mrg5m-dc8_merge_'+dates[n]+'_RL.sav'
           
            restore, dc8fi

            tmp_c3h8_obs_was = dc8.Propane_WAS_BLAKE/1e3
            tmp_Propane_obs_toga = dc8.Propane_TOGA_APEL/1E3
            tmp_Propane_obs_iwas = dc8.Propane_NOAAiWAS_GILMAN
            c3h8_obs_was = [c3h8_obs_was, tmp_c3h8_obs_was]
            Propane_obs_toga = [Propane_obs_toga, tmp_Propane_obs_toga]
            Propane_obs_iwas = [Propane_obs_iwas, tmp_Propane_obs_iwas]
            undefine,dc8    
        endfor
        c3h8_obs_was = c3h8_obs_was[1:*]
        Propane_obs_toga = Propane_obs_toga[1:*]
        Propane_obs_iwas = Propane_obs_iwas[1:*]
    endfor
        
        if n_elements(co_obs_dacom) ne n_elements(co_gc_gfas) or $
            n_elements(co_obs_dacom) ne n_elements(co_gc_threegfas) $
            then stop 
;;setting up koh of different VOCs
        temperature_obs = temp_obs
        temperature_mod = temp_gc_gfas
        prs_mod = prs_gc_gfas
        
;; number density of air: molec/m-3
;; used to get number density of trace gas 
        scalefactor_obs = avo*(prs_obs*100)/(8.31*temperature_obs)
        scalefactor_mod = avo*(prs_mod*100)/(8.31*temperature_mod)

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

;; test how many datapoints avaiable for each measurements
;; CO, both dacom and lgr have high temporal resolution
;; NO2, CANOE>CL>ACES
;; O3, ROZE>CL
;; HCN, CIMS ~ PTR ~ TOGA all have very nice tepmoral resolution.
;; CH3CN, PTR>TOGA>IWAS
;; CH2O, PTR=ISAF>CAMS>TOGA
;; ALD2, TOGA>NOAAPTR
;; Acetone, TOGA> IWAS
;; Propane, TOGA> IWAS
;; MEK, TOGA>WAS>IWAS
;; Benzene, NOAAPTR>TOGA>WAS>IWAS
;; Toluene, NOAAPTR>TOGA>WAS>IWAS
;; Xylenes, NOAAPTR>TOGA>WAS>IWAS
;; Ethane, CAMS>WAS>IWAS

;ind = where(o3_obs_cl lt 0, ct)
;if ct gt 0 then o3_obs_cl[ind] = !values.f_nan
;print, ct
;remove_ind = where(finite(o3_obs_cl),ct)
;if ct gt 0 then o3_obs_cl = o3_obs_cl[remove_ind]

;help, o3_obs_cl
;help, o3_obs_roze
;continue

; calculate the percentage of west and east
;print, 'this is dp for whole'
;help, lon_obs
;print, 'this is dp for west'
;ind_test = where(lon_obs le -105, ct_lon)
;help, ind_test
;; ==============================
;; set up detection limit for PTR
;; ==============================
        ;if keyword_set(lod_ptr) then begin
            LoD = 50.0/1000
            ind = where(CH2O_obs_noaaptr le LoD and CH2O_obs_noaaptr gt 0, ct)
            if ct gt 0 then begin
                ;CH2O_obs_noaaptr[ind] = !VALUES.F_NAN
                CH2O_obs_noaaptr[ind] = CH2O_obs_toga[ind]
            endif

            LoD = 10.0/1000
            ind = where(CH3CHO_obs_noaaptr le LoD and CH3CHO_obs_noaaptr gt 0, ct)
            if ct gt 0 then begin
                ;CH3CHO_obs_noaaptr[ind] = !VALUES.F_NAN
                CH3CHO_obs_noaaptr[ind] = CH3CHO_obs_toga[ind]
            endif

            LoD = 20.0/1000
            ind = where(Acetone_obs_noaaptr le LoD and Acetone_obs_noaaptr gt 0, ct)
            if ct gt 0 then begin
                ;Acetone_obs_noaaptr[ind] = !VALUES.F_NAN
                Acetone_obs_noaaptr[ind] = Acetone_obs_toga[ind]
            endif

            LoD = 10.0/1000
            ind = where(HCOOH_obs_noaaptr le LoD and HCOOH_obs_noaaptr gt 0, ct)
            if ct gt 0 then begin
                HCOOH_obs_noaaptr[ind] = HCOOH_obs_cims[ind]
            endif

            LoD = 10.0/1000
            ind = where(MEK_obs_noaaptr le LoD and MEK_obs_noaaptr gt 0, ct)
            if ct gt 0 then begin
                ;MEK_obs_noaaptr[ind] = !VALUES.F_NAN
                MEK_obs_noaaptr[ind] = MEK_obs_toga[ind]
            endif

            LoD = 30.0/1000  ; change 15.0 to 10.0
            ind = where(Benzene_obs_noaaptr le LoD and Benzene_obs_noaaptr gt 0, ct)
            if ct gt 0 then begin
                ;Benzene_obs_noaaptr[ind] = !VALUES.F_NAN
                Benzene_obs_noaaptr[ind] = Benzene_obs_toga[ind]
            endif

            LoD = 50.0/1000 ; change 20.0 to 10.0
            ind = where(Toluene_obs_noaaptr le LoD and Toluene_obs_noaaptr gt 0, ct)
            if ct gt 0 then begin
                ;Toluene_obs_noaaptr[ind] = !VALUES.F_NAN
                Toluene_obs_noaaptr[ind] = Toluene_obs_toga[ind]
            endif

            LoD = 50.0/1000 
            ind = where(Xylenes_obs_noaaptr le LoD and Xylenes_obs_noaaptr gt 0, ct)
            if ct gt 0 then begin
                ;Xylenes_obs_noaaptr[ind] = !VALUES.F_NAN
                Xylenes_obs_noaaptr[ind] = Xylenes_obs_toga[ind]
            endif
            
            LoD = 50.0/1000 
            ind = where(C8Aromatics_obs_noaaptr le LoD and C8Aromatics_obs_noaaptr gt 0, ct)
            if ct gt 0 then begin
                ;Xylenes_obs_noaaptr[ind] = !VALUES.F_NAN
                C8Aromatics_obs_noaaptr[ind] = C8Aromatics_obs_toga[ind]
            endif
        ;endif

;; ========================
;; Prepare for the filters 
;; ========================
    ;; jupyter time
        ;doy_obs = jday_obs + utc_obs / 3600. / 24. 

    ;; lixu, 12/18/2019, 0208/2021   
    ;; plume filter setting, use in the plotting part
        no2_thresh = 4. ; ppb
        co_thresh  = 150. ; ppb
    ;; urban influence
        x224TrimePentane_thresh = 20/1000
        C2Cl4_thresh = 2/1000
        HFC134a_thresh = 125/1000
        HCFC22_thresh = 275/1000
    
    ;; biomass burning threshhold 
        hcn_thresh   = 275./1000 ;pptv to ppbv
        ch3cn_thresh = 156/1000. ;;ppt to ppb, 131
        co_thresh =  100
    ;; strat filter
        strat_thresh = 1.25 ; O3/CO ratio
        
    ;; CO threshold of fresh and aged smoke 
        ;; young
        methylfuran_thresh = 0.7/1000 
        ;; median
        acrolein_thresh = 7.4/1000
        ;; old
        acrylonitrile_thresh = 2.9/1000
        
    ;; get hhmm of each day
        utc_obs = 24.*(doy_obs mod 1.)
        utc_gc_gfas  = 24.*(doy_gc_gfas  mod 1.)

        
    ;; for vertical profile
    ;;    yrange = [0,7]
    ;; lxu, 06202021
        acn_co_tresh = 2.01 ;ppb/ppm
;; ========keep the same for regression, profiles and time series=================
        if keyword_set(addup) then rows = 1
        if keyword_set(addup) then cols = 1

        if keyword_set(all) or keyword_set(SI) or keyword_set(obs_only) then rows = 3
        if keyword_set(all) or keyword_set(SI) or keyword_set(obs_only) then cols = 3

        if keyword_set(each) then rows = 1
        if keyword_set(each) then cols = 1
    
        if keyword_set(onebyone) then rows = 1
        if keyword_set(onebyone) then cols = 1
        
        if keyword_set(chemistry) then rows = 2
        if keyword_set(chemistry) then cols = 2
        
        if keyword_set(OVOCs) then rows = 3
        if keyword_set(OVOCs) then cols = 3

        if keyword_set(NMHCs) then rows = 2
        if keyword_set(NMHCs) then cols = 2
        
        
        multipanel, rows = rows, cols = cols

        if keyword_set(all)  then s_end = 17
        if keyword_set(each) then s_end = 0
        if keyword_set(addup) then s_end = 0
        if keyword_set(chemistry) then s_end = 3
        if keyword_set(nmhcs) then s_end = 3

        if keyword_set(ovocs) or keyword_set(SI) or keyword_set(obs_only)  then s_end = 8

        if keyword_set(measurements) then s_end = 12

        for s = 0, s_end do begin
        
            if keyword_set(obs_only) then begin
                case s of                     
                    ;;c2h6(AWAS, discrete data)
                    ;;hcho
                    0:begin
                        ovoc_ptr = ch2o_obs_cams
                        ovoc_other = ch2o_obs_isaf
                        gvoc_gfas = hcho_gc_gfas
                        gvoc_threegfas=hcho_gc_threegfas
                        gvoc_nobb=hcho_gc_nobb
                        
                        title = 'Formaldehyde'
                        xrange = [0,5]
                        if keyword_set(wus) then xrange = [0,5]  
                        if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) then xrange = [0,5] 
                        ;if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) and keyword_set(wus) then xrange = [0,5]    

                        ;; OH reactivity
                        ;; Units
                        ;; scalefactor is used to convert ppm to molec/m3. In 1e-15, 1e-6 from m3, 1e-9 from ppb
                        koh_voc_obs = tmp_koh_ch2o_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_ch2o_mod*scalefactor_mod*1e-15
                        if keyword_set(OHrplot) then xrange = [0,1]
                    end
                    ;;ald2
                    1:begin
                        ovoc_ptr = CH3CHO_obs_noaaptr
                        ovoc_other = CH3CHO_obs_toga
                        gvoc_gfas = ald2_gc_gfas
                        gvoc_threegfas = ald2_gc_threegfas
                        gvoc_nobb = ald2_gc_nobb
                        
                        title = 'Acetaldehyde'
                        xrange = [0,2]
                        if keyword_set(wus) then xrange = [0,2]  
                        ;if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) then xrange = [0,0.5]    
                        if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) and keyword_set(wus) then xrange = [0,0.8]   
                        
                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_ald2_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_ald2_mod*scalefactor_mod*1e-15
                        if keyword_set(OHrplot) then xrange = [0,0.3]

                    end
                    ;;acet
                    2:begin
                        ovoc_ptr = AcetonePropanal_obs_noaaptr*0.78
                        ovoc_other = Acetone_obs_toga
                        ;ovoc = Acetone_obs_toga # TOGA could be even higher
                        gvoc_gfas = acet_gc_gfas
                        gvoc_threegfas = acet_gc_threegfas
                        gvoc_nobb = acet_gc_nobb
                        title = 'Acetone'
                        xrange = [0,4]
                        ;if keyword_set(wus) then xrange = [0,3]  
                        if keyword_set(wus) then xrange = [0,4.2]  

                        if keyword_set(filter2_nobb) or keyword_set(filter3_nobb)  or keyword_set(filter5_nobb) then xrange = [0,2.5]    
                        ;if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) and keyword_set(wus) then xrange = [0,4]                            
                    
                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_acet_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_acet_mod*scalefactor_mod*1e-15
                        IF keyword_set(OHrplot) then xrange = [0,0.03]
 

                    end
                    ;;HCOOH
                    3:begin
                        ovoc_ptr = HCOOH_obs_noaaptr
                        ovoc_other = HCOOH_obs_cims
                        gvoc_gfas = hcooh_gc_gfas
                        gvoc_threegfas = hcooh_gc_threegfas
                        gvoc_nobb = hcooh_gc_nobb

                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_hcooh_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_hcooh_mod*scalefactor_mod*1e-15
            
                        title = 'HCOOH from PTR'
                        xrange = [0,15]
                        if keyword_set(wus) then xrange = [0,6]
                        if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) and keyword_set(wus) then xrange = [0,6]    

                        if keyword_set(OHrplot) then xrange = [0,0.15]

                    end
                    ;;mek
                    4:begin
                        ovoc_ptr = C4Carbonyls_obs_noaaptr*0.8
                        ovoc_other = MEK_obs_toga
                        gvoc_gfas = mek_gc_gfas
                        gvoc_threegfas = mek_gc_threegfas
                        gvoc_nobb = mek_gc_nobb
                        
                        title = 'MEK'
                        xrange = [0,2]
                        if keyword_set(wus) then xrange = [0,0.2]  

                        ;if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) then xrange = [0,0.3]   
                        if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) and keyword_set(wus) then xrange = [0,0.4]    
                        
                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_mek_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_mek_mod*scalefactor_mod*1e-15
                        if keyword_set(OHrplot) then xrange = [0,0.02]

                        
                    end
                    ;;eoh
                    5:begin
                        ovoc_ptr = C2H5OH_obs_noaaptr
                        ovoc_other = C2H5OH_obs_toga
                        gvoc_gfas = c2h5oh_gc_gfas
                        gvoc_threegfas = c2h5oh_gc_threegfas
                        gvoc_nobb = c2h5oh_gc_nobb

                        
                        title = 'Ethanol'
                        xrange = [0,8]
                        if keyword_set(wus) then xrange = [0,6]  
                        if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) then xrange = [0,2]   

                        ; waiting to be corrected for EOH
                        koh_voc_obs = tmp_koh_co_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_co_mod*scalefactor_mod*1e-15
                        if keyword_set(OHrplot) then xrange = [0,0.05]

                    end
                    ;;benz
                    6:begin
                        ovoc_ptr = Benzene_obs_noaaptr
                        ovoc_other = Benzene_obs_toga
                        gvoc_gfas = benz_gc_gfas
                        gvoc_threegfas = benz_gc_threegfas
                        gvoc_nobb = benz_gc_nobb
                        
                        title = 'Benzene'
                        xrange = [0,0.1]
                        if keyword_set(wus) then xrange = [0,0.2]  

                        if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) then xrange = [0,0.05]    
                        ;if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) and keyword_set(wus) then xrange = [0,0.1]   
                        
                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_benz_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_benz_mod*scalefactor_mod*1e-15
                        if keyword_set(OHrplot) then xrange = [0,0.003]

                    end
                    ;;tolu
                    7:begin
                        ovoc_ptr = Toluene_obs_noaaptr
                        ovoc_other = Toluene_obs_toga
                        gvoc_gfas = c7h8_gc_gfas
                        gvoc_threegfas = c7h8_gc_threegfas
                        gvoc_nobb = c7h8_gc_nobb
                        
                        title = 'Toluene'
                        xrange = [0,0.1]
                        if keyword_set(wus) then xrange = [0,0.2] 

                        if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) then xrange = [0,0.015]    
                        ;if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) and keyword_set(wus) then xrange = [0,0.1]    

                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_tolu_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_tolu_mod*scalefactor_mod*1e-15
                        if keyword_set(OHrplot) then xrange = [0,0.01]

                    end  
                    ;;xyle, either was or ptr would be fine, ptr measures C8 aromatics
                    8:begin
                        ;ovoc_was = mpXylene_obs_was + oXylene_obs_was
                        ;ovoc_iwas = mpXylene_obs_iwas + oXylene_obs_iwas
                        ovoc_ptr = C8Aromatics_obs_noaaptr*0.65
                        ovoc_other = mpXylene_obs_toga + oXylene_obs_toga
                        gvoc_gfas = c8h10_gc_gfas
                        gvoc_threegfas = c8h10_gc_threegfas
                        gvoc_nobb = c8h10_gc_nobb

                        title = 'Xylene'
                        xrange = [0, 0.05]
                        ;if keyword_set(wus) then xrange = [0,0.05]
                        if keyword_set(wus) then xrange = [0,0.04]

                        if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) then xrange = [0,0.002]    
                        ;if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) and keyword_set(wus) then xrange = [0,0.05]    

                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_xyle_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_xyle_mod*scalefactor_mod*1e-15
                        if keyword_set(OHrplot) then xrange = [0,0.02]
                    end

                endcase
            endif
            
    ;;all
            if keyword_set(all) then begin
                case s of                     
                    ;;c2h6(AWAS, discrete data)
                    ;;hcho
                    0:begin
                        ovoc_ptr = ch2o_obs_cams
                        ovoc_other = ch2o_obs_isaf
                        gvoc_gfas = hcho_gc_gfas
                        gvoc_threegfas=hcho_gc_threegfas
                        gvoc_nobb=hcho_gc_nobb
                        
                        title = 'Formaldehyde'
                        xrange = [0,5]
                        if keyword_set(wus) then xrange = [0,5]  
                        if keyword_set(filter2_nobb) and keyword_set(filter3_nobb) and keyword_set(wus) then xrange = [0,5] 
                        ;if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) and keyword_set(wus) then xrange = [0,5]    

                        ;; OH reactivity
                        ;; Units
                        ;; scalefactor is used to convert ppm to molec/m3. In 1e-15, 1e-6 from m3, 1e-9 from ppb
                        koh_voc_obs = tmp_koh_ch2o_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_ch2o_mod*scalefactor_mod*1e-15
                        if keyword_set(OHrplot) then xrange = [0,1]
                    end
                    ;;ald2
                    1:begin
                        ovoc_ptr = CH3CHO_obs_noaaptr
                        ovoc_other = CH3CHO_obs_toga
                        gvoc_gfas = ald2_gc_gfas
                        gvoc_threegfas = ald2_gc_threegfas
                        gvoc_nobb = ald2_gc_nobb
                        
                        title = 'Acetaldehyde'
                        xrange = [0,2]
                        if keyword_set(wus) then xrange = [0,2]  
                        ;if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) then xrange = [0,0.5]    
                        if keyword_set(filter2_nobb) and keyword_set(filter3_nobb) and keyword_set(wus) then xrange = [0,0.6]   
                        
                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_ald2_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_ald2_mod*scalefactor_mod*1e-15
                        if keyword_set(OHrplot) then xrange = [0,0.3]

                    end
                    ;;acet
                    2:begin
                        ovoc_ptr       = Acetone_obs_noaaptr
                        ovoc_other     = Acetone_obs_toga
                        ;ovoc = Acetone_obs_toga # TOGA could be even higher
                        gvoc_gfas = acet_gc_gfas
                        gvoc_threegfas = acet_gc_threegfas
                        gvoc_nobb = acet_gc_nobb
                        title = 'Acetone'
                        xrange = [0,4]
                        ;if keyword_set(wus) then xrange = [0,3]  
                        if keyword_set(wus) then xrange = [0,4]  

                        if keyword_set(filter2_nobb) or keyword_set(filter3_nobb)  or keyword_set(filter5_nobb) then xrange = [0,2.5]    
                        ;if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) and keyword_set(wus) then xrange = [0,4]                            
                    
                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_acet_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_acet_mod*scalefactor_mod*1e-15
                        IF keyword_set(OHrplot) then xrange = [0,0.03]
 

                    end
                    ;;c3h8
                    3:begin
                        ovoc_ptr =  Propane_obs_toga
                        ovoc_other = Propane_obs_toga ;Propane_obs_iwas, iwas has less datapoints
                        gvoc_gfas = c3h8_gc_gfas
                        gvoc_threegfas = c3h8_gc_threegfas
                        gvoc_nobb = c3h8_gc_nobb
                        
                        title = 'Propane'
                        xrange = [0,3]
                        if keyword_set(wus) then xrange = [0,0.35]

                        ;if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) then xrange = [0,0.6]    
                        if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) and keyword_set(wus) then xrange = [0,0.12]    
    
                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_c3h8_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_c3h8_mod*scalefactor_mod*1e-15
                        IF keyword_set(OHrplot) then xrange = [0,0.01]

                    end
                    ;;mek
                    4:begin
                        ovoc_ptr        = MEK_obs_noaaptr
                        ovoc_other      = MEK_obs_toga
                        gvoc_gfas = mek_gc_gfas
                        gvoc_threegfas = mek_gc_threegfas
                        gvoc_nobb = mek_gc_nobb
                        
                        title = 'MEK'
                        xrange = [0,2]
                        if keyword_set(wus) then xrange = [0,0.4]  

                        ;if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) then xrange = [0,0.3]   
                        if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) and keyword_set(wus) then xrange = [0,0.3]    
                        
                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_mek_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_mek_mod*scalefactor_mod*1e-15
                        if keyword_set(OHrplot) then xrange = [0,0.02]  
                    end
                    ;;eoh
                    5:begin
                        ovoc_ptr = C2H5OH_obs_noaaptr
                        ovoc_other = C2H5OH_obs_toga
                        gvoc_gfas = c2h5oh_gc_gfas
                        gvoc_threegfas = c2h5oh_gc_threegfas
                        gvoc_nobb = c2h5oh_gc_nobb
                        title = 'Ethanol'
                        xrange = [0,8]
                        if keyword_set(wus) then xrange = [0,6]  
                        if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) then xrange = [0,2]   

                        ; waiting to be corrected for EOH
                        koh_voc_obs = tmp_koh_co_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_co_mod*scalefactor_mod*1e-15
                        if keyword_set(OHrplot) then xrange = [0,0.05]

                    end
                    ;;benz
                    6:begin
                        ovoc_ptr = Benzene_obs_noaaptr
                        ovoc_other = Benzene_obs_toga
                        gvoc_gfas = benz_gc_gfas
                        gvoc_threegfas = benz_gc_threegfas
                        gvoc_nobb = benz_gc_nobb
                        
                        title = 'Benzene'
                        xrange = [0,0.1]
                        if keyword_set(wus) then xrange = [0,0.2]  

                        if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) then xrange = [0,0.05]    
                        ;if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) and keyword_set(wus) then xrange = [0,0.1]   
                        
                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_benz_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_benz_mod*scalefactor_mod*1e-15
                        if keyword_set(OHrplot) then xrange = [0,0.003]

                    end
                    ;;tolu
                    7:begin
                        ovoc_ptr = Toluene_obs_noaaptr
                        ovoc_other = Toluene_obs_toga
                        gvoc_gfas = c7h8_gc_gfas
                        gvoc_threegfas = c7h8_gc_threegfas
                        gvoc_nobb = c7h8_gc_nobb
                        
                        title = 'Toluene'
                        xrange = [0,0.1]
                        if keyword_set(wus) then xrange = [0,0.2] 

                        if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) then xrange = [0,0.012]    
                        ;if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) and keyword_set(wus) then xrange = [0,0.1]    

                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_tolu_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_tolu_mod*scalefactor_mod*1e-15
                        if keyword_set(OHrplot) then xrange = [0,0.01]

                    end  
                    ;;xyle, either was or ptr would be fine, ptr measures C8 aromatics
                    8:begin
                        ;ovoc_was = mpXylene_obs_was + oXylene_obs_was
                        ;ovoc_iwas = mpXylene_obs_iwas + oXylene_obs_iwas
                        ovoc_ptr       = Xylenes_obs_noaaptr
                        ovoc_other     = Xylenes_obs_toga
                        gvoc_gfas = c8h10_gc_gfas
                        gvoc_threegfas = c8h10_gc_threegfas
                        gvoc_nobb = c8h10_gc_nobb

                        title = 'Xylene'
                        xrange = [0, 0.05]
                        ;if keyword_set(wus) then xrange = [0,0.05]
                        if keyword_set(wus) then xrange = [0,0.03]
                        
                        if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) then xrange = [0,0.0025]    
                        ;if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) and keyword_set(wus) then xrange = [0,0.05]    

                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_xyle_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_xyle_mod*scalefactor_mod*1e-15
                        if keyword_set(OHrplot) then xrange = [0,0.02]
                    end
                    ;;CO
                    9:begin
                        ovoc_ptr = co_obs_dacom
                        ovoc_other = co_obs_lgr
                        gvoc_gfas = co_gc_gfas
                        gvoc_threegfas = co_gc_threegfas
                        gvoc_nobb = co_gc_nobb
                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_co_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_co_mod*scalefactor_mod*1e-15
            
                        title = 'CO'
                        
                        xrange = [0,200]
                        
                        if keyword_set(wus) then xrange = [0,300]
                        if keyword_set(OHrplot) then xrange = [0,5]
                        ;if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) and keyword_set(wus) then xrange = [0,200]    
                    end
                    ;;ethane
                    10:begin
                        ovoc_ptr = c2h6_obs_cams; we also have was measurements
                        ovoc_other = c2h6_obs_was; we also have was measurements
                        gvoc_gfas = c2h6_gc_gfas
                        gvoc_threegfas = c2h6_gc_threegfas
                        gvoc_nobb = c2h6_gc_nobb
                        
                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_c2h6_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_c2h6_mod*scalefactor_mod*1e-15
            
                        title = 'C2H6'
                        xrange = [0,4]
                        if keyword_set(wus) then xrange = [0,3]
                        ;if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) and keyword_set(wus) then xrange = [0,4]    

                    end
                    
                    ;;ALK4
                    11:begin
                        ovoc_ptr =  alk4_obs_iwas
                        ovoc_other = alk4_obs_was
                        gvoc_gfas = alk4_gc_gfas
                        gvoc_threegfas = alk4_gc_threegfas
                        gvoc_nobb = alk4_gc_nobb
                        
                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_alk4_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_alk4_mod*scalefactor_mod*1e-15
            
                        title = 'ALK4'
                        xrange = [0,2]
                        if keyword_set(wus) then xrange = [0,4]
                        ;if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) and keyword_set(wus) then xrange = [0,1.5]    

                        
                    end
                    
                    ;;RCHO
                    12:begin                    
                        ovoc_ptr = rcho_obs_toga
                        ovoc_other = rcho_obs_toga
                        gvoc_gfas = rcho_gc_gfas
                        gvoc_threegfas = rcho_gc_threegfas
                        gvoc_nobb = rcho_gc_nobb

                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_rcho_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_rcho_mod*scalefactor_mod*1e-15
                        
                        title = 'RCHO'
                        xrange = [0,10]
                        if keyword_set(wus) then xrange = [0,0.3]
                        if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) and keyword_set(wus) then xrange = [0,0.4]    

                        if keyword_set(OHrplot) then xrange = [0,0.15]


                    end
                    
                    ;;PRPE
                    13:begin
                        ovoc_ptr = prpe_obs_was
                        ovoc_other = prpe_obs_iwas
                        gvoc_gfas = prpe_gc_gfas
                        gvoc_threegfas = prpe_gc_threegfas
                        gvoc_nobb = prpe_gc_nobb

                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_prpe_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_prpe_mod*scalefactor_mod*1e-15
            
                        title = 'PRPE'
                        xrange = [0,1]
                        if keyword_set(wus) then xrange = [0,1]
                        ;if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) and keyword_set(wus) then xrange = [0,10]    

                    end
                    ;;PRPE
                    14:begin
                        ovoc_ptr = prpe_obs_was
                        ovoc_other = prpe_obs_iwas
                        gvoc_gfas = prpe_gc_gfas
                        gvoc_threegfas = prpe_gc_threegfas
                        gvoc_nobb = prpe_gc_nobb

                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_prpe_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_prpe_mod*scalefactor_mod*1e-15
            
                        title = 'PRPE'
                        xrange = [0,10]
                        if keyword_set(wus) then xrange = [0,10]
                        ;if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) and keyword_set(wus) then xrange = [0,0.5]    

                    end
                    ;;Formic acid
                    15:begin
                        ovoc_ptr = HCOOH_obs_cims
                        ovoc_other = HCOOH_obs_noaaptr
                        gvoc_gfas = hcooh_gc_gfas
                        gvoc_threegfas = hcooh_gc_threegfas
                        gvoc_nobb = hcooh_gc_nobb

                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_hcooh_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_hcooh_mod*scalefactor_mod*1e-15
            
                        title = 'HCOOH'
                        xrange = [0,15]
                        if keyword_set(wus) then xrange = [0,6.2]
                        if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) and keyword_set(wus) then xrange = [0,6]    

                        if keyword_set(OHrplot) then xrange = [0,0.15]
                    end
;; test!!!!!
                    ;;Acetic acid
                    16:begin
                        ovoc_ptr = GlycolaldehydeCH3COOH_obs_noaaptr
                        ovoc_other = GlycolaldehydeCH3COOH_obs_noaaptr
                        gvoc_gfas = acta_gc_gfas
                        gvoc_threegfas = acta_gc_threegfas
                        gvoc_nobb = acta_gc_nobb

                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_acta_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_acta_mod*scalefactor_mod*1e-15
            
                        title = 'Acetic'
                        xrange = [0,10]
                        if keyword_set(wus) then xrange = [0,2.5]
                        if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) and keyword_set(wus) then xrange = [0,2.0]    

                        
                        if keyword_set(OHrplot) then xrange = [0,0.15]

                    end
                    ;;RCHO
                    17:begin                    
                        ovoc_ptr = rcho_obs_toga
                        ovoc_other = rcho_obs_toga
                        gvoc_gfas = rcho_gc_gfas
                        gvoc_threegfas = rcho_gc_threegfas
                        gvoc_nobb = rcho_gc_nobb

                        ;; OH reactivity
                        koh_voc_obs = tmp_koh_rcho_obs*scalefactor_obs*1e-15
                        koh_voc_mod = tmp_koh_rcho_mod*scalefactor_mod*1e-15
                        
                        title = 'RCHO'
                        xrange = [0,10]
                        if keyword_set(wus) then xrange = [0,0.3]
                        if keyword_set(filter2_nobb) or keyword_set(filter3_nobb) or keyword_set(filter5_nobb) and keyword_set(wus) then xrange = [0,0.4]    

                        if keyword_set(OHrplot) then xrange = [0,0.15]

                    end
                    
                endcase
            endif
            print, 'processing ' + title
            
;///////////////////////////////////////////////////////////////////////////////////////////
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;++                                      STEP 3                                           ++
;++                     Setting filters:NA/bb emission/fresh/aged                         ++
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;///////////////////////////////////////////////////////////////////////////////////////////

    
        ;for vertical profile
            zvar_obs  = alt_obs
            zvar_gc_gfas   = alt_gc_gfas
            zvar_gc_threegfas   = alt_gc_threegfas
            zvar_gc_nobb  = alt_gc_nobb

    ;test,lxu
    if keyword_set(prs) then begin
            zvar_obs  = prs_obs
            zvar_gc_gfas   = prs_obs
            zvar_gc_threegfas   = prs_obs
            zvar_gc_nobb   = alt_gc_nobb

            
    endif
    ;for time series, reset the time array every loop (diff)
            xtime_obs  = doy_obs
            xtime_gc   = doy_gc_gfas

    ;for voc .vs. co
            oco_tmp=co_obs_dacom
            gco_gfas=co_gc_gfas
            gco_threegfas=co_gc_threegfas
            gco_nobb =co_gc_nobb

            

    ;for filters
            no2_filter=no2_obs_canoe
            o3_filter=o3_obs_roze
            co_filter=co_obs_dacom
            ch3cn_filter=CH3CN_obs_noaaptr
            hcn_filter = HCN_obs_noaaptr
            acn_co_filter = CH3CN_obs_noaaptr/co_obs_dacom*1000 ; convert it into ppb/ppm
            
    ;;tracks
            lat_gc = lat_gc_gfas
            lon_gc = lon_gc_gfas
    
    ;;setting tmp value
            tmp_lat_obs = lat_obs
            tmp_lon_obs = lon_obs
            tmp_lat_gc  = lat_gc
            tmp_lon_gc  = lon_gc

    ;; For removing missing values!!!!!!!    
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
    ;;Values out of western US and eastern US
            ;if keyword_set(WestandEast) then begin
                if keyword_set(wus) then begin
                    ind_lon = where(lon_obs ge -105, ct_lon)
                    if ct_lon gt 0 then begin
                        ovoc_ptr[ind_lon] = !VALUES.F_NAN
                        ovoc_other[ind_lon] = !VALUES.F_NAN

                        gvoc_gfas[ind_lon] = !VALUES.F_NAN
                        gvoc_threegfas[ind_lon] = !VALUES.F_NAN
                        gvoc_nobb[ind_lon] = !VALUES.F_NAN

                    endif
                endif
                if keyword_set(sus) then begin
                    ind_lon = where(lon_obs le -105, ct_lon)
                    if ct_lon gt 0 then begin
                        ovoc_ptr[ind_lon] = !VALUES.F_NAN
                        ovoc_other[ind_lon] = !VALUES.F_NAN
                        gvoc_gfas[ind_lon] = !VALUES.F_NAN
                        gvoc_threegfas[ind_lon] = !VALUES.F_NAN
                        gvoc_nobb[ind_lon] = !VALUES.F_NAN

                    endif
                endif
            ;endif


;;=================
;;plume filters
;;=================
            if plume eq 1 then begin
    ;; notice that we don't use doy to remove data
                ;;filter1:filter STE, urban plumes
                if keyword_set(filter1_other) then begin
                    ind = where(no2_filter ge no2_thresh or $
                        o3_filter/co_filter ge strat_thresh, ct)
                    if ct gt 0 then begin
                        ovoc_ptr[ind] = !VALUES.F_NAN
                        ovoc_other[ind] = !VALUES.F_NAN
                        gvoc_gfas[ind] = !VALUES.F_NAN
                        gvoc_threegfas[ind] = !VALUES.F_NAN
                        gvoc_nobb[ind] = !VALUES.F_NAN
                    endif
                    ;print, 'this is the data points to be deleted for flitler1', ct
                endif
                ;;filter2:use 25 percentile of acetonitrile
                if keyword_set(filter2_bb) then begin
                    ind = where(ch3cn_filter le ch3cn_thresh,ct)
                    if ct gt 0 then begin
                        ovoc_ptr[ind] = !VALUES.F_NAN
                        ovoc_other[ind] = !VALUES.F_NAN
                        gvoc_gfas[ind] = !VALUES.F_NAN
                        gvoc_threegfas[ind] = !VALUES.F_NAN
                        gvoc_nobb[ind] = !VALUES.F_NAN
                    endif
                    ;print, 'this is the data points to be deleted for flitler2', ct
                endif
                if keyword_set(filter2_nobb) then begin
                    ind = where(ch3cn_filter ge ch3cn_thresh,ct)
                    if ct gt 0 then begin
                        ovoc_ptr[ind] = !VALUES.F_NAN
                        ovoc_other[ind] = !VALUES.F_NAN
                        gvoc_gfas[ind] = !VALUES.F_NAN
                        gvoc_threegfas[ind] = !VALUES.F_NAN
                        gvoc_nobb[ind] = !VALUES.F_NAN
                    endif
                    ;print, 'this is the data points to be deleted for flitler2', ct
                endif
                ;;fitler3: use acn/co ratio (2.01ppb/ppm) lt the threshold: nonbb (doesn't work well)
                if keyword_set(filter3_bb) then begin
                    ind = where(acn_co_filter le acn_co_tresh,ct)
                    if ct gt 0 then begin
                        ovoc_ptr[ind] = !VALUES.F_NAN
                        ovoc_other[ind] = !VALUES.F_NAN
                        gvoc_gfas[ind] = !VALUES.F_NAN
                        gvoc_threegfas[ind] = !VALUES.F_NAN
                        gvoc_nobb[ind] = !VALUES.F_NAN
                    endif
                    ;print, 'this is the data points to be deleted for flitler3', ct
                endif
                if keyword_set(filter3_nobb) then begin
                    ind = where(acn_co_filter ge acn_co_tresh,ct)
                    if ct gt 0 then begin
                        ovoc_ptr[ind] = !VALUES.F_NAN
                        ovoc_other[ind] = !VALUES.F_NAN
                        gvoc_gfas[ind] = !VALUES.F_NAN
                        gvoc_threegfas[ind] = !VALUES.F_NAN
                        gvoc_nobb[ind] = !VALUES.F_NAN
                    endif
                    ;print, 'this is the data points to be deleted for flitler3', ct
                endif
                if keyword_set(filter4_bb) then begin
                    ind = where(hcn_filter le hcn_thresh,ct)
                    if ct gt 0 then begin
                        ovoc_ptr[ind] = !VALUES.F_NAN
                        ovoc_other[ind] = !VALUES.F_NAN
                        gvoc_gfas[ind] = !VALUES.F_NAN
                        gvoc_threegfas[ind] = !VALUES.F_NAN
                        gvoc_nobb[ind] = !VALUES.F_NAN
                    endif
                    ovoc[ind] = !values.f_nan
                endif
                if keyword_set(filter4_nobb) then begin
                    ind = where(hcn_filter ge hcn_thresh,ct)
                    if ct gt 0 then begin
                        ovoc_ptr[ind] = !VALUES.F_NAN
                        ovoc_other[ind] = !VALUES.F_NAN
                        gvoc_gfas[ind] = !VALUES.F_NAN
                        gvoc_threegfas[ind] = !VALUES.F_NAN
                        gvoc_nobb[ind] = !VALUES.F_NAN
                    endif
                    ;print, 'this is the data points to be deleted for flitler4', ct
                endif
                if keyword_set(filter5_nobb) then begin
                    ind = where(co_filter ge co_thresh,ct)
                    if ct gt 0 then begin
                        ovoc_ptr[ind] = !VALUES.F_NAN
                        ovoc_other[ind] = !VALUES.F_NAN
                        gvoc_gfas[ind] = !VALUES.F_NAN                    
                        gvoc_threegfas[ind] = !VALUES.F_NAN                    
                        gvoc_nobb[ind] = !VALUES.F_NAN                    
                    endif
                    ;print, 'this is the data points to be deleted for flitler5', ct
                endif                     
           endif
           

;///////////////////////////////////////////////////////////////////////////////////////////
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;++                                      STEP 4                                           ++
;++                               Removing missing value                                  ++
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;///////////////////////////////////////////////////////////////////////////////////////////
            if keyword_set(codeployed) then remove_ind = where(finite(ovoc_ptr) and finite(ovoc_other)) else remove_ind = where(finite(ovoc_ptr))
            if ct gt 0 then begin
                ovoc_ptr = ovoc_ptr[remove_ind]
                ovoc_other = ovoc_other[remove_ind]
                gvoc_gfas = gvoc_gfas[remove_ind]
                gvoc_threegfas = gvoc_threegfas[remove_ind]
                gvoc_nobb = gvoc_nobb[remove_ind]

                zvar_obs = zvar_obs[remove_ind]
                zvar_gc_gfas = zvar_gc_gfas[remove_ind]
                zvar_gc_threegfas = zvar_gc_threegfas[remove_ind]
                zvar_gc_nobb = zvar_gc_nobb[remove_ind]
                
                koh_voc_obs = koh_voc_obs[remove_ind]
                koh_voc_mod = koh_voc_mod[remove_ind]
            endif
           
           ;remove_ind_noaaptr = where(finite(ovoc_noaaptr), ct1) ; ALD2 (similar, less data), EOH (PTR), BENZ (PTR), TOLU (PTR), XYLE (PTR)
           ;remove_ind_cams = where(finite(ovoc_cams), ct2) ;CH2O(higher than ISAF, less data)
           ;remove_ind_isaf = where(finite(ovoc_isaf), ct3) ;CH2O(less than CAMS, more data)
           ;remove_ind_iwas = where(finite(ovoc_iwas), ct4) 
           ;remove_ind_toga = where(finite(ovoc_toga), ct5) ;ALD2 (similar, more data), ACET & C3H8 &MEK(TOGA, but toga is a disaster, may try IWAS or WAS), EOH, BENZ, TOLU, XYLE
           ;remove_ind_uioptr = where(finite(ovoc_uioptr), ct6) 
           ;remove_ind_was = where(finite(ovoc_was), ct7) ; MEK (TOGA is messy, try WAS or IWAS), BENZ, TOLU, XYLE (but much less temporal resolution)
           ;print, 'This is the avaialbe data left...'
           ;print, ct1, ct2, ct3, ct4, ct5, ct6, ct7
           ;continue


;///////////////////////////////////////////////////////////////////////////////////////////
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;++                                      STEP 5                                           ++
;++                                  DO THE PLOTTING                                      ++
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;///////////////////////////////////////////////////////////////////////////////////////////    


    ;; Vertical profiles
            zvar_obs   = zvar_obs
            zvar_gc    = zvar_gc_gfas

            if keyword_set(prs) then begin
                ;if keyword_set(OHrplot) then begin
   
                    oVert_ptr = bin_vert(Data=ovoc_ptr, ZData = zvar_obs, $
                                        ZEdges=PEdges,/Press)
                    oVert_other = bin_vert(Data=ovoc_other, ZData = zvar_obs, $
                                        ZEdges=PEdges,/Press)
                                        
                                        
                    gVert_gfas  = bin_vert(Data=gvoc_gfas, ZData=zvar_gc, $
                                    ZEdges=PEdges,/Press)
                    gVert_threegfas  = bin_vert(Data=gvoc_threegfas, ZData=zvar_gc, $
                                    ZEdges=PEdges,/Press)
                    gVert_nobb  = bin_vert(Data=gvoc_nobb, ZData=zvar_gc, $
                                    ZEdges=PEdges,/Press)
                ;endelse
            endif
            
; not sure why it contains diff but will ignore it for now.            
;if s eq 12 then begin
;    print, total(ovoc_ptr - ovoc_other,/nan)
;    print, total(oVert_ptr.DataMed - oVert_other.DataMed,/nan)
;    stop
;endif
            ;; charsize
            if keyword_set(all) or keyword_set(ovocs) or keyword_set(SI) then charsize = 1

            ;; err_thick
            if keyword_set(all) or keyword_set(ovocs) and not keyword_set(checkGFAS) or keyword_set(SI) then err_thick = 4
            if keyword_set(each)  then err_thick = 4
            
            if keyword_set(all) or keyword_set(ovocs) or keyword_set(SI) or keyword_set(obs_only) then begin
                charsize = 2.5
                charthick = 4
                thick1 = 20
                thick2 = 12
            endif
            
            if keyword_set(each) then begin
                charsize = 1.8
                charthick = 4
                thick1 = 16
                thick2 = 12
            endif
            if keyword_set(chemistry) then begin
                charsize = 1.8
                charthick = 4
                thick1 = 16
                thick2 = 12
            endif
    ;; Plotting part            
            if keyword_set(average) then begin
                obs_ptr = oVert_ptr.DataMean
                obs_other = oVert_other.DataMean
                
                gfas  = gVert_gfas.DataMean
                threegfas  = gVert_threegfas.DataMean
                nobb  = gVert_nobb.DataMean
            endif
            if keyword_set(med) then begin
                obs_ptr = oVert_ptr.DataMed
                obs_other = oVert_other.DataMed
                
                gfas  = gVert_gfas.DataMed             
                threegfas  = gVert_threegfas.DataMed             
                nobb  = gVert_nobb.DataMed             
            endif  



            ;; OUTPUTS
            ;  alt
            z_obs_ptr = oVert_ptr.zMean
            z_obs_other = oVert_other.zMean
            
            z_gfas = gVert_gfas.zMean
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

            ;; gfas
            gfas_q10 = gVert_gfas.DataQ10
            gfas_q90 = gVert_gfas.DataQ90
            gfas_q25 = gVert_gfas.DataQ25
            gfas_q75 = gVert_gfas.DataQ75
            gfas_numpts = gVert_gfas.NumPts
            
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

            if plume eq 0 then yrange = [1000,400]
            if plume eq 1 then yrange = [1000,400]
            

            
            ;; PTR
            if plume eq 0 then begin
                numpt_thresh = 10
            endif else begin
                numpt_thresh = 10
            endelse
            
            ind_num = where(obs_ptr_numpts lt numpt_thresh , ct)
            if ct gt 0 then begin                    
                obs_ptr[ind_num] = !VALUES.F_NAN
                obs_other[ind_num] = !VALUES.F_NAN
            endif


            
            ind =  where(finite(obs_ptr),ct)
            if ct gt 0 then begin
                obs_ptr = obs_ptr[ind]
                obs_other = obs_other[ind]
                gfas  = gfas[ind]
                threegfas  = threegfas[ind]
                nobb  = nobb[ind]

                z_obs_ptr = z_obs_ptr[ind]
                z_obs_other = z_obs_other[ind]
                
                z_gfas = z_gfas[ind]
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
                
                gfas_q10 = gfas_q10[ind]
                gfas_q90 = gfas_q90[ind]
                gfas_q25 = gfas_q25[ind]
                gfas_q75 = gfas_q75[ind]
                gfas_numpts = gfas_numpts[ind]
                
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

            if (keyword_set(obs_major) or plume eq 1) and s eq 15 then obs_ptr = obs_other ;; set for HCOOH

            ;if keyword_set(OHrplot) then begin
                if keyword_set(all) or keyword_set(each) or keyword_set(SI) or keyword_set(obs_only) then begin
                    if s eq 0 or s eq 3 or s eq 6 or s eq 9 or s eq 12 or s eq 15 then begin                
                        if keyword_set(errorbar) then begin
                            if s eq 0 and keyword_set(obs_only) then begin ;; for formaldehyde
                                ;; coyote library by using error bar to represent percentile lines
                                ;; sources: http://www.idlcoyote.com/idldoc/cg/cgplot.html
                                cgplot,obs_ptr,z_obs_ptr,col=12,xrange=xrange,yrange=yrange,$;, title=title,,$; title=title,$
                                ;xtitle=xtitle,$;,ytitle='Altitude[km]',$     
                                thick=thick1,charsize=charsize,charthick=charthick, $
                                ERR_XLow=(obs_ptr-obs_ptr_q25), ERR_XHigh=(obs_ptr_q75-obs_ptr), $
                                ERR_Color='black', $
                                ERR_THICK = ERR_THICK, $
                                /ERR_CLIP
                                ;YTICKFORMAT="(A1)", XTICKFORMAT="(A1)"
                                ;,ystyle=4,$
                            endif else begin
                        
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
                            endelse
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
                    if (s ne 3) and (s ne 0) then oplot, obs_other, z_obs_other, col=3,thick=thick1,line=0
                    if (s eq 0) then oplot, obs_other, z_obs_other, col=4,thick=thick1,line=0
                    if (s eq 3) then oplot, obs_other, z_obs_other, col=9,thick=thick1,line=0
                endif
                
                if not keyword_set(obs_only) then begin
                    ;; only set sepcies implemented in the model
                    oplot,gfas,z_gfas,col=4,thick=thick1,line=0

                    if not keyword_set(filter2_nobb) then begin
                        oplot,nobb,z_nobb,col=7,thick=thick1,line=1
                        oplot,threegfas,z_threegfas,col=12,thick=thick1,line=1
                    endif
                endif
                
                

            if keyword_set(all) or keyword_set(ovocs) or keyword_set(SI) or keyword_set(obs_only) then begin
                charsize =1.1
                charthick =3
            endif 
            
            

            ;; change the postition later.
            ;; xyouts for NumPoints for each layers            
            for num = 0, n_elements(obs_ptr_numpts)-1 do begin
                if z_obs_ptr[num] ge yrange[0] then continue            
                xyouts, 1.015*max(xrange),z_obs_ptr[num], STRTRIM(string(obs_ptr_numpts[num]), 1), /data, col=1, charsize=charsize, $
                    charthick=charthick, alignment=0.      
            endfor
            
            ;if dd eq 0 then begin
                spec_obs_total   = [0, obs_ptr] 
                spec_gfas_total  = [1, gfas]
                spec_pres        = [2, z_obs_ptr]
                print, 'currently saving out data...'
                if plume eq 0 then write_csv, './tmp/FIREX-AQ_full_' + title + '.csv' ,  spec_pres, spec_obs_total, spec_gfas_total
                if plume eq 1 then write_csv, './tmp/FIREX-AQ_nobb_' + title + '.csv' ,  spec_pres, spec_obs_total, spec_gfas_total
            ;endif
            
        endfor ;; s loop
    close_device
    ;DEVICE, /CLOSE
    print,'done!'
end
