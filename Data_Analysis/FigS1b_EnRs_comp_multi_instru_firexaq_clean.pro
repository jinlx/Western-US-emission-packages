;; Lixu, 09/30/2021, this script is used to calculate the enhancement ratio of VOCs vs CO (in V1)
;; In Version 2, we are trying to add hcho VS no to show the ozone scheme.
;; This version used emission passes with sorting order
@org_boot
pro EnRs_comp_multi_instru_firexaq_clean,$
        plume=plume,$
        emipass=emipass, filter2_nobb=filter2_nobb,filter3_nobb=filter3_nobb,$
        all=all,$
        test=test, save=save,$
        nested=nested,fbf=fbf,$
        sus=sus, wus=wus, $
        raw_60s=raw_60s, avg_5min=avg_5min, TOGA_merge=TOGA_merge
        
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
    if keyword_set(test) then fi = './test_emipass'  
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
                    
    if keyword_set(TOGA_merge) then     dates_all=  ['20190722', '20190724', '20190725', '20190729', '20190730',$
                                                        '20190802', '20190803', '20190806', '20190807', '20190808', $
                                                        '20190812', '20190813', '20190815', '20190816', '20190819', $
                                                        '20190821', '20190823', '20190826', '20190830', '20190831', $ 
                                                        '20190903']

    avo= 6.022e23  ;; molec/mol
    
;; emission passes time stamp from Wade
    dates_emipass = ['20190725','20190725','20190729','20190730','20190802','20190803','20190803',$
                        '20190803', '20190806',$
                        '20190806','20190808',$
                        '20190812','20190813',$
                        '20190816','20190830']
    
    interval = 15
                        
    start_time_emipass = [82080, 85620, 84060, 9600+86400, 84000, 4380+86400,80520,$
                        2460+86400, 84120,$
                        84480, 6540+86400,$
                        1080+86400, 5520+86400,$
                        2520+86400, 61860] - interval*60
    
    end_time_emipass = [82080, 85620, 84060, 9600+86400, 84000, 4380+86400, 80520,$
                        2460+86400, 84120,$
                        84480, 6540+86400,$
                        1080+86400, 5520+86400,$
                        2520+86400, 61860] + interval*60

    if keyword_set(emipass) then begin
        dates_all = dates_emipass
        enter_time = start_time_emipass
        exit_time = end_time_emipass
        help, dates_all
    endif

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
            if keyword_set(avg_5min) then dc8fi = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                                                    'mrg2sav/FIREX-AQ/firexaq_avg/firexaq-mrg5m-dc8_merge_'+dates[n]+'_R1.sav'
            if keyword_set(raw_60s) then dc8fi = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                                                    'mrg2sav/FIREX-AQ/firexaq_raw/firexaq-mrg60-dc8_merge_'+dates[n]+'_R1.sav'
                                                    
            if keyword_set(TOGA_merge) then dc8fi = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                                                    'mrg2sav/FIREX-AQ/TOGA_merges_raw/firexaq-mrgTOGA-dc8_merge_'+dates[n]+'_R1.sav'
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
                        ovoc_ptr        = ch2o_obs_cams
                        ovoc_other      = ch2o_obs_isaf
                        title = 'Formaldehyde'
                        if plume eq 0 then begin
                            xrange = [0, 120]
                            yrange = [0, 120]
                        endif
                        if plume eq 1 then begin
                            xrange = [0, 15]
                            yrange = [0, 15]
                        endif
                    end
                    ;;ald2
                    1:begin
                        ovoc_ptr       = CH3CHO_obs_noaaptr
                        ovoc_other     = CH3CHO_obs_toga
                        title = 'Acetaldehyde'
                        if plume eq 0 then begin
                            xrange = [0, 60]
                            yrange = [0, 60]
                        endif
                        if plume eq 1 then begin
                            xrange = [0, 1.5]
                            yrange = [0, 1.5]
                        endif
                    end
                    ;;acet
                    2:begin
                        ovoc_ptr       = Acetone_obs_noaaptr
                        ovoc_other     = Acetone_obs_toga
                        ;ovoc_other     = Acetone_obs_iwas
                        title = 'Acetone'
                        if plume eq 0 then begin
                            xrange = [0, 30]
                            yrange = [0, 30]
                        endif
                        if plume eq 1 then begin
                            xrange = [0, 8]
                            yrange = [0, 8]
                        endif
                    end
                    ;;formic acid
                    3:begin
                        ovoc_ptr       = HCOOH_obs_noaaptr
                        ovoc_other     = HCOOH_obs_cims
                        title = 'Formic acid'
                        if plume eq 0 then begin
                            xrange = [0, 20]
                            yrange = [0, 20]
                        endif
                        if plume eq 1 then begin
                            xrange = [0, 10]
                            yrange = [0, 10]
                        endif
                    end
                    ;;mek
                    4:begin
                        ovoc_ptr        = MEK_obs_noaaptr
                        ovoc_other      = MEK_obs_toga
                        ;ovoc_other      = MEK_obs_iwas
                        title = 'MEK'
                        if plume eq 0 then begin
                            xrange = [0, 5]
                            yrange = [0, 5]
                        endif
                        if plume eq 1 then begin
                            xrange = [0, 2.0]
                            yrange = [0, 2.0]
                        endif
                    end
                    ;;eoh
                    5:begin                    
                        ovoc_ptr        = C2H5OH_obs_noaaptr
                        ovoc_other      = C2H5OH_obs_toga

                        title = 'Ethanol'
                        if plume eq 0 then begin
                            xrange = [0, 5.0]
                            yrange = [0, 5.0]
                        endif
                        if plume eq 1 then begin
                            xrange = [0, 1]
                            yrange = [0, 1]
                        endif
                    end
                    ;;benz
                    6:begin
                        ovoc_ptr        = Benzene_obs_noaaptr
                        ovoc_other      = Benzene_obs_toga
                        ;ovoc_other      = Benzene_obs_iwas
                        title = 'Benzene'
                        if plume eq 0 then begin
                            xrange = [0, 15]
                            yrange = [0, 15]
                        endif
                        if plume eq 1 then begin
                            xrange = [0, 0.5]
                            yrange = [0, 0.5]
                        endif
                    end
                    ;;tolu
                    7:begin
                        ovoc_ptr        = Toluene_obs_noaaptr
                        ovoc_other      = Toluene_obs_toga
                        ;ovoc_other      = Toluene_obs_iwas
                        title = 'Toluene'
                        if plume eq 0 then begin
                            xrange = [0, 10]
                            yrange = [0, 10]
                        endif
                        if plume eq 1 then begin
                            xrange = [0, 0.2]
                            yrange = [0, 0.2]
                        endif
                    end  
                    ;;xyle
                    8:begin
                        ovoc_ptr       = Xylenes_obs_noaaptr
                        ovoc_other     = Xylenes_obs_toga
                        ;ovoc_ptr       = C8Aromatics_obs_noaaptr
                        ;ovoc_other     = C8Aromatics_obs_toga
                        ;ovoc_other     = mpXylene_obs_iwas + oXylene_obs_iwas
                        title = 'Xylene'
                        if plume eq 0 then begin
                            xrange = [0, 3.0]
                            yrange = [0, 3.0]
                        endif
                        if plume eq 1 then begin
                            xrange = [0, 0.05]
                            yrange = [0, 0.05]
                        endif
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
    ;for voc .vs. co
            oco_tmp=co_obs_dacom
    ;for filters
            no2_filter=no2_obs_cl
            o3_filter=o3_obs_cl
            co_filter=co_obs_dacom
            ch3cn_filter=CH3CN_obs_noaaptr
            hcn_filter = HCN_obs_noaaptr
            acn_co_filter = CH3CN_obs_noaaptr/co_obs_dacom*1000 ; convert it into ppb/ppm
            
    ;cloud
            tmp_rhum = rh_obs  
;;=================
;;plume filters
;;=================
            if keyword_set(wus) then begin
                remove_ind = where(lon_obs ge -105, ct)
                if ct gt 0 then begin
                    ovoc_ptr[ind] = !VALUES.F_NAN
                    ovoc_other[ind] = !VALUES.F_NAN
                endif
            endif
            if keyword_set(sus) then begin
                remove_ind = where(lon_obs le -105, ct)
                if ct gt 0 then begin
                    ovoc_ptr[ind] = !VALUES.F_NAN
                    ovoc_other[ind] = !VALUES.F_NAN
                endif
            endif

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
                charsize=2.5, charthick = 4, thick = 10, xrange=xrange, yrange=yrange, /log;,xrange=[0,200];,yrange=yrange,xrange=xrange;,$
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
















