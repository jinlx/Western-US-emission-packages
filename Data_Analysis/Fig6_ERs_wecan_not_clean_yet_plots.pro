;; modified by lixu, 10/17/2019, used for compare the profiles from wecan and gc 
;; original script mk_flight_v2, lhu
;; modicication: parameters, set missing value for observation data, cancel filter. 
;@tvmap
;; include all VOCs we would study, addup VOCs and some old setting
;; use R4 observation data, change tag names (ethanol_toga, acetdaldehyde_toga) 09/26/2020
;; all

;; this is used for compare ERs for obs and mod for each VOCs

;;.compile ERs_emipass_single_inclusion_v5
;; cloud1, /emipass
;; ERs_emipass_comp_multi, plume=0, /all,/test,/cloud1,/save,/emipass,inventories='gfed4',/nested,/onebyone
;; ERs_emipass_comp_multi, plume=0, /all,/test,/cloud1,/save,/emipass,inventories='finn',/nested,/onebyone
;; ERs_emipass_comp_multi, plume=0, /all,/test,/cloud1,/save,/emipass,inventories='gfas',/nested,/onebyone
;; ERs_emipass_comp_multi, plume=0, /all,/test,/cloud1,/save,/emipass,inventories='qfed',/nested,/onebyone

;; cloud1, whole campaign, only keep BB by use CH3CN and CO, cloud1
;; ERs_emipass_comp_multi, plume=1, /filter1, /filter2, /all,/test,/cloud1,/save,inventories='gfed4',/nested,/onebyone
;; ERs_emipass_comp_multi, plume=1, /filter1, /filter2, /all,/test,/cloud1,/save,inventories='finn',/nested,/onebyone
;; ERs_emipass_comp_multi, plume=1, /filter1, /filter2, /all,/test,/cloud1,/save,inventories='gfas',/nested,/onebyone
;; ERs_emipass_comp_multi, plume=1, /filter1, /filter2, /all,/test,/cloud1,/save,inventories='qfed',/nested,/onebyone

;; cloud1, whole campaign, only keep BB by use CH3CN and CO, cloud1, using fresh plume
;; ERs_emipass_comp_multi, plume=1, /filter2, /all,/test,/cloud1,/save,inventories='gfed4',/nested,/fresh,/onebyone
;; ERs_emipass_comp_multi, plume=1, /filter2, /all,/test,/cloud1,/save,inventories='finn',/nested,/fresh,/onebyone
;; ERs_emipass_comp_multi, plume=1, /filter2, /all,/test,/cloud1,/save,inventories='gfas',/nested,/fresh,/onebyone
;; ERs_emipass_comp_multi, plume=1, /filter2, /all,/test,/cloud1,/save,inventories='qfed',/nested,/fresh,/onebyone

;;addup
;; cloud1, /emipass
;; ERs_emipass_comp_multi, plume=0, /addup,/test,/cloud1,/save,/emipass,inventories='gfed4',/nested,/onebyone
;; ERs_emipass_comp_multi, plume=0, /addup,/test,/cloud1,/save,/emipass,inventories='finn',/nested,/onebyone
;; ERs_emipass_comp_multi, plume=0, /addup,/test,/cloud1,/save,/emipass,inventories='gfas',/nested,/onebyone
;; ERs_emipass_comp_multi, plume=0, /addup,/test,/cloud1,/save,/emipass,inventories='qfed',/nested,/onebyone

;; In version 3, we gonna try using expand the time a little bit.
;; In version 4, we add nonbb 
;; In version 5, we add more species from AWAS
;; In version 6
;; we add least fire influenced for different methods 
;; we add OVOCs/benzene (ongoing, need to tune the format)
;; we add OVOCs/toluene
;; v8: add chemistry keyword for O3, NO, NO2, PAN;; add NMHC and OVOCs keywords
;; v9: add CH3CN/CO thresh
@org_boot
pro ERs_wecan_not_clean_yet_plots,plume=plume, $
        young=young,intermed=intermed,nonsmoke=nonsmoke,aged = aged, GT4=GT4,LT4=LT4,emipass = emipass,$
        addup=addup,onebyone=onebyone,all=all, lumped=lumped, $
        COthresh_fresh=COthresh_fresh, COthresh_aged = COthresh_aged, $
        chemistry = chemistry, nmhcs=nmhcs, ovocs=ovocs,$
        test=test,$
        n_criteria=n_criteria,save=save,$
        inventories = inventories,$
        nested=nested,fbf=fbf,$
        metrics=metrics,scfactor=scfactor,$
        xeqco = xeqco, xeqbenz=xeqbenz, xeqtolu = xeqtolu, xeqisop=xeqisop, xeqacet=xeqacet, $
        filter2_nobb=filter2_nobb,filter2_bb=filter2_bb,$
        filter3_nobb=filter3_nobb,filter3_bb=filter3_bb,$
        newthresh_nobb = newthresh_nobb , newthresh_bb = newthresh_bb, $
        mod_only = mod_only

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
               
    ;delete 20180823/RF14 overboard
    ;if keyword_set(cloud1) then   dates_all=  ['20180724','20180726','20180730','20180731','20180802','20180803',$
    ;                                        '20180806','20180808','20180809','20180813','20180815','20180816','20180820',$
    ;                                        '20180826','20180828','20180906','20180910','20180913']
    
    ;delete 20180823/RF14 overboard; 20180803/RF06, 20180813/RF10, 20180816/RF12, 20180820/RF13
    if keyword_set(cloud1) then   dates_all=  ['20180724','20180726','20180730','20180731','20180802',$
                                            '20180806','20180808','20180809','20180815',$
                                            '20180826','20180828','20180906','20180910','20180913']

    if keyword_set(cloud2) then  dates_all=  ['20180724','20180802',$
                                            '20180806','20180808','20180809','20180815','20180820',$
                                            '20180826','20180828','20180906','20180910','20180913']

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
                     
                     
;; environments time stamp from Kate
;; young plume
; Read the file in and assign it to the sed_data variable;
; assign the header information to variables for use later.
    test_file = '/glade/work/lixujin/timestamp/young_timestamp.csv'
    sed_data = READ_CSV(test_file, HEADER=SedHeader, $
        N_TABLE_HEADER=1, TABLE_HEADER=SedTableHeader)
   
    dates_young = sed_data.FIELD07
    start_time_young = sed_data.FIELD03
    end_time_young = sed_data.FIELD04
    for i = 0, n_elements(dates_young) - 1 do begin
        dates_young[i] = dates_all[dates_young[i]-1]
    endfor


;; median plume
    test_file = '/glade/work/lixujin/timestamp/med_timestamp.csv'
    sed_data = READ_CSV(test_file, HEADER=SedHeader, $
        N_TABLE_HEADER=1, TABLE_HEADER=SedTableHeader)
   
    dates_med = sed_data.FIELD07
    start_time_med = sed_data.FIELD03
    end_time_med = sed_data.FIELD04
    for i = 0, n_elements(dates_med) - 1 do begin
        dates_med[i] = dates_all[dates_med[i]-1]
    endfor


;; age plume: greater than 4 days
    test_file = '/glade/work/lixujin/timestamp/oldGT4_timestamp.csv'
    sed_data = READ_CSV(test_file, HEADER=SedHeader, $
        N_TABLE_HEADER=1, TABLE_HEADER=SedTableHeader)
   
    dates_GT4days = sed_data.FIELD07
    start_time_GT4days = sed_data.FIELD03
    end_time_GT4days = sed_data.FIELD04
    for i = 0, n_elements(dates_GT4days) - 1 do begin
        dates_GT4days[i] = dates_all[dates_GT4days[i]-1]
    endfor


;; age plume: less than 4 days
    test_file = '/glade/work/lixujin/timestamp/oldLT4_timestamp.csv'
    sed_data = READ_CSV(test_file, HEADER=SedHeader, $
        N_TABLE_HEADER=1, TABLE_HEADER=SedTableHeader)
   
    dates_LT4days = sed_data.FIELD07
    start_time_LT4days = sed_data.FIELD03
    end_time_LT4days = sed_data.FIELD04
    for i = 0, n_elements(dates_LT4days) - 1 do begin
        dates_LT4days[i] = dates_all[dates_LT4days[i]-1]
    endfor

;; aged plumes, combine LT4 and GT4
    test_file = '/glade/work/lixujin/timestamp/oldGT2_timestamp.csv'
    sed_data = READ_CSV(test_file, HEADER=SedHeader, $
        N_TABLE_HEADER=1, TABLE_HEADER=SedTableHeader)
   
    dates_aged = sed_data.FIELD07
    start_time_aged = sed_data.FIELD03
    end_time_aged = sed_data.FIELD04
    for i = 0, n_elements(dates_aged) - 1 do begin
        dates_aged[i] = dates_all[dates_aged[i]-1]
    endfor
    
;; non smoke
    test_file = '/glade/work/lixujin/timestamp/nonsmoke_timestamp.csv'
    sed_data = READ_CSV(test_file, HEADER=SedHeader, $
        N_TABLE_HEADER=1, TABLE_HEADER=SedTableHeader)
   
    dates_nonsmoke = sed_data.FIELD07
    start_time_nonsmoke = sed_data.FIELD03
    end_time_nonsmoke = sed_data.FIELD04
    for i = 0, n_elements(dates_nonsmoke) - 1 do begin
        dates_nonsmoke[i] = dates_all[dates_nonsmoke[i]-1]
    endfor
    

    if keyword_set(emipass) then begin
        dates_all = dates_emipass
        enter_time = start_time_emipass
        exit_time = end_time_emipass
        help, dates_all
    endif
    if keyword_set(young) then begin
        dates_all = dates_young
        enter_time = start_time_young
        exit_time = end_time_young
        help, dates_all

    endif
    if keyword_set(intermed) then begin
        dates_all = dates_med
        enter_time = start_time_med
        exit_time = end_time_med
        help, dates_all

    endif
    if keyword_set(GT4) then begin
        dates_all = dates_GT4days
        enter_time = start_time_GT4days
        exit_time = end_time_GT4days
        help, dates_all

    endif
    if keyword_set(LT4) then begin
        dates_all = dates_LT4days
        dates_all = dates_LT4days
        enter_time = start_time_LT4days
        exit_time = end_time_LT4days
    endif    
    
    if keyword_set(aged) then begin
        dates_all = dates_aged
        dates_all = dates_aged
        enter_time = start_time_aged
        exit_time = end_time_aged
    endif  
    
    if keyword_set(nonsmoke) then begin
        dates_all = dates_nonsmoke
        dates_all = dates_nonsmoke
        enter_time = start_time_nonsmoke
        exit_time = end_time_nonsmoke
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
;; ==============
;; GEOS-Chem DATA
;; ==============
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
        mtpo_gc_threegfas  = [0]
        
        c2h4_gc_threegfas   = [0]
        c2h6_gc_threegfas   = [0]
        
    ;; add lumped species
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
            
            ;c130fi = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
            ;    'mrg2sav/R4_merges_avg/wecan-mrg10m-c130_merge_'+dates[n]+'_R4.sav'
                
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
                    'planelog2sav/output_'+'gfas_4x5_opt_7' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'
            
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
            tmp_c2h6_gc_gfed4   = gc.c2h6*1e9/2   ; 2C 

    ;; concert hhmm(utc) into ss(utc)      
            ;print,tmp_utc_gc_gfed4
            hh = floor(tmp_utc_gc_gfed4/100)
            ind = where(hh gt 12, ct)
            if ct gt 0 then hh[ind] = hh[ind]
            ind = where(hh lt 12, ct)
            if ct gt 0 then hh[ind] = hh[ind] + 24
            
            mm = tmp_utc_gc_gfed4- floor(tmp_utc_gc_gfed4/100)*100
            tmp_utc_gc_gfed4 = float(hh)*60*60+float(mm)*60
            if keyword_set(emipass) or keyword_set(young) or $
                keyword_set(intermed) or keyword_set(GT4) or keyword_set(LT4) or $
                keyword_set(nonsmoke) or keyword_set(aged) then begin
                ind_enter_time = where(abs(tmp_utc_gc_gfed4 - enter_time[n]) le 30,ct1)
                ind_exit_time = where(abs(tmp_utc_gc_gfed4 - exit_time[n]) le 30, ct2)    
                if ct1 gt 0 and ct2 gt 0 then begin
                    tmp_co_gc_gfed4   = tmp_co_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_O3_gc_gfed4   = tmp_O3_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_pan_gc_gfed4  = tmp_pan_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_hcho_gc_gfed4 = tmp_hcho_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_acet_gc_gfed4 = tmp_acet_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_benz_gc_gfed4 = tmp_benz_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]
            ;        ch3oh_gc= [ch3oh_gc, gc.MOH*1e9]
                    tmp_ald2_gc_gfed4 = tmp_ald2_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]

                    tmp_no_gc_gfed4   = tmp_no_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_no2_gc_gfed4  = tmp_no2_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_so2_gc_gfed4  = tmp_so2_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]

                    ;tmp_na = tmp_na[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_oh_gc_gfed4  = tmp_oh_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]

                    tmp_date_gc_gfed4 = tmp_date_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_utc_gc_gfed4  = tmp_utc_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_doy_gc_gfed4  = tmp_doy_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_lat_gc_gfed4  = tmp_lat_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_lon_gc_gfed4  = tmp_lon_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_alt_gc_gfed4  = tmp_alt_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_prs_gc_gfed4  = tmp_prs_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]
             ;; add vocs
                    tmp_c3h8_gc_gfed4   = tmp_c3h8_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c3h6_gc_gfed4   = tmp_c3h6_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c2h5oh_gc_gfed4 = tmp_c2h5oh_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c5h8_gc_gfed4   = tmp_c5h8_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c7h8_gc_gfed4   = tmp_c7h8_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_dms_gc_gfed4    = tmp_dms_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_mek_gc_gfed4    = tmp_mek_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c8h10_gc_gfed4  = tmp_c8h10_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]

                    tmp_acta_gc_gfed4   = tmp_acta_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_macr_mvk_gc_gfed4 = tmp_macr_mvk_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_hcooh_gc_gfed4  = tmp_hcooh_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]

                    tmp_mtpa_gc_gfed4  = tmp_mtpa_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_limo_gc_gfed4  = tmp_limo_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_mtpo_gc_gfed4  = tmp_mtpo_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]

                ;;add lumped species
                    tmp_alk4_gc_gfed4   = tmp_alk4_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_prpe_gc_gfed4   = tmp_prpe_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_rcho_gc_gfed4   = tmp_rcho_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]

                    ;tmp_c2h4_gc_gfed4   = tmp_c2h4_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c2h6_gc_gfed4   = tmp_c2h6_gc_gfed4[ind_enter_time[0]:ind_exit_time[-1]]

                endif
                if ct1 le 0 or ct2 le 0 then stop
            endif
            help,tmp_co_gc_gfed4

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
                    'planelog2sav/output_'+'gfas_tmp' + '/mrg5m_wecan_c130_'+dates[n]+'.sav'
            endif
            if keyword_set(fbf) then begin
                gcfi_finn   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                    'planelog2sav/output_'+'finn_4x5' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'
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
            tmp_c2h6_gc_finn   = gc.c2h6*1e9/2 ; 2C
            
            ;; concert hhmm(utc) into ss(utc)      
            ;print,tmp_utc_gc_finn
            hh = floor(tmp_utc_gc_finn/100)
            ind = where(hh gt 12, ct)
            if ct gt 0 then hh[ind] = hh[ind]
            ind = where(hh lt 12, ct)
            if ct gt 0 then hh[ind] = hh[ind] + 24
            
            mm = tmp_utc_gc_finn- floor(tmp_utc_gc_finn/100)*100
            tmp_utc_gc_finn = float(hh)*60*60+float(mm)*60
                        

            
            if keyword_set(emipass) or keyword_set(young) or $
                keyword_set(intermed) or keyword_set(GT4) or keyword_set(LT4) or $
                keyword_set(nonsmoke) or keyword_set(aged) then begin
                ind_enter_time = where(abs(tmp_utc_gc_finn - enter_time[n]) le 30,ct1)
                ind_exit_time = where(abs(tmp_utc_gc_finn - exit_time[n]) le 30, ct2)    
                if ct1 gt 0 and ct2 gt 0 then begin
                    tmp_co_gc_finn   = tmp_co_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_O3_gc_finn   = tmp_O3_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_pan_gc_finn  = tmp_pan_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_hcho_gc_finn = tmp_hcho_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_acet_gc_finn = tmp_acet_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_benz_gc_finn = tmp_benz_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]
            ;        ch3oh_gc= [ch3oh_gc, gc.MOH*1e9]
                    tmp_ald2_gc_finn = tmp_ald2_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]

                    tmp_no_gc_finn   = tmp_no_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_no2_gc_finn  = tmp_no2_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_so2_gc_finn  = tmp_so2_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]

                    ;tmp_na = tmp_na[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_oh_gc_finn  = tmp_oh_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]

                    tmp_date_gc_finn = tmp_date_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_utc_gc_finn  = tmp_utc_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_doy_gc_finn  = tmp_doy_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_lat_gc_finn  = tmp_lat_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_lon_gc_finn  = tmp_lon_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_alt_gc_finn  = tmp_alt_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_prs_gc_finn  = tmp_prs_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]
             ;; add vocs
                    tmp_c3h8_gc_finn   = tmp_c3h8_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c3h6_gc_finn   = tmp_c3h6_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c2h5oh_gc_finn = tmp_c2h5oh_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c5h8_gc_finn   = tmp_c5h8_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c7h8_gc_finn   = tmp_c7h8_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_dms_gc_finn    = tmp_dms_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_mek_gc_finn    = tmp_mek_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c8h10_gc_finn  = tmp_c8h10_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]

                    tmp_acta_gc_finn   = tmp_acta_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_macr_mvk_gc_finn = tmp_macr_mvk_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_hcooh_gc_finn  = tmp_hcooh_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]

                    tmp_mtpa_gc_finn  = tmp_mtpa_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_limo_gc_finn  = tmp_limo_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_mtpo_gc_finn  = tmp_mtpo_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]

                ;;add lumped species
                    tmp_alk4_gc_finn   = tmp_alk4_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_prpe_gc_finn   = tmp_prpe_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_rcho_gc_finn   = tmp_rcho_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]

                    ;tmp_c2h4_gc_finn   = tmp_c2h4_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c2h6_gc_finn   = tmp_c2h6_gc_finn[ind_enter_time[0]:ind_exit_time[-1]]
                    
                endif
                if ct1 le 0 or ct2 le 0 then stop
            endif
            help,tmp_co_gc_finn
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
                ;gcfi_gfas   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                ;    'planelog2sav/output_'+'gfas_tmp' + '/mrg10m_wecan_c130_'+dates[n]+'.sav'   
                gcfi_gfas   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                    'planelog2sav/output_'+'gfas_tmp' + '/mrg5m_wecan_c130_'+dates[n]+'.sav'   
            endif
            if keyword_set(fbf) then begin
                gcfi_gfas   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                    'planelog2sav/output_'+'gfas_4x5' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'
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
            tmp_c2h6_gc_gfas   = gc.c2h6*1e9/2 ;2C

    ;; concert hhmm(utc) into ss(utc)      
            ;print,tmp_utc_gc_gfas
            hh = floor(tmp_utc_gc_gfas/100)
            ind = where(hh gt 12, ct)
            if ct gt 0 then hh[ind] = hh[ind]
            ind = where(hh lt 12, ct)
            if ct gt 0 then hh[ind] = hh[ind] + 24
            
            mm = tmp_utc_gc_gfas- floor(tmp_utc_gc_gfas/100)*100
            tmp_utc_gc_gfas = float(hh)*60*60+float(mm)*60
            

            if keyword_set(emipass) or keyword_set(young) or $
                keyword_set(intermed) or keyword_set(GT4) or keyword_set(LT4) or $
                keyword_set(nonsmoke) or keyword_set(aged) then begin
                ind_enter_time = where(abs(tmp_utc_gc_gfas - enter_time[n]) le 30,ct1)
                ind_exit_time = where(abs(tmp_utc_gc_gfas - exit_time[n]) le 30, ct2)                
                if ct1 gt 0 and ct2 gt 0 then begin
                    ;print, ind_enter_time[0], ind_enter_time[-1]
                    tmp_co_gc_gfas = tmp_co_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_O3_gc_gfas   = tmp_O3_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_pan_gc_gfas  = tmp_pan_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_hcho_gc_gfas = tmp_hcho_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_acet_gc_gfas = tmp_acet_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_benz_gc_gfas = tmp_benz_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]
            ;        ch3oh_gc= [ch3oh_gc, gc.MOH*1e9]
                    tmp_ald2_gc_gfas = tmp_ald2_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]

                    tmp_no_gc_gfas   = tmp_no_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_no2_gc_gfas  = tmp_no2_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_so2_gc_gfas  = tmp_so2_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]

                    ;tmp_na = tmp_na[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_oh_gc_gfas  = tmp_oh_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]

                    tmp_date_gc_gfas = tmp_date_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_utc_gc_gfas  = tmp_utc_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_doy_gc_gfas  = tmp_doy_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_lat_gc_gfas  = tmp_lat_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_lon_gc_gfas  = tmp_lon_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_alt_gc_gfas  = tmp_alt_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_prs_gc_gfas  = tmp_prs_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]
             ;; add vocs
                    tmp_c3h8_gc_gfas   = tmp_c3h8_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c3h6_gc_gfas   = tmp_c3h6_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c2h5oh_gc_gfas = tmp_c2h5oh_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c5h8_gc_gfas   = tmp_c5h8_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c7h8_gc_gfas   = tmp_c7h8_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_dms_gc_gfas    = tmp_dms_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_mek_gc_gfas    = tmp_mek_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c8h10_gc_gfas  = tmp_c8h10_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]

                    tmp_acta_gc_gfas   = tmp_acta_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_macr_mvk_gc_gfas = tmp_macr_mvk_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_hcooh_gc_gfas  = tmp_hcooh_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]

                    tmp_mtpa_gc_gfas  = tmp_mtpa_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_limo_gc_gfas  = tmp_limo_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_mtpo_gc_gfas  = tmp_mtpo_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]

                ;;add lumped species
                    tmp_alk4_gc_gfas   = tmp_alk4_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_prpe_gc_gfas   = tmp_prpe_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_rcho_gc_gfas   = tmp_rcho_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]

                    ;tmp_c2h4_gc_gfas   = tmp_c2h4_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c2h6_gc_gfas   = tmp_c2h6_gc_gfas[ind_enter_time[0]:ind_exit_time[-1]]
                    
                    
                endif
                if ct1 le 0 or ct2 le 0 then stop
            endif
            help,tmp_co_gc_gfas

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
                    'planelog2sav/output_'+'qfed_4x5' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'
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
            tmp_c2h6_gc_qfed   = gc.c2h6*1e9/2 ; 2C
            
    ;; concert hhmm(utc) into ss(utc)      
            ;print,tmp_utc_gc_qfed
            hh = floor(tmp_utc_gc_qfed/100)
            ind = where(hh gt 12, ct)
            if ct gt 0 then hh[ind] = hh[ind]
            ind = where(hh lt 12, ct)
            if ct gt 0 then hh[ind] = hh[ind] + 24
            mm = tmp_utc_gc_qfed- floor(tmp_utc_gc_qfed/100)*100
            tmp_utc_gc_qfed = float(hh)*60*60+float(mm)*60
            
            if keyword_set(emipass) or keyword_set(young) or $
                keyword_set(intermed) or keyword_set(GT4) or keyword_set(LT4) or $
                keyword_set(nonsmoke) or keyword_set(aged) then begin
                ind_enter_time = where(abs(tmp_utc_gc_qfed - enter_time[n]) le 30,ct1)
                ind_exit_time = where(abs(tmp_utc_gc_qfed - exit_time[n]) le 30, ct2)    
                if ct1 gt 0 and ct2 gt 0 then begin
                    tmp_co_gc_qfed   = tmp_co_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_O3_gc_qfed   = tmp_O3_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_pan_gc_qfed  = tmp_pan_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_hcho_gc_qfed = tmp_hcho_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_acet_gc_qfed = tmp_acet_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_benz_gc_qfed = tmp_benz_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]
            ;        ch3oh_gc= [ch3oh_gc, gc.MOH*1e9]
                    tmp_ald2_gc_qfed = tmp_ald2_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]

                    tmp_no_gc_qfed   = tmp_no_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_no2_gc_qfed  = tmp_no2_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_so2_gc_qfed  = tmp_so2_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]

                    ;tmp_na = tmp_na[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_oh_gc_qfed  = tmp_oh_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]

                    tmp_date_gc_qfed = tmp_date_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_utc_gc_qfed  = tmp_utc_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_doy_gc_qfed  = tmp_doy_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_lat_gc_qfed  = tmp_lat_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_lon_gc_qfed  = tmp_lon_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_alt_gc_qfed  = tmp_alt_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_prs_gc_qfed  = tmp_prs_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]
             ;; add vocs
                    tmp_c3h8_gc_qfed   = tmp_c3h8_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c3h6_gc_qfed   = tmp_c3h6_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c2h5oh_gc_qfed = tmp_c2h5oh_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c5h8_gc_qfed   = tmp_c5h8_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c7h8_gc_qfed   = tmp_c7h8_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_dms_gc_qfed    = tmp_dms_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_mek_gc_qfed    = tmp_mek_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c8h10_gc_qfed  = tmp_c8h10_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]

                    tmp_acta_gc_qfed   = tmp_acta_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_macr_mvk_gc_qfed = tmp_macr_mvk_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_hcooh_gc_qfed  = tmp_hcooh_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]

                    tmp_mtpa_gc_qfed  = tmp_mtpa_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_limo_gc_qfed  = tmp_limo_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_mtpo_gc_qfed  = tmp_mtpo_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]

                ;;add lumped species
                    tmp_alk4_gc_qfed   = tmp_alk4_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_prpe_gc_qfed   = tmp_prpe_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_rcho_gc_qfed   = tmp_rcho_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]

                    ;tmp_c2h4_gc_qfed   = tmp_c2h4_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c2h6_gc_qfed   = tmp_c2h6_gc_qfed[ind_enter_time[0]:ind_exit_time[-1]]
                    
                endif
                if ct1 le 0 and ct2 le 0 then stop
            endif
            help,tmp_co_gc_qfed

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
                ;gcfi_threegfas   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                ;    'planelog2sav/output_'+'threegfas_tmp' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'
                ;gcfi_threegfas   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                ;    'planelog2sav/output_'+'3Xgfas_tmp' + '/mrg10m_wecan_c130_'+dates[n]+'.sav'
                gcfi_threegfas   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                    'planelog2sav/output_'+'3Xgfas_tmp' + '/mrg5m_wecan_c130_'+dates[n]+'.sav'
                
            endif
            if keyword_set(fbf) then begin
                gcfi_threegfas   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                    'planelog2sav/output_'+'threegfas_4x5' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'
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
            tmp_c2h6_gc_threegfas   = gc.c2h6*1e9/2 ; 2C
            
    ;; concert hhmm(utc) into ss(utc)      
            ;print,tmp_utc_gc_threegfas
            hh = floor(tmp_utc_gc_threegfas/100)
            ind = where(hh gt 12, ct)
            if ct gt 0 then hh[ind] = hh[ind]
            ind = where(hh lt 12, ct)
            if ct gt 0 then hh[ind] = hh[ind] + 24
            mm = tmp_utc_gc_threegfas- floor(tmp_utc_gc_threegfas/100)*100
            tmp_utc_gc_threegfas = float(hh)*60*60+float(mm)*60
            
            if keyword_set(emipass) or keyword_set(young) or $
                keyword_set(intermed) or keyword_set(GT4) or keyword_set(LT4) or $
                keyword_set(nonsmoke) or keyword_set(aged) then begin
                ind_enter_time = where(abs(tmp_utc_gc_threegfas - enter_time[n]) le 30,ct1)
                ind_exit_time = where(abs(tmp_utc_gc_threegfas - exit_time[n]) le 30, ct2)    
                if ct1 gt 0 and ct2 gt 0 then begin
                    tmp_co_gc_threegfas   = tmp_co_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_O3_gc_threegfas   = tmp_O3_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_pan_gc_threegfas  = tmp_pan_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_hcho_gc_threegfas = tmp_hcho_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_acet_gc_threegfas = tmp_acet_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_benz_gc_threegfas = tmp_benz_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]
            ;        ch3oh_gc= [ch3oh_gc, gc.MOH*1e9]
                    tmp_ald2_gc_threegfas = tmp_ald2_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]

                    tmp_no_gc_threegfas   = tmp_no_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_no2_gc_threegfas  = tmp_no2_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_so2_gc_threegfas  = tmp_so2_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]

                    ;tmp_na = tmp_na[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_oh_gc_threegfas  = tmp_oh_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]

                    tmp_date_gc_threegfas = tmp_date_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_utc_gc_threegfas  = tmp_utc_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_doy_gc_threegfas  = tmp_doy_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_lat_gc_threegfas  = tmp_lat_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_lon_gc_threegfas  = tmp_lon_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_alt_gc_threegfas  = tmp_alt_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_prs_gc_threegfas  = tmp_prs_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]
             ;; add vocs
                    tmp_c3h8_gc_threegfas   = tmp_c3h8_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c3h6_gc_threegfas   = tmp_c3h6_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c2h5oh_gc_threegfas = tmp_c2h5oh_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c5h8_gc_threegfas   = tmp_c5h8_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c7h8_gc_threegfas   = tmp_c7h8_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_dms_gc_threegfas    = tmp_dms_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_mek_gc_threegfas    = tmp_mek_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c8h10_gc_threegfas  = tmp_c8h10_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]

                    tmp_acta_gc_threegfas   = tmp_acta_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_macr_mvk_gc_threegfas = tmp_macr_mvk_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_hcooh_gc_threegfas  = tmp_hcooh_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]

                    tmp_mtpa_gc_threegfas  = tmp_mtpa_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_limo_gc_threegfas  = tmp_limo_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_mtpo_gc_threegfas  = tmp_mtpo_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]

                ;;add lumped species
                    tmp_alk4_gc_threegfas   = tmp_alk4_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_prpe_gc_threegfas   = tmp_prpe_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_rcho_gc_threegfas   = tmp_rcho_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]

                    ;tmp_c2h4_gc_threegfas   = tmp_c2h4_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c2h6_gc_threegfas   = tmp_c2h6_gc_threegfas[ind_enter_time[0]:ind_exit_time[-1]]
                    
                endif
                if ct1 le 0 and ct2 le 0 then stop
            endif
            help,tmp_co_gc_threegfas

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
;; use finn for now!!!!!!!!!!!
            if keyword_set(nested) then begin
                ;gcfi_nobb   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                ;    'planelog2sav/output_'+'nobb' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'
                gcfi_nobb   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                    'planelog2sav/output_'+'nobb_tmp' + '/mrg5m_wecan_c130_'+dates[n]+'.sav'
            endif
            if keyword_set(fbf) then begin
                gcfi_nobb   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                    'planelog2sav/output_'+'gfas_4x5' + '/mrg1m_wecan_c130_'+dates[n]+'.sav'
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
            tmp_alk4_gc_nobb   = gc.alk4*1e9/4.3 ; 4.3C
            tmp_prpe_gc_nobb   = gc.prpe*1e9/3   ; 3C 
            tmp_rcho_gc_nobb   = gc.rcho*1e9

            ;tmp_c2h4_gc_nobb   = gc.c2h4*1e9
            tmp_c2h6_gc_nobb   = gc.c2h6*1e9/2 ; 2C
        
    ;; concert hhmm(utc) into ss(utc)      
            ;print,tmp_utc_gc_nobb
            hh = floor(tmp_utc_gc_nobb/100)
            ind = where(hh gt 12, ct)
            if ct gt 0 then hh[ind] = hh[ind]
            ind = where(hh lt 12, ct)
            if ct gt 0 then hh[ind] = hh[ind] + 24
            
            mm = tmp_utc_gc_nobb- floor(tmp_utc_gc_nobb/100)*100
            tmp_utc_gc_nobb = float(hh)*60*60+float(mm)*60
                        
            
            if keyword_set(emipass) or keyword_set(young) or $
                keyword_set(intermed) or keyword_set(GT4) or keyword_set(LT4) or $
                keyword_set(nonsmoke) or keyword_set(aged) then begin
                ind_enter_time = where(abs(tmp_utc_obs - enter_time[n]) le 30,ct1)
                ind_exit_time = where(abs(tmp_utc_obs - exit_time[n]) le 30, ct2)    
                if ct1 gt 0 and ct2 gt 0 then begin
                    tmp_co_gc_nobb   = tmp_co_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_O3_gc_nobb   = tmp_O3_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_pan_gc_nobb  = tmp_pan_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_hcho_gc_nobb = tmp_hcho_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_acet_gc_nobb = tmp_acet_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_benz_gc_nobb = tmp_benz_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]
            ;        ch3oh_gc= [ch3oh_gc, gc.MOH*1e9]
                    tmp_ald2_gc_nobb = tmp_ald2_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]

                    tmp_no_gc_nobb   = tmp_no_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_no2_gc_nobb  = tmp_no2_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_so2_gc_nobb  = tmp_so2_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]

                    ;tmp_na = tmp_na[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_oh_gc_nobb  = tmp_oh_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]

                    tmp_date_gc_nobb = tmp_date_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_utc_gc_nobb  = tmp_utc_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_doy_gc_nobb  = tmp_doy_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_lat_gc_nobb  = tmp_lat_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_lon_gc_nobb  = tmp_lon_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_alt_gc_nobb  = tmp_alt_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_prs_gc_nobb  = tmp_prs_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]
             ;; add vocs
                    tmp_c3h8_gc_nobb   = tmp_c3h8_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c3h6_gc_nobb   = tmp_c3h6_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c2h5oh_gc_nobb = tmp_c2h5oh_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c5h8_gc_nobb   = tmp_c5h8_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c7h8_gc_nobb   = tmp_c7h8_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_dms_gc_nobb    = tmp_dms_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_mek_gc_nobb    = tmp_mek_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c8h10_gc_nobb  = tmp_c8h10_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]

                    tmp_acta_gc_nobb   = tmp_acta_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_macr_mvk_gc_nobb = tmp_macr_mvk_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_hcooh_gc_nobb  = tmp_hcooh_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]

                    tmp_mtpa_gc_nobb  = tmp_mtpa_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_limo_gc_nobb  = tmp_limo_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_mtpo_gc_nobb  = tmp_mtpo_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]

                ;;add lumped species
                    tmp_alk4_gc_nobb   = tmp_alk4_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_prpe_gc_nobb   = tmp_prpe_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_rcho_gc_nobb   = tmp_rcho_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]

                    ;tmp_c2h4_gc_nobb   = tmp_c2h4_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]
                    tmp_c2h6_gc_nobb   = tmp_c2h6_gc_nobb[ind_enter_time[0]:ind_exit_time[-1]]
                    
                endif
                if ct1 le 0 or ct2 le 0 then stop
            endif

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

        ;c2h4_gc_gfed4 = c2h4_gc_gfed4[1:*]
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

        ;c2h4_gc_finn = c2h4_gc_finn[1:*]
        c2h6_gc_finn = c2h6_gc_finn[1:*]
    ;; GFAS
        co_gc_gfas   = co_gc_gfas[1:*]
        o3_gc_gfsa   = o3_gc_gfas[1:*]
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

        ;c2h4_gc_gfas = c2h4_gc_gfas[1:*]
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

        ;c2h4_gc_qfed = c2h4_gc_qfed[1:*]
        c2h6_gc_qfed = c2h6_gc_qfed[1:*]
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

        ;c2h4_gc_nobb = c2h4_gc_nobb[1:*]
        c2h6_gc_nobb = c2h6_gc_nobb[1:*]
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
        
        ;c2h4_gc_threegfas = c2h4_gc_threegfas[1:*]
        c2h6_gc_threegfas = c2h6_gc_threegfas[1:*]


            
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

        ;;test
        if n_elements(co_obs) ne n_elements(co_gc_gfed4) or $
            n_elements(co_obs) ne n_elements(co_gc_finn) or $
            n_elements(co_obs) ne n_elements(co_gc_gfas) or $
            n_elements(co_obs) ne n_elements(co_gc_qfed) or $
            n_elements(co_obs) ne n_elements(co_gc_nobb) then stop 
;; ------------------------------------------------

    ;; jupyter time
        ;doy_obs = jday_obs + utc_obs / 3600. / 24. 

    ;; lixu, 12/18/2019
    ;; plume filter setting, use in the plotting part
        no2_thresh = 4. ; ppb
        co_thresh  = 150. ; ppb
    ;; biomass burning threshhold 
        hcn_thresh   = 500./1000 ;pptv to ppbv
        ch3cn_thresh = 225/1000. ;;ppt to ppb
    ;; strat filter
        strat_thresh=1.25 ; O3/CO ratio
    ;; CO threshold of fresh and aged smoke 
        smoke_thresh = 1000 ;ppb
        ;smoke_thresh = 400
    ;; get hhmm of each day
        ;utc_obs = 24.*(doy_obs mod 1.)
        ;utc_gc_gfas  = 24.*(doy_gc_gfas  mod 1.)
        ;utc_gc_gfed4  = 24.*(doy_gc_gfed4  mod 1.)
        ;utc_gc_qfed  = 24.*(doy_gc_qfed  mod 1.)
        ;utc_gc_finn  = 24.*(doy_gc_finn  mod 1.)

    ;; for vertical profile
        ;yrange = [0,7]
    ;; lxu, 06202021
        acn_co_tresh = 2.01
;; ========keep the same for regression, profiles and time series=================
        if keyword_set(primary) then rows = 3
        if keyword_set(primary) then cols = 3

        if keyword_set(secondary) then rows = 2
        if keyword_set(secondary) then cols = 2

        if keyword_set(addup) then rows = 1
        if keyword_set(addup) then cols = 1
        
        if keyword_set(all) then rows = 3
        if keyword_set(all) then cols = 3
        
        if keyword_set(onebyone) then rows = 1 
        if keyword_set(onebyone) then cols = 1 
    
        if keyword_set(lumped) then rows = 1
        if keyword_set(lumped) then cols = 3
        
        if keyword_set(chemistry) then rows = 2
        if keyword_set(chemistry) then cols = 2
        if keyword_set(nmhcs) then rows = 2
        if keyword_set(nmhcs) then cols = 2

        if keyword_set(ovocs) then rows = 3
        if keyword_set(ovocs) then cols = 3
        
        multipanel, rows = rows, cols = cols
        
        n_species = rows*cols
        
        if keyword_set(all) then n_species = 6
        if keyword_set(addup) then n_species = 1
        if keyword_set(chemistry) then n_species = 4
        if keyword_set(nmhcs) then n_species = 4

        if keyword_set(ovocs) then n_species = 9
        if keyword_set(lumped) then n_species = 5
        
        ; initialize save-out data
        spec_obs_total   = []
        spec_gfed4_total = []
        spec_gfas_total  = []
        spec_qfed_total  = []
        spec_obs_std_total   = []
        spec_gfed4_std_total = []
        spec_gfas_std_total  = []
        spec_qfed_std_total  = []
        corr_obs_total   = []
        corr_gfed4_total = []
        corr_gfas_total  = []
        corr_qfed_total  = []
        title_all = []
        
        for s = 0, n_species - 1 do begin

;; all/9 VOCs here
            if keyword_set(all) then begin
                case s of                     
                    ;;c2h6(AWAS, discrete data)
                    ;;hcho
                    0:begin
                        ovoc = ch2o_obs_ptr
                        gvoc_gfas = hcho_gc_gfas
                        gvoc_gfed4 = hcho_gc_gfed4
                        gvoc_qfed  = hcho_gc_qfed
                        gvoc_finn  = hcho_gc_finn
                        gvoc_threegfas  = hcho_gc_threegfas
                        gvoc_nobb       = hcho_gc_nobb
                        title = 'Formaldehyde'
                    end
                    ;;ald2
                    1:begin
                        ovoc = ald2_obs_ptr
                        gvoc_gfas = ald2_gc_gfas
                        gvoc_gfed4= ald2_gc_gfed4
                        gvoc_qfed = ald2_gc_qfed
                        gvoc_finn = ald2_gc_finn
                        gvoc_threegfas = ald2_gc_threegfas
                        gvoc_nobb      = ald2_gc_nobb
                        title = 'Acetaldehyde'
                    end
                    ;;acet
                    2:begin
                        ovoc = acet_obs_ptr
                        gvoc_gfas = acet_gc_gfas
                        gvoc_gfed4 = acet_gc_gfed4
                        gvoc_qfed = acet_gc_qfed
                        gvoc_finn = acet_gc_finn
                        gvoc_threegfas = acet_gc_threegfas
                        gvoc_nobb      = acet_gc_nobb
                        title = 'Acetone'  
                    end

                    ;;benz
                    3:begin
                        ovoc = benz_obs_ptr
                        gvoc_gfas = benz_gc_gfas
                        gvoc_gfed4 = benz_gc_gfed4
                        gvoc_qfed = benz_gc_qfed
                        gvoc_finn = benz_gc_finn
                        gvoc_threegfas = benz_gc_threegfas
                        gvoc_nobb      = benz_gc_nobb
                        title = 'Benzene'
                    end
                    ;;tolu
                    4:begin
                        ovoc = tolu_obs_ptr
                        gvoc_gfas = c7h8_gc_gfas
                        gvoc_gfed4 = c7h8_gc_gfed4
                        gvoc_qfed = c7h8_gc_qfed
                        gvoc_finn = c7h8_gc_finn
                        gvoc_threegfas = c7h8_gc_threegfas
                        gvoc_nobb      = c7h8_gc_nobb
                        title = 'Toluene'
                    end  
                    ;;xyle
                    5:begin
                        ovoc = Xylenes_obs_ptr
                        gvoc_gfas = c8h10_gc_gfas
                        gvoc_gfed4 = c8h10_gc_gfed4
                        gvoc_qfed  = c8h10_gc_qfed
                        gvoc_finn  = c8h10_gc_finn
                        gvoc_threegfas = c8h10_gc_threegfas
                        gvoc_nobb      = c8h10_gc_nobb
                        title = 'Xylene'
                    end

                    ;;CO
                    6:begin
                        ovoc = co_obs
                        gvoc_gfed4 = co_gc_gfed4
                        gvoc_finn = co_gc_finn
                        gvoc_gfas = co_gc_gfas
                        gvoc_qfed = co_gc_qfed
                        gvoc_threegfas = co_gc_threegfas
                        gvoc_nobb = co_gc_nobb
                        title = 'CO'
                    end
                    
                endcase
            endif

            ;if s eq 19 then begin
            ;    print, gvoc_gfed4
            ;    print, co_gc_gfed4
            ;    stop
            ;endif
            ;;formic acid, follow up Catie's project
            ;;MACR_MVK, MOH, ISOP, monotherpenes follow up Wade's project

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
            zvar_gc_nobb        = alt_gc_nobb

    ;for time series, reset the time array every loop (diff)
            ;xtime_obs  = doy_obs
            ;xtime_gc   = doy_gc_gfas
            
    ;for voc .vs. co
            oco_tmp=co_obs
            gco_gfas=co_gc_gfas
            gco_gfed4=co_gc_gfed4
            gco_qfed=co_gc_qfed     
            gco_finn=co_gc_finn
            gco_threegfas=co_gc_threegfas
            gco_nobb     =co_gc_nobb

    ;;add for OVOC VS. ISOP
            oisop_ptr = isop_obs_ptr
            oisop_toga = isop_obs_toga
            oisop_tmp = isop_obs_ptr
            
            gisop_gfas = c5h8_gc_gfas
            gisop_gfed4 = c5h8_gc_gfed4
            gisop_qfed = c5h8_gc_qfed
            gisop_finn = c5h8_gc_finn
            gisop_threegfas = c5h8_gc_threegfas
            gisop_nobb      = c5h8_gc_nobb
            
    ;;add for OVOC vs acetone
            oacet_tmp = acet_obs_ptr
            gacet_gfed4 = acet_gc_gfed4
            gacet_finn = acet_gc_finn
            gacet_gfas = acet_gc_gfas
            gacet_qfed = acet_gc_qfed
            gacet_threegfas = acet_gc_threegfas
            gacet_nobb = acet_gc_nobb
    
    ;;add for OVOC vs benz
            obenz_tmp = benz_obs
            gbenz_gfed4 = benz_gc_gfed4
            gbenz_finn = benz_gc_finn
            gbenz_gfas = benz_gc_gfas
            gbenz_qfed = benz_gc_qfed
            gbenz_threegfas = benz_gc_threegfas
            gbenz_nobb = benz_gc_nobb

    ;;add for OVOC vs tolu
            otolu_tmp = tolu_obs_ptr
            gtolu_gfed4 = c7h8_gc_gfed4
            gtolu_finn = c7h8_gc_finn
            gtolu_gfas = c7h8_gc_gfas
            gtolu_qfed = c7h8_gc_qfed
            gtolu_threegfas = c7h8_gc_threegfas
            gtolu_nobb = c7h8_gc_nobb
            
    
    ;for filters
            no2_filter=no2_obs
            o3_filter=o3_obs
            co_filter=co_obs
            ch3cn_filter=ch3cn_obs
            hcn_filter = hcn_obs
            acn_co_filter = ch3cn_obs/co_obs*1000 ; convert it into ppb/ppm
            
    ;for new criteria
            isop_filter = isop_obs_ptr
            mek_filter = mek_obs_ptr 
            benz_filter = benz_obs_ptr
            tolu_filter = tolu_obs_ptr
            xyle_filter = c8h10_obs_toga
            mono_filter = tricyclene_obs_toga + aPinene_obs_TOGA + Camphene_obs_TOGA + bPineneMyrcene_obs_toga
            ch2o_filter = ch2o_obs_ptr
            c2h5oh_filter = c2h5oh_obs
            
    ;;tracks
            lat_gc = lat_gc_gfas
            lon_gc = lon_gc_gfas
    ;;setting tmp value
            tmp_isop_ptr = oisop_ptr
            tmp_isop_toga = oisop_toga
            tmp_lat_obs = lat_obs
            tmp_lon_obs = lon_obs
            tmp_lat_gc  = lat_gc
            tmp_lon_gc  = lon_gc            
    ;cloud
            tmp_rhum = rhum  


        ;;remove missing value in co_obs_picarro/oco_tmp
            if keyword_set(xeqco) then i = where(oco_tmp le 0, ct)
            if keyword_set(xeqisop) then i = where(oisop_tmp le 0, ct)
            if keyword_set(xeqacet) then i = where(oacet_tmp le 0, ct)
            if keyword_set(xeqbenz) then i = where(obenz_tmp le 0, ct)
            if keyword_set(xeqtolu) then i = where(otolu_tmp le 0, ct)
            
            if not keyword_set(mod_only) then begin
                if ct gt 0 then begin
                    ovoc[i] = !VALUES.F_NAN
                end
                ind = where(ovoc le 0, count)
                if count gt 0 then ovoc[ind] = !VALUES.F_NAN
            endif
            
;;=================
;;plume filters
;;=================
            if plume eq 1 then begin
    ;; notice that we don't use noy to remove data
                ;;filter1:filter STE, urban plumes
                if keyword_set(filter1_other) then begin
                    ind = where(no2_filter ge no2_thresh or $
                        o3_filter/co_filter ge strat_thresh, ct)
                    if ct gt 0 then ovoc[ind] = !values.f_nan
                    print, 'this is the data points to be deleted for flitler1', ct
                endif
                ;;filter2:use 25 percentile of acetonitrile
                if keyword_set(filter2_bb) then begin
                    ind = where(ch3cn_filter le ch3cn_thresh,ct)
                    if ct gt 0 then ovoc[ind] = !values.f_nan
                    print, 'this is the data points to be deleted for flitler2', ct
                endif
                if keyword_set(filter2_nobb) then begin
                    ind = where(ch3cn_filter ge ch3cn_thresh,ct)
                    if ct gt 0 then ovoc[ind] = !values.f_nan
                    print, 'this is the data points to be deleted for flitler2', ct
                endif
                ;;fitler3: use acn/co ratio (2.01ppb/ppm) lt the threshold: nonbb (doesn't work well)
                if keyword_set(filter3_bb) then begin
                    ind = where(acn_co_filter le acn_co_tresh,ct)
                    if ct gt 0 then ovoc[ind] = !values.f_nan
                    print, 'this is the data points to be deleted for flitler3', ct
                endif
                if keyword_set(filter3_nobb) then begin
                    ind = where(acn_co_filter ge acn_co_tresh,ct)
                    if ct gt 0 then ovoc[ind] = !values.f_nan
                    print, 'this is the data points to be deleted for flitler3', ct
                endif
                if keyword_set(filter4_bb) then begin
                    ind = where(hcn_filter le hcn_thresh,ct)
                    if ct gt 0 then ovoc[ind] = !values.f_nan
                    print, 'this is the data points to be deleted for flitler4', ct
                endif
                if keyword_set(filter4_nobb) then begin
                    ind = where(hcn_filter ge hcn_thresh,ct)
                    if ct gt 0 then ovoc[ind] = !values.f_nan
                    print, 'this is the data points to be deleted for flitler4', ct
                endif
                if keyword_set(filter5_nobb) then begin
                    ind = where(co_filter ge co_thresh,ct)
                    if ct gt 0 then ovoc[ind] = !values.f_nan
                    print, 'this is the data points to be deleted for flitler5', ct
                endif                
                
                
           endif

            remove_ind = where(finite(ovoc))

;///////////////////////////////////////////////////////////////////////////////////////////
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;++                                      STEP 4                                           ++
;++                               Removing missing value                                  ++
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;///////////////////////////////////////////////////////////////////////////////////////////
            control_ind = 0
            if not keyword_set(lumped) then control_ind = 1
            if control_ind eq 1 then begin
                remove_ind1 = where(finite(ovoc))
                if keyword_set(xeqco) then remove_ind2 = where(finite(oco_tmp))
                if keyword_set(xeqisop) then remove_ind2 = where(finite(oisop_tmp))
                if keyword_set(xeqacet) then remove_ind2 = where(finite(oacet_tmp))
                if keyword_set(xeqbenz) then remove_ind2 = where(finite(obenz_tmp))
                if keyword_set(xeqtolu) then remove_ind2 = where(finite(otolu_tmp))
                com_ind = setintersection(remove_ind1,remove_ind2)
             
                if (remove_ind1[0] ne -1) and (remove_ind2[0] ne -1) then begin
                    ;ovoc_ptr = ovoc_ptr[remove_ind]
                    ;ovoc_toga = ovoc_toga[remove_ind] 
                    ovoc=ovoc[com_ind]
                    gvoc_gfed4 = gvoc_gfed4[com_ind]
                    gvoc_finn = gvoc_finn[com_ind]
                    gvoc_gfas = gvoc_gfas[com_ind]
                    gvoc_qfed = gvoc_qfed[com_ind]
                    gvoc_threegfas = gvoc_threegfas[com_ind]
                    gvoc_nobb = gvoc_nobb[com_ind]

                    ;for vertical profile
                    zvar_obs=zvar_obs[com_ind]
                    zvar_gc_gfed4=zvar_gc_gfed4[com_ind]
                    zvar_gc_finn = zvar_gc_finn[com_ind]
                    zvar_gc_gfas=zvar_gc_gfas[com_ind]
                    zvar_gc_qfed=zvar_gc_qfed[com_ind]
                    zvar_gc_threegfas = zvar_gc_threegfas[com_ind]
                    zvar_gc_nobb = zvar_gc_nobb[com_ind]

                    ;for time series
                    ;xtime_obs=xtime_obs[remove_ind] 
                    ;xtime_gc = xtime_gc[remove_ind]

                    ;for filters
                    no2_filter=no2_filter[com_ind]
                    o3_filter=o3_filter[com_ind]
                    co_filter = co_filter[com_ind]
                    ch3cn_filter=ch3cn_filter[com_ind]
                    hcn_filter = hcn_filter[com_ind]

                    isop_filter = isop_filter[com_ind]
                    mek_filter = mek_filter[com_ind]
                    benz_filter = benz_filter[com_ind]
                    tolu_filter = tolu_filter[com_ind]
                    xyle_filter = xyle_filter[com_ind]
                    mono_filter = mono_filter[com_ind]
                    ch2o_filter = ch2o_filter[com_ind]
                    c2h5oh_filter = c2h5oh_filter[com_ind]

                    ;for voc vs. co
                    oco_tmp=oco_tmp[com_ind]
                    gco_gfed4=gco_gfed4[com_ind]
                    gco_finn=gco_finn[com_ind]
                    gco_gfas=gco_gfas[com_ind]
                    gco_qfed=gco_qfed[com_ind]
                    gco_threegfas=gco_threegfas[com_ind]
                    gco_nobb=gco_nobb[com_ind]

                    ;for voc vs. isop
                    tmp_isop_ptr = tmp_isop_ptr[com_ind]
                    tmp_isop_toga = tmp_isop_toga[com_ind]
                    oisop_tmp = isop_obs_ptr[com_ind]
                    
                    gisop_gfed4 = gisop_gfed4[com_ind]
                    gisop_finn = gisop_finn[com_ind]
                    gisop_gfas = gisop_gfas[com_ind]
                    gisop_qfed = gisop_qfed[com_ind]
                    gisop_threegfas = gisop_threegfas[com_ind]
                    gisop_nobb = gisop_nobb[com_ind]

                    ;for voc vs. acet
                    oacet_tmp = acet_obs[com_ind]
                    gacet_gfed4 = gacet_gfed4[com_ind]
                    gacet_finn = gacet_finn[com_ind]
                    gacet_gfas = gacet_gfas[com_ind]
                    gacet_qfed = gacet_qfed[com_ind]
                    gacet_threegfas = gacet_threegfas[com_ind]
                    gacet_nobb = gacet_nobb[com_ind]

                    ;for voc vs. benz
                    obenz_tmp = benz_obs[com_ind]
                    gbenz_gfed4 = gbenz_gfed4[com_ind]
                    gbenz_finn = gbenz_finn[com_ind]
                    gbenz_gfas = gbenz_gfas[com_ind]
                    gbenz_qfed = gbenz_qfed[com_ind]
                    gbenz_threegfas = gbenz_threegfas[com_ind]
                    gbenz_nobb = gbenz_nobb[com_ind]

                    ;for voc vs. tolu
                    otolu_tmp = tolu_obs_ptr[com_ind]
                    gtolu_gfed4 = gtolu_gfed4[com_ind]
                    gtolu_finn = gtolu_finn[com_ind]
                    gtolu_gfas = gtolu_gfas[com_ind]
                    gtolu_qfed = gtolu_qfed[com_ind]
                    gtolu_threegfas = gtolu_threegfas[com_ind]
                    gtolu_nobb = gtolu_nobb[com_ind]

                    ;;tracks
                    tmp_lat_obs = tmp_lat_obs[com_ind]
                    tmp_lon_obs = tmp_lon_obs[com_ind]
                    tmp_lat_gc = tmp_lat_gc[com_ind]
                    tmp_lon_gc = tmp_lon_gc[com_ind]
                endif
                ;;set NA 
                if com_ind[0] eq -1 then gvoc_gfas = !values.f_nan
                if com_ind[0] eq -1 then gvoc_gfed4 = !values.f_nan
                if com_ind[0] eq -1 then gvoc_qfed = !values.f_nan
                if com_ind[0] eq -1 then gvoc_finn = !values.f_nan
                if com_ind[0] eq -1 then gvoc_threegfas = !values.f_nan
                if com_ind[0] eq -1 then gvoc_nobb = !values.f_nan

                ;if remove_ind[0] eq -1 then ovoc_ptr = !values.f_nan
                ;if remove_ind[0] eq -1 then ovoc_toga = !values.f_nan
                if com_ind[0] eq -1 then ovoc = !values.f_nan
            
            endif
            
            
            num_array = n_elements(c2h5oh_filter)
            tmp_ind = FLTARR(num_array)
            
;; anthropogenic setting, 1
            remove_ind1 = where(c2h5oh_filter/benz_filter gt 0.5, ct1)
            remove_ind2 = where(ch2o_filter/benz_filter gt 17, ct2)
            com_ind1 = setintersection(remove_ind1,remove_ind2)            
            if (remove_ind1[0] ne -1) and (remove_ind2[0] ne -1) then begin
                tmp_ind[com_ind1] = 1
            endif

;; biomasss burning setting, 2
            remove_ind1 = where(c2h5oh_filter/benz_filter le 0.5, ct1)
            remove_ind2 = where(ch2o_filter/benz_filter le 17.5, ct2)
            com_ind2 = setintersection(remove_ind1,remove_ind2)
            if (remove_ind1[0] ne -1) and (remove_ind2[0] ne -1) then begin
                tmp_ind[com_ind2] = 2
            endif 
;///////////////////////////////////////////////////////////////////////////////////////////
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;++                                      STEP 5                                           ++
;++                                  DO THE PLOTTING                                      ++
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;///////////////////////////////////////////////////////////////////////////////////////////

        ;; unit scale factor
            if keyword_set(xeqco) then unitscale = 1000
            if keyword_set(xeqbenz) then unitscale = 1
            if keyword_set(xeqtolu) then unitscale = 1
            if keyword_set(xeqisop) then unitscale = 1
            if keyword_set(xeqacet) then unitscale = 1
            
            
            y1=ovoc
            ;y1=ovoc_ptr
            ;y2=ovoc_toga
            y2=gvoc_gfed4
            y3=gvoc_finn
            y4=gvoc_gfas
            y5=gvoc_qfed
            y6=gvoc_threegfas
            y7 = gvoc_nobb
            
            
            xa1=oco_tmp
            xa2=gco_gfed4
            xa3=gco_finn
            xa4=gco_gfas
            xa5=gco_qfed
            xa6=gco_threegfas
            xa7= gco_nobb
            
            xb1 = obenz_tmp
            xb2 = gbenz_gfed4
            xb3 = gbenz_finn
            xb4 = gbenz_gfas
            xb5 = gbenz_qfed
            xb6 = gbenz_threegfas
            xb7 = gbenz_nobb
            
            xc1 = otolu_tmp
            xc2 = gtolu_gfed4
            xc3 = gtolu_finn
            xc4 = gtolu_gfas
            xc5 = gtolu_qfed
            xc6 = gtolu_threegfas
            xc7 = gtolu_nobb
            
            xd1 = oisop_tmp
            xd2 = gisop_gfed4
            xd3 = gisop_finn
            xd4 = gisop_gfas
            xd5 = gisop_qfed
            xd6 = gisop_threegfas
            xd7 = gisop_nobb
            
            xe1 = oacet_tmp
            xe2 = gacet_gfed4
            xe3 = gacet_finn
            xe4 = gacet_gfas
            xe5 = gacet_qfed
            xe6 = gacet_threegfas
            xe7 = gacet_nobb
            if keyword_set(newthresh_nobb) then begin
                ind = where(tmp_ind eq 1.0, ct)
                print, ct
                if ct gt 0 then begin
                    y1 = y1[ind]
                    y2 = y2[ind]
                    y3 = y3[ind]
                    y4 = y4[ind]
                    y5 = y5[ind]
                    y6 = y6[ind]
                    y7 = y7[ind]

                    xa1=xa1[ind]
                    xa2=xa2[ind]
                    xa3=xa3[ind]
                    xa4=xa4[ind]
                    xa5=xa5[ind]
                    xa6= xa6[ind]
                    xa7= xa7[ind]

                    xb1=xb1[ind]
                    xb2=xb2[ind]
                    xb3=xb3[ind]
                    xb4=xb4[ind]
                    xb5=xb5[ind]
                    xb6= xb6[ind]
                    xb7= xb7[ind]

                    xc1=xc1[ind]
                    xc2=xc2[ind]
                    xc3=xc3[ind]
                    xc4=xc4[ind]
                    xc5=xc5[ind]
                    xc6= xc6[ind]
                    xc7= xc7[ind]
                    
                    xd1=xd1[ind]
                    xd2=xd2[ind]
                    xd3=xd3[ind]
                    xd4=xd4[ind]
                    xd5=xd5[ind]
                    xd6= xd6[ind]
                    xd7= xd7[ind]

                    xe1=xe1[ind]
                    xe2=xe2[ind]
                    xe3=xe3[ind]
                    xe4=xe4[ind]
                    xe5=xe5[ind]
                    xe6= xe6[ind]
                    xe7= xe7[ind]
                endif
            endif

            if keyword_set(newthresh_bb) then begin
                ind = where(tmp_ind eq 2, ct)
                if ct gt 0 then begin
                    y1 = y1[ind]
                    y2 = y2[ind]
                    y3 = y3[ind]
                    y4 = y4[ind]
                    y5 = y5[ind]
                    y6 = y6[ind]
                    y7 = y7[ind]

                    xa1=xa1[ind]
                    xa2=xa2[ind]
                    xa3=xa3[ind]
                    xa4=xa4[ind]
                    xa5=xa5[ind]
                    xa6= xa6[ind]
                    xa7= xa7[ind]

                    xb1=xb1[ind]
                    xb2=xb2[ind]
                    xb3=xb3[ind]
                    xb4=xb4[ind]
                    xb5=xb5[ind]
                    xb6= xb6[ind]
                    xb7= xb7[ind]

                    xc1=xc1[ind]
                    xc2=xc2[ind]
                    xc3=xc3[ind]
                    xc4=xc4[ind]
                    xc5=xc5[ind]
                    xc6= xc6[ind]
                    xc7= xc7[ind]
                    
                    xd1=xd1[ind]
                    xd2=xd2[ind]
                    xd3=xd3[ind]
                    xd4=xd4[ind]
                    xd5=xd5[ind]
                    xd6= xd6[ind]
                    xd7= xd7[ind]

                    xe1=xe1[ind]
                    xe2=xe2[ind]
                    xe3=xe3[ind]
                    xe4=xe4[ind]
                    xe5=xe5[ind]
                    xe6= xe6[ind]
                    xe7= xe7[ind]
                endif
            endif

            if keyword_set(xeqco) then begin
                x1 = xa1
                x2 = xa2
                x3 = xa3
                x4 = xa4
                x5 = xa5
                x6 = xa6
                x7 = xa7
            end
            
            if keyword_set(xeqbenz) then begin
                x1 = xb1
                x2 = xb2
                x3 = xb3
                x4 = xb4
                x5 = xb5
                x6 = xb6
                x7 = xb7                
            end
            
            if keyword_set(xeqtolu) then begin
                x1 = xc1
                x2 = xc2
                x3 = xc3
                x4 = xc4
                x5 = xc5
                x6 = xc6
                x7 = xc7
            end
            
            if keyword_set(xeqisop) then begin
                x1 = xd1
                x2 = xd2
                x3 = xd3
                x4 = xd4
                x5 = xd5
                x6 = xd6
                x7 = xd7
            end
            
            if keyword_set(xeqacet) then begin
                x1 = xe1
                x2 = xe2
                x3 = xe3
                x4 = xe4
                x5 = xe5
                x6 = xe6
                x7 = xe7
            end
            
            ind = where(finite(y1))            
            if ind[0] eq -1 then continue
                         
;;calculate for ER, intercept and R
            print,'this is ',title
            print,'slope, slope_se, intercept, R'
    ;;obs:ptr
            stat1 = org_boot(X=x1,Y=y1,ntrials=1000l) ;M (slope), MSE, MCI_boot, B, BSE, BCI_Boot, R, RCI_Boot, N            
            print,':OBS',' slope',stat1.M*unitscale, $
                ' SE of slope', stat1.MSE*unitscale, $
                ' CI of slope', stat1.MCI_boot*unitscale, $
                ' intercept', stat1.B, $
                ' SE of intercept', stat1.BSE, $
                ' CI of intercept', stat1.BCI_Boot, $
                ' correlations', stat1.R, $
                ' CI of correlations', stat1.RCI_Boot, $
                ' # of dp', stat1.N
    ;;GFED4
            stat2 = org_boot(X=x2,Y=y2,ntrials=1000l) ;M (slope), MSE, MCI_boot, B, BSE, BCI_Boot, R, RCI_Boot, N            
            print,':GFED4',' slope',stat2.M*unitscale,$
                ' SE of slope', stat2.MSE*unitscale,$
                ' CI of slope', stat2.MCI_boot*unitscale, $
                ' intercept', stat2.B, $
                ' SE of intercept', stat2.BSE, $
                ' CI of intercept', stat2.BCI_Boot, $
                ' correlations', stat2.R, $
                ' CI of correlations', stat2.RCI_Boot, $
                ' # of dp', stat2.N
 
    ;;FINN
            stat3 = org_boot(X=x3,Y=y3,ntrials=1000l) ;M (slope), MSE, MCI_boot, B, BSE, BCI_Boot, R, RCI_Boot, N            
            print,':FINN',' slope',stat3.M*unitscale,$
                ' SE of slope', stat3.MSE*unitscale,$
                ' CI of slope', stat3.MCI_boot*unitscale, $
                ' intercept', stat3.B, $
                ' SE of intercept', stat3.BSE, $
                ' CI of intercept', stat3.BCI_Boot, $
                ' correlations', stat3.R, $
                ' CI of correlations', stat3.RCI_Boot, $
                ' # of dp', stat3.N

    ;GFAS
            stat4 = org_boot(X=x4,Y=y4,ntrials=1000l) ;M (slope), MSE, MCI_boot, B, BSE, BCI_Boot, R, RCI_Boot, N            
            print,':GFAS',' slope',stat4.M*unitscale,$
                ' SE of slope', stat4.MSE*unitscale,$
                ' CI of slope', stat4.MCI_boot*unitscale, $
                ' intercept', stat4.B, $
                ' SE of intercept', stat4.BSE, $
                ' CI of intercept', stat4.BCI_Boot, $
                ' correlations', stat4.R, $
                ' CI of correlations', stat4.RCI_Boot, $
                ' # of dp', stat4.N
    ;;QFED
            stat5 = org_boot(X=x5,Y=y5,ntrials=1000l) ;M (slope), MSE, MCI_boot, B, BSE, BCI_Boot, R, RCI_Boot, N            
            print,':QFED',' slope',stat5.M*unitscale,$
                ' SE of slope', stat5.MSE*unitscale,$
                ' CI of slope', stat5.MCI_boot*unitscale, $
                ' intercept', stat5.B, $
                ' SE of intercept', stat5.BSE, $
                ' CI of intercept', stat5.BCI_Boot, $
                ' correlations', stat5.R, $
                ' CI of correlations', stat5.RCI_Boot, $
                ' # of dp', stat5.N
            
    ;;THREEGFAS
            stat6 = org_boot(X=x6,Y=y6,ntrials=1000l) ;M (slope), MSE, MCI_boot, B, BSE, BCI_Boot, R, RCI_Boot, N            
            print,':3XGFAS',' slope',stat6.M*unitscale,$
                ' SE of slope', stat6.MSE*unitscale,$
                ' CI of slope', stat6.MCI_boot*unitscale, $
                ' intercept', stat6.B, $
                ' SE of intercept', stat6.BSE, $
                ' CI of intercept', stat6.BCI_Boot, $
                ' correlations', stat6.R, $
                ' CI of correlations', stat6.RCI_Boot, $
                ' # of dp', stat6.N
            
    ;;NOBB
            stat7 = org_boot(X=x7,Y=y7,ntrials=1000l) ;M (slope), MSE, MCI_boot, B, BSE, BCI_Boot, R, RCI_Boot, N            
            print,':NOBB',' slope',stat7.M*unitscale,$
                ' SE of slope', stat7.MSE*unitscale,$
                ' CI of slope', stat7.MCI_boot*unitscale, $
                ' intercept', stat7.B, $
                ' SE of intercept', stat7.BSE, $
                ' CI of intercept', stat7.BCI_Boot, $
                ' correlations', stat7.R, $
                ' CI of correlations', stat7.RCI_Boot, $
                ' # of dp', stat7.N 

;; assign values
            slope_obs = (stat1.MCI_boot[0]+stat1.MCI_boot[1])/2.
            slope_gfed4 = (stat2.MCI_boot[0]+stat2.MCI_boot[1])/2.
            slope_finn = (stat3.MCI_boot[0]+stat3.MCI_boot[1])/2.
            slope_gfas = (stat4.MCI_boot[0]+stat4.MCI_boot[1])/2.
            slope_qfed =(stat5.MCI_boot[0]+stat5.MCI_boot[1])/2.
            slope_threegfas = (stat6.MCI_boot[0]+stat6.MCI_boot[1])/2.
            slope_nobb  = (stat7.MCI_boot[0]+stat7.MCI_boot[1])/2.

            slope_se_obs = (stat1.MCI_boot[1]-stat1.MCI_boot[0])/2.
            slope_se_gfed4 = (stat2.MCI_boot[1]-stat2.MCI_boot[0])/2.
            slope_se_finn = (stat3.MCI_boot[1]-stat3.MCI_boot[0])/2.
            slope_se_gfas = (stat4.MCI_boot[1]-stat4.MCI_boot[0])/2.
            slope_se_qfed = (stat5.MCI_boot[1]-stat5.MCI_boot[0])/2.
            slope_se_threegfas = (stat6.MCI_boot[1]-stat6.MCI_boot[0])/2.
            slope_se_nobb = (stat7.MCI_boot[1]-stat7.MCI_boot[0])/2.

            intercept_obs = (stat1.BCI_boot[0]+stat1.BCI_boot[1])/2.    
            intercept_gfed4 = (stat2.BCI_boot[0]+stat2.BCI_boot[1])/2.
            intercept_finn = (stat3.BCI_boot[0]+stat3.BCI_boot[1])/2.
            intercept_gfas = (stat4.BCI_boot[0]+stat4.BCI_boot[1])/2.
            intercept_qfed = (stat5.BCI_boot[0]+stat5.BCI_boot[1])/2.
            intercept_threegfas = (stat6.BCI_boot[0]+stat6.BCI_boot[1])/2.
            intercept_nobb = (stat7.BCI_boot[0]+stat7.BCI_boot[1])/2.

            corr_obs = (stat1.RCI_boot[0]+stat1.RCI_boot[1])/2.
            corr_gfed4 = (stat2.RCI_boot[0]+stat2.RCI_boot[1])/2.
            corr_finn = (stat3.RCI_boot[0]+stat3.RCI_boot[1])/2.
            corr_gfas = (stat4.RCI_boot[0]+stat4.RCI_boot[1])/2.
            corr_qfed = (stat5.RCI_boot[0]+stat5.RCI_boot[1])/2.
            corr_threegfas = (stat6.RCI_boot[0]+stat6.RCI_boot[1])/2.
            corr_nobb = (stat7.RCI_boot[0]+stat7.RCI_boot[1])/2.

            spec_obs_total   = [spec_obs_total, slope_obs*unitscale] 
            spec_gfed4_total = [spec_gfed4_total, slope_gfed4*unitscale]
            spec_gfas_total  = [spec_gfas_total, slope_gfas*unitscale]
            spec_qfed_total  = [spec_qfed_total, slope_qfed*unitscale] 

            spec_obs_std_total   = [spec_obs_std_total, slope_se_obs*unitscale] 
            spec_gfed4_std_total = [spec_gfed4_std_total, slope_se_gfed4*unitscale]
            spec_gfas_std_total  = [spec_gfas_std_total, slope_se_gfas*unitscale]
            spec_qfed_std_total  = [spec_qfed_std_total, slope_se_qfed*unitscale] 

            corr_obs_total  = [corr_obs_total, corr_obs]
            corr_gfed4_total = [corr_gfed4_total, corr_gfed4]
            corr_gfas_total  = [corr_gfas_total, corr_gfas]
            corr_qfed_total  = [corr_qfed_total, corr_qfed]
            title_all = [title_all, title]
            print, '               '
            help, x1
            
;; Plotting part
            ;if control_ind eq 1 then xrange=[0,10]
            ;ovoc_ptr .VS. oco_picarro
            SCATTERPLOT, X1, Y1,_Extra = _Extra, $
                charsize=2.5, charthick = 4, thick = 10,xrange=xrange;,xrange=[0,200];,yrange=yrange,xrange=xrange;,$
                ;,XTICKFORMAT="(A1)",YTICKFORMAT="(A1)", 
            oplot,!x.crange,intercept_obs+!x.crange*slope_obs,col=1,line=0, thick=10

;; add math symbol, waiting for test!!
;; http://www.idlcoyote.com/ps_tips/greeksym.php
        thisletter = "262B
        greekLetter = '!9' + String(thisLetter) + '!X'

;; Format setting
        if keyword_set(xeqco) then begin
            if keyword_set(all) and (s eq 1 or s eq 3) $
                then begin
                format1 = '(F5.2,F4.1)' ;; XX.XX  formaldehyde(for all), formic acid(), and acetic acid
                format2 = '(F4.2)' ;; X.XX, 
            endif else begin
                format1 = '(F4.2,F4.1)'
                format2 = '(F4.2)'
            endelse
            
            if keyword_set(all) and (keyword_set(emipass) or keyword_set(young) or keyword_set(intermed) or keyword_set(aged)) $
                then begin
                format1 = '(F5.2,F4.1)'
                format2 = '(F4.2)'
                format3 = '(F5.2)'
            endif
            if (keyword_set(all) and (keyword_set(emipass) or keyword_set(young))) and (s eq 1 or s eq 2) $
                then format1 = '(F5.2,F4.1)'

            if (keyword_set(all) and keyword_set(intermed)) and (s eq 1 or s eq  3) $
                then format1 = '(F5.2,F4.1)'

            if (keyword_set(all) and keyword_set(aged)) and (s eq 1 or s eq 3 or s eq 7) $
                then format1 = '(F5.2,F4.1)'

            if (keyword_set(all) and keyword_set(aged)) and (s eq 7) $
                then format2 = '(F5.2)'

            if keyword_set(all) and keyword_set(nonsmoke) and (s eq 0 or s eq 2 or s eq 3 or s eq 4 or s eq 8) then begin
                format1 = '(F5.2,F4.1)'
                format2 = '(F5.2)'
            endif

            if (keyword_set(all) and (keyword_set(emipass) or keyword_set(young))) and (s eq 1 or s eq 2 or s eq 3 or s eq 4) $
                then format1 = '(F5.2,F4.1)'
                
                
            if keyword_set(chemistry) then begin
                format1 = '(F6.1,F4.1)'
                format2 = '(F5.2)'
                format3 = '(F5.2)'
            endif
        endif
        
        if keyword_set(xeqbenz) or keyword_set(xeqtolu) then begin
            if keyword_set(all) and (s eq 1 or s eq 3) $
                then begin
                format1 = '(F5.2,F4.1)'
                format2 = '(F4.2)'
            endif else begin
                format1 = '(F4.2,F4.1)'
                format2 = '(F4.2)'
            endelse

            if keyword_set(all) and (keyword_set(emipass) or keyword_set(young) or keyword_set(intermed) or keyword_set(aged)) $
                then begin
                format1 = '(F4.2,F4.1)'
                format2 = '(F4.2)'
            endif

            if (keyword_set(all) and (keyword_set(emipass) or keyword_set(young))) and (s eq 1 or s eq 2) $
                then format1 = '(F5.2,F4.1)'

            if (keyword_set(all) and keyword_set(intermed)) and (s eq 1 or s eq  3) $
                then format1 = '(F5.2,F4.1)'

            if (keyword_set(all) and keyword_set(aged)) and (s eq 1 or s eq 3 or s eq 7) $
                then format1 = '(F5.2,F4.1)'

            if (keyword_set(all) and keyword_set(aged)) and (s eq 7) $
                then format2 = '(F5.2)'

            if keyword_set(all) and keyword_set(nonsmoke) and (s eq 0 or s eq 2 or s eq 3 or s eq 4 or s eq 8) then begin
                format1 = '(F5.2,F4.1)'
                format2 = '(F5.2)'
            endif

            if keyword_set(all) and (keyword_set(emipass) or keyword_set(young)) and (s eq 1 or s eq 2 or s eq 3 or s eq 4) $
                then format1 = '(F5.2,F4.1)'
            ;; start from here
            if keyword_set(all) and (s eq 2 or s eq 4 or s eq 8) then format1 = '(F5.2,F4.1)'

            if keyword_set(all) and keyword_set(young) and (s eq 1 or s eq 3 or s eq 8) then begin
                format1 = '(F6.2,F4.1)'
                format2 = '(F5.2)'
            endif
            if keyword_set(all) and (keyword_set(intermed) or keyword_set(aged)) and ((s eq 2) or (s eq 3)) then format1 = '(F5.2,F4.1)'
            if keyword_set(all) and keyword_set(aged) and ((s eq 8)) then format1 = '(F5.2,F4.1)'
            if keyword_set(all) and keyword_set(nonsmoke) and ((s eq 1) or (s eq 3)) then format1 = '(F6.2,F4.1)'
            
        endif
            
        

        ;if keyword_set(~chemistry) then format3 = '(F4.2,F4.1)'
    

    
;; setting 
        charsize = 1.5
        charthick = 4.2
        
        thesymbol =   cgSymbol('+-')

;; Test!!! Change charsize for 3x3
        if keyword_set(all) or keyword_set(ovocs) or keyword_set(chemistry) then $
            xyouts,!x.crange(0) + (!x.crange(1)-!x.crange(0))/25,!y.crange(1)*0.90,$
               /data,col=1,'Obs:slope='+string(strtrim(slope_obs*unitscale, 1), format=format1) +thesymbol + string(strtrim(slope_se_obs*unitscale, 1), format='(F4.2)') +', R='+$
               string(strtrim(corr_obs, 1),format=format2),charsize = charsize, charthick = charthick
        if keyword_set(addup) then $
            xyouts,!x.crange(0) + (!x.crange(1)-!x.crange(0))/15,!y.crange(1)*0.90,$
               /data,col=1,'Obs:slope='+string(strtrim(slope_obs*unitscale, 1), format='(F5.1,F4.1)') +thesymbol + string(strtrim(slope_se_obs*unitscale, 1), format=format3) +', R='+$
               string(strtrim(correlate(x1,y1), 1),format='(F4.2)'),charsize = 1.5, charthick = 4

        if keyword_set(lumped) then $
            xyouts,!x.crange(0) + (!x.crange(1)-!x.crange(0))/25,!y.crange(1)*0.90,$
               /data,col=1,'Obs:slope='+string(strtrim(slope_obs*unitscale, 1), format=format1) +thesymbol + string(strtrim(slope_se_obs*unitscale, 1), format=format3) +', R='+$
               string(strtrim(correlate(x1,y1), 1),format=format2),charsize = charsize, charthick = charthick

;; Test for different BB inventoreis for all
        if inventories eq 'gfed4' then begin
            print, 'this GFED4'
            ;gvoc_gfed4 .VS. gvoc_gfed4
            SCATTERPLOT, X2, Y2,_Extra = _Extra,charsize=2,col=2,overplot=1;,$
                ;,XTICKFORMAT="(A1)",YTICKFORMAT="(A1)"
            oplot,!x.crange,intercept_gfed4+!x.crange*slope_gfed4,col=2,line=3, thick=8
            
            SCATTERPLOT, X7, Y7,_Extra = _Extra,charsize=2,col=3,overplot=1;,$
                ;,XTICKFORMAT="(A1)",YTICKFORMAT="(A1)"
            oplot,!x.crange,intercept_nobb+!x.crange*slope_nobb,col=3,line=3, thick=8
            
            
            if keyword_set(all) then begin
                xyouts,!x.crange(0) + (!x.crange(1)-!x.crange(0))/25,!y.crange(1)*0.8,$
                   /data,col=2,'Mod:slope='+string(strtrim(slope_gfed4*unitscale, 1), format=format1) +thesymbol + string(strtrim(slope_se_gfed4*unitscale, 1), format=format3) +', R='+$
                   string(strtrim(correlate(x2,y2), 1),format=format2),charsize = charsize, charthick = charthick

                xyouts,!x.crange(0) + (!x.crange(1)-!x.crange(0))/25,!y.crange(1)*0.7,$
                   /data,col=3,'NOBB:slope='+string(strtrim(slope_nobb*unitscale, 1), format=format1) +thesymbol + string(strtrim(slope_se_nobb*unitscale, 1), format=format3) +', R='+$
                   string(strtrim(correlate(X7, Y7), 1),format=format2),charsize = charsize, charthick = charthick
            endif
            
            if keyword_set(addup) then $
            xyouts,!x.crange(0) + (!x.crange(1)-!x.crange(0))/15,!y.crange(1)*0.85,$
               /data,col=2,'Mod:slope='+string(strtrim(slope_gfed4*unitscale, 1), format='(F5.1,F4.1)') +thesymbol+ string(strtrim(slope_se_gfed4*unitscale, 1), format='(F3.1,F4.1)') +', R='+$
               string(strtrim(correlate(x2,y2), 1),format='(F4.2)'),charsize = 1.5, charthick = 4
               
               
            if keyword_set(lumped) then begin
                xyouts,!x.crange(0) + (!x.crange(1)-!x.crange(0))/25,!y.crange(1)*0.88,$
                   /data,col=2,'Mod:slope='+string(strtrim(slope_gfed4*unitscale, 1), format=format1) +thesymbol + string(strtrim(slope_se_gfed4*unitscale, 1), format=format3) +', R='+$
                   string(strtrim(correlate(x2,y2), 1),format=format2),charsize = charsize, charthick = charthick
                   
                xyouts,!x.crange(0) + (!x.crange(1)-!x.crange(0))/25,!y.crange(1)*0.86,$
                   /data,col=3,'NOBB:slope='+string(strtrim(slope_nobb*unitscale, 1), format=format1) +thesymbol + string(strtrim(slope_se_nobb*unitscale, 1), format=format3) +', R='+$
                   string(strtrim(correlate(X7, Y7), 1),format=format2),charsize = charsize, charthick = charthick        
            endif
        endif
        
        if inventories eq 'finn' then begin
            print, 'this FINN'

            ;gvoc_finn .VS. gvoc_finn
            SCATTERPLOT, X3, Y3,_Extra = _Extra,charsize=2,col=3,overplot=1;,$
                ;,XTICKFORMAT="(A1)",YTICKFORMAT="(A1)"
            oplot,!x.crange,intercept_finn+!x.crange*slope_finn,col=3,line=3, thick=8

            SCATTERPLOT,X7, Y7,_Extra = _Extra,charsize=2,col=3,overplot=1;,$
                ;,XTICKFORMAT="(A1)",YTICKFORMAT="(A1)"
            oplot,!x.crange,intercept_nobb+!x.crange*slope_nobb,col=3,line=3, thick=8

            if keyword_set(all) or keyword_set(ovocs) then $
            xyouts,!x.crange(0) + (!x.crange(1)-!x.crange(0))/25,!y.crange(1)*0.80,$
               /data,col=3,'Mod:slope='+string(strtrim(slope_finn*unitscale, 1), format=format1) +thesymbol+ string(strtrim(slope_se_finn*unitscale, 1), format=format3) +', R='+$
               string(strtrim(correlate(x3,y3), 1),format=format2),charsize = charsize, charthick = charthick          
               
            xyouts,!x.crange(0) + (!x.crange(1)-!x.crange(0))/25,!y.crange(1)*0.7,$
               /data,col=3,'NOBB:slope='+string(strtrim(slope_nobb*unitscale, 1), format=format1) +thesymbol+ string(strtrim(slope_se_nobb*unitscale, 1), format=format3) +', R='+$
               string(strtrim(correlate(X7, Y7), 1),format=format2),charsize = charsize, charthick = charthick    
               
            if keyword_set(addup) then $
            xyouts,!x.crange(0) + (!x.crange(1)-!x.crange(0))/15,!y.crange(1)*0.85,$
               /data,col=2,'Mod:slope='+string(strtrim(slope_finn*unitscale, 1), format='(F5.1,F4.1)') +thesymbol + string(strtrim(slope_se_finn*unitscale, 1), format='(F4.1,F4.1)') +', R='+$
               string(strtrim(correlate(x3,y3), 1),format='(F4.2)'),charsize = 1.5, charthick = 4   
               
            if keyword_set(lumped) then begin
                xyouts,!x.crange(0) + (!x.crange(1)-!x.crange(0))/25,!y.crange(1)*0.88,$
                   /data,col=2,'Mod:slope='+string(strtrim(slope_gfed4*unitscale, 1), format=format1) +thesymbol + string(strtrim(slope_se_gfed4*unitscale, 1), format=format1) +', R='+$
                   string(strtrim(correlate(x2,y2), 1),format=format2),charsize = charsize, charthick = charthick
               
                xyouts,!x.crange(0) + (!x.crange(1)-!x.crange(0))/25,!y.crange(1)*0.86,$
                   /data,col=3,'NOBB:slope='+string(strtrim(slope_nobb*unitscale, 1), format=format1) +thesymbol+ string(strtrim(slope_se_nobb*unitscale, 1), format=format1) +', R='+$
                   string(strtrim(correlate(X7, Y7), 1),format=format2),charsize = charsize, charthick = charthick     
            endif
        endif
        
        if inventories eq 'gfas' then begin
            print, 'this GFAS'
            ;gvoc_gfas .VS. gvoc_gfas
            SCATTERPLOT, X4, Y4,_Extra = _Extra,charsize=2,col=4,overplot=1;,$
                ;,XTICKFORMAT="(A1)",YTICKFORMAT="(A1)"
            ;org_corr,x4,y4,correlate(x4,y4),n_elements(x4),slope,intercept
            oplot,!x.crange,intercept_gfas+!x.crange*slope_gfas,col=4,line=3, thick=8
            
            ;SCATTERPLOT, X7, Y7,_Extra = _Extra,charsize=2,col=3,overplot=1;,$
                ;,XTICKFORMAT="(A1)",YTICKFORMAT="(A1)"
            ;oplot,!x.crange,intercept_nobb+!x.crange*slope_nobb,col=3,line=3, thick=8
            
            ;SCATTERPLOT, X6, Y6,_Extra = _Extra,charsize=2,col=12,overplot=1;,$
                ;,XTICKFORMAT="(A1)",YTICKFORMAT="(A1)"
            ;oplot,!x.crange,intercept_threegfas+!x.crange*slope_threegfas,col=12,line=3, thick=8
            
        if keyword_set(all) or keyword_set(ovocs) or keyword_set(chemistry) then $
            xyouts,!x.crange(0) + (!x.crange(1)-!x.crange(0))/25,!y.crange(1)*0.82,$
               /data,col=4,'GFAS:slope='+string(strtrim(slope_gfas*unitscale, 1), format=format1) +thesymbol + string(strtrim(slope_se_gfas*unitscale, 1),format='(F4.2)') +', R='+$
               string(strtrim(corr_gfas, 1),format=format2),charsize = charsize, charthick = charthick 
               
            ;xyouts,!x.crange(0) + (!x.crange(1)-!x.crange(0))/25,!y.crange(1)*0.74,$
            ;   /data,col=12,'3XGFAS:slope='+string(strtrim(slope_threegfas*unitscale, 1), format=format1) +thesymbol + string(strtrim(slope_se_threegfas*unitscale, 1), format=format3) +', R='+$
            ;   string(strtrim(correlate(X6, Y6), 1),format=format2),charsize = charsize, charthick = charthick   
               
            ;xyouts,!x.crange(0) + (!x.crange(1)-!x.crange(0))/25,!y.crange(1)*0.66,$
            ;   /data,col=3,'NOBB:slope='+string(strtrim(slope_nobb*unitscale, 1), format=format1) +thesymbol+ string(strtrim(slope_se_nobb*unitscale, 1), format=format3) +', R='+$
            ;   string(strtrim(correlate(X7, Y7), 1),format=format2),charsize = charsize, charthick = charthick    
            
            
            if keyword_set(addup) then $
            xyouts,!x.crange(0) + (!x.crange(1)-!x.crange(0))/15,!y.crange(1)*0.85,$
               /data,col=2,'Mod:slope='+string(strtrim(slope_gfas*unitscale, 1), format='(F5.1,F4.1)') +thesymbol + string(strtrim(slope_se_gfas*unitscale, 1), format='(F3.1,F4.1)') +', R='+$
               string(strtrim(correlate(x4,y4), 1),format='(F4.2)'),charsize = 1.5, charthick = 4
               
            if keyword_set(lumped) then begin
                xyouts,!x.crange(0) + (!x.crange(1)-!x.crange(0))/25,!y.crange(1)*0.88,$
                   /data,col=2,'Mod:slope='+string(strtrim(slope_gfas*unitscale, 1), format=format1) +thesymbol + string(strtrim(slope_se_gfas*unitscale, 1), format=format3) +', R='+$
                   string(strtrim(correlate(x2,y2), 1),format=format2),charsize = charsize, charthick = charthick
                   
                xyouts,!x.crange(0) + (!x.crange(1)-!x.crange(0))/25,!y.crange(1)*0.86,$
                   /data,col=3,'NOBB:slope='+string(strtrim(slope_nobb*unitscale, 1), format=format1) +thesymbol+ string(strtrim(slope_se_nobb*unitscale, 1), format=format3) +', R='+$
                   string(strtrim(correlate(X7, Y7), 1),format=format2),charsize = charsize, charthick = charthick     
            endif
               
        endif


        if inventories eq 'qfed' then begin
            print, 'this QFED'
            ;gvoc_qfed .VS. gvoc_qfed
            SCATTERPLOT, X5, Y5,_Extra = _Extra,charsize=2,col=5,overplot=1;,$
                ;,XTICKFORMAT="(A1)",YTICKFORMAT="(A1)"
            ;org_corr,x5,y5,correlate(x5,y5),n_elements(x5),slope,intercept
            oplot,!x.crange,intercept_qfed+!x.crange*slope_qfed,col=5,line=3, thick=8
            
            SCATTERPLOT,X7, Y7,_Extra = _Extra,charsize=2,col=3,overplot=1;,$
                ;,XTICKFORMAT="(A1)",YTICKFORMAT="(A1)"
            oplot,!x.crange,intercept_nobb+!x.crange*slope_nobb,col=3,line=3, thick=8
            
            if keyword_set(all) or keyword_set(ovocs) then $  
            xyouts,!x.crange(0) + (!x.crange(1)-!x.crange(0))/25,!y.crange(1)*0.80,$
               /data,col=5,'Mod:slope='+string(strtrim(slope_qfed*unitscale, 1), format=format1) +thesymbol + string(strtrim(slope_se_qfed*unitscale, 1), format=format3) +', R='+$
               string(strtrim(correlate(x5,y5), 1),format=format2),charsize = charsize, charthick = charthick              
            xyouts,!x.crange(0) + (!x.crange(1)-!x.crange(0))/25,!y.crange(1)*0.7,$
               /data,col=3,'NOBB:slope='+string(strtrim(slope_nobb*unitscale, 1), format=format1) +thesymbol + string(strtrim(slope_se_nobb*unitscale, 1), format=format3) +', R='+$
               string(strtrim(correlate(X7, Y7), 1),format=format2),charsize = charsize, charthick = charthick    
            
            if keyword_set(addup) then $
            xyouts,!x.crange(0) + (!x.crange(1)-!x.crange(0))/15,!y.crange(1)*0.85,$
               /data,col=2,'Mod:slope='+string(strtrim(slope_qfed*unitscale, 1), format='(F5.1,F4.1)') +thesymbol + string(strtrim(slope_se_qfed*unitscale, 1), format='(F3.1,F4.1)') +', R='+$
               string(strtrim(correlate(x5,y5), 1),format='(F5.2)'),charsize = 1.5, charthick = 4  
               
            if keyword_set(lumped) then begin
                xyouts,!x.crange(0) + (!x.crange(1)-!x.crange(0))/25,!y.crange(1)*0.88,$
                   /data,col=2,'Mod:slope='+string(strtrim(slope_gfed4*unitscale, 1), format=format1) +thesymbol + string(strtrim(slope_se_gfed4*unitscale, 1), format=format3) +', R='+$
                   string(strtrim(correlate(x2,y2), 1),format=format2),charsize = charsize, charthick = charthick
                
                xyouts,!x.crange(0) + (!x.crange(1)-!x.crange(0))/25,!y.crange(1)*0.86,$
                   /data,col=3,'NOBB:slope='+string(strtrim(slope_nobb*unitscale, 1), format=format1) +thesymbol+ string(strtrim(slope_se_nobb*unitscale, 1), format=format3) +', R='+$
                   string(strtrim(correlate(X7, Y7), 1),format=format2),charsize = charsize, charthick = charthick     
            endif
        endif


;; add bar
        xpos1 = 0.78-0.01-0.6-0.05-0.03+0.5+0.05
        ypos1 = [0.95,0.9,0.85,0.8,0.75,0.7]-0.1


        endfor ;; s loop
        
        help, title_all, spec_obs_total,     spec_gfed4_total,     spec_gfas_total,     spec_qfed_total    
        help, title_all, spec_obs_std_total, spec_gfed4_std_total, spec_gfas_std_total, spec_qfed_std_total    
        help, title_all, corr_obs_total,     corr_gfed4_total,     corr_gfas_total,     corr_qfed_total 


        
        if keyword_set(all) then begin
        ;; Lixu, pressure
            xyouts,0.03,0.46,'VOCs(ppb)',/normal, col = 1, charsize=1.8, $
            charthick=4, alignment=0., ORIENTATION = 90
            
            if keyword_set(xeqco) then $
            xyouts,0.48,0.01,'CO(ppb)',/normal, col = 1, charsize=1.8, $
            charthick=4, alignment=0., ORIENTATION = 0
            
            if keyword_set(xeqbenz) then $
            xyouts,0.48,0.01,'Benzene(ppb)',/normal, col = 1, charsize=1.8, $
            charthick=4, alignment=0., ORIENTATION = 0
            
            if keyword_set(xeqtolu) then $
            xyouts,0.48,0.01,'Toluene(ppb)',/normal, col = 1, charsize=1.8, $
            charthick=4, alignment=0., ORIENTATION = 0
        endif
    endfor ;; dd loop
close_device
;DEVICE, /CLOSE
print,'done!'
end
