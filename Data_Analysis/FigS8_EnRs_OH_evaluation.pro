;; modified by lixu, 10/17/2019, used for compare the profiles from wecan and gc 
;; original script mk_flight_v2, lhu
;; modicication: parameters, set missing value for observation data, cancel filter. 
;@bin_vert
pro EnRs_OH_evaluation,$
    plume=plume, test=test,save=save

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
    ;!P.CHARTHICK=2.5
    !P.CHARTHICK=1
    ;P.THICK=2.5
    !P.THICK=1
    ;!X.THICK=4
    !X.THICK=2
    ;!Y.THICK=4
    !Y.THICK=4
    PEdges = [0.5,1,1.5,2.,2.5,3.,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12]
    avo= 6.022e23  ;; molec/mol

;;setting
    if keyword_set(test) then fi = 'test_emipass'  
    if keyword_set(save) then $
    open_device, /ps, /color, bits=8, $
       filename='./ps/'+ fi + '.ps'

;///////////////////////////////////////////////////////////////////////////////////////////
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;++                                      STEP 2                                           ++
;++                       Read merge file and plane log files                             ++
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;///////////////////////////////////////////////////////////////////////////////////////////                                  
    avo= 6.022e23  ;; molec/mol
    ;; ==============
    ;; WECAN DATA
    ;; ==============
    dates_all=  ['20180724','20180726','20180730','20180731','20180802',$
                 '20180803','20180806','20180808','20180809','20180813','20180815','20180816','20180820',$
                 '20180823','20180826','20180828','20180906','20180910','20180913']    
    benz_obs_ptr_wecan = [0]
    tolu_obs_ptr_wecan = [0]
    benz_obs_toga_wecan = [0]
    tolu_obs_toga_wecan = [0]
    benz_gc_gfas_wecan = [0]
    tolu_gc_gfas_wecan = [0]
    
    for dd = 0, n_elements(dates_all) do begin
        ;;only consider whole flight
        if dd ne 0 then break
        if dd eq 0 then dates = dates_all else dates=dates_all[dd-1]
        for n=0,n_elements(dates)-1 do begin
            str_date=dates[n]
            print,'processing '+str_date
            ;; observation
            c130fi = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                'mrg2sav/WECAN/R4_merges_avg/wecan-mrg5m-c130_merge_'+dates[n]+'_R4.sav'
            restore, c130fi
            tmp_benz_obs_ptr_wecan  = c130.Benzene_MixingRatio_PTR
            tmp_tolu_obs_ptr_wecan  = c130.Toluene_MixingRatio_PTR
            tmp_benz_obs_toga_wecan = c130.Benzene_TOGA/1e3
            tmp_tolu_obs_toga_wecan = c130.Toluene_TOGA/1e3
            
            benz_obs_ptr_wecan = [benz_obs_ptr_wecan, tmp_benz_obs_ptr_wecan]
            tolu_obs_ptr_wecan = [tolu_obs_ptr_wecan, tmp_tolu_obs_ptr_wecan]
            benz_obs_toga_wecan = [benz_obs_toga_wecan, tmp_benz_obs_toga_wecan]
            tolu_obs_toga_wecan = [tolu_obs_toga_wecan, tmp_tolu_obs_toga_wecan]
            
            undefine,c130
            ;; GEOS-Chem + GFAS
            gcfi_gfas   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                'planelog2sav/output_'+'gfas_tmp' + '/mrg5m_wecan_c130_'+dates[n]+'.sav'  
            restore, gcfi_gfas  
            tmp_benz_gc_gfas_wecan   = gc.benz*1e9/6; 6C benzene in ppb
            tmp_tolu_gc_gfas_wecan   = gc.tolu*1e9/7; 7C toluene in ppb
            
            benz_gc_gfas_wecan = [benz_gc_gfas_wecan, tmp_benz_gc_gfas_wecan]
            tolu_gc_gfas_wecan = [tolu_gc_gfas_wecan, tmp_tolu_gc_gfas_wecan]
            undefine,gc

        endfor ;; loop for dates
    endfor
;; ======================
;; FIREX-AQ observations
;; ======================
    dates_all=  ['20190722', '20190724', '20190725', '20190729', '20190730',$
                    '20190802', '20190803', '20190806', '20190807', '20190808', $
                    '20190812', '20190813', '20190815', '20190816', '20190819', $
                    '20190821', '20190823', '20190826', '20190830', '20190831', $ 
                    '20190903']
    Benzene_obs_noaaptr_firexaq = [0]
    Toluene_obs_noaaptr_firexaq = [0]
    Benzene_obs_toga_firexaq    = [0]
    Toluene_obs_toga_firexaq    = [0]
    
    Benzene_gc_gfas_firexaq = [0]
    Toluene_gc_gfas_firexaq = [0]
    
    lat_obs_firexaq  = [0]  
    lon_obs_firexaq  = [0]
    
    for dd = 0, n_elements(dates_all) do begin
        ;;only consider whole flight
        if dd ne 0 then break
        if dd eq 0 then dates = dates_all else dates=dates_all[dd-1]
        for n=0,n_elements(dates)-1 do begin
            str_date=dates[n]
            print,'processing '+str_date
            ; observation
            dc8fi = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                'mrg2sav/FIREX-AQ/firexaq_avg/firexaq-mrg5m-dc8_merge_'+dates[n]+'_R1.sav'
            restore, dc8fi
            tmp_Benzene_obs_noaaptr_firexaq = dc8.Benzene_PTR_WARNEKE_WISTHALER
            tmp_Toluene_obs_noaaptr_firexaq = dc8.Toluene_PTR_WARNEKE_WISTHALER
            tmp_Benzene_obs_toga_firexaq    = dc8.Benzene_TOGA_APEL/1E3
            tmp_Toluene_obs_toga_firexaq    = dc8.Toluene_TOGA_APEL/1E3
            tmp_lat_obs_firexaq = dc8.Latitude_YANG
            tmp_lon_obs_firexaq = dc8.Longitude_YANG

            Benzene_obs_noaaptr_firexaq = [Benzene_obs_noaaptr_firexaq, tmp_Benzene_obs_noaaptr_firexaq]
            Toluene_obs_noaaptr_firexaq = [Toluene_obs_noaaptr_firexaq, tmp_Toluene_obs_noaaptr_firexaq]
            Benzene_obs_toga_firexaq    = [Benzene_obs_toga_firexaq, tmp_Benzene_obs_toga_firexaq]
            Toluene_obs_toga_firexaq    = [Toluene_obs_toga_firexaq, tmp_Toluene_obs_toga_firexaq]
            lat_obs_firexaq  = [lat_obs_firexaq, tmp_lat_obs_firexaq]  
            lon_obs_firexaq  = [lon_obs_firexaq, tmp_lon_obs_firexaq]
            undefine,dc8fi
            
            ; GEOS-Chem + GFAS
            gcfi_gfas   = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
                'planelog2sav/firexaq_test' + '/mrg5m_firexaq_dc8_'+dates[n]+'.sav' 
            restore, gcfi_gfas  
            tmp_Benzene_gc_gfas_firexaq = gc.benz*1e9/6; 6C benzene in ppb
            tmp_Toluene_gc_gfas_firexaq = gc.tolu*1e9/7; 7C toluene in ppb
            Benzene_gc_gfas_firexaq = [Benzene_gc_gfas_firexaq, tmp_Benzene_gc_gfas_firexaq]
            Toluene_gc_gfas_firexaq = [Toluene_gc_gfas_firexaq, tmp_Toluene_gc_gfas_firexaq]
            undefine,gc            
        endfor ;; loop for dates
    endfor
    ;; Remove placeholder
    benz_obs_ptr_wecan = benz_obs_ptr_wecan[1:*]
    tolu_obs_ptr_wecan = tolu_obs_ptr_wecan[1:*]
    benz_obs_toga_wecan = benz_obs_toga_wecan[1:*]
    tolu_obs_toga_wecan = tolu_obs_toga_wecan[1:*]
    benz_gc_gfas_wecan = benz_gc_gfas_wecan[1:*]
    tolu_gc_gfas_wecan = tolu_gc_gfas_wecan[1:*]
    
    Benzene_obs_noaaptr_firexaq = Benzene_obs_noaaptr_firexaq[1:*]
    Toluene_obs_noaaptr_firexaq = Toluene_obs_noaaptr_firexaq[1:*]
    Benzene_obs_toga_firexaq = Benzene_obs_toga_firexaq[1:*]
    Toluene_obs_toga_firexaq = Toluene_obs_toga_firexaq[1:*]
    Benzene_gc_gfas_firexaq = Benzene_gc_gfas_firexaq[1:*]
    Toluene_gc_gfas_firexaq = Toluene_gc_gfas_firexaq[1:*]
    lat_obs_firexaq = lat_obs_firexaq[1:*]
    lon_obs_firexaq = lon_obs_firexaq[1:*]
    
    if n_elements(benz_obs_ptr_wecan) ne n_elements(tolu_obs_ptr_wecan) or $
        n_elements(benz_obs_ptr_wecan) or n_elements(benz_obs_toga_wecan) or $
        n_elements(benz_obs_ptr_wecan) or n_elements(tolu_obs_toga_wecan) or $
        n_elements(benz_obs_ptr_wecan) or n_elements(benz_gc_gfas_wecan) or $
        n_elements(benz_obs_ptr_wecan) or n_elements(tolu_gc_gfas_wecan) then stop
    
    if n_elements(Benzene_obs_noaaptr_firexaq) ne n_elements(Toluene_obs_noaaptr_firexaq) or $
        n_elements(Benzene_obs_noaaptr_firexaq) or n_elements(Benzene_obs_toga_firexaq) or $
        n_elements(Benzene_obs_noaaptr_firexaq) or n_elements(Toluene_obs_toga_firexaq) or $
        n_elements(Benzene_obs_noaaptr_firexaq) or n_elements(Benzene_gc_gfas_firexaq) or $
        n_elements(Benzene_obs_noaaptr_firexaq) or n_elements(Toluene_gc_gfas_firexaq) or $
        n_elements(Benzene_obs_noaaptr_firexaq) or n_elements(lat_obs_firexaq) or $
        n_elements(Benzene_obs_noaaptr_firexaq) or n_elements(lon_obs_firexaq) then stop
    
    
    ;; LoD
    LoD_benz = 30.0/1000 
    LoD_tolu = 30.0/1000 ; change it into 30 ppt as well to get more data points otherwise the FIREX-AQ will be untrustful due to limited datapoints. The LoD is loosely defined
    ind = where(benz_obs_ptr_wecan le LoD_benz and benz_obs_ptr_wecan gt 0, ct)
    if ct gt 0 then begin
        benz_obs_ptr_wecan[ind] = !VALUES.F_NAN
        ;benz_obs_ptr_wecan[ind] = benz_obs_toga_wecan[ind]
    endif
    ind = where(tolu_obs_ptr_wecan le LoD_tolu and tolu_obs_ptr_wecan gt 0, ct)
    if ct gt 0 then begin
        tolu_obs_ptr_wecan[ind] = !VALUES.F_NAN
        ;tolu_obs_ptr_wecan[ind] = tolu_obs_toga_wecan[ind]
    endif
    ind = where(Benzene_obs_noaaptr_firexaq le LoD_benz and Benzene_obs_noaaptr_firexaq gt 0, ct)
    if ct gt 0 then begin
        Benzene_obs_noaaptr_firexaq[ind] = !VALUES.F_NAN
        ;Benzene_obs_noaaptr_firexaq[ind] = Benzene_obs_toga_firexaq[ind]
    endif
    ind = where(Toluene_obs_noaaptr_firexaq le LoD_tolu and Toluene_obs_noaaptr_firexaq gt 0, ct)
    if ct gt 0 then begin
        Toluene_obs_noaaptr_firexaq[ind] = !VALUES.F_NAN
        ;Toluene_obs_noaaptr_firexaq[ind] = Toluene_obs_toga_firexaq[ind]
    endif
    
    ; low ends model is not doing okay
    ind = where(benz_gc_gfas_wecan le LoD_benz or tolu_gc_gfas_wecan le LoD_tolu, ct)
    if ct gt 0 then begin
        benz_obs_ptr_wecan[ind] = !VALUES.F_NAN
        tolu_obs_ptr_wecan[ind] = !VALUES.F_NAN
    endif
    ind = where(Benzene_gc_gfas_firexaq le LoD_benz or Toluene_gc_gfas_firexaq le LoD_tolu, ct)

    if ct gt 0 then begin
        Benzene_obs_noaaptr_firexaq[ind] = !VALUES.F_NAN
        Toluene_obs_noaaptr_firexaq[ind] = !VALUES.F_NAN
    endif
    

    ;; filter out missing values or eastern values
    ind = where(benz_obs_ptr_wecan le 0 or tolu_obs_ptr_wecan le 0, ct) 
    if ct gt 0 then benz_obs_ptr_wecan[ind] = !VALUES.F_NAN
    ind = where(Benzene_obs_noaaptr_firexaq le 0 or Toluene_obs_noaaptr_firexaq le 0 or lon_obs_firexaq ge -105, ct) 
    if ct gt 0 then Benzene_obs_noaaptr_firexaq[ind] = !VALUES.F_NAN
    
    ;; remove missing values
    remove_ind = where(finite(benz_obs_ptr_wecan) and finite(tolu_obs_ptr_wecan), ct)
    if ct gt 0 then begin
        benz_obs_ptr_wecan = benz_obs_ptr_wecan[remove_ind]
        tolu_obs_ptr_wecan = tolu_obs_ptr_wecan[remove_ind]
        benz_gc_gfas_wecan = benz_gc_gfas_wecan[remove_ind]
        tolu_gc_gfas_wecan = tolu_gc_gfas_wecan[remove_ind]
    endif

    remove_ind = where(finite(Benzene_obs_noaaptr_firexaq) and finite(Toluene_obs_noaaptr_firexaq), ct)
    if ct gt 0 then begin
        Benzene_obs_noaaptr_firexaq = Benzene_obs_noaaptr_firexaq[remove_ind]
        Toluene_obs_noaaptr_firexaq = Toluene_obs_noaaptr_firexaq[remove_ind]
        Benzene_gc_gfas_firexaq = Benzene_gc_gfas_firexaq[remove_ind]
        Toluene_gc_gfas_firexaq = Toluene_gc_gfas_firexaq[remove_ind]
    endif

;///////////////////////////////////////////////////////////////////////////////////////////
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;++                                      STEP 3                                           ++
;++                                  Stat. Analysis                                       ++
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;///////////////////////////////////////////////////////////////////////////////////////////
    unitscale = 1

    ;; make it log scale
    stat_wecan_obs = org_boot(X=alog(benz_obs_ptr_wecan),Y=alog(tolu_obs_ptr_wecan),ntrials=1000l) ;M (slope), MSE, MCI_boot, B, BSE, BCI_Boot, R, RCI_Boot, N            
    print,':WE-CAN obs',' slope',stat_wecan_obs.M*unitscale, $
            ' SE of slope', stat_wecan_obs.MSE*unitscale, $
            ' CI of slope', stat_wecan_obs.MCI_boot*unitscale, $
            ' intercept', stat_wecan_obs.B, $
            ' SE of intercept', stat_wecan_obs.BSE, $
            ' CI of intercept', stat_wecan_obs.BCI_Boot, $
            ' correlations', stat_wecan_obs.R, $
            ' CI of correlations', stat_wecan_obs.RCI_Boot, $
            ' # of dp', stat_wecan_obs.N
        
    stat_wecan_mod = org_boot(X=alog(benz_gc_gfas_wecan),Y=alog(tolu_gc_gfas_wecan),ntrials=1000l) ;M (slope), MSE, MCI_boot, B, BSE, BCI_Boot, R, RCI_Boot, N            
    print,':WE-CAN mod',' slope',stat_wecan_mod.M*unitscale,$
            ' SE of slope', stat_wecan_mod.MSE*unitscale,$
            ' CI of slope', stat_wecan_mod.MCI_boot*unitscale, $
            ' intercept', stat_wecan_mod.B, $
            ' SE of intercept', stat_wecan_mod.BSE, $
            ' CI of intercept', stat_wecan_mod.BCI_Boot, $
            ' correlations', stat_wecan_mod.R, $
            ' CI of correlations', stat_wecan_mod.RCI_Boot, $
            ' # of dp', stat_wecan_mod.N
            
    stat_firexaq_obs = org_boot(X=alog(Benzene_obs_noaaptr_firexaq),Y=alog(Toluene_obs_noaaptr_firexaq),ntrials=10000l) ;M (slope), MSE, MCI_boot, B, BSE, BCI_Boot, R, RCI_Boot, N            
    ;stat_firexaq_obs = org_boot(X=(Benzene_obs_noaaptr_firexaq),Y=(Toluene_obs_noaaptr_firexaq),ntrials=1000l) ;M (slope), MSE, MCI_boot, B, BSE, BCI_Boot, R, RCI_Boot, N            
    print,':FIREX-AQ obs',' slope',stat_firexaq_obs.M*unitscale,$
            ' SE of slope', stat_firexaq_obs.MSE*unitscale,$
            ' CI of slope', stat_firexaq_obs.MCI_boot*unitscale, $
            ' intercept', stat_firexaq_obs.B, $
            ' SE of intercept', stat_firexaq_obs.BSE, $
            ' CI of intercept', stat_firexaq_obs.BCI_Boot, $
            ' correlations', stat_firexaq_obs.R, $
            ' CI of correlations', stat_firexaq_obs.RCI_Boot, $
            ' # of dp', stat_firexaq_obs.N
            
    stat_firexaq_mod = org_boot(X=alog(Benzene_gc_gfas_firexaq),Y=alog(Toluene_gc_gfas_firexaq),ntrials=10000l) ;M (slope), MSE, MCI_boot, B, BSE, BCI_Boot, R, RCI_Boot, N           
    ;stat_firexaq_mod = org_boot(X=(Benzene_gc_gfas_firexaq),Y=(Toluene_gc_gfas_firexaq),ntrials=1000l) ;M (slope), MSE, MCI_boot, B, BSE, BCI_Boot, R, RCI_Boot, N            
    print,':FIREX-AQ mod',' slope',stat_firexaq_mod.M*unitscale,$
            ' SE of slope', stat_firexaq_mod.MSE*unitscale,$
            ' CI of slope', stat_firexaq_mod.MCI_boot*unitscale, $
            ' intercept', stat_firexaq_mod.B, $
            ' SE of intercept', stat_firexaq_mod.BSE, $
            ' CI of intercept', stat_firexaq_mod.BCI_Boot, $
            ' correlations', stat_firexaq_mod.R, $
            ' CI of correlations', stat_firexaq_mod.RCI_Boot, $
            ' # of dp', stat_firexaq_mod.N
        
;; assign values
    slope_wecan_obs = (stat_wecan_obs.MCI_boot[0]+stat_wecan_obs.MCI_boot[1])/2.
    slope_wecan_mod = (stat_wecan_mod.MCI_boot[0]+stat_wecan_mod.MCI_boot[1])/2.
    slope_firexaq_obs = (stat_firexaq_obs.MCI_boot[0]+stat_firexaq_obs.MCI_boot[1])/2.
    slope_firexaq_mod = (stat_firexaq_mod.MCI_boot[0]+stat_firexaq_mod.MCI_boot[1])/2.
    
    slope_se_wecan_obs = (stat_wecan_obs.MCI_boot[1]-stat_wecan_obs.MCI_boot[0])/2.
    slope_se_wecan_mod = (stat_wecan_mod.MCI_boot[1]-stat_wecan_mod.MCI_boot[0])/2.
    slope_se_firexaq_obs = (stat_firexaq_obs.MCI_boot[1]-stat_firexaq_obs.MCI_boot[0])/2.
    slope_se_firexaq_mod = (stat_firexaq_mod.MCI_boot[1]-stat_firexaq_mod.MCI_boot[0])/2.
                
    intercept_wecan_obs = (stat_wecan_obs.BCI_boot[0]+stat_wecan_obs.BCI_boot[1])/2.    
    intercept_wecan_mod = (stat_wecan_mod.BCI_boot[0]+stat_wecan_mod.BCI_boot[1])/2.    
    intercept_firexaq_obs = (stat_firexaq_obs.BCI_boot[0]+stat_firexaq_obs.BCI_boot[1])/2.                
    intercept_firexaq_mod = (stat_firexaq_mod.BCI_boot[0]+stat_firexaq_mod.BCI_boot[1])/2.                
                
    corr_wecan_obs = (stat_wecan_obs.RCI_boot[0]+stat_wecan_obs.RCI_boot[1])/2.
    corr_wecan_mod = (stat_wecan_mod.RCI_boot[0]+stat_wecan_mod.RCI_boot[1])/2.
    corr_firexaq_obs = (stat_firexaq_obs.RCI_boot[0]+stat_firexaq_obs.RCI_boot[1])/2.
    corr_firexaq_mod = (stat_firexaq_mod.RCI_boot[0]+stat_firexaq_mod.RCI_boot[1])/2.

;///////////////////////////////////////////////////////////////////////////////////////////
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;++                                      STEP 5                                           ++
;++                                  DO THE PLOTTING                                      ++
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;///////////////////////////////////////////////////////////////////////////////////////////
            
    ;; add math symbol, waiting for test!!
    ;; http://www.idlcoyote.com/ps_tips/greeksym.php
    thisletter = "262B
    greekLetter = '!9' + String(thisLetter) + '!X'
    thesymbol =   cgSymbol('+-')

    ;; plot setting
    format1 = '(F4.1,F4.1)'
    format2 = '(F4.1)'
    
    charsize  = 2
    charthick = 6
    thick = 10
    
    !x.crange = [1e-2,1e1]
    range =  [1e-2,1e1]
    ;; WECAN
    SCATTERPLOT, benz_obs_ptr_wecan, tolu_obs_ptr_wecan,_Extra = _Extra, $
            charsize=charsize, charthick = charthick, thick = thick, /xlog, /ylog, xrange = range, yrange = range;,xrange=[0,200];,yrange=yrange,xrange=xrange;,$
    oplot,!x.crange,intercept_wecan_obs+!x.crange*slope_wecan_obs,col=1,line=0, thick=thick
    
    ;xyouts,!x.crange(0) + (!x.crange(1)-!x.crange(0))/25,!y.crange(1)*0.90,$
    ;        /data,col=1,'Obs:slope='+string(strtrim(slope_wecan_obs*unitscale, 1), format=format1) +thesymbol + string(strtrim(slope_se_wecan_obs*unitscale, 1), format='(F4.2)') +', R='+$
    ;        string(strtrim(corr_wecan_obs, 1),format=format2),charsize = charsize, charthick = charthick
               
    SCATTERPLOT, benz_gc_gfas_wecan, tolu_gc_gfas_wecan,_Extra = _Extra,charsize=charsize,col=4,overplot=1, /xlog, /ylog;,$;,$
    oplot,!x.crange,intercept_wecan_mod+!x.crange*slope_wecan_mod,col=4,line=3, thick=thick
    
    ;xyouts,!x.crange(0) + (!x.crange(1)-!x.crange(0))/25,!y.crange(1)*0.82,$
    ;       /data,col=4,'GFAS:slope='+string(strtrim(slope_wecan_mod*unitscale, 1), format=format1) +thesymbol + string(strtrim(slope_se_wecan_mod*unitscale, 1),format='(F4.2)') +', R='+$
    ;       string(strtrim(corr_wecan_mod, 1),format=format2),charsize = charsize, charthick = charthick 
        
        
    ;; FIREXAQ
    SCATTERPLOT, Benzene_obs_noaaptr_firexaq, Toluene_obs_noaaptr_firexaq, _Extra = _Extra, $
            charsize=charsize, charthick = charthick, thick = thick, /xlog, /ylog, xrange = range, yrange = range;,xrange=[0,200];,yrange=yrange,xrange=xrange;,$
    oplot,!x.crange,intercept_firexaq_obs+!x.crange*slope_firexaq_obs,col=1,line=0, thick=thick
    
    ;xyouts,!x.crange(0) + (!x.crange(1)-!x.crange(0))/25,!y.crange(1)*0.90,$
    ;        /data,col=1,'Obs:slope='+string(strtrim(slope_firexaq_obs*unitscale, 1), format=format1) +thesymbol + string(strtrim(slope_se_firexaq_obs*unitscale, 1), format='(F4.2)') +', R='+$
    ;        string(strtrim(corr_firexaq_obs, 1),format=format2),charsize = charsize, charthick = charthick
        
    SCATTERPLOT, Benzene_gc_gfas_firexaq, Toluene_gc_gfas_firexaq,_Extra = _Extra,charsize=charsize,col=4,overplot=1, /xlog, /ylog
    oplot,!x.crange,intercept_firexaq_mod+!x.crange*slope_firexaq_mod,col=4,line=3, thick=thick
    
    ;xyouts,!x.crange(0) + (!x.crange(1)-!x.crange(0))/25,!y.crange(1)*0.82,$
    ;       /data,col=4,'GFAS:slope='+string(strtrim(slope_firexaq_mod*unitscale, 1), format=format1) +thesymbol + string(strtrim(slope_se_firexaq_mod*unitscale, 1),format='(F4.2)') +', R='+$
    ;       string(strtrim(corr_firexaq_mod, 1),format=format2),charsize = charsize, charthick = charthick 

close_device
;DEVICE, /CLOSE
print,'done!'
end
