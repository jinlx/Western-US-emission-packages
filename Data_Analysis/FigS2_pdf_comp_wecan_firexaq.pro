pro pdf_comp_wecan_firexaq, save=save
    
;; ==============
;; save the file
;; ==============
    if keyword_set(save) then $
    open_device, /ps, /color, /landscape, $
               fi = './ps/wecan_firexaq_pdf.ps'
    DEVICE,  /LANDSCAPE

;; =============
;; READ WE-CAN
;; =============
    dates_wecan =  ['20180724','20180726','20180730','20180731','20180802',$
                 '20180803','20180806','20180808','20180809','20180813','20180815','20180816','20180820',$
                 '20180823','20180826','20180828','20180906','20180910','20180913']
    for dd = 0, n_elements(dates_wecan) do begin
        ;;only consider whole flight
        if dd ne 0 then break
        LATITUDE_wecan = [0]
        HCN_TOGA_wecan = [0]
        HCN_UWCIMS_wecan = [0]
        CH3CN_TOGA_wecan = [0]
        CH3CN_PTR_wecan = [0]
    
    ; smart loop for each fligiht and total flights from lu
        if dd eq 0 then dates = dates_wecan else dates=dates_wecan[dd-1]
        dates = STRTRIM(string(dates),1)
        for n=0,n_elements(dates)-1 do begin
    ;; GC file for all data
            str_date=dates[n]
            print,'processing '+str_date
            c130fi = '/glade/scratch/lixujin/wecan_obs/mrg01data/wecan-mrg01-c130_merge_'+dates[n]+'_R0.sav'
            ;c130fi = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
            ;    'mrg2sav/R4_merges_avg/wecan-mrg5m-c130_merge_'+dates[n]+'_R4.sav'
            ;; get 1 sec data!
            restore, c130fi
            
            tmp_LATITUDE_wecan = c130.LATITUDE
            tmp_HCN_TOGA_wecan = c130.HCN_TOGA
            tmp_HCN_UWCIMS_wecan = c130.HCN_UWCIMS
            tmp_CH3CN_TOGA_wecan = c130.CH3CN_TOGA
            tmp_CH3CN_PTR_wecan = c130.Acetonitrile_MixingRatio_PTR
        
    ;; prepare data for 20180826, lat is less than zero and it is not simulated by the obs,fixed
    ;; add lat, add HCN
            i = where(tmp_LATITUDE_wecan le 0,ct) 
            if ct gt 0 then begin
                tmp_LATITUDE_wecan[i] = !VALUES.F_NAN
                remove_ind=where(finite(tmp_LATITUDE_wecan),ct)
                print, ct
                if remove_ind[0] ne -1 then begin                    
                    tmp_HCN_TOGA_wecan = tmp_HCN_TOGA_wecan[remove_ind]
                    tmp_HCN_UWCIMS_wecan = tmp_HCN_UWCIMS_wecan[remove_ind]
                    tmp_CH3CN_TOGA_wecan = tmp_CH3CN_TOGA_wecan[remove_ind]
                    tmp_CH3CN_PTR_wecan = tmp_CH3CN_PTR_wecan[remove_ind]
                endif
            endif
            HCN_TOGA_wecan = [HCN_TOGA_wecan, tmp_HCN_TOGA_wecan]
            HCN_UWCIMS_wecan = [HCN_UWCIMS_wecan, tmp_HCN_UWCIMS_wecan]
            CH3CN_TOGA_wecan = [CH3CN_TOGA_wecan, tmp_CH3CN_TOGA_wecan]
            CH3CN_PTR_wecan = [CH3CN_PTR_wecan, tmp_CH3CN_PTR_wecan]
            
            undefine,c130
        endfor ;; loop for dates, read files
        
        HCN_TOGA_wecan = HCN_TOGA_wecan[1:*]
        HCN_UWCIMS_wecan = HCN_UWCIMS_wecan[1:*]
        CH3CN_TOGA_wecan = CH3CN_TOGA_wecan[1:*]        
        CH3CN_PTR_wecan = CH3CN_PTR_wecan[1:*]

    endfor

;; =============
;; READ FIREX-AQ
;; =============
    dates_firexaq =  ['20190722', '20190724', '20190725', '20190729', '20190730',$
                    '20190802', '20190803', '20190806', '20190807', '20190808', $
                    '20190812', '20190813', '20190815', '20190816', '20190819', $
                    '20190821', '20190823', '20190826', '20190830', '20190831', $ 
                    '20190903']
    for dd = 0, n_elements(dates_firexaq) do begin
        ;;only consider whole flight
        if dd ne 0 then break
        LATITUDE_firexaq = [0]
        LONGITUDE_firexaq = [0]
        
        HCN_obs_toga_firexaq = [0]
        HCN_obs_noaaptr_firexaq = [0]
        HCN_obs_cims_firexaq = [0]

        CH3CN_obs_toga_firexaq = [0]
        CH3CN_obs_iwas_firexaq = [0]
        CH3CN_obs_noaaptr_firexaq = [0]
    ; smart loop for each fligiht and total flights from lu
        if dd eq 0 then dates = dates_firexaq else dates=dates_firexaq[dd-1]
        dates = STRTRIM(string(dates),1)
        for n=0,n_elements(dates)-1 do begin
    ;; GC file for all data
            str_date=dates[n]
            print,'processing '+str_date
            
            dc8fi = '/glade/scratch/lixujin/firexaq_obs/firexaq-mrg01-dc8_merge_'+dates[n]+'_R1.sav'
            ;dc8fi = '/glade/u/home/lixujin/IDL/SEAC4RS-IDL-master/' + $
            ;    'mrg2sav/firexaq_avg/firexaq-mrg5m-dc8_merge_'+dates[n]+'_RL.sav'
            
            
            restore, dc8fi

            tmp_LATITUDE_firexaq   = dc8.Latitude_YANG
            tmp_LONGITUDE_firexaq  = dc8.Longitude_YANG

            tmp_HCN_obs_toga_firexaq = dc8.HCN_TOGA_APEL/1E3
            ;tmp_HCN_obs_noaaptr_firexaq = dc8.HCN_NOAAPTR_ppbv_WARNEKE
            tmp_HCN_obs_cims_firexaq = dc8.HCN_NOAACIMS_VERES/1e3

            tmp_CH3CN_obs_toga_firexaq = dc8.CH3CN_TOGA_APEL/1E3
            tmp_CH3CN_obs_iwas_firexaq = dc8.CH3CN_NOAAiWAS_GILMAN
            ;tmp_CH3CN_obs_noaaptr_firexaq = dc8.CH3CN_NOAAPTR_ppbv_WARNEKE
            
            
            
            LATITUDE_firexaq   = [LATITUDE_firexaq, tmp_LATITUDE_firexaq]
            LONGITUDE_firexaq  =  [LONGITUDE_firexaq, tmp_LONGITUDE_firexaq]
            HCN_obs_toga_firexaq = [HCN_obs_toga_firexaq, tmp_HCN_obs_toga_firexaq]
            ;HCN_obs_noaaptr_firexaq = [HCN_obs_noaaptr_firexaq, tmp_HCN_obs_noaaptr_firexaq]
            HCN_obs_cims_firexaq = [HCN_obs_cims_firexaq, tmp_HCN_obs_cims_firexaq]

            CH3CN_obs_toga_firexaq = [CH3CN_obs_toga_firexaq, tmp_CH3CN_obs_toga_firexaq]
            CH3CN_obs_iwas_firexaq = [CH3CN_obs_iwas_firexaq, tmp_CH3CN_obs_iwas_firexaq]
            ;CH3CN_obs_noaaptr_firexaq = [CH3CN_obs_noaaptr_firexaq, tmp_CH3CN_obs_noaaptr_firexaq]
            undefine,dc8
        endfor ;; loop for dates, read files
        
        
        LATITUDE_firexaq = LATITUDE_firexaq[1:*]
        LONGITUDE_firexaq = LONGITUDE_firexaq[1:*]
        
        HCN_obs_toga_firexaq = HCN_obs_toga_firexaq[1:*]
        ;HCN_obs_noaaptr_firexaq = HCN_obs_noaaptr_firexaq[1:*]
        HCN_obs_cims_firexaq = HCN_obs_cims_firexaq[1:*]
        
        CH3CN_obs_toga_firexaq = CH3CN_obs_toga_firexaq[1:*]
        CH3CN_obs_iwas_firexaq = CH3CN_obs_iwas_firexaq[1:*]
        ;CH3CN_obs_noaaptr_firexaq = CH3CN_obs_noaaptr_firexaq[1:*]
    endfor
    ;; only select western US for FIREXAQ
    ind_lon = where(LONGITUDE_firexaq le -105, ct)
    
    HCN_obs_toga_firexaq = HCN_obs_toga_firexaq[ind_lon]
    HCN_obs_cims_firexaq = HCN_obs_cims_firexaq[ind_lon]
    CH3CN_obs_toga_firexaq = CH3CN_obs_toga_firexaq[ind_lon]
    CH3CN_obs_iwas_firexaq = CH3CN_obs_iwas_firexaq[ind_lon]
    
    ;; data points GT 0
    ind = where(HCN_TOGA_wecan > 0 and $
                HCN_UWCIMS_wecan > 0 and $
                CH3CN_TOGA_wecan > 0 and $
                CH3CN_PTR_wecan > 0)
    HCN_TOGA_wecan = HCN_TOGA_wecan[ind] ; ppt
    HCN_UWCIMS_wecan = HCN_UWCIMS_wecan[ind] ; ppt
    CH3CN_TOGA_wecan = CH3CN_TOGA_wecan[ind]; ppt
    CH3CN_PTR_wecan = CH3CN_PTR_wecan[ind]*1000 ; ppb2ppt
    
    ind = where(HCN_obs_toga_firexaq > 0 and $
                HCN_obs_cims_firexaq > 0 and $
                CH3CN_obs_toga_firexaq > 0 )
    
    HCN_obs_toga_firexaq = HCN_obs_toga_firexaq[ind]*1000 ; ppb2ppt
    HCN_obs_cims_firexaq = HCN_obs_cims_firexaq[ind]*1000; ppb2ppt
    CH3CN_obs_toga_firexaq = CH3CN_obs_toga_firexaq[ind]*1000 ; ppb2ppt
    CH3CN_obs_iwas_firexaq = CH3CN_obs_iwas_firexaq[ind]*1000 ; ppb2ppt
    
    ;help, HCN_TOGA_wecan, HCN_UWCIMS_wecan
    ;help, CH3CN_TOGA_wecan, CH3CN_PTR_wecan
    ;help, HCN_obs_toga_firexaq, HCN_obs_cims_firexaq
    ;help, CH3CN_obs_toga_firexaq, CH3CN_obs_iwas_firexaq
; ----------------------------------------------------------------------
; For histograrm plot
; we choose HCN_UWCIMS_wecan, HCN_obs_cims_firexaq, CH3CN_PTR_wecan, CH3CN_obs_toga_firexaq
; we use ppt as unit
; ----------------------------------------------------------------------
;print,HCN_UWCIMS_wecan[0], HCN_TOGA_wecan[0], HCN_obs_toga_firexaq[0], HCN_obs_cims_firexaq[0]
;print, CH3CN_TOGA_wecan[0], CH3CN_PTR_wecan[0], CH3CN_obs_toga_firexaq[0], CH3CN_obs_iwas_firexaq[0]


    data_1 = CH3CN_PTR_wecan
    data_2 = CH3CN_obs_toga_firexaq
;; set percentile
    percentile_wecan = cgPercentiles(data_1, Percentiles=0.5)
    percentile_firexaq = cgPercentiles(data_2, Percentiles=0.5)
    print, percentile_wecan, percentile_firexaq
    
;; set binsize
    binsize = 50
    xrange = [0, 1000]
    ;yrange = [0, 30000]
    yrange = [0, 50]
    charsize = 2
    charthick = 5
    thick = 5
    
    ;cgHistoplot, CH3CN_PTR_wecan, Binsize=binsize, /Fill, XRange=[0, 1000], yrange = [0, 30000], POLYCOLOR='yellow'
    ;cgHistoplot, CH3CN_obs_toga_firexaq, Binsize=binsize, /OPlot,POLYCOLOR='olive',/FILL;, /Fill, XRange=[-25, -4]
    
    
    ;;; first
    ;cgHistoplot, data_2, POLYCOLOR='olive', /FILL, $
    ;   XRange=[0, 1000], yrange = [0, 30000], BINSIZE=binsize
       
    ;cgHistoplot, data_1, POLYCOLOR='yellow', /FILL, /OPLOT, BINSIZE=binsize
    
    ;firstPlot = cgSnapshot() 

    ;;; second
    ;cgHistoplot, data_1, POLYCOLOR='yellow', /FILL, $
    ;   XRange=[0, 1000], yrange = [0, 30000], BINSIZE=binsize
       
    ;cgHistoplot, data_2, /OPLOT, POLYCOLOR='olive', /FILL, BINSIZE=binsize
    
    ;secondPlot = cgSnapshot() 
    
    ;;; merge the first one and the second one
    ;cgBlendimage, firstPlot, secondPlot, ALPHA=0.75 
    
    y1 = HISTOGRAM(data_1, binsize=binsize, Locations=xbin1)
    y2 = HISTOGRAM(data_2, binsize=binsize, Locations=xbin2)

    cgPlot, xbin1, y1/total(y1)*100, COLOR='black', xrange=xrange, yrange=yrange, ytitle='Frequency (%)' ,xtitle='CH3CN (ppt)', thick = thick, charsize=charsize, charthick=charthick
    cgPlot, xbin2, y2/total(y2)*100, /OVERPLOT, COLOR='red', thick=thick
    
   ; Draw a red line at the true speed of light.
    cgPlot, [percentile_wecan, percentile_wecan], !Y.CRange,  Color='black', Thick=thick, /OVERPLOT, line = 3
    cgPlot, [percentile_firexaq, percentile_firexaq], !Y.CRange, Color='red', Thick=thick, /OVERPLOT, line = 3
    
    close_device
    ;DEVICE, /CLOSE
    print,'done!'
end
