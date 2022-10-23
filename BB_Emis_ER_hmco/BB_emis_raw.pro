
pro BB_emis_raw, unit=unit, region=region, fireseason=fireseason, fireyear=fireyear, year=year
;; set reigon
    if n_elements(region) eq 0 then region=''
    case strlowcase(region) of 
        'southeast': limit=[25,-100,40,-75]
        'se'       : limit=[25,-100,40,-75]
        'wus'     : limit=[36,-127,49.5,-105]
        ;lixu, 10/10/2019
        'na'       : limit=[10,-140,70,-70]    
        else       : limit=[25,-127,52,-65]
    endcase
    
    if keyword_set(fireseason) then months  = ['06', '07','08','09']
    if keyword_set(fireyear) then months  = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
    
    ;; CO along with 14 VOCs
    compounds = ['CO', 'C2H6', 'C3H8', 'ALK4', 'PRPE', $
                'CH2O', 'ALD2', 'ACET', 'MEK', $
                'BENZ', 'TOLU', 'XYLE', $
                'FA', 'ACTA', 'RCHO']

    ratio_mass2Cmass = [12/28.0, 24/28.0, 36/44.0, 48/58.12, 36/42.09, $
                            12/30.0, 24/44.06, 36/58.09, 48/72.11, $
                            72/78.0, 84/92.15, 96/106.18, $
                            12/46.03, 24/60.06, 36/58.09]
;; =================================================================
;; GFED4: emissions += DM_emissions * contribution * EF_CO[source]
;; =================================================================
    print, 'This is GFED4!'
    ;; Presetting for EFs of CO and VOCs
    ;; see emissions_diagnositics excel file
    ; raw data of the EF factor of each compound
    EF_CO = [63,127,88,93,210,102]
    EF_C2H6 = [0.66,1.79,0.63,0.71,0.71,0.91]
    EF_C3H8 = [0.1,0.44,0.22,0.126,0.126,0.28]
    EF_ALK4 = [0.133,0.385,0.369,0.267,0.267,0.333]
    EF_PRPE = [0.79 + 0.133, 1.13 + 0.385, 0.61 + 0.369, 0.64 + 0.267, 3.05 + 0.267, 0.68 + 0.333]
    EF_CH2O = [0.73,1.86,2.09,1.73,1.4,2.08]
    EF_ALD2 = [0.57,0.77,0.77,1.55,3.27,1.24]
    EF_ACET = [0.16,0.75,0.54,0.63,1.25,0.45]    
    EF_MEK  = [0.181,0.22,0.13, 0.5, 0.5, 0.9]
    EF_BENZ       = [0.2,1.11,0.27,0.39,3.19,0.15]
    EF_TOLU       = [0.08,0.48,0.19,0.26,1.55,0.19]
    EF_XYLE       = [0.014,0.18,0.13,0.11,0.11,0.114]
    EF_FA         = [0.21, 0.57, 0.28, 0.79, 0.38, 1]
    EF_ACTA       = [3.55, 4.41, 2.13, 3.05, 8.97, 5.59]
    EF_RCHO       = [0.25, 0.25, 0.25, 0.25, 0.25, 0.25] ; USE PERMAR's REPORT
    
    ;; make a in-order compilation group by each fuel type
    EF_SAVA = [EF_CO[0], EF_C2H6[0], EF_C3H8[0], EF_ALK4[0], EF_PRPE[0], $
                EF_CH2O[0], EF_ALD2[0], EF_ACET[0], EF_MEK[0], $
                EF_BENZ[0], EF_TOLU[0], EF_XYLE[0], $
                EF_FA[0], EF_ACTA[0], EF_RCHO[0]]
                
    EF_BORF = [EF_CO[1], EF_C2H6[1], EF_C3H8[1], EF_ALK4[1], EF_PRPE[1], $
                EF_CH2O[1], EF_ALD2[1], EF_ACET[1], EF_MEK[1], $
                EF_BENZ[1], EF_TOLU[1], EF_XYLE[1], $
                EF_FA[1], EF_ACTA[1], EF_RCHO[1]]

    EF_TEMP = [EF_CO[2], EF_C2H6[2], EF_C3H8[2], EF_ALK4[2], EF_PRPE[2], $
                EF_CH2O[2], EF_ALD2[2], EF_ACET[2], EF_MEK[2], $
                EF_BENZ[2], EF_TOLU[2], EF_XYLE[2], $
                EF_FA[2], EF_ACTA[2], EF_RCHO[2]]

    EF_DEFO = [EF_CO[3], EF_C2H6[3], EF_C3H8[3], EF_ALK4[3], EF_PRPE[3], $
                EF_CH2O[3], EF_ALD2[3], EF_ACET[3], EF_MEK[3], $
                EF_BENZ[3], EF_TOLU[3], EF_XYLE[3], $
                EF_FA[3], EF_ACTA[3], EF_RCHO[3]]
    EF_PEAT = [EF_CO[4], EF_C2H6[4], EF_C3H8[4], EF_ALK4[4], EF_PRPE[4], $
                EF_CH2O[4], EF_ALD2[4], EF_ACET[4], EF_MEK[4], $
                EF_BENZ[4], EF_TOLU[4], EF_XYLE[4], $
                EF_FA[4], EF_ACTA[4], EF_RCHO[4]]
                
    EF_AGRI = [EF_CO[5], EF_C2H6[5], EF_C3H8[5], EF_ALK4[5], EF_PRPE[5], $
                EF_CH2O[5], EF_ALD2[5], EF_ACET[5], EF_MEK[5], $
                EF_BENZ[5], EF_TOLU[5], EF_XYLE[5], $
                EF_FA[5], EF_ACTA[5], EF_RCHO[5]]
    ;===============================
    ;Get index for different region
    ;===============================
    ; read area for 025x025 resolution
    fi = '/glade/work/lixujin/rundirs/tools/area_025x025.nc'
    ncdf_read, tmp,fi=fi, VARIABLES = 'AREA'
    area = tmp.area
    lat = tmp.lat
    lon = tmp.lon

    ;GLB
    index_lat_s_glb = 0
    index_lat_e_glb = -1
    index_lon_s_glb = 0
    index_lon_e_glb = -1

    ;NA
    index_lat1 = where(lat gt 9.75 and lat lt 70, ct1)
    index_lon1 = where(lon gt -140 and lon lt -70, ct2)        
    if ct1 gt 0 and ct2 gt 2 then begin
        index_lat_s_na = index_lat1[0]
        index_lat_e_na = n_elements(index_lat1)-1+ index_lat1[0]
        index_lon_s_na = index_lon1[0]
        index_lon_e_na = n_elements(index_lon1)-1+ index_lon1[0]
    endif
    
    ;Western US
    index_lat1 = where(lat gt 36 and lat lt 49.5, ct1)
    index_lon1 = where(lon gt -127 and lon lt -105, ct2)        
    if ct1 gt 0 and ct2 gt 2 then begin
        index_lat_s_wus = index_lat1[0]
        index_lat_e_wus = n_elements(index_lat1)-1+ index_lat1[0]
        index_lon_s_wus = index_lon1[0]
        index_lon_e_wus = n_elements(index_lon1)-1+ index_lon1[0]
    endif

    ;;CONUS
    indext_lat1 = where(lat gt 25 and lat lt 49, ct1)
    index_lon1  = where(lon gt -125 and lon lt -67, ct2)
    if ct1 gt 0 and ct2 gt 2 then begin
        index_lat_s_conus = index_lat1[0]
        index_lat_e_conus = n_elements(index_lat1)-1+ index_lat1[0]
        index_lon_s_conus = index_lon1[0]
        index_lon_e_conus = n_elements(index_lon1)-1+ index_lon1[0]
    endif 
    
    if region eq 'wus' then begin
        index_lat_s = index_lat_s_wus
        index_lat_e = index_lat_e_wus
        index_lon_s = index_lon_s_wus
        index_lon_e = index_lon_e_wus
        area_target = area[index_lon_s:index_lon_e,index_lat_s:index_lat_e]
        
        lat_target = lat[index_lat_s:index_lat_e]
        lon_target = lon[index_lon_s:index_lon_e]
    endif
    if region eq 'glb' then begin
        index_lat_s = index_lat_s_glb
        index_lat_e = index_lat_e_glb
        index_lon_s = index_lon_s_glb
        index_lon_e = index_lon_e_glb
        area_target = area[index_lon_s:index_lon_e,index_lat_s:index_lat_e]
        lat_target = lat[index_lat_s:index_lat_e]
        lon_target = lon[index_lon_s:index_lon_e]
    endif
    ;; ==========================
    ;; read DM data 
    ;; ==========================
    dir = '/glade/p/univ/uumm0001/geos-chem/ExtData/HEMCO/GFED4/v2020-02/' + year + '/GFED4_gen.025x025.' + year
    ;; get total DM for fire seasons
    for mm = 0, n_elements(months) - 1 do begin
        ;; get number of days
        if months[mm] eq '01' or months[mm] eq '03' or months[mm] eq '05' or months[mm] eq '07' or months[mm] eq '08' or months[mm] eq '10' or months[mm] eq '12' then n_days = 31
        if months[mm] eq '04' or months[mm] eq '06' or months[mm] eq '09' or months[mm] eq '11' then n_days = 30
        if months[mm] eq '02' then n_days = 28
        ;; read each month data
        print,'Month:', months[mm]
        ncdf_read, tmp,fi=dir + months[mm]+'.nc', /all
        ; read the DM data: kg/m2/s
        ; get total DM for each fuel type: kg
        if mm eq 0 then begin
            DM_SAVA = total(tmp.DM_SAVA[index_lon_s:index_lon_e,index_lat_s:index_lat_e] * area_target * n_days * 24 * 60 * 60)
            DM_BORF = total(tmp.DM_BORF[index_lon_s:index_lon_e,index_lat_s:index_lat_e] * area_target * n_days * 24 * 60 * 60)
            DM_TEMP = total(tmp.DM_TEMP[index_lon_s:index_lon_e,index_lat_s:index_lat_e] * area_target * n_days * 24 * 60 * 60)
            DM_DEFO = total(tmp.DM_DEFO[index_lon_s:index_lon_e,index_lat_s:index_lat_e] * area_target * n_days * 24 * 60 * 60)
            DM_PEAT = total(tmp.DM_PEAT[index_lon_s:index_lon_e,index_lat_s:index_lat_e] * area_target * n_days * 24 * 60 * 60)
            DM_AGRI = total(tmp.DM_AGRI[index_lon_s:index_lon_e,index_lat_s:index_lat_e] * area_target * n_days * 24 * 60 * 60)
        endif else begin
            DM_SAVA = total(tmp.DM_SAVA[index_lon_s:index_lon_e,index_lat_s:index_lat_e] * area_target * n_days * 24 * 60 * 60) + DM_SAVA
            DM_BORF = total(tmp.DM_BORF[index_lon_s:index_lon_e,index_lat_s:index_lat_e] * area_target * n_days * 24 * 60 * 60) + DM_BORF
            DM_TEMP = total(tmp.DM_TEMP[index_lon_s:index_lon_e,index_lat_s:index_lat_e] * area_target * n_days * 24 * 60 * 60) + DM_TEMP
            DM_DEFO = total(tmp.DM_DEFO[index_lon_s:index_lon_e,index_lat_s:index_lat_e] * area_target * n_days * 24 * 60 * 60) + DM_DEFO
            DM_PEAT = total(tmp.DM_PEAT[index_lon_s:index_lon_e,index_lat_s:index_lat_e] * area_target * n_days * 24 * 60 * 60) + DM_PEAT
            DM_AGRI = total(tmp.DM_AGRI[index_lon_s:index_lon_e,index_lat_s:index_lat_e] * area_target * n_days * 24 * 60 * 60) + DM_AGRI
        endelse
    endfor
    
    ;; get the total emissions of each compound for fire seasons
    Emis_SAVA = MAKE_ARRAY(n_elements(compounds), /FLOAT, VALUE = 0)
    Emis_BORF = MAKE_ARRAY(n_elements(compounds), /FLOAT, VALUE = 0)
    Emis_TEMP = MAKE_ARRAY(n_elements(compounds), /FLOAT, VALUE = 0)
    Emis_DEFO = MAKE_ARRAY(n_elements(compounds), /FLOAT, VALUE = 0)
    Emis_PEAT = MAKE_ARRAY(n_elements(compounds), /FLOAT, VALUE = 0)
    Emis_AGRI = MAKE_ARRAY(n_elements(compounds), /FLOAT, VALUE = 0)
    for comp=0, n_elements(compounds) - 1 do begin
        Emis_SAVA[comp] = EF_SAVA[comp] * DM_SAVA ; g/kg DM * kg DM -> g 
        Emis_BORF[comp] = EF_BORF[comp] * DM_BORF
        Emis_TEMP[comp] = EF_TEMP[comp] * DM_TEMP
        Emis_DEFO[comp] = EF_DEFO[comp] * DM_DEFO
        Emis_PEAT[comp] = EF_PEAT[comp] * DM_PEAT
        Emis_AGRI[comp] = EF_AGRI[comp] * DM_AGRI
    endfor    
    
    ; conver the unit from g to Gg or GgC
    if unit eq 'mass' then begin
        Emis_total_GFED4 = (Emis_SAVA + Emis_BORF + Emis_TEMP + Emis_DEFO + Emis_PEAT + Emis_AGRI)/1e9
    endif 
    if unit eq 'carbon' then begin
        Emis_total_GFED4 = (Emis_SAVA + Emis_BORF + Emis_TEMP + Emis_DEFO + Emis_PEAT + Emis_AGRI)/1e9*ratio_mass2Cmass
    endif 
    
;; =====
;; GFAS
;; =====
    print, 'This is GFAS!'
    dir = '/glade/u/home/lixujin/project/GEOS-CHEM/ExtData/HEMCO/GFAS/v2018-09/' + year + '/'

    variables = ['cofire', 'c2h6fire', 'c3h8fire', 'hialkanesfire', $
                'c3h6fire', 'hialkenesfire', $
                'ch3ohfire', 'c2h5ohfire', $
                'c6h6fire', 'c7h8fire','c8h10fire', $
                'ch2ofire', 'c2h4ofire', $
                'c3h6ofire']
    ;===============================
    ;Get index for different region
    ;===============================
    ; read in area data
    ncdf_read,tmp,fi='/glade/work/lixujin/rundirs/tools/area_01x01.nc', VARIABLES=['area']
    area = tmp.area
    lon = tmp.lon
    lat = tmp.lat

    ;; get lat and lon for NA and western US
    ;NA
    index_lat1 = where(lat gt 9.75 and lat lt 70, ct1)
    index_lon1 = where(lon gt -140 and lon lt -70, ct2)        
    if ct1 gt 0 and ct2 gt 2 then begin
        index_lat_s_na = index_lat1[0]
        index_lat_e_na = n_elements(index_lat1)-1+ index_lat1[0]
        index_lon_s_na = index_lon1[0]
        index_lon_e_na = n_elements(index_lon1)-1+ index_lon1[0]
    endif
    
    ;Western US
    index_lat1 = where(lat gt 36 and lat lt 49.5, ct1)
    index_lon1 = where(lon gt -127 and lon lt -105, ct2)        
    if ct1 gt 0 and ct2 gt 2 then begin
        index_lat_s_wus = index_lat1[0]
        index_lat_e_wus = n_elements(index_lat1)-1+ index_lat1[0]
        index_lon_s_wus = index_lon1[0]
        index_lon_e_wus = n_elements(index_lon1)-1+ index_lon1[0]
    endif
    
    ;;CONUS
    indext_lat1 = where(lat gt 25 and lat lt 49, ct1)
    index_lon1  = where(lon gt -125 and lon lt -67, ct2)
    if ct1 gt 0 and ct2 gt 2 then begin
        index_lat_s_conus = index_lat1[0]
        index_lat_e_conus = n_elements(index_lat1)-1+ index_lat1[0]
        index_lon_s_conus = index_lon1[0]
        index_lon_e_conus = n_elements(index_lon1)-1+ index_lon1[0]
    endif 

    if region eq 'wus' then begin
        index_lat_s = index_lat_s_wus
        index_lat_e = index_lat_e_wus
        index_lon_s = index_lon_s_wus
        index_lon_e = index_lon_e_wus
        area_target = area[index_lon_s:index_lon_e,index_lat_s:index_lat_e]
        
        lat_target = lat[index_lat_s:index_lat_e]
        lon_target = lon[index_lon_s:index_lon_e]
    endif
    if region eq 'glb' then begin
        index_lat_s = index_lat_s_glb
        index_lat_e = index_lat_e_glb
        index_lon_s = index_lon_s_glb
        index_lon_e = index_lon_e_glb
        area_target = area[index_lon_s:index_lon_e,index_lat_s:index_lat_e]
        lat_target = lat[index_lat_s:index_lat_e]
        lon_target = lon[index_lon_s:index_lon_e]
    endif

    ;; read in daily data but make it into montly resolution
    for mm = 0, n_elements(months)-1 do begin
        if months[mm] eq '01' or months[mm] eq '03' or months[mm] eq '05' or months[mm] eq '07' or months[mm] eq '08' or months[mm] eq '10' or months[mm] eq '12' then n_days = 31
        if months[mm] eq '04' or months[mm] eq '06' or months[mm] eq '09' or months[mm] eq '11' then n_days = 30
        if months[mm] eq '02' then n_days = 28
        print,'Month:', months[mm]
        ncdf_read,tmp,fi=dir+'GFAS_' + year + months[mm] + '_monthly_grid.nc'  , VARIABLES=variables
        if mm eq 0 then begin
            ;; convert kg/m2/s in to Gg
            co   = tmp.cofire[index_lon_s:index_lon_e,index_lat_s:index_lat_e]*60*60*24*n_days * area_target * 1e-6
            c2h6 = tmp.c2h6fire[index_lon_s:index_lon_e,index_lat_s:index_lat_e]*60*60*24*n_days * area_target * 1e-6
            c3h8 = tmp.c3h8fire[index_lon_s:index_lon_e,index_lat_s:index_lat_e]*60*60*24*n_days * area_target * 1e-6
            alk4 = tmp.hialkanesfire[index_lon_s:index_lon_e,index_lat_s:index_lat_e]*60*60*24*n_days *area_target * 1e-6
            c3h6 = tmp.c3h6fire[index_lon_s:index_lon_e,index_lat_s:index_lat_e]*60*60*24*n_days * area_target * 1e-6
            prpe = tmp.hialkenesfire[index_lon_s:index_lon_e,index_lat_s:index_lat_e]*60*60*24*n_days * area_target * 1e-6
            benz = tmp.c6h6fire[index_lon_s:index_lon_e,index_lat_s:index_lat_e]*60*60*24*n_days * area_target * 1e-6
            tolu = tmp.c7h8fire[index_lon_s:index_lon_e,index_lat_s:index_lat_e]*60*60*24*n_days * area_target * 1e-6
            xyle = tmp.c8h10fire[index_lon_s:index_lon_e,index_lat_s:index_lat_e]*60*60*24*n_days * area_target * 1e-6
            hcho = tmp.ch2ofire[index_lon_s:index_lon_e,index_lat_s:index_lat_e]*60*60*24*n_days * area_target * 1e-6
            ald2 = tmp.c2h4ofire[index_lon_s:index_lon_e,index_lat_s:index_lat_e]*60*60*24*n_days * area_target * 1e-6
            acet = tmp.c3h6ofire[index_lon_s:index_lon_e,index_lat_s:index_lat_e]*60*60*24*n_days * area_target * 1e-6
        endif else begin
            co   = co + tmp.cofire[index_lon_s:index_lon_e,index_lat_s:index_lat_e]*60*60*24*n_days * area_target * 1e-6
            c2h6 = c2h6 + tmp.c2h6fire[index_lon_s:index_lon_e,index_lat_s:index_lat_e]*60*60*24*n_days * area_target * 1e-6
            c3h8 = c3h8 + tmp.c3h8fire[index_lon_s:index_lon_e,index_lat_s:index_lat_e]*60*60*24*n_days * area_target * 1e-6
            alk4 = alk4 + tmp.hialkanesfire[index_lon_s:index_lon_e,index_lat_s:index_lat_e]*60*60*24*n_days * area_target * 1e-6
            c3h6 = c3h6 + tmp.c3h6fire[index_lon_s:index_lon_e,index_lat_s:index_lat_e]*60*60*24*n_days * area_target * 1e-6
            prpe = prpe + tmp.hialkenesfire[index_lon_s:index_lon_e,index_lat_s:index_lat_e]*60*60*24*n_days * area_target * 1e-6
            benz = benz + tmp.c6h6fire[index_lon_s:index_lon_e,index_lat_s:index_lat_e]*60*60*24*n_days * area_target * 1e-6
            tolu = tolu + tmp.c7h8fire[index_lon_s:index_lon_e,index_lat_s:index_lat_e]*60*60*24*n_days * area_target * 1e-6
            xyle = xyle + tmp.c8h10fire[index_lon_s:index_lon_e,index_lat_s:index_lat_e]*60*60*24*n_days * area_target * 1e-6
            hcho = hcho + tmp.ch2ofire[index_lon_s:index_lon_e,index_lat_s:index_lat_e]*60*60*24*n_days * area_target * 1e-6
            ald2 = ald2 + tmp.c2h4ofire[index_lon_s:index_lon_e,index_lat_s:index_lat_e]*60*60*24*n_days * area_target * 1e-6
            acet = acet + tmp.c3h6ofire[index_lon_s:index_lon_e,index_lat_s:index_lat_e]*60*60*24*n_days * area_target * 1e-6
        endelse
    endfor
    
    Emissions = [total(co), total(c2h6), total(c3h8), total(alk4), total(c3h6) + total(prpe), $
                 total(hcho), total(ald2), total(acet), total(co)*0.73*(1e-3)/28.01*72.11, $
                 total(benz), total(tolu), total(xyle), $
                 total(co)*46/28*9.5*1E-3, total(co)*11.48*0.75*(1e-3)/28.01*60.06, total(co)*1.01*(1e-3)*58.09/28]
                    
    if unit eq 'mass' then begin
        Emis_total_GFAS = Emissions
    endif 
    if unit eq 'carbon' then begin
        Emis_total_GFAS = Emissions*ratio_mass2Cmass
    endif
    
;====
;QFED
;====
    print, 'This is QFED!'

    dir = '/glade/p/univ/uumm0001/geos-chem/ExtData/HEMCO/QFED/v2018-07/' + year 
    ;===============================
    ;Get index for different region
    ;===============================
    
    ; read in area data
    ncdf_read,tmp,fi='/glade/work/lixujin/rundirs/tools/area_01x01.nc', VARIABLES=['area']
    area = tmp.area
    lon = tmp.lon
    lat = tmp.lat

    ;; get lat and lon for NA and western US
    ;NA
    index_lat1 = where(lat gt 9.75 and lat lt 70, ct1)
    index_lon1 = where(lon gt -140 and lon lt -70, ct2)        
    if ct1 gt 0 and ct2 gt 2 then begin
        index_lat_s_na = index_lat1[0]
        index_lat_e_na = n_elements(index_lat1)-1+ index_lat1[0]
        index_lon_s_na = index_lon1[0]
        index_lon_e_na = n_elements(index_lon1)-1+ index_lon1[0]
    endif
    
    ;Western US
    index_lat1 = where(lat gt 36 and lat lt 49.5, ct1)
    index_lon1 = where(lon gt -127 and lon lt -105, ct2)        
    if ct1 gt 0 and ct2 gt 2 then begin
        index_lat_s_wus = index_lat1[0]
        index_lat_e_wus = n_elements(index_lat1)-1+ index_lat1[0]
        index_lon_s_wus = index_lon1[0]
        index_lon_e_wus = n_elements(index_lon1)-1+ index_lon1[0]
    endif
    
    ;;CONUS
    indext_lat1 = where(lat gt 25 and lat lt 49, ct1)
    index_lon1  = where(lon gt -125 and lon lt -67, ct2)
    if ct1 gt 0 and ct2 gt 2 then begin
        index_lat_s_conus = index_lat1[0]
        index_lat_e_conus = n_elements(index_lat1)-1+ index_lat1[0]
        index_lon_s_conus = index_lon1[0]
        index_lon_e_conus = n_elements(index_lon1)-1+ index_lon1[0]
    endif 

    if region eq 'wus' then begin
        index_lat_s = index_lat_s_wus
        index_lat_e = index_lat_e_wus
        index_lon_s = index_lon_s_wus
        index_lon_e = index_lon_e_wus
        area_target = area[index_lon_s:index_lon_e,index_lat_s:index_lat_e]
        
        lat_target = lat[index_lat_s:index_lat_e]
        lon_target = lon[index_lon_s:index_lon_e]
    endif
    if region eq 'glb' then begin
        index_lat_s = index_lat_s_glb
        index_lat_e = index_lat_e_glb
        index_lon_s = index_lon_s_glb
        index_lon_e = index_lon_e_glb
        area_target = area[index_lon_s:index_lon_e,index_lat_s:index_lat_e]
        lat_target = lat[index_lat_s:index_lat_e]
        lon_target = lon[index_lon_s:index_lon_e]
    endif

    
    ;; read in daily data but make it into montly resolution
    for mm = 0, n_elements(months)-1 do begin
        if months[mm] eq '01' or months[mm] eq '03' or months[mm] eq '05' or months[mm] eq '07' or months[mm] eq '08' or months[mm] eq '10' or months[mm] eq '12' then begin
            n_days = 31
            days = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', $
                        '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31']
        endif
        if months[mm] eq '04' or months[mm] eq '06' or months[mm] eq '09' or months[mm] eq '11' then begin 
            n_days = 30
            days = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', $
                        '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30']
        endif
        if months[mm] eq '02' then begin
            n_days = 28
            days = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', $
                        '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28']
        endif
        print,'Month:', months[mm]
        for dd = 0, n_elements(days)-1 do begin
            ;; unit: kg s-1 m-2
            ncdf_read, tmp_co, fi=dir + '/' + months[mm] + '/' + 'qfed2.emis_' + 'co' + '.006.' + year + months[mm] + days[dd] + '.nc4', VARIABLES=['biomass']
            if mm eq 0 and dd eq 0 then begin
                co  = tmp_co.biomass
            endif else begin
                co = co + tmp_co.biomass
            endelse

            ncdf_read, tmp_c2h6, fi=dir + '/' + months[mm] + '/' + 'qfed2.emis_' + 'c2h6' + '.006.' + year + months[mm] + days[dd] + '.nc4', VARIABLES=['biomass']
            if mm eq 0 and dd eq 0 then begin
                c2h6  = tmp_c2h6.biomass
            endif else begin
                c2h6 = c2h6 + tmp_c2h6.biomass
            endelse

            ncdf_read, tmp_c3h8, fi=dir + '/' + months[mm] + '/' + 'qfed2.emis_' + 'c3h8' + '.006.' + year + months[mm] + days[dd] + '.nc4', VARIABLES=['biomass']
            if mm eq 0 and dd eq 0 then begin
                c3h8  = tmp_c3h8.biomass
            endif else begin
                c3h8 = c3h8 + tmp_c3h8.biomass
            endelse

            ncdf_read, tmp_alk4, fi=dir + '/' + months[mm] + '/' + 'qfed2.emis_' + 'alk4' + '.006.' + year + months[mm] + days[dd] + '.nc4', VARIABLES=['biomass']
            if mm eq 0 and dd eq 0 then begin
                alk4  = tmp_alk4.biomass
            endif else begin
                alk4 = alk4 + tmp_alk4.biomass
            endelse

            ncdf_read, tmp_c3h6, fi=dir + '/' + months[mm] + '/' + 'qfed2.emis_' + 'c3h6' + '.006.' + year + months[mm] + days[dd] + '.nc4', VARIABLES=['biomass']
            if mm eq 0 and dd eq 0 then begin
                c3h6  = tmp_c3h6.biomass
            endif else begin
                c3h6 = c3h6 + tmp_c3h6.biomass
            endelse

            ncdf_read, tmp_ch2o, fi=dir + '/' + months[mm] + '/' + 'qfed2.emis_' + 'ch2o' + '.006.' + year + months[mm] + days[dd] + '.nc4', VARIABLES=['biomass']
            if mm eq 0 and dd eq 0 then begin
                ch2o  = tmp_ch2o.biomass
            endif else begin
                ch2o = ch2o + tmp_ch2o.biomass
            endelse

            ncdf_read, tmp_ald2, fi=dir + '/' + months[mm] + '/' + 'qfed2.emis_' + 'ald2' + '.006.' + year + months[mm] + days[dd] + '.nc4', VARIABLES=['biomass']
            if mm eq 0 and dd eq 0 then begin
                ald2  = tmp_ald2.biomass
            endif else begin
                ald2 = ald2 + tmp_ald2.biomass
            endelse

            ncdf_read, tmp_acet, fi=dir + '/' + months[mm] + '/' + 'qfed2.emis_' + 'acet' + '.006.' + year + months[mm] + days[dd] + '.nc4', VARIABLES=['biomass']
            if mm eq 0 and dd eq 0 then begin
                acet  = tmp_acet.biomass
            endif else begin
                acet = acet + tmp_acet.biomass
            endelse
        endfor
    endfor
    ;; convert kg/m2/s to Gg
    co_target   = co[index_lon_s:index_lon_e,index_lat_s:index_lat_e]   * area_target * 60 * 60 * 24 / 1E6
    c2h6_target = c2h6[index_lon_s:index_lon_e,index_lat_s:index_lat_e] * area_target * 60 * 60 * 24 / 1E6
    c3h8_target = c3h8[index_lon_s:index_lon_e,index_lat_s:index_lat_e] * area_target * 60 * 60 * 24 / 1E6
    alk4_target = alk4[index_lon_s:index_lon_e,index_lat_s:index_lat_e] * area_target * 60 * 60 * 24 / 1E6
    c3h6_target = c3h6[index_lon_s:index_lon_e,index_lat_s:index_lat_e] * area_target * 60 * 60 * 24 / 1E6
    ch2o_target = ch2o[index_lon_s:index_lon_e,index_lat_s:index_lat_e] * area_target * 60 * 60 * 24 / 1E6
    ald2_target = ald2[index_lon_s:index_lon_e,index_lat_s:index_lat_e] * area_target * 60 * 60 * 24 / 1E6
    acet_target = acet[index_lon_s:index_lon_e,index_lat_s:index_lat_e] * area_target * 60 * 60 * 24 / 1E6
    
    Emissions = [total(co_target), total(c2h6_target), total(c3h8_target), total(alk4_target), total(c3h6_target), $
                    total(ch2o_target), total(ald2_target), total(acet_target), total(co_target)*0.73*(1e-3)/28.01*72.11, $
                    total(co_target)*1.764*(1e-3)/28.01*78, total(co_target)*1.230*(1e-3)/28.01*92.15, total(co_target)*0.222*(1E-3)/28.01*106.18, $
                    total(co_target)*46/28*9.5*1E-3, total(co_target)*11.48*0.75*(1e-3)/28.01*60.06, total(co_target)*1.01*(1e-3)*58.09/28]
                
    if unit eq 'mass' then begin
        Emis_total_QFED = Emissions

        write_csv, 'BB_raw_emissions_mass.csv', compounds, Emis_total_GFED4, Emis_total_GFAS, Emis_total_QFED 
    endif 
    if unit eq 'carbon' then begin
        Emis_total_QFED = Emissions*ratio_mass2Cmass

        write_csv, 'BB_raw_emissions_carbon.csv', compounds, Emis_total_GFED4, Emis_total_GFAS, Emis_total_QFED
    endif 
end
