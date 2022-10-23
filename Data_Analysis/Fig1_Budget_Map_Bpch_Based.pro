;; This script is used for plotting spatial distributions of total VOCs from different sources, used as figure 1 in WE-CAN emission paper
;; lxu, 0426/2021, updated the 1) ground sites and 2) flight tracks on anthropogenic emission map and biomass burning map.
;; To reproduce the figure 1 in the emission paper, you need to rename the file and type command below
;; The input is from bpch files and it's obsolete after v12.5.0. Thus, we encourage to use netcdf hemco output in the future.
;; comp_src_vocs_budget_nested_singlefile_v4,inventories='gfas',region='west',/test,/save
pro comp_src_vocs_budget_nested_singlefile_v4,inventories=inventories,region=region,test=test,save=save

myct,/WhGrYlRd

;; plotting set
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

;; file directories
    fi = '/glade/u/home/lixujin/GEOS-CHEM/rundirs/tools/area.nc'
    ncdf_read,tmp,fi=fi,/all
    area=tmp.DXYP__DXYP
    lat = tmp.LAT
    lon = tmp.LON

; Reform longitude
    index = where(lon GT 180.,ct)
    lon[index] = lon[index]-360.

    dir_tmp = '/glade/scratch/lixujin/rundirs/GC12.5.0/emissions_only/'
    case strlowcase(inventories) of 
        'gfed4'    : subdir='geosfp_025x03125_tropchem_na_gfed4'
        'finn'     : subdir='geosfp_025x03125_tropchem_na_gfas'
        'qfed'     : subdir='geosfp_025x03125_tropchem_na_qfed'
        else       : subdir='geosfp_025x03125_tropchem_na_gfas' ;gfas
    endcase


    if n_elements(region) eq 0 then region=''
    case strlowcase(region) of 
        'world'    : limit = [0,-160,90,160]
        'southeast': limit=[25,-100,40,-75]
        'se'       : limit=[25,-100,40,-75]
        'west'     : limit=[36,-127,49.5,-102]
        ;lixu, 10/10/2019
        'na'       : limit=[10,-140,70,-70]    
        else       : limit=[25,-127,52,-65]
    endcase

    filename = './ps/'+'buget_'+inventories+ region + '_nested.ps'

    if keyword_set(test) then filename = './ps/test.ps'

    if keyword_set(save) then open_device, /ps, /color, bits=8, $
                               filename = filename
                               
;; read flight location informaiton.
    ; Set up arrays for plotting tracks

    lon_flight = [0]
    lat_flight = [0]

; Get c130 flight tracks
    ;if Keyword_Set(c130) then begin
    NewFiles = [MFindFile('/glade/u/home/lixujin/PROJECT/data/we-can/R4_merges/wecan-mrg60-c130_merge_*.ict')]
    ;NewFiles = [MFindFile('/glade/u/home/lixujin/PROJECT/data/we-can/R2_merges/wecan-mrg01-c130_merge_*.ict')]
    Date = strarr(n_elements(NewFiles))

    for i = 0L, n_elements( NewFiles )-1L do begin
    ; Get flight date
        loc = strpos(NewFiles[i],'2018')
        date[i]=string((byte(NewFiles[i]))[loc:loc+7])

        ;; WE-CAN data were collected by C-130, but data format as DC8's
        Get_c130, NewFiles[i],dc8_Long = lon_tmp, dc8_Lat =lat_tmp
        if lat_tmp[0] eq 0 then print, 'issue in this data:',date[i]
        ;; issue in 20180826
        if date[i] eq '20180826' then continue
        lat_flight = [lat_flight, lat_tmp]
        lon_flight = [lon_flight, lon_tmp]     
    endfor
    lat_flight = lat_flight[1:*]
    lon_flight = lon_flight[1:*]
    
    index=where(lon_flight GT 180., ct)
    if index[0] NE -1 then lon_flight[index] = lon_flight[index] - 360.
    

    
;; read ground sites location information
    lon_grd = [-122.3086, -116.3479, -105.0052, -121.2685, -119.7732, -119.8077, -121.8405, -113.9851, -121.687]
    lat_grd = [47.56824,   43.6007,   39.77949,  37.95074,  36.78538,  39.52508,  39.76168,  46.8585,   43.979]
    
;==============
;define basics
;==============
    avo=6.022e+23

    dir = dir_tmp + subdir

    ;; biomass burning emission:14
    variables = ['BIOBSRCE_C2H6','BIOBSRCE_C3H8','BIOBSRCE_ALK4','BIOBSRCE_PRPE',$
            'BIOBSRCE_CH2O','BIOBSRCE_ALD2', $
            'BIOBSRCE_BENZ','BIOBSRCE_TOLU','BIOBSRCE_XYLE',$
            'BIOBSRCE_ACET', 'BIOBSRCE_CO', $ ; add MEK, RCHO, HCOOH, and ACTA later
    ;; anth:14, see the excel
            'ANTHSRCE_C2H6','ANTHSRCE_C3H8','ANTHSRCE_ALK4','ANTHSRCE_PRPE', $
            'ANTHSRCE_CH2O','ANTHSRCE_ALD2','ANTHSRCE_RCHO','ANTHSRCE_ACET','ANTHSRCE_MEK', $
            'ANTHSRCE_BENZ','ANTHSRCE_TOLU','ANTHSRCE_XYLE', $
            'ANTHSRCE_EOH','ANTHSRCE_MACR', $
    ;; biogenic emission:23
            'BIOGSRCE_ISOP','BIOGSRCE_ACET','BIOGSRCE_PRPE','BIOGSRCE_MONX',$
            'BIOGSRCE_MBOX','BIOGSRCE_APIN','BIOGSRCE_BPIN','BIOGSRCE_LIMO',$
            'BIOGSRCE_SABI','BIOGSRCE_MYRC','BIOGSRCE_CARE','BIOGSRCE_OCIM',$
            'BIOGSRCE_HCOOH','BIOGSRC_ACTA','BIOGSRCE_ALD2','BIOGSRCE_OMON',$
            'BIOGSRCE_MOHX','BIOGSRCE_EOH','BIOGSRCE_FARN','BIOGSRCE_BCAR',$
            'BIOGSRCE_OSQT','BIOGSRCE_CHBr3','BIOGSRCE_CH2Br2']
            
            ;; SSBr2, ALD2sen, EOHsene, ALD2osr are not included
            ;;added in MEGAN
            ;;'BIOGSRCE__C2H4','BIOGSRCE__MTPA','BIOGSRCE__MTPO','BIOGSRCE__SESQ']
   
    
    
    sources = ['biogenic emission', 'anthropogenic emission', 'biomass burning emission']                

    months = ['07','08','09']

    days1 = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14',$
            '15','16','17','18','19','20','21','22','23','24','25','26','27','28',$
            '29','30','31']
    days2= ['01','02','03','04','05','06','07','08','09','10','11','12','13','14',$
            '15','16','17','18','19','20','21','22','23','24','25','26','27','28',$
            '29','30']

    days_start= ['19','20','21','22','23','24','25','26','27','28',$
                '29','30','31']

    ;for mm = 0, n_elements(months)-1 do begin
    ;;;test
    ;    ;if mm ne 0 then continue
    ;    print,'month:', months[mm]
    ;    if mm eq 0 then days = days_start
    ;    if mm eq 1 then days = days1
    ;    if mm eq 2 then days = days2
    ;    for dd = 0, n_elements(days)-1 do begin
            fi = dir + '/trac_avg.geosfp_025x03125_tropchem_na.201806010000_diagns.nc'
            result = file_test(fi)
    ;        ;if fi doesn't exist, then ignore the file
            if result eq 0 then stop
            
            ncdf_read,tmp,fi=fi,variables=variables
    ;        if dd eq 0 and mm eq 0 then begin
            ;;lumped species: ALK4, MEK
            ;biomass burning emission:14;anth:14;biogenic emission:23
                ;; Read in all VOCs of biomass burning emission
                
;; VOC reactivity
        ;; assume T = 298K
;; BB: 13 = 3 + 1 + 9
                bb_c2h6 = total(tmp.BIOBSRCE_C2H6,3)
                ;help,bb_c2h6, okk
                bb_c3h8 = total(tmp.BIOBSRCE_C3H8,3)
                ;help,bb_c3h8, okk
                bb_alk4 = total(tmp.BIOBSRCE_ALK4,3)
                ;help,bb_alk4, okk
                bb_prpe = total(tmp.BIOBSRCE_PRPE,3)
                ;help,bb_prpe, okk
                bb_ch2o = total(tmp.BIOBSRCE_CH2O,3)
                ;help,bb_ch2o, okk
                bb_ald2 = total(tmp.BIOBSRCE_ALD2,3)
                ;help,bb_ald2, okk
                bb_benz = total(tmp.BIOBSRCE_BENZ,3)
                ;help,bb_benz, okk
                bb_tolu = total(tmp.BIOBSRCE_TOLU,3)
                ;help,bb_tolu, okk
                bb_xyle = total(tmp.BIOBSRCE_XYLE,3)
                ;help,bb_xyle, okk
                bb_acet = total(tmp.BIOBSRCE_ACET,3)
                ;help,bb_acet, okk
                bb_mek   =  total(tmp.BIOBSRCE_CO,3) * 0.73/1000*4       ; atoms C/cm2/s
                bb_hcooh =  total(tmp.BIOBSRCE_CO,3) * 9.5/1000        ; atoms C/cm2/s
                bb_acta  =  total(tmp.BIOBSRCE_CO,3) * 8.61/1000*2       ; atoms C/cm2/s
                bb_rcho  =  total(tmp.BIOBSRCE_CO,3) * 1.01/1000*3      ; atoms C/cm2/s


;; biogenic emission:12 + 1 (MONX includes 8 compounds;) 
        ;; Read in all VOCs of biogenic emission
                biog_isop = total(tmp.BIOGSRCE_ISOP,3)
                ;help,biog_isop, okk
                biog_mbox = total(tmp.BIOGSRCE_MBOX,3)
                ;help,biog_mbox, okk
                biog_ald2 = total(tmp.BIOGSRCE_ALD2,3)
                ;help,biog_ald2, okk    

;; these biogenic emissiosn have been added up in MONX
                ;biog_apin = tmp.BIOGSRCE__APIN
                ;help,biog_apin, okk
                ;biog_bpin = tmp.BIOGSRCE__BPIN
                ;help,biog_bpin, okk
                ;biog_limo = tmp.BIOGSRCE__LIMO
                ;help,biog_limo, okk
                ;biog_sabi = tmp.BIOGSRCE__SABI
                ;help,biog_sabi, okk
                ;biog_myrc = tmp.BIOGSRCE__MYRC
                ;help,biog_myrc, okk
                ;biog_care = tmp.BIOGSRCE__CARE
                ;help,biog_care, okk
                ;biog_OMON = tmp.BIOGSRCE__OMON
                ;help,biog_omon, okk

                
;; don't repeatedly add monoterpenes
;; total monoterpenes could be a-pinene, b-pinene, Limonene, sabinene, carene, erpinene, terpinolene, myrcene, ocimene, other monoterpenes
;; why variables name are different between output and extention codes? Where does it change the names of variables?
;; Could be a question for GCST...
                ; MONX:Total monoterpenes(APIN () + BPIN + LIMO + SABI + MYRC + CARE + OCIM + OMON)
                biog_monx = total(tmp.BIOGSRCE_MONX,3)
                ;help,biog_monx, okk
                
                ; biog_mohx = tmp.BIOGSRCE_MOH
                ; help,biog_mohx, not in the output. Does model have it? Yes, but GC doesn't have MOH in the version before V12.9.X
                
                biog_EOH = total(tmp.BIOGSRCE_EOH,3)
                ;help,biog_eoh, okk 

                ;biog_hcooh = tmp.BIOGSRCE__HCOOH
                ;help,biog_hcooh, not in the output. Does model have it? Yes, not sure if model links to EMIS_FAXX very well... 
                ;biog_acta = tmp.BIOGSRCE__ACTA
                ;help,biog_acta, not in the output. Does model have it? Yes, not sure if model links EMIS_AAXX very well... 
;; could be questions for either Lu or Sree!!!!!!!!
; scientific question: are hcooh and acetic acid important from biogenic sources? Could be a scientifc question fro venessa !!!!!!
                biog_acet = total(tmp.BIOGSRCE_ACET,3)
                ;help,biog_acet, okk
                biog_prpe = total(tmp.BIOGSRCE_PRPE,3)
                ;help,biog_prpe, okk
        
                ; biog_c2h4 = tmp.BIOGSRCE__C2H4
                ; help,biog_c2h4, not in the output. Does model have it? Yes, but GC doesn't have C2H4 in the chemistry schemes
; scientific question: are C2H4 important from biogenic sources? Could be a scientifc question fro venessa !!!!!!

                biog_farn = total(tmp.BIOGSRCE_FARN,3)
                ;help,biog_farn, okk
                biog_bcar = total(tmp.BIOGSRCE_BCAR,3)
                ;help,biog_bcar, okk
                biog_osqt = total(tmp.BIOGSRCE_OSQT,3)        
                ;help,biog_osqt, okk
    
                ;; GC output doesn't have OTHR. Weird lumped monoterpenes.
                biog_chbr3 = total(tmp.BIOGSRCE_CHBR3,3) 
                ;help,biog_chbr3 ,okk
                biog_ch2br2 = total(tmp.BIOGSRCE_CH2BR2,3)   
                ;help,biog_ch2br2 ,okk
            
        ;; Read in all VOCs of anthropogenic emission
                anth_c2h6 = total(tmp.ANTHSRCE_C2H6,4)
                ;help,anth_c2h6   
                anth_c3h8 = total(tmp.ANTHSRCE_C3H8,4)
                ;help,anth_c3h8, okk
                anth_alk4 = total(tmp.ANTHSRCE_ALK4,4)
                ;help,anth_alk4,okk
                anth_prpe = total(tmp.ANTHSRCE_PRPE,4)
                ;help,anth_prpe,okk
                anth_ch2o = total(tmp.ANTHSRCE_CH2O,4)
                ;help,anth_ch2o, okk
                anth_ald2 = total(tmp.ANTHSRCE_ALD2,4)
                ;help,anth_ald2
                anth_rcho = total(tmp.ANTHSRCE_RCHO,4)
                ;help,anth_rcho, okk
                anth_acet = total(tmp.ANTHSRCE_ACET,4)
                ;help,anth_acet, okk
                anth_mek = total(tmp.ANTHSRCE_MEK,4)
                ;help,anth_mek, okk
                anth_benz = total(tmp.ANTHSRCE_BENZ,4)
                ;help,anth_benz,okk
                anth_tolu = total(tmp.ANTHSRCE_TOLU,4)
                ;help,anth_tolu,okk
                anth_xyle = total(tmp.ANTHSRCE_XYLE,4)
                ;help,anth_xyle,okk
                anth_eoh  = total(tmp.ANTHSRCE_EOH,4)
                ;help,anth_eoh,okk
                anth_MACR = total(tmp.ANTHSRCE_MACR,4)
                ;help,anth_MACR, okk
    
    anth = anth_c2h6 + anth_c3h8 + anth_alk4 + anth_prpe + $
            anth_ch2o + anth_ald2 + anth_rcho + $
            anth_acet + anth_mek + $
            anth_benz + anth_tolu + anth_xyle + $
            anth_eoh + anth_macr

    biog = biog_isop + biog_mbox + biog_ald2 + biog_monx + $
                biog_EOH + biog_acet + biog_prpe + $
                biog_farn + biog_bcar + biog_osqt + biog_chbr3 + $
                biog_ch2br2 
                
    bb = bb_c2h6 + bb_c3h8 + bb_alk4 + bb_prpe + $
            bb_ch2o + bb_ald2 + $ 
            bb_benz  + bb_tolu + bb_xyle + $
            bb_acet + $
            bb_mek + bb_hcooh + bb_acta + bb_rcho
    
    ; remove the lev dimension
    anth = total(anth,1)
    
; --------------------------------------
; western region:[36,-127,49.5,-105]
; --------------------------------------
    index_lat = where(lat gt 36 and lat lt 49.5, ct1)
    index_lon = where(lon gt -127 and lon lt -102, ct2)
        if ct1 gt 0 and ct2 gt 2 then begin
            index_lat_s = index_lat[0]
            index_lat_e = n_elements(index_lat)-1+ index_lat[0]
            index_lon_s = index_lon[0]
            index_lon_e = n_elements(index_lon)-1+ index_lon[0]
        endif
    ;print,lat[index_lat_s:index_lat_e],lon[index_lon_s:index_lon_e]

    ;convert [atoms C/cm2/s] into [atmos C/cm2], daily data
    anth_tt  = (anth)*60*60*24
    biog_tt = (biog)*60*60*24
    bb_tt   = (bb)*60*60*24



    ;;convert m2 to cm2
    area = area *100.*100. ;cm2
    

    ;; transpose data
    anth_tt = transpose(anth_tt)
    biog_tt = transpose(biog_tt)
    bb_tt = transpose(bb_tt)
    
    ;convert it into mass：[Gg C]
    anth_mass = (anth_tt)*area/avo*12/1e9
    biog_mass = (biog_tt)*area/avo*12/1e9
    bb_mass   = (bb_tt)*area/avo*12/1e9
    
            
            
    anth_ttprint, 'this is for western region'
    print, 'anthropogenic emission for west US[Gg C]', total(anth_mass[index_lon_s:index_lon_e,index_lat_s:index_lat_e])
    print, 'biogenic emission for west US[Gg C]',total(biog_mass[index_lon_s:index_lon_e,index_lat_s:index_lat_e])
    print, 'biomass burning emisison for west US[Gg C]',total(bb_mass[index_lon_s:index_lon_e,index_lat_s:index_lat_e])
    
    print, 'this is for NA'
    print,'anthropogenic emission for NA[Gg C]', total(anth_mass)
    print,'biogenic emission for NA[Gg]', total(biog_mass)
    print,'biomass burning emission for NA[Gg C]', total(bb_mass)
 
    print, total(transpose(anth_c2h6[index_lon_s:index_lon_e,index_lat_s:index_lat_e])*60*60*24*area[index_lon_s:index_lon_e,index_lat_s:index_lat_e]/avo*12/1e9)
    print, total(transpose(anth_c3h8[index_lon_s:index_lon_e,index_lat_s:index_lat_e])*60*60*24*area[index_lon_s:index_lon_e,index_lat_s:index_lat_e]/avo*12/1e9)
    print, total(transpose(anth_alk4[index_lon_s:index_lon_e,index_lat_s:index_lat_e])*60*60*24*area[index_lon_s:index_lon_e,index_lat_s:index_lat_e]/avo*12/1e9)
    print, total(transpose(anth_prpe[index_lon_s:index_lon_e,index_lat_s:index_lat_e])*60*60*24*area[index_lon_s:index_lon_e,index_lat_s:index_lat_e]/avo*12/1e9)
    print, total(transpose(anth_ch2o[index_lon_s:index_lon_e,index_lat_s:index_lat_e])*60*60*24*area[index_lon_s:index_lon_e,index_lat_s:index_lat_e]/avo*12/1e9)
    print, total(transpose(anth_ald2[index_lon_s:index_lon_e,index_lat_s:index_lat_e])*60*60*24*area[index_lon_s:index_lon_e,index_lat_s:index_lat_e]/avo*12/1e9)
    print, total(transpose(anth_acet[index_lon_s:index_lon_e,index_lat_s:index_lat_e])*60*60*24*area[index_lon_s:index_lon_e,index_lat_s:index_lat_e]/avo*12/1e9)
    print, total(transpose(anth_rcho[index_lon_s:index_lon_e,index_lat_s:index_lat_e])*60*60*24*area[index_lon_s:index_lon_e,index_lat_s:index_lat_e]/avo*12/1e9)

    print, total(transpose(anth_mek[index_lon_s:index_lon_e,index_lat_s:index_lat_e])*60*60*24*area[index_lon_s:index_lon_e,index_lat_s:index_lat_e]/avo*12/1e9)
    print, total(transpose(anth_benz[index_lon_s:index_lon_e,index_lat_s:index_lat_e])*60*60*24*area[index_lon_s:index_lon_e,index_lat_s:index_lat_e]/avo*12/1e9)
    print, total(transpose(anth_tolu[index_lon_s:index_lon_e,index_lat_s:index_lat_e])*60*60*24*area[index_lon_s:index_lon_e,index_lat_s:index_lat_e]/avo*12/1e9)
    print, total(transpose(anth_xyle[index_lon_s:index_lon_e,index_lat_s:index_lat_e])*60*60*24*area[index_lon_s:index_lon_e,index_lat_s:index_lat_e]/avo*12/1e9)
    print, total(transpose(anth_eoh[index_lon_s:index_lon_e,index_lat_s:index_lat_e])*60*60*24*area[index_lon_s:index_lon_e,index_lat_s:index_lat_e]/avo*12/1e9)
    print, total(transpose(anth_macr[index_lon_s:index_lon_e,index_lat_s:index_lat_e])*60*60*24*area[index_lon_s:index_lon_e,index_lat_s:index_lat_e]/avo*12/1e9)
 
 
;; ====================
;; plotting part
;; ====================
    s_end = 2
    for s = 0, 2 do begin
        if s eq 0 or s eq 1 then clevels = [0,1,2,5,10,20,35,50,70,90]
        if s eq 2 then clevels = [0,5,10,20,50,100,180,250,350]
        case s of
            0:begin
                distribution = transpose(bb)
            end
            1:begin
                distribution = transpose(anth)
            end
            2:begin
                distribution = transpose(biog)
            end
        endcase
        tvmap, distribution/1e12, lon, lat, $
            ;/cbar, $
            c_levels=clevels,                      $
            /continents, /isotropic,/us,                             $;, /grid
            /sample, $;title = 'Ozone at ' +  pres_str + ' | ' + title +  ' (UTC)', $
             title = title, $;charsize=1.3, charthick = 2,                $
            /noadvance, $;, /nogxlabel,/nogylabel,                                      $
            TCSFAC = 1,    CSFAC = 1.2,                                            $
            ;CBUNIT = '[1e12 * atoms C/cm2/s]',$
            limit = limit
            
;; oplot should be followed by colorbar
        if s eq 0 then oplot,lon_flight,lat_flight,color=1,thick=4
        if s eq 1 then oplot, lon_grd, lat_grd, color=1, PSym = 7, thick = 12, SYMSIZE = 2.5
        

        colorbar,  charthick=4.,$
            c_levels = clevels
    endfor
        
    if keyword_set(save) then close_device
    ;DEVICE, /CLOSE
    print,'done!'
end
