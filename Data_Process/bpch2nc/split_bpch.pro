;;lxu
;; IDL script may have some limits then we need to look into python and shell scripts
@bpch2nc
@bpch_sep
pro split_bpch,nested=nested,fbf=fbf,nc=nc
   ; Splits into daily bpch files (change as necesary)
   ; Output directory
   ; Set up
    dir = './'   
    if keyword_set(nested) then tmp = 'trac_avg.geosfp_025x03125_tropchem_na.2018'
    if keyword_set(fbf) then tmp = 'trac_avg.geosfp_4x5_benchmark.2018'


    mms= ['07','08','09']
    
    big_mm = ['01','02','03','04','05','06','07','08',$
                '09','10','11','12','13','14','15','16',$
                '17','18','19','20','21','22','23','24',$
                '25','26','27','28','29','30','31']
    small_mm = ['01','02','03','04','05','06','07','08',$
                '09','10','11','12','13','14','15','16',$
                '17','18','19','20','21','22','23','24',$
                '25','26','27','28','29','30']
                
    n_mms = INDGEN(3)+7
    n_big_mm = INDGEN(31)+1
    n_small_mm = INDGEN(30)+1

    ;big_mm = string(INDGEN(31)+1,format=i02)
    ;big_mm = STRTRIM(big_mm, 1)
    ;small_mm = string(INDGEN(30)+1)
    ;small_mm = STRTRIM(small_mm, 1)
    
    for mm = 0, n_elements(mms)-1 do begin
        if mms[mm] eq '07' or mms[mm] eq '08' then dds=big_mm
        if mms[mm] eq '07' or mms[mm] eq '08' then n_dds=n_big_mm
        if mms[mm] eq '09' then dds=small_mm
        if mms[mm] eq '09' then n_dds=n_small_mm
        for dd =0,n_elements(dds)-1 do begin
            infile = tmp + mms[mm] + dds[dd] +'0000'
            outfile = tmp + mms[mm] + dds[dd] + '1111'
            ;to see if the file exists
            tmp1 = file_test(dir + infile)
            tmp2 = file_test(dir + outfile)
            ;if infile doesn't exist then skip it
            if tmp1 eq 0 then continue
            ;if outfile exist then skip it
            ;if tmp2 eq 1 then continue
            ;tau0 time
            time=10000L*2018+n_mms[mm]*100L+n_dds[dd]*1L
            print,time
            bpch_sep,dir+infile,dir+outfile,diagn=['BIOBSRCE','ANTHSRCE','BIOGSRCE','IJ-CHK-$','IJ-AVG-$']
            ;bpch_sep,dir+infile,dir+'trac_avg.geosfp_025x03125_tropchem_na.%DATE%.%TIME%'+'1111',diagn=['BIOBSRCE','ANTHSRCE','BIOGSRCE','IJ-CHK-$','IJ-AVG-$']
            
            ;bpch_sep,dir+infile,dir+outfile+'1',diagn=['IJ-AVG-$']
            ctm_cleanup
            CLOSE,/ALL
        endfor
        if keyword_set(nc) then begin
            for dd =0,n_elements(dds)-1 do begin
                infile = tmp + mms[mm] + dds[dd] +'0000'
                outfile = tmp + mms[mm] + dds[dd] + '1111'
                ;to see if the file exists
                tmp1 = file_test(dir + infile)
                ;tmp2 = file_test(dir + 'trac_avg.geosfp_025x03125_tropchem_na.%DATE%.%TIME%.nc')
                ;if file doesn't exist then skip it
                if tmp1 eq 0 then continue
                ;bpch2nc, dir+outfile, dir+'trac_avg.geosfp_025x03125_tropchem_na.%DATE%.%TIME%.nc'
            
                bpch2nc, dir+outfile, dir + 'trac_avg.geosfp_4x5_tropchem.%DATE%.%TIME%.nc'
                ctm_cleanup
                CLOSE,/ALL
            endfor
        endif
    endfor
print,'done'
end
