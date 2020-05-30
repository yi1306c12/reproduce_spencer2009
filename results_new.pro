pro results_new, run, str1, str2=str2, nophc=nophc

common path, home_path, single_path, avg_path, ps_path, latadj_path, dv_path, pca_path, bhv_path, $
             img_path, ers_path, ica_path
common params, n_ids
common net_params, N_all, N_e, N_ir, N_if, dt, inv_dt, update, t_init, t_pre, t_stim, t_post, t_all
common vars, W, W_noise

if not(keyword_set(str2)) then  str2=''

conds = run + str1 + ['010','020','030','040','050','060','070','080','090','100'] + str2
conds = [run,conds]
n_conds = n_elements(conds)

Vm = dblarr(N_all,t_all/inv_dt)
fft_Vm = complexarr(N_all,t_stim/inv_dt)
pow_Vm = fltarr(N_all,t_stim/inv_dt)
spk = intarr(N_all,t_all/inv_dt)  ; spike data
fft_spk = complexarr(N_all,t_stim/inv_dt)
pow_spk = fltarr(N_all,t_stim/inv_dt)

pwt_e = fltarr(n_ids+n_conds, t_stim/inv_dt) & pwt_e[2,*]=0.
pwt_i = fltarr(n_ids+n_conds, t_stim/inv_dt) & pwt_i[2,*]=1.
pwt_a = fltarr(n_ids+n_conds, t_stim/inv_dt) & pwt_a[2,*]=2.
spk_e = fltarr(n_ids+n_conds, t_stim/inv_dt)
spk_i = fltarr(n_ids+n_conds, t_stim/inv_dt)
spk_a = fltarr(n_ids+n_conds, t_stim/inv_dt)
pwe_e = fltarr(n_ids+n_conds, t_stim/inv_dt) & pwe_e[2,*]=3.
pwe_i = fltarr(n_ids+n_conds, t_stim/inv_dt) & pwe_i[2,*]=4.
pwe_a = fltarr(n_ids+n_conds, t_stim/inv_dt) & pwe_a[2,*]=5.
phl_e = fltarr(n_ids+n_conds, t_stim/inv_dt) & phl_e[2,*]=6.
phl_i = fltarr(n_ids+n_conds, t_stim/inv_dt) & phl_i[2,*]=7.
phl_a = fltarr(n_ids+n_conds, t_stim/inv_dt) & phl_a[2,*]=8.
phc_e_ir = fltarr(n_ids+n_conds, t_stim/inv_dt) & phc_e_ir[2,*]=9.
phc_e_if = fltarr(n_ids+n_conds, t_stim/inv_dt) & phc_e_if[2,*]=10.
phc_ir_if = fltarr(n_ids+n_conds, t_stim/inv_dt) & phc_ir_if[2,*]=11.

stim_ix = lindgen(t_stim/inv_dt) + (t_init + t_pre)/inv_dt
;hwin = hanning(t_stim/inv_dt)
freq_axis = (1000./dt)*findgen(t_stim)/(t_stim)

openw, dat_fu, ers_path + run + str1 + str2 + '_Vm.wdat', /get_lun
openw, pwt_fu, dv_path + run + str1 + str2 + '_pwt.txt', /get_lun
printf, pwt_fu, conds[0], '...', conds[n_conds-1]
openw, pwt2_fu, dv_path + run + str1 + str2 + '_pwt2.txt', /get_lun
printf, pwt2_fu, conds[0], '...', conds[n_conds-1]
openw, pwe_fu, dv_path + run + str1 + str2 + '_pwe.txt', /get_lun
printf, pwe_fu, conds[0], '...', conds[n_conds-1]
openw, phl_fu, dv_path + run + str1 + str2 + '_phl.txt', /get_lun
printf, phl_fu, conds[0], '...', conds[n_conds-1]
openw, phc_fu, dv_path + run + str1 + str2 + '_phc.txt', /get_lun
printf, phc_fu, conds[0], '...', conds[n_conds-1]

for i=0,n_conds-1 do begin
  openr, fu, single_path + conds[i] + '_Vm.dat', /get_lun
  readu, fu, Vm  &  free_lun, fu
  openr, fu, single_path + conds[i] + '_spk.dat', /get_lun
  readu, fu, spk  &  free_lun, fu

  ; Power Spectra
  for n=0,N_all-1 do begin
    fft_Vm[n,*] = fft(Vm[n,stim_ix], -1)
    pow_Vm[n,*] = abs(fft_Vm[n,*])^2
    ;fft_spk[n,*] = fft(spk[n,stim_ix], -1)
    ;pow_spk[n,*] = abs(fft_spk[n,*])^2
  endfor

  pwt_e[n_ids+i,*] = total(pow_Vm[0:N_e-1,*],1)/N_e
  pwt_i[n_ids+i,*] = total(pow_Vm[N_e:N_e+N_ir-1,*],1)/N_ir
  pwt_a[n_ids+i,*] = total(pow_Vm[N_e+N_ir:*,*],1)/N_if

  spk_e[n_ids+i,*] = total(spk[0:N_e-1,stim_ix],1)/N_e
  spk_i[n_ids+i,*] = total(spk[N_e:N_e+N_ir-1,stim_ix],1)/N_ir
  spk_a[n_ids+i,*] = total(spk[N_e+N_ir:*,stim_ix],1)/N_if

  avg_Vm_e = total(Vm[0:N_e-1,stim_ix],1)/N_e
  avg_Vm_i = total(Vm[N_e:N_e+N_ir-1,stim_ix],1)/N_ir
  avg_Vm_a = total(Vm[N_e+N_ir:*,stim_ix],1)/N_if

  pwe_e[n_ids+i,*] = abs(fft(avg_Vm_e, -1))^2
  pwe_i[n_ids+i,*] = abs(fft(avg_Vm_i, -1))^2
  pwe_a[n_ids+i,*] = abs(fft(avg_Vm_a, -1))^2

  phl_e[n_ids+i,*] = abs(total((fft_Vm[0:N_e-1,*]/abs(fft_Vm[0:N_e-1,*])),1)/N_e)
  phl_i[n_ids+i,*] = abs(total((fft_Vm[N_e:N_e+N_ir-1,*]/abs(fft_Vm[N_e:N_e+N_ir-1,*])),1)/N_ir)
  phl_a[n_ids+i,*] = abs(total((fft_Vm[N_e+N_ir:*,*]/abs(fft_Vm[N_e+N_ir:*,*])),1)/N_if)

  if not keyword_set(nophc) then begin
    diff=fltarr(t_stim/inv_dt,N_e*N_ir)  &  count=long(0)
    for n=0,N_e-1 do $
      for m=N_e,N_e+N_ir-1 do begin
        diff[*,count] = atan(fft_Vm[n,*],/phase) - atan(fft_Vm[m,*],/phase)
        count++
      endfor
    phc_e_ir[n_ids+i,*] = abs(total(exp(complex(0,1)*diff[*,*]),2))/(N_e*N_ir)

    diff=fltarr(t_stim/inv_dt,N_e*N_if)  &  count=long(0)
    for n=0,N_e-1 do $
      for m=N_e+N_ir,N_all-1 do begin
        diff[*,count] = atan(fft_Vm[n,*],/phase) - atan(fft_Vm[m,*],/phase)
        count++
      endfor
    phc_e_if[n_ids+i,*] = abs(total(exp(complex(0,1)*diff[*,*]),2))/(N_e*N_if)

    diff=fltarr(t_stim/inv_dt,N_ir*N_if)  &  count=long(0)
    for n=N_e,N_e+N_ir-1 do $
      for m=N_e+N_ir,N_all-1 do begin
        diff[*,count] = atan(fft_Vm[n,*],/phase) - atan(fft_Vm[m,*],/phase)
        count++
      endfor
    phc_ir_if[n_ids+i,*] = abs(total(exp(complex(0,1)*diff[*,*]),2))/(N_ir*N_if)

    phc_e_ir_max = max(phc_e_ir[n_ids+i,1:49], phc_e_ir_max_i)
    phc_e_if_max = max(phc_e_if[n_ids+i,1:49], phc_e_if_max_i)
    phc_ir_if_max = max(phc_ir_if[n_ids+i,1:49], phc_ir_if_max_i)

    printf, phc_fu, phc_e_ir_max, phc_e_if_max, phc_ir_if_max, $
            round(freq_axis[phc_e_ir_max_i+1]), round(freq_axis[phc_e_if_max_i+1]), round(freq_axis[phc_ir_if_max_i+1])
  endif  ; no_phc

  pwt_e_max = max(pwt_e[n_ids+i,1:99], pwt_e_max_i)
  pwt_i_max = max(pwt_i[n_ids+i,1:99], pwt_i_max_i)
  pwt_a_max = max(pwt_a[n_ids+i,1:99], pwt_a_max_i)

  pwt_e_mean = mean(pwt_e[n_ids+i,1:99])
  pwt_i_mean = mean(pwt_i[n_ids+i,1:99])
  pwt_a_mean = mean(pwt_a[n_ids+i,1:99])

  spk_e_mean = mean(spk_e[n_ids+i,*])*1000.  ; converting to spikes/s
  spk_i_mean = mean(spk_i[n_ids+i,*])*1000.
  spk_a_mean = mean(spk_a[n_ids+i,*])*1000.

  pwe_e_max = max(pwe_e[n_ids+i,1:99], pwe_e_max_i)
  pwe_i_max = max(pwe_i[n_ids+i,1:99], pwe_i_max_i)
  pwe_a_max = max(pwe_a[n_ids+i,1:99], pwe_a_max_i)

  phl_e_max = max(phl_e[n_ids+i,1:99], phl_e_max_i)
  phl_i_max = max(phl_i[n_ids+i,1:99], phl_i_max_i)
  phl_a_max = max(phl_a[n_ids+i,1:99], phl_a_max_i)

  printf, pwt_fu, pwt_e_max, pwt_i_max, pwt_a_max, $
          round(freq_axis[pwt_e_max_i+1]), round(freq_axis[pwt_i_max_i+1]), round(freq_axis[pwt_a_max_i+1])
  printf, pwt2_fu, pwt_e_mean, pwt_i_mean, pwt_a_mean, spk_e_mean, spk_i_mean, spk_a_mean
  printf, pwe_fu, pwe_e_max, pwe_i_max, pwe_a_max, $
          round(freq_axis[pwe_e_max_i+1]), round(freq_axis[pwe_i_max_i+1]), round(freq_axis[pwe_a_max_i+1])
  printf, phl_fu, phl_e_max, phl_i_max, phl_a_max, $
          round(freq_axis[phl_e_max_i+1]), round(freq_axis[phl_i_max_i+1]), round(freq_axis[phl_a_max_i+1])
endfor

writeu, dat_fu, [[[pwt_e]], [[pwt_i]], [[pwt_a]], [[pwe_e]], [[pwe_i]], [[pwe_a]], $
                 [[phl_e]], [[phl_i]], [[phl_a]], [[phc_e_ir]], [[phc_e_if]], [[phc_ir_if]]]
free_lun, dat_fu
free_lun, pwt_fu, pwt2_fu, pwe_fu, phl_fu, phc_fu

print, 'results_new: done'
return
end
