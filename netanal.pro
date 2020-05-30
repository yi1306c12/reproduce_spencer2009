pro netanal, str, stim=stim, print=print, results=results

common path, home_path, single_path, avg_path, ps_path, latadj_path, dv_path, pca_path, bhv_path, $
             img_path, ers_path, ica_path
common params, n_ids
common net_params, N_all, N_e, N_ir, N_if, dt, inv_dt, update, t_init, t_pre, t_stim, t_post, t_all
common vars, W, W_noise

noise = intarr(N_all,t_all/inv_dt)
if keyword_set(stim) then  stim = dblarr(N_all,t_all/inv_dt)
spikes = intarr(N_all,t_all/inv_dt)
Vm = dblarr(N_all,t_all/inv_dt)

openr, fu, single_path + str + '_noise.dat', /get_lun
readu, fu, noise  &  free_lun, fu

if keyword_set(stim) then begin
  openr, fu, single_path + str + '_stim.dat', /get_lun
  readu, fu, stim  &  free_lun, fu
endif

openr, fu, single_path + str + '_spk.dat', /get_lun
readu, fu, spikes  &  free_lun, fu

openr, fu, single_path + str + '_Vm.dat', /get_lun
readu, fu, Vm  &  free_lun, fu

; Power Spectra
stim_ix = lindgen(t_stim/inv_dt) + (t_init + t_pre)/inv_dt
;hwin = hanning(t_stim/inv_dt)
;hwin = hanning(t_stim/inv_dt, alpha=.54)
hwin = dblarr(t_stim/inv_dt) + 1.

pow_spk = fltarr(N_all,t_stim/inv_dt)
pow_Vm = fltarr(N_all,t_stim/inv_dt)
fft_Vm = complexarr(N_all,t_stim/inv_dt)
pow_noise = fltarr(N_all,t_stim/inv_dt)
if keyword_set(stim) then  pow_stim = fltarr(N_all,t_stim/inv_dt)

for n=0,N_all-1 do begin
  pow_spk[n,*] = abs(fft(spikes[n,stim_ix]*hwin, -1))^2
  pow_Vm[n,*] = abs(fft(Vm[n,stim_ix]*hwin, -1))^2
  fft_Vm[n,*] = fft(Vm[n,stim_ix]*hwin, -1)
  pow_noise[n,*] = abs(fft(noise[n,stim_ix]*hwin, -1))^2
endfor

e_pwt_spk = total(pow_spk[0:N_e-1,*],1)/N_e
i_pwt_spk = total(pow_spk[N_e:N_e+N_ir-1,*],1)/N_ir
a_pwt_spk = total(pow_spk[N_e+N_ir:*,*],1)/N_if
n_pwt_spk = total(pow_noise[*,*],1)/N_all

e_pwt_Vm = total(pow_Vm[0:N_e-1,*],1)/N_e
i_pwt_Vm = total(pow_Vm[N_e:N_e+N_ir-1,*],1)/N_ir
a_pwt_Vm = total(pow_Vm[N_e+N_ir:*,*],1)/N_if

e_avg_spk = total(spikes[0:N_e-1,stim_ix],1)/N_e
i_avg_spk = total(spikes[N_e:N_e+N_ir-1,stim_ix],1)/N_ir
a_avg_spk = total(spikes[N_e+N_ir:*,stim_ix],1)/N_if
n_avg_spk = total(noise[*,stim_ix],1)/N_all

e_avg_Vm = total(Vm[0:N_e-1,stim_ix],1)/N_e
i_avg_Vm = total(Vm[N_e:N_e+N_ir-1,stim_ix],1)/N_ir
a_avg_Vm = total(Vm[N_e+N_ir:*,stim_ix],1)/N_if
if keyword_set(stim) then  stim_avg = total(stim[*,stim_ix],1)/N_all

e_pwe_spk = abs(fft(e_avg_spk*hwin, -1))^2
i_pwe_spk = abs(fft(i_avg_spk*hwin, -1))^2
a_pwe_spk = abs(fft(a_avg_spk*hwin, -1))^2
n_pwe_spk = abs(fft(n_avg_spk*hwin, -1))^2

e_pwe_Vm = abs(fft(e_avg_Vm*hwin, -1))^2
i_pwe_Vm = abs(fft(i_avg_Vm*hwin, -1))^2
a_pwe_Vm = abs(fft(a_avg_Vm*hwin, -1))^2
if keyword_set(stim) then  stim_pwe = abs(fft(stim_avg*hwin, -1))^2

e_phl_Vm = abs(total((fft_Vm[0:N_e-1,*]/abs(fft_Vm[0:N_e-1,*])),1)/N_e)
i_phl_Vm = abs(total((fft_Vm[N_e:N_e+N_ir-1,*]/abs(fft_Vm[N_e:N_e+N_ir-1,*])),1)/N_ir)
a_phl_Vm = abs(total((fft_Vm[N_e+N_ir:*,*]/abs(fft_Vm[N_e+N_ir:*,*])),1)/N_if)

if keyword_set(print) then begin
  set_plot, 'PS'
  !P.font=0
  device, /landscape
  ;device, /portrait, /inches, xsize=8.0, xoffset=.18, ysize=9.0, yoffset=1.0
  device, filename=img_path + str + '_netanal.ps'
  !P.charsize = 1.3
endif else begin
  window, 0, xsize=1100*1.1, ysize=850*1.1, title='netanal: ' + str
  !P.charsize = 1.8
endelse

!P.multi[1:2] = [4,7]
;!P.thick = 3.

!X.thick = 2.
!X.style=1
!X.margin = [9,2]
!X.ticklen = .06
!Y.margin = [3,3]
!Y.thick = 2.
!Y.style=1
!Y.charsize = 1.
!Y.minor = 2
ylog=0

x_axis = findgen(t_stim/inv_dt) + (t_init+t_pre)/inv_dt
stim_range = [t_init, t_init+t_stim]/inv_dt
freq_axis = (1000./dt)*findgen(t_stim)/(t_stim)
freq_range = [freq_axis[1],200.]
y_Vm = [0,0]

plot, x_axis, e_avg_spk, xrange=stim_range, title='PC Spiking: Average', xtitle='ms'
plot, x_axis, i_avg_spk, xrange=stim_range, title='RSI Spiking: Average', xtitle='ms'
plot, x_axis, a_avg_spk, xrange=stim_range, title='FSI Spiking: Average', xtitle='ms'
plot, x_axis, n_avg_spk, xrange=stim_range, yrange=[0.,0.5], title='Noise Spiking: Average', xtitle='!4ms!x'

plot, freq_axis, e_pwt_spk, xrange=freq_range, title='PC Spiking: Total Power', xtitle='!4Hz!x', ylog=ylog
plot, freq_axis, i_pwt_spk, xrange=freq_range, title='RSI Spiking: Total Power', xtitle='!4Hz!x', ylog=ylog
plot, freq_axis, a_pwt_spk, xrange=freq_range, title='FSI Spiking: Total Power', xtitle='!4Hz!x', ylog=ylog
plot, freq_axis, n_pwt_spk, xrange=freq_range, yrange=[0.,0.0006], title='Noise Spiking: Total Power', xtitle='!4Hz!x', ylog=ylog

plot, freq_axis, e_pwe_spk, xrange=freq_range, title='PC Spiking: Evoked Power', xtitle='!4Hz!x', ylog=ylog
plot, freq_axis, i_pwe_spk, xrange=freq_range, title='RSI Spiking: Evoked Power', xtitle='!4Hz!x', ylog=ylog
plot, freq_axis, a_pwe_spk, xrange=freq_range, title='FSI Spiking: Evoked Power', xtitle='!4Hz!x', ylog=ylog
plot, freq_axis, n_pwe_spk, xrange=freq_range, title='Noise Spiking: Evoked Power', xtitle='!4Hz!x', ylog=ylog

Vm_yrange = [-65,-50]
plot, x_axis, e_avg_Vm, xrange=stim_range, yrange=Vm_yrange, yminor=5, title='PC Vm: Average', xtitle='!4ms!x'
plot, x_axis, i_avg_Vm, xrange=stim_range, yrange=Vm_yrange, yminor=5, title='RSI Vm: Average', xtitle='!4ms!x'
plot, x_axis, a_avg_Vm, xrange=stim_range, yrange=Vm_yrange, yminor=5, title='FSI Vm: Average', xtitle='!4ms!x'
if keyword_set(stim) then $
  ;plot, x_axis, stim_avg, xrange=stim_range, yrange=[-1,1], title='Stimulus: Average', xtitle='!4ms!x' $
  plot, x_axis, stim_avg, xrange=stim_range, title='Stimulus: Average', xtitle='!4ms!x' $
else $
  plot, x_axis, /nodata, xstyle=5, ystyle=5

y_Vm = [0,0]
plot, freq_axis, e_pwt_Vm, xrange=freq_range, yrange=y_Vm, title='PC Vm: Total Power', xtitle='!4Hz!x', ylog=ylog
plot, freq_axis, i_pwt_Vm, xrange=freq_range, yrange=y_Vm, title='RSI Vm: Total Power', xtitle='!4Hz!x', ylog=ylog
plot, freq_axis, a_pwt_Vm, xrange=freq_range, yrange=y_Vm, title='FSI Vm: Total Power', xtitle='!4Hz!x', ylog=ylog
plot, x_axis, /nodata, xstyle=5, ystyle=5, yrange=[0,1]
xyouts, 0, .5, str

y_Vm = [0,0]
plot, freq_axis, e_pwe_Vm, xrange=freq_range, yrange=y_Vm, title='PC Vm: Evoked Power', xtitle='!4Hz!x', ylog=ylog
plot, freq_axis, i_pwe_Vm, xrange=freq_range, yrange=y_Vm, title='RSI Vm: Evoked Power', xtitle='!4Hz!x', ylog=ylog
plot, freq_axis, a_pwe_Vm, xrange=freq_range, yrange=y_Vm, title='FSI Vm: Evoked Power', xtitle='!4Hz!x', ylog=ylog
if keyword_set(stim) then $
  plot, freq_axis, stim_pwe, xrange=freq_range, title='Stimulus: Evoked Power', xtitle='!4Hz!x', ylog=ylog $
else $
  plot, x_axis, /nodata, xstyle=5, ystyle=5

plot, freq_axis, e_phl_Vm, xrange=freq_range, yrange=[0.,1.], title='PC Vm: PLF', xtitle='!4Hz!x'
plot, freq_axis, i_phl_Vm, xrange=freq_range, yrange=[0.,1.], title='RSI Vm: PLF', xtitle='!4Hz!x'
plot, freq_axis, a_phl_Vm, xrange=freq_range, yrange=[0,1.], title='FSI Vm: PLF', xtitle='!4Hz!x'
plot, x_axis, /nodata, xstyle=5, ystyle=5

if keyword_set(print) then begin
  device, /close_file
  set_plot, 'WIN'
endif

!P.multi=0
!P.charsize=0
!P.thick=0
!X.range=0
!X.margin=[10,4]
!X.ticklen = 0
!X.thick=0
!Y.margin=[4,2]
!Y.charsize=0
!Y.thick=0
!Y.minor = 0

return
end
