; Izhikevich neuron model, modified from IF4
; 3/12/09: added GABA decay time constant keyword (gabad)
; 4/6/09: changed back from Izh version
; 4/7/09: changed "gabad" keyword to "decay_gaba" (conflict with "gaba")
; 4/18/09: added "con_i_all" keyword (passed to init_net.pro)
pro ctx, str, rseed=rseed, ampa_scl=ampa_scl, nmda_e=nmda_e, nmda_rs=nmda_rs, nmda_fs=nmda_fs, gaba=gaba, $
              con_t=con_t, con_p=con_p, con_r=con_r, con_if=con_if, con_ir=con_ir, con_i_all=con_i_all, $
              i_del=i_del, noise_wt=noise_wt, W_scl=W_scl, decay_gaba=decay_gaba

common path, home_path, single_path, avg_path, ps_path, latadj_path, dv_path, pca_path, bhv_path, $
             image_path, ers_path, ica_path
common params, n_ids
common net_params, N_all, N_e, N_ir, N_if, dt, inv_dt, update, t_init, t_pre, t_stim, t_post, t_all
common vars, W, W_noise
common rand, seed

; Get seed for random number generator
if keyword_set(rseed) then begin
  if (size(rseed,/type) EQ 7) then begin  seed=rseed
    seed = lonarr(36)
    openr, fu, single_path + rseed + '_rand.dat', /get_lun
    readu, fu, seed  &  free_lun, fu
  endif  else  seed=long(rseed)
endif else  dum=randomu(seed,1)

openw, fu, single_path + str + '_rand.dat', /get_lun
writeu, fu, seed  &  free_lun, fu

if not(keyword_set(ampa_scl)) then  ampa_scl=1. $
else  print, 'ctx: ampa=', ampa_scl

if not(keyword_set(nmda_e)) then  nmda_e=0. $
else  print, 'ctx: nmda_e=', nmda_e
if not(keyword_set(nmda_rs)) then  nmda_rs=0. $
else  print, 'ctx: nmda_rs=', nmda_rs
if not(keyword_set(nmda_fs)) then  nmda_fs=0. $
else  print, 'ctx: nmda_fs=', nmda_fs
nmda_scl = $
[fltarr(N_e)+(0.45*(1. - nmda_e)), fltarr(N_ir)+(0.1*(1. - nmda_rs)), fltarr(N_if)+(0.1*(1. - nmda_fs))]

if not(keyword_set(gaba)) then  gaba=0
if not(keyword_set(con_t)) then  con_t=0
if not(keyword_set(con_p)) then  con_p=0
if not(keyword_set(con_r)) then  con_r=0
if not(keyword_set(con_if)) then  con_if=0
if not(keyword_set(con_ir)) then  con_ir=0
if not(keyword_set(con_i_all)) then  con_i_all=0
if not(keyword_set(i_del)) then  i_del=0
if not(keyword_set(noise_wt)) then  noise_wt=0
if not(keyword_set(W_scl)) then  W_scl=0
if not(keyword_set(decay_gaba)) then  decay_gaba=5.

spikes = fltarr(N_all,inv_dt)
current_spike = inv_dt-1
v = dblarr(N_all)
temp = intarr(N_all)
noise = fltarr(N_all,update)

g_r_ampa=fltarr(N_all)  &  g_r_nmda=fltarr(N_all)  &  g_r_gaba=fltarr(N_all)
g_d_ampa=fltarr(N_all)  &  g_d_nmda=fltarr(N_all)  &  g_d_gaba=fltarr(N_all)
g_ampa=0.               &  g_nmda=0.               &  g_gaba=0.
inv_tau_r_ampa=-dt/.5   &  inv_tau_r_nmda=-dt/2.   &  inv_tau_r_gaba=-dt/.5
;inv_tau_d_ampa=-dt/2.   &  inv_tau_d_nmda=-dt/100. &  inv_tau_d_gaba=-dt/5.
inv_tau_d_ampa=-dt/2.   &  inv_tau_d_nmda=-dt/100. &  inv_tau_d_gaba=-dt/decay_gaba
norm_ampa=.472689       &  norm_nmda=.904821       &  norm_gaba=.697016

inv_tau_mem=-1./[fltarr(N_e)+20.,fltarr(N_ir)+20.,fltarr(N_if)+10.]  ; used in non-Izh version
;thresh = [fltarr(N_e)-52.,fltarr(N_ir)-55.,fltarr(N_if)-52.]   ; non-Izh
thresh = [fltarr(N_e)-52.,fltarr(N_ir)-52.,fltarr(N_if)-52.]   ; non-Izh
;thresh = [fltarr(N_e)+30.,fltarr(N_ir)+30.,fltarr(N_if)+30.]   ; Izh

init_net, str, gaba=gaba, con_t=con_t, con_p=con_p, con_r=con_r, con_if=con_if, con_ir=con_ir, con_i_all=con_i_all, $
               i_del=i_del, noise_wt=noise_wt, W_scl=W_scl
noise_gen, str

openw, spike_fu, single_path + str + '_spk.dat', /get_lun
openw, Vm_fu, single_path + str + '_v.dat', /get_lun
openw, Isyn_fu, single_path + str + '_I_syn.dat', /get_lun
openw, g_ampa_fu, single_path + str + '_g_ampa.dat', /get_lun
openw, g_r_ampa_fu, single_path + str + '_g_r_ampa.dat', /get_lun
openw, g_d_ampa_fu, single_path + str + '_g_d_ampa.dat', /get_lun
openw, g_nmda_fu, single_path + str + '_g_nmda.dat', /get_lun
openw, g_r_nmda_fu, single_path + str + '_g_r_nmda.dat', /get_lun
openw, g_d_nmda_fu, single_path + str + '_g_d_nmda.dat', /get_lun
openw, g_gaba_fu, single_path + str + '_g_gaba.dat', /get_lun
openw, g_r_gaba_fu, single_path + str + '_g_r_gaba.dat', /get_lun
openw, g_d_gaba_fu, single_path + str + '_g_d_gaba.dat', /get_lun
openr, noise_fu, single_path + str + '_noise.dat', /get_lun

;spikes[*,*]=0  &  v[*]=double(-65.)  &  g_ampa[*]=0.  &  g_nmda[*]=0.  &  g_gaba[*]=0.  ; Izh
spikes[*,*]=0  &  v[*]=double(-70.)  &  g_ampa[*]=0.  &  g_nmda[*]=0.  &  g_gaba[*]=0.  ; non-Izh

;; Izh parameters
;a = [0.02 + fltarr(N_e), 0.02 + fltarr(N_ir), 0.1 + fltarr(N_if)]
;b = [0.2 + fltarr(N_e), 0.25 + fltarr(N_ir), 0.2 + fltarr(N_if)]
;c = [-65. + fltarr(N_e), -65. + fltarr(N_ir), -65. + fltarr(N_if)]
;d = [8. + fltarr(N_e), 2. + fltarr(N_ir), 2. + fltarr(N_if)]
;u = b * v

for t=long(0),t_all-1 do begin
  fired = where(v GE thresh, count)
  if (count GT 0) then begin
    v[fired] = double(-59.)   ; non-Izh
    ;v[fired] = c[fired]  ; Izh
    ;u[fired] += d[fired]  ; Izh
    spikes[fired,current_spike] = 1.
  endif

  if ((t mod update) EQ 0) then begin
    readu, noise_fu, temp
    noise[*,0] = temp
  endif

  glu_act = (W > 0.)#spikes[*,0] + W_noise*noise[*,t mod update]

  g_r_ampa += inv_tau_r_ampa*g_r_ampa + glu_act
  g_d_ampa += inv_tau_d_ampa*g_d_ampa + glu_act
  g_ampa = (g_d_ampa - g_r_ampa)/norm_ampa

  g_r_nmda += inv_tau_r_nmda*g_r_nmda + glu_act
  g_d_nmda += inv_tau_d_nmda*g_d_nmda + glu_act
  g_nmda = (g_d_nmda - g_r_nmda)/norm_nmda

  g_r_gaba += inv_tau_r_gaba*g_r_gaba - (W < 0.)#spikes[*,0]
  g_d_gaba += inv_tau_d_gaba*g_d_gaba - (W < 0.)#spikes[*,0]
  g_gaba = (g_d_gaba - g_r_gaba)/norm_gaba

  I_syn = -(ampa_scl*g_ampa*v + nmda_scl*g_nmda*v/(1. + exp(-0.062*v)/3.57) + g_gaba*(v+70.))

  if ((t mod update) EQ 0) then begin
    writeu, spike_fu, fix(total(spikes[*,*],2) GT 0)
    writeu, Vm_fu, v
    writeu, Isyn_fu, I_syn
    writeu, g_ampa_fu, g_ampa
    writeu, g_r_ampa_fu, g_r_ampa
    writeu, g_d_ampa_fu, g_d_ampa
    writeu, g_nmda_fu, g_nmda
    writeu, g_r_nmda_fu, g_r_nmda
    writeu, g_d_nmda_fu, g_d_nmda
    writeu, g_gaba_fu, g_gaba
    writeu, g_r_gaba_fu, g_r_gaba
    writeu, g_d_gaba_fu, g_d_gaba
  endif

  v += dt*(inv_tau_mem*(v + 70.) + I_syn)   ; non-Izh
  ;v += dt*(0.04 * v^2 + 5.*v + 140. - u + I_syn)
  ;u += dt*(a * ((b * v) - u))

  spikes=shift(spikes,0,-1)  &  spikes[*,current_spike]=0.
end ; t

free_lun, spike_fu, Vm_fu, Isyn_fu, noise_fu, g_ampa_fu, g_r_ampa_fu, g_d_ampa_fu, g_nmda_fu, g_r_nmda_fu, g_d_nmda_fu, g_gaba_fu, g_r_gaba_fu, g_d_gaba_fu

if (total(finite(v)) NE N_all) then  print, 'error: v=NaN'

print, 'ctx: ', str, ' done'
return
end
