pro noise_gen, str

common path, home_path, single_path, avg_path, ps_path, latadj_path, dv_path, pca_path, bhv_path, $
             image_path, ers_path, ica_path
common net_params, N_all, N_e, N_ir, N_if, dt, inv_dt, update, t_init, t_pre, t_stim, t_post, t_all
common vars, W, W_noise
common rand, seed

freq = 100.  ; in Hz
noise = intarr(N_all,t_all/inv_dt)

for i=0,N_all-1 do begin
  a = randomu(seed, t_all/inv_dt, gamma=1.)/freq
  spike_times = total(long(1000.*a),/cumulative)
  ix = where(spike_times LT t_all/inv_dt, count)
  noise[i,spike_times[ix]] = 1
endfor

openw, fu, single_path + str + '_noise.dat', /get_lun
writeu, fu, noise  &  free_lun, fu

print, 'noise_gen: ', str, ' done'
return
end