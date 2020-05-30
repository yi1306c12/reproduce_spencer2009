; added wt_scl parameter to make same as IF5 init_net (except for constant values)
; 4/18/09: added "con_i_all" parameter
pro init_net, str, weights=weights, verbose=verbose, gaba=gaba, con_t=con_t, con_p=con_p, con_r=con_r, $
              con_if=con_if, con_ir=con_ir, con_i_all=con_i_all, i_del=i_del, noise_wt=noise_wt, W_scl=W_scl

common path, home_path, single_path, avg_path, ps_path, latadj_path, dv_path, pca_path, bhv_path, $
             image_path, ers_path, ica_path
common net_params, N_all, N_e, N_ir, N_if, dt, inv_dt, update, t_init, t_pre, t_stim, t_post, t_all
common vars, W, W_noise
common rand, seed

W = fltarr(N_all,N_all)
W_noise = fltarr(N_all)

if not(keyword_set(noise_wt)) then  noise_wt=0.016

if keyword_set(weights) then begin
  openr, fu, single_path + weights + '_weights.dat', /get_lun
  readu, fu, W, W_noise  &  free_lun, fu
endif else begin

  if not(keyword_set(gaba)) then  gaba=0. $
  else  print, 'init_net: gaba = ', gaba

  ; Set up W by determining units to which unit i makes connections
  ; E->E
  i=0 & e_e_con=fix(N_e*.10)
  while (i LE N_e-1) do begin
    ix = fix(randomu(seed,e_e_con) * N_e)
    if (n_elements(uniq(ix,sort(ix))) EQ e_e_con) then begin
      W[ix,i] = fltarr(e_e_con) + 1.
      i+=1
    endif
  endwhile
  print, 'e-e done'

  ; E->Ir
  i=0 & e_ir_con=N_ir*.40
  while (i LE N_e-1) do begin
    a=randomu(seed,N_ir)  &  ix=uniq(a,sort(a))
    if (n_elements(ix) EQ N_ir) then begin
      W[N_e+ix[0:e_ir_con-1],i] = fltarr(e_ir_con) + 1.
      i+=1
    endif
  endwhile
  print, 'e-ir done'

  ; E->If
  i=0 & e_if_con=N_if*.40
  while (i LE N_e-1) do begin
    ix = fix(randomu(seed,e_if_con) * N_if)
    if (n_elements(uniq(ix,sort(ix))) EQ e_if_con) then begin
      W[ix+N_e+N_ir,i] = fltarr(e_if_con) + 1.
      i+=1
    endif
  endwhile
  print, 'e-if done'

  ; Ir->E
  i=N_e & ir_e_con=N_e*.50
  while (i LE N_e+N_ir-1) do begin
    a=randomu(seed,N_e)  &  ix=uniq(a,sort(a))
    if (n_elements(ix) EQ N_e) then begin
      W[ix[0:ir_e_con-1],i] = fltarr(ir_e_con) - 1.
      i+=1
    endif
  endwhile
  print, 'ir-e done'

  ; Ir->Ir
  i=N_e & ir_ir_con=N_ir*.15
  while (i LE N_e+N_ir-1) do begin
    a=randomu(seed,N_ir)  &  ix=uniq(a,sort(a))
    if (n_elements(ix) EQ N_ir) then begin
      W[N_e+ix[0:ir_ir_con-1],i] = fltarr(ir_ir_con) - 1.
      i+=1
    endif
  endwhile
  print, 'ir-ir done'

  ; Ir->If
  i=N_e & ir_if_con=N_if*.50
  while (i LE N_e+N_ir-1) do begin
    a = fix(randomu(seed,N_if) * 20000)
    ix = uniq(a,sort(a))

    if (n_elements(ix) EQ N_if) then begin
      W[ix[0:ir_if_con-1]+N_e+N_ir,i] = fltarr(ir_if_con) - 1.
      i+=1
    endif
  endwhile
  print, 'ir-if done'

  ; If->E
  i=N_e+N_ir & if_e_con=N_e*.50
  while (i LE N_all-1) do begin
    a=randomu(seed,N_e)  &  ix=uniq(a,sort(a))
    if (n_elements(ix) EQ N_e) then begin
      W[ix[0:if_e_con-1],i] = fltarr(if_e_con) - 1.
      i+=1
    endif
  endwhile
  print, 'if-e done'

  ; If->Ir
  i=N_e+N_ir & if_ir_con=N_ir*.35
  while (i LE N_all-1) do begin
    a=randomu(seed,N_ir)  &  ix=uniq(a,sort(a))
    if (n_elements(ix) EQ N_ir) then begin
      W[N_e+ix[0:if_ir_con-1],i] = fltarr(if_ir_con) - 1.
      i+=1
    endif
  endwhile
  print, 'if-ir done'

  ; If->If
  i=N_e+N_ir & if_if_con=N_if*.60
  while (i LE N_all-1) do begin
    a=randomu(seed,N_if)  &  ix=uniq(a,sort(a))
    if (n_elements(ix) EQ N_if) then begin
      W[N_e+N_ir+ix[0:if_if_con-1],i] = fltarr(if_if_con) - 1.
      i+=1
    endif
  endwhile
  print, 'if-if done'

  W[0:N_e-1, 0:N_e-1] *= 0.01         ; E->E
  W[N_e:N_e+N_ir-1, 0:N_e-1] *= 0.008  ; E->Ir
  W[N_e+N_ir:*, 0:N_e-1] *= 0.019      ; E->If

  W[0:N_e-1, N_e:N_e+N_ir-1] *= 0.008         ; Ir->E
  W[N_e:N_e+N_ir-1, N_e:N_e+N_ir-1] *= 0.008  ; Ir->Ir
  W[N_e+N_ir:*, N_e:N_e+N_ir-1] *= 0.008      ; Ir->If

  W[0:N_e-1, N_e+N_ir:*] *= 0.019*(1. - gaba)        ; If->E
  W[N_e:N_e+N_ir-1, N_e+N_ir:*] *= 0.01*(1. - gaba)  ; If->Ir
  W[N_e+N_ir:*, N_e+N_ir:*] *= 0.01*(1. - gaba)      ; If->If

  if keyword_set(W_scl) then begin
    W *= W_scl
    print, 'global weight scaling = ', W_scl
  endif else $
    W *= 0.0825  ; global weight scaling

  W_noise = fltarr(N_all) + noise_wt   ; external noise input weights
endelse

print, 'Initial # weights: ', total(W NE 0.)
print, '        # excitatory weights: ', total(W GT 0.)
print, '        # inhibitory weights: ', total(W LT 0.)
print, '        # E->E weights: ', total(W[0:N_e-1,0:N_e-1] GT 0.)

if keyword_set(con_t) then begin
  print, '       con_t = ', con_t
  weight_ix = where(W[*] NE 0., n_weights)
  W2 = W[weight_ix]   ; nonzero connections
  del_n = con_t*n_weights

  ok=0
  while (not(OK)) do begin
    a=randomu(seed2,n_weights,/long)  &  ix=uniq(a,sort(a))
    if (n_elements(ix) EQ n_weights) then  ok=1
  endwhile
  delete_ix = ix[0:del_n-1]
  W2[delete_ix] = 0.
  W[weight_ix] = W2

  n_weights_rem = total(W NE 0.)
  print, '       # remaining weights = ', n_weights_rem,'  ', float(n_weights_rem)/n_weights
endif

if keyword_set(con_p) then begin  ; connections to pyramidal cells
  print, '       con_p = ', con_p
  W_p = W[0:N_e-1,*]
  weight_ix = where(W_p[*] NE 0., n_weights)
  W2 = W_p[weight_ix]
  del_n = con_p*n_weights

  ok=0
  while (not(OK)) do begin
    a=randomu(seed2,n_weights,/long)  &  ix=uniq(a,sort(a))
    if (n_elements(ix) EQ n_weights) then  ok=1
  endwhile
  delete_ix = ix[0:del_n-1]
  W2[delete_ix] = 0.
  W_p[weight_ix] = W2
  W[0:N_e-1,*] = W_p

  n_weights_rem = total(W NE 0.)
  print, '       # remaining weights = ', n_weights_rem,'  ', float(n_weights_rem)/n_weights
endif

if keyword_set(con_r) then begin  ; recurrent excitatory connections (btw pyramidal cells)
  print, '       con_r = ', con_r
  W_p = W[0:N_e-1,0:N_e-1]
  weight_ix = where(W_p[*] NE 0., n_weights)
  W2 = W_p[weight_ix]
  del_n = con_r*n_weights

  ok=0
  while (not(OK)) do begin
    a=randomu(seed2,n_weights,/long)  &  ix=uniq(a,sort(a))
    if (n_elements(ix) EQ n_weights) then  ok=1
  endwhile
  delete_ix = ix[0:del_n-1]
  W2[delete_ix] = 0.
  W_p[weight_ix] = W2
  W[0:N_e-1,0:N_e-1] = W_p

  n_weights_rem = total(W NE 0.)
  print, '       # remaining weights = ', n_weights_rem,'  ', float(n_weights_rem)/n_weights
endif

if keyword_set(con_if) then begin  ; If-cell connections
  print, '       con_if = ', con_if

  ; projections from If-cells to all other cells
  W_p = W[*,N_e+N_ir:*]
  weight_ix = where(W_p[*] NE 0., n_weights)
  W2 = W_p[weight_ix]
  del_n = con_if*n_weights

  ok=0
  while (not(OK)) do begin
    a=randomu(seed2,n_weights,/long)  &  ix=uniq(a,sort(a))
    if (n_elements(ix) EQ n_weights) then  ok=1
  endwhile
  delete_ix = ix[0:del_n-1]
  W2[delete_ix] = 0.
  W_p[weight_ix] = W2
  W[*,N_e+N_ir:*] = W_p

  ; projections to If-cells from E- and Ir-cells
  W_p = W[N_e+N_ir:*,0:N_e+N_ir-1]
  weight_ix = where(W_p[*] NE 0., n_weights)
  W2 = W_p[weight_ix]
  del_n = con_if*n_weights

  ok=0
  while (not(OK)) do begin
    a=randomu(seed2,n_weights,/long)  &  ix=uniq(a,sort(a))
    if (n_elements(ix) EQ n_weights) then  ok=1
  endwhile
  delete_ix = ix[0:del_n-1]
  W2[delete_ix] = 0.
  W_p[weight_ix] = W2
  W[N_e+N_ir:*,0:N_e+N_ir-1] = W_p

  n_weights_rem = total(W NE 0.)
  print, '       # remaining weights = ', n_weights_rem,'  ', float(n_weights_rem)/n_weights
endif

if keyword_set(con_ir) then begin  ; Ir-cell connections
  print, '       con_ir = ', con_ir

  ; projections from Ir-cells to all other cells
  W_p = W[*,N_e:N_e+N_ir-1]
  weight_ix = where(W_p[*] NE 0., n_weights)
  W2 = W_p[weight_ix]
  del_n = con_ir*n_weights

  ok=0
  while (not(OK)) do begin
    a=randomu(seed2,n_weights,/long)  &  ix=uniq(a,sort(a))
    if (n_elements(ix) EQ n_weights) then  ok=1
  endwhile
  delete_ix = ix[0:del_n-1]
  W2[delete_ix] = 0.
  W_p[weight_ix] = W2
  W[*,N_e:N_e+N_ir-1] = W_p

  ; projections from E-cells to Ir-cells
  W_p = W[N_e:N_e+N_ir-1,0:N_e-1]
  weight_ix = where(W_p[*] NE 0., n_weights)
  W2 = W_p[weight_ix]
  del_n = con_ir*n_weights

  ok=0
  while (not(OK)) do begin
    a=randomu(seed2,n_weights,/long)  &  ix=uniq(a,sort(a))
    if (n_elements(ix) EQ n_weights) then  ok=1
  endwhile
  delete_ix = ix[0:del_n-1]
  W2[delete_ix] = 0.
  W_p[weight_ix] = W2
  W[N_e:N_e+N_ir-1,0:N_e-1] = W_p

  ; projections from If-cells to Ir-cells
  W_p = W[N_e:N_e+N_ir-1,N_e+N_ir:*]
  weight_ix = where(W_p[*] NE 0., n_weights)
  W2 = W_p[weight_ix]
  del_n = con_ir*n_weights

  ok=0
  while (not(OK)) do begin
    a=randomu(seed2,n_weights,/long)  &  ix=uniq(a,sort(a))
    if (n_elements(ix) EQ n_weights) then  ok=1
  endwhile
  delete_ix = ix[0:del_n-1]
  W2[delete_ix] = 0.
  W_p[weight_ix] = W2
  W[N_e:N_e+N_ir-1,N_e+N_ir:*] = W_p

  n_weights_rem = total(W NE 0.)
  print, '       # remaining weights = ', n_weights_rem,'  ', float(n_weights_rem)/n_weights
endif

if keyword_set(con_i_all) then begin  ; Ir- and If-cell connections
  print, '       con_ir = ', con_ir

  ; projections from I-cells to all other cells
  W_p = W[*,N_e:*]
  weight_ix = where(W_p[*] NE 0., n_weights)
  W2 = W_p[weight_ix]
  del_n = con_i_all*n_weights

  ok=0
  while (not(OK)) do begin
    a=randomu(seed2,n_weights,/long)  &  ix=uniq(a,sort(a))
    if (n_elements(ix) EQ n_weights) then  ok=1
  endwhile
  delete_ix = ix[0:del_n-1]
  W2[delete_ix] = 0.
  W_p[weight_ix] = W2
  W[*,N_e:*] = W_p

  ; projections from E-cells to I-cells
  W_p = W[N_e:*,0:N_e-1]
  weight_ix = where(W_p[*] NE 0., n_weights)
  W2 = W_p[weight_ix]
  del_n = con_i_all*n_weights

  ok=0
  while (not(OK)) do begin
    a=randomu(seed2,n_weights,/long)  &  ix=uniq(a,sort(a))
    if (n_elements(ix) EQ n_weights) then  ok=1
  endwhile
  delete_ix = ix[0:del_n-1]
  W2[delete_ix] = 0.
  W_p[weight_ix] = W2
  W[N_e:*,0:N_e-1] = W_p

  n_weights_rem = total(W NE 0.)
  print, '       # remaining weights = ', n_weights_rem,'  ', float(n_weights_rem)/n_weights
endif

if keyword_set(i_del) then begin
  del_n = i_del*(N_ir+N_if)
  print, 'init_net: i_del units = ', del_n
  delete_ix = N_e + long(randomu(seed2,del_n) * (N_ir+N_if))
  W[*,delete_ix] = 0.
  W[delete_ix,*] = 0.
endif

print, 'Final # weights: ', total(W NE 0.)
print, '      # excitatory weights: ', total(W GT 0.)
print, '      # inhibitory weights: ', total(W LT 0.)
print, '      # E->E weights: ', total(W[0:N_e-1,0:N_e-1] GT 0.)

if keyword_set(verbose) then begin
  print, 'mean weights: ', total(W)/N_all, total(W[0:N_e-1,*])/N_e, total(W[N_e:N_e+N_ir-1,*])/N_ir, $
                           total(W[N_e+N_ir:*,*])/N_if

  W2 = W GT 0.
  print, 'mean E->E connections per E-cell: ', total(W2[0:N_e-1,0:N_e-1])/N_e
  print, 'mean E->Ir connections per E-cell: ', total(W2[N_e:N_e+N_ir-1,0:N_e-1])/N_e
  print, 'mean E->If connections per E-cell: ', total(W2[N_e+N_ir:*,0:N_e-1])/N_e

  W2 = W LT 0.
  print, 'mean Ir->E connections per Ir-cell: ', total(W2[0:N_e-1,N_e:N_e+N_ir-1])/N_ir
  print, 'mean Ir->Ir connections per Ir-cell: ', total(W2[N_e:N_e+N_ir-1,N_e:N_e+N_ir-1])/N_ir
  print, 'mean Ir->If connections per Ir-cell: ', total(W2[N_e+N_ir:*,N_e:N_e+N_ir-1])/N_ir

  print, 'mean If->E connections per If-cell: ', total(W2[0:N_e-1,N_e+N_ir:*])/N_if
  print, 'mean If->Ir connections per If-cell: ', total(W2[N_e:N_e+N_ir-1,N_e+N_ir:*])/N_if
  print, 'mean If->If connections per If-cell: ', total(W2[N_e+N_ir:*,N_e+N_ir:*])/N_if

  scale=1
  window, 0, xsize=N_all*scale, ysize=N_all*scale, title='init_net: ' + str
  tvscl, congrid(bytscl(W), N_all*scale, N_all*scale)
endif

openw, fu, single_path + str + '_param.txt', /get_lun
printf, fu, N_all, N_e, N_ir, N_if, dt, inv_dt, update, t_init, t_pre, t_stim, t_post, t_all
free_lun, fu

openw, fu, single_path + str + '_weights.dat', /get_lun
writeu, fu, W, W_noise
free_lun, fu

print, 'init_net: ', str, ' done'
return
end
