pro werstf, fname, xrange=xrange, freq_range=freq_range, zrange=zrange, thresh=thresh, $
            neg=neg, baseline=baseline, save=save, print=print

common path, home_path, single_path, avg_path, ps_path, latadj_path, dv_path, $
             pca_path, bhv_path, image_path, ers_path, ica_path
common params, n_ids, n_points, n_chans, period, epoch_range, filt_width, base_range
common params, n_ids
common wers_params, n_wers_pts, wers_range, wers_base_range, n_scales, scales

;on_error, 2

if not(keyword_set(xrange)) then  xrange=[0,100]
if not(keyword_set(freq_range)) then  freq_range=[2,100]

data = read_wers(ers_path + fname)
data = reverse(data,2)

; log10 scaling of power measures
data[n_ids:*,*,[0,1,2,3,4,5]] = alog10(data[n_ids:*,*,[0,1,2,3,4,5]])

if keyword_set(baseline) then  data=wers_baseline(data,baseline)
if not(keyword_set(thresh)) then  neg=0
if keyword_set(thresh) then  data=wers_thresh(data,thresh,freq_range,neg=neg)

orient = 'landscape'
if keyword_set(print) then begin
  set_plot, 'PS'
  !P.font=0
  if (orient EQ 'landscape') then $
    device, /landscape $
  else $
    device, /portrait, /inches, xsize=8.0, xoffset=.18, ysize=9.0, yoffset=1.0
  device, /helvetica, font_size=font_size, filename=ps_path + fname + '.ps'
  device, /color
endif else begin
  xsize=1100  &  ysize=700
  window, xsize=xsize, ysize=ysize, retain=2, title='Wavelet ERS'
endelse

cols=4  &  rows=3
!P.multi[1:2] = [cols,rows]
!P.charsize = 2.
!X.margin = [3,2]
!X.ticklen = 0.05
!X.minor = 2
!Y.margin = [3,2]
!Y.ticklen = 0.04
n_plots = cols*rows

;xtickname = [' ',' ',' ',' ',' ',' ']
;ytickname = [' ',' ',' ',' ',' ']

legend_box=200  &  text_box1=201  &  text_box2=202  &  empty_box=203
layout = $
;[text_box1, 0,  1,  2, $
[text_box1, 3,  4,  5, $
 empty_box, 6,  7,  8, $
 empty_box, 9, 10, 11]

; Set plot titles
titles = $
;[   '', 'PWT: PCs', 'PWT: RSIs', 'PWT: FSIs', $
[    '', 'Power: PCs', 'Power: RSIs', 'Power: FSIs', $
    '', 'Synchrony: PCs', 'Synchrony: RSIs', 'Synchrony: FSIs', $
    '', 'Synchrony: PC-RSI', 'Synchrony: PC-FSI', 'Synchrony: RSI-FSI']

x_disp_pts = ms_to_pts(xrange,/wers)
n_x_disp_pts = x_disp_pts[1] - x_disp_pts[0] + 1
x_disp_ix = n_ids + indgen(n_x_disp_pts) + x_disp_pts[0]
x_axis = indgen(n_x_disp_pts)*period + xrange[0]

freqs = 1000./scales
freq_ix = where((freqs GE freq_range[0]) AND (freqs LE freq_range[1]), n_freq_ix)
y_axis = freqs[freq_ix]
;print, 'n_scales = ', n_freq_ix

n_levels = 253          ; for Rainbow 18 & Purple-Red+Stripes colormaps
;temp = data[x_disp_ix,freq_ix,*]
;if not(keyword_set(zrange)) then begin
;  maxv = max(temp[*,*,*], min=minv)
;endif else begin
;  minv = zrange[0]
;  maxv = zrange[1]
;endelse

;level_inc = (maxv-minv)/(n_levels-1)
;levels1 = (findgen(n_levels) * level_inc) + minv
;;levels = [minv-level_inc, levels1, maxv+level_inc, maxv+(2*level_inc)]  ; Rainbow 18
;levels = [minv-(3*level_inc), minv-(2*level_inc), minv-(1*level_inc), levels1]  ; Purple-Red+Stripes

for i=0,n_plots-1 do begin
  if (layout[i] EQ empty_box) then  !P.multi[0]=!P.multi[0]-1  $
  else  if (layout[i] EQ text_box1) then  text_plot, x_axis, fname  $
  else if (layout[i] EQ legend_box) then begin
    temp = congrid(levels[3:*], n_x_disp_pts, n_freq_ix)
    contour, temp[*,freq_ix], x_axis, reverse(y_axis), xstyle=9, ystyle=1, /fill, levels=levels, $
             xtickname=['',' ','',' ','',' '], xtitle='ms', ytitle='Hz', title=' '
    axis, xaxis=1, xticks=2, xcharsize=0.8, $
          xtickn=[strtrim(string(minv,'(f6.2)'),2), ' ', strtrim(string(maxv,'(f6.2)'),2)]
  endif else begin
    chan_ix = where(data[2,0,*] EQ layout[i])
    temp = reform(data[x_disp_ix,*,chan_ix])
    ;if (layout[i] EQ 1) then  temp = alog10(temp)
    ;contour, reverse(temp[*,freq_ix],2), x_axis, reverse(y_axis), levels=levels, xstyle=1, ystyle=1, $
    ;         xtickname=xtickname, ytickname=ytickname, title=titles[i], /fill
    
    maxv = max(temp[*,freq_ix], min=minv)
    if (chan_ix GT 5) then begin
      minv=0.  &  maxv=1.
    endif
    print, titles[i], ': ', minv, maxv
    level_inc = (maxv-minv)/(n_levels-1)
    levels1 = (findgen(n_levels) * level_inc) + minv
    levels = [minv-(3*level_inc), minv-(2*level_inc), minv-(1*level_inc), levels1]  ; Purple-Red+Stripes
    
    contour, reverse(temp[*,freq_ix],2), x_axis, reverse(y_axis), levels=levels, xstyle=1, ystyle=1, $
             xtickname=xtickname, ytickname=ytickname, title=titles[i], /fill
  endelse
endfor

; Legend
if keyword_set(print) then begin
  ;s = reform(levels[3:*], 1, n_levels)
  ;s = congrid(bytscl(s, min=minv, max=maxv), 2, 253)
  ;tv, transpose(s), 8., 6.3, /inches, xsize=1.2, ysize=.15
  ;xyouts, .932, .875, string(maxv,'(f6.2)'), /normal, charsize=1.
  ;xyouts, .837, .875, string(minv,'(f6.2)'), /normal, charsize=1.
endif else begin
  s_length = xsize/(!P.multi[1]+1)
  s_width = (ysize/!P.multi[2])*.05
  ;s = reform(levels[0:n_levels-2], 1, n_levels-1)
  s = reform(levels[3:*], 1, n_levels)
  s = congrid(bytscl(s, min=minv, max=maxv), s_width, s_length)

  x_offset=10  &  y_offset=30
  tv, transpose(s), x_offset, ysize - (s_width+y_offset+50)
  xyouts, x_offset, ysize - (s_width+y_offset+70), $
          strtrim(string(minv,'(f8.2)'),1), /device, align=0.
  xyouts, s_length+x_offset, ysize - (s_width+y_offset+70), string(maxv,'(f8.2)'), /device, align=1.
endelse

!P.multi = 0
!P.charsize = 1.
!X.charsize = 1.
!X.ticklen = 0
!X.minor = 0
!Y.charsize = 1.
!Y.ticklen = 0

if keyword_set(print) then begin
  device, /close_file
  set_plot, 'WIN'
endif

if keyword_set(save) then  saveimage, image_path + fname + '.bmp', /bmp

return
end
