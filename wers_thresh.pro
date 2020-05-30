function wers_thresh, data, thresh, freq_range, neg=neg

common path, home_path, single_path, avg_path, ps_path, latadj_path, dv_path, $
             pca_path, bhv_path, image_path, ers_path, ica_path
common params, n_ids, n_points, n_chans, period, epoch_range, filt_width, base_range
common wers_params, n_wers_pts, wers_range, wers_base_range, n_scales, scales

data2 = data
data2[n_ids:*,*,*] = 0.

chan_ix = [indgen(61),64]  ; monopolar channels, excluding Nas (not used) (note A1/A2 are mastoids)
n_chan_ix = n_elements(chan_ix)

;scales = dindgen(n_scales) * 0.125
;scales = 2d0^(scales)*(2*period)
scale_range = 1000./freq_range
scale_ix = where((scales GE scale_range[1]) AND (scales LE scale_range[0]), n_scale_ix)

pts_ix = n_ids + indgen(n_wers_pts) + ms_to_pts(wers_range[0],/wers)

for i=0,n_chan_ix-1 do $
  for j=0,n_scale_ix-1 do begin
    cycle = round(scales[scale_ix[j]]/period)
    ;cycle = round(0.5 * scales[scale_ix[j]]/period)
    ;cycle=1

    if (n_elements(thresh) EQ 1) then $
      if (keyword_set(neg)) then $
        for k=0,n_wers_pts-cycle do begin
          sel = data[pts_ix[k]:pts_ix[k]+cycle-1,scale_ix[j],chan_ix[i]] LE -(thresh)
          if (total(sel) EQ cycle) then $
            data2[pts_ix[k]:pts_ix[k]+cycle-1,scale_ix[j],chan_ix[i]] = data[pts_ix[k]:pts_ix[k]+cycle-1,scale_ix[j],chan_ix[i]]
        endfor $
      else $
        for k=0,n_wers_pts-cycle do begin
          sel = data[pts_ix[k]:pts_ix[k]+cycle-1,scale_ix[j],chan_ix[i]] GE thresh
          if (total(sel) EQ cycle) then $
            data2[pts_ix[k]:pts_ix[k]+cycle-1,scale_ix[j],chan_ix[i]] = data[pts_ix[k]:pts_ix[k]+cycle-1,scale_ix[j],chan_ix[i]]
        endfor $
    else $
      for k=0,n_wers_pts-cycle do begin
        sel_lo = data[pts_ix[k]:pts_ix[k]+cycle-1,scale_ix[j],chan_ix[i]] LE thresh[0]
        sel_hi = data[pts_ix[k]:pts_ix[k]+cycle-1,scale_ix[j],chan_ix[i]] GE thresh[1]
        if ((total(sel_lo) EQ cycle) OR (total(sel_hi) EQ cycle)) then $
          data2[pts_ix[k]:pts_ix[k]+cycle-1,scale_ix[j],chan_ix[i]] = data[pts_ix[k]:pts_ix[k]+cycle-1,scale_ix[j],chan_ix[i]]
      endfor
  endfor ;j

print, 'wers_thresh: done'
return, data2
end
