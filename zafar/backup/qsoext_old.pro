function alambda,x,p
c1 = p(0)
c2 = p(1)
c3 = p(2)
c4 = p(3)
c5 = p(4)
gamma = p(5)
R_V = p(6)
x0    = p(7)
c=299792458.
z = 0.92
lambdax = c/(10^x)/1e-10
lambdax = lambdax/(1.+z)
lambdax = 10000./lambdax
 xcutuv = 10000.0/2700.0
 xspluv = 10000.0/[2700.0,2600.0]
index = where(lambdax le 100.)
lamneed = lambdax(index)
curve = lamneed*0.
 iuv = where(lamneed ge xcutuv, N_UV, complement = iopir, Ncomp = Nopir)
 IF (N_UV GT 0) THEN xuv = [xspluv,lamneed[iuv]] ELSE  xuv = xspluv

    yuv = c1  + c2*xuv
    yuv = yuv + c3*xuv^2/((xuv^2-x0^2)^2 +(xuv*gamma)^2)
    yuv = yuv + c4*(((xuv>c5)-c5)^2)
    yuv = yuv + R_V
    yspluv  = yuv[0:1]                  ; save spline points

 IF (N_UV GT 0) THEN curve[iuv] = yuv[2:*]      ; remove spline points
 
; Compute optical portion of A(lambda)/E(B-V) curve
; using cubic spline anchored in UV, optical, and IR

 xsplopir = [0,10000.0/[26500.0,12200.0,6000.0,5470.0,4670.0,4110.0]]
 ysplir   = [0.0,0.26469,0.82925]*R_V/3.1 
 ysplop   = [poly(R_V, [-4.22809e-01, 1.00270, 2.13572e-04] ), $
             poly(R_V, [-5.13540e-02, 1.00216, -7.35778e-05] ), $
             poly(R_V, [ 7.00127e-01, 1.00184, -3.32598e-05] ), $
             poly(R_V, [ 1.19456, 1.01707, -5.46959e-03, 7.97809e-04, $ 
                     -4.45636e-05] ) ]
  
 ysplopir = [ysplir,ysplop]

 if (Nopir GT 0) then $
          curve[iopir] = CSPLINE([xsplopir,xspluv],[ysplopir,yspluv],lamneed[iopir])

return,curve  
 end
;-----------------------------------------------
function tzfm,x,p
flux = p(8) + p(9)*x - 0.4*p(10)*alambda(x,p[0:7])/p(6)
return, flux
end
;---------------------------------------------------

pro qsoext

c=299792458.
z = 0.92

readcol,'phot.txt',ra,dec,wlph,tf,tferr

;for plotting purposes in F_lambda, conversion from log(nu_Fu) to F_lambda and lambda to log(nu)
tw = wlph*1e-10
tw = c/tw
tw = alog10(tw)
mage = tf + tferr
mag = 10^(tf - tw)
mage = 10^(mage - tw)
mag = 2.99792458E-05*mag/wlph^2
mage = 2.99792458E-05*mage/wlph^2 - mag

plot, tw,tf,psym=4,yr=[9,12]
oploterror,tw,tf,tferr,psym=4

;Including spectrum 
tab = mrdfits('spec-4195-55452-0817.fits',1,h)
wl = tab.loglam
wl = 10^wl
spec = tab.flux
spec = spec*1e-17
sig = tab.ivar
;Scalling spectrum to the i’ band
ii = where(wl GE 7070 and wl LE 8300)
norm = mean(spec(ii))
normi = mag[5]
spec = spec/norm*normi
;converting variance to 1 sigma
sig = sqrt(1/sig)
sig = sig*1e-17
funrederr = sig+spec

;converting spectrum and errors to the log(nu_Fnu) and log(nu)
funred = 3.33564095E+04*spec*wl^2
funrederr = 3.33564095E+04*funrederr*wl^2
wl = wl*1e-10
wl = c/wl
funrederr = funrederr*wl
funred = funred*wl
wl = alog10(wl)
funred = alog10(funred)
funrederr = alog10(funrederr) - funred

plot, wl,funred,yr=[9.3,11.5],/ysty,xr=[14,15.5],/xsty
oplot, tw,tf,psym=4,color=250

;Taking non-emission regions of the spectra
ind = where(wl LT 14.73 and wl GT 14.56)
wll = wl(ind)
funredl = funred(ind)
funrederrl = funrederr(ind)
ind = where(wl LT 14.85 and wl GT 14.77)
wlu = wl(ind)
funredu = funred(ind)
funrederru = funrederr(ind)

wl = [wlu,wll]
funred = [funredu,funredl]
funrederr = [funrederru,funrederrl]

;combining data to make full SED in a way that photometry is not included where spectrum is available
allx = [tw[2],[wl],tw[6],tw[9]]
ally = [tf[2],[funred],tf[6],tf[9]]
allerr = [tf[2],[funrederr],tferr[6],tferr[9]]
ally(where(ally)) = 10.34
allerr(where(allerr)) = 0.15

;QSO template. Here trying to see the slope of the template and scalling it to the observed K band data.
readcol, 'compoM.data',wl1,flux1
wl1 = wl1*(1.+z)
flux1 = flux1*2.2e-14
ii = where(wl1 GE 19000 and wl1 LE 40000)
normm = mean(flux1(ii))
normk = mag[9]
flux1 = flux1/normm*normk
index = where(wl1 ge 1200. and wl1 le 40000.)
wl1 = wl1(index)
flux1 = flux1(index)
flux1 = 3.33564095E+04* flux1*wl1^2

wl1 = wl1*1e-10
wl1 = c/wl1
flux1 = flux1*wl1
wl1 = alog10(wl1)
flux1 = alog10(flux1)

;fit the SED with a slope which is clear of emission lines
; first match the slope and normalisation to the template. 
normk = mean(flux1(ii))
normsed = tf[9]-(tw[9]*0.7531)

pi = replicate({ fixed :0,limited:[0,0], limits:[0.D,0.D]},11)
pi(0).fixed = 0 	        ;c1
pi(1).fixed = 0			;c2
pi(2).fixed = 1			;c3
pi(3).fixed = 0  		;c4
pi(4).fixed = 0			;c5
pi(5).fixed = 1			;gamma
pi(6).fixed = 0			;R_V
pi(7).fixed = 1			;x0
pi(8).fixed = 1			;norm
pi(9).fixed = 1			;slope
pi(10).fixed = 0		;A_V
start = [-4.959,2.264,0.0,0.46,5.9,1.0,2.74,4.6,normsed,0.75,0.1] 
;Gordon et al. 2003 SMC parameters
;start = [-4.959,2.264,0.389,0.461,5.9,1.0,2.74,4.6,0.6,-0.1,0.75,0.67] 
pi(0).limited = 1
pi(0).limits = [-5.0,2.0]
pi(1).limited = 1
pi(1).limits = [0.0,3.0]
pi(3).limited = 1
pi(3).limits = [0.15,1.45]
pi(4).limited = 1
pi(4).limits = [5.0,6.5]
pi(6).limited = 1
pi(6).limits = [1.5,6.0]
result = mpfitfun('tzfm',allx,ally,allerr,start,parinfo=pi,perror=perror,bestnorm=b)
newx = findgen(440)/100+12.5
inp = where(newx ge 1.5 and newx le 15.84)
newxo = newx(inp)


plot,wl1,flux1,xr=[13.5,16.0],yr=[9.3,11.5],/xsty,/ysty
flux1 = [flux1 - 0.4*result(10)*alambda(wl1,result[0:7])/result(6)]
;oplot, newxo, [2.20+0.57*newxo],color=250
oplot, newxo, [result[8]+result[9]*newxo],color=250,linest=2
oplot, wl,funred,color=250,linest=3
oplot, wl1, flux1,linest=3
fit = tzfm(wl1,result[0:10])
oplot,tw, tf, psym=4,thick=4,color=250
oplot,wl1,fit,linest=2,color=200

;------------------------------------------------
;simulation MC for power-law
mean=fltarr(n_elements(allx))
sigma=fltarr(n_elements(allx))
array=fltarr(n_elements(allx))
new_array=fltarr(n_elements(allx))
result_array=fltarr(11,1000)

FOR j=0,999 DO BEGIN
	FOR i=0,n_elements(allx)-1 DO BEGIN	
		mean=ally
		sigma=allerr
		array[i]=RANDOMN(i+j*10,1)
	ENDFOR
new_array=array*sigma+mean
result_array(*,j)=mpfitfun('tzfm',allx,new_array,allerr,start,parinfo=pi,perror=perror,dof=dof)
ENDFOR
p0 = result_array(0,*)
p0_mean=total(p0)/n_elements(p0)
;p0_resid = p0-p0_mean
;p0_var = total(p0_resid^2)/(n_elements(p0)-1.0)
p0_var = stddev(p0)

p1 = result_array(1,*)
p1_mean=total(p1)/n_elements(p1)
p1_var = stddev(p1)

p2 = result_array(2,*)
p2_mean=total(p2)/n_elements(p2)
p2_var = stddev(p2)

p3 = result_array(3,*)
p3_mean=total(p3)/n_elements(p3)
p3_var = stddev(p3)

p4 = result_array(4,*)
p4_mean=total(p4)/n_elements(p4)
p4_var = stddev(p4)

p5 = result_array(5,*)
p5_mean=total(p5)/n_elements(p5)
p5_var = stddev(p5)

p6 = result_array(6,*)
p6_mean=total(p6)/n_elements(p6)
p6_var = stddev(p6)
psort6 = p6[sort(p6)]
p6_min = p6_mean-psort6[159]
p6_plus = psort6[839]-p6_mean 
;plothist,p6,bin=0.1

p7 = result_array(7,*)
p7_mean=total(p7)/n_elements(p7)
p7_var = stddev(p7)

p8 = result_array(8,*)
p8_mean=total(p8)/n_elements(p8)
p8_var = stddev(p8)

p9 = result_array(9,*)
p9_mean=total(p9)/n_elements(p9)
p9_var = stddev(p9)
psort9 = p9[sort(p9)]
p9_min = p9_mean-psort9[159]
p9_plus = psort9[839]-p9_mean 

p10 = result_array(10,*)
p10_mean=total(p10)/n_elements(p10)
p10_var = stddev(p10)
psort10 = p10[sort(p10)]
p10_min = p10_mean-psort10[159]
p10_plus = psort10[839]-p10_mean 



;printing results
print,'c1',result[0],'  +/-',p0_var
print,'c2',result[1],'  +/-',p1_var
print,'c3',result[2],'  +/-',p2_var
print,'c4',result[3],'  +/-',p3_var
print,'c5',result[4],'  +/-',p4_var
print,'gamma',result[5],'  +/-',p5_var
print,'R_V',result[6],'  +/-',p6_plus,p6_min,p6_var
print,'x0',result[7],'  +/-',p7_var
;print,'norm',result[8],'  +/-',p8_var
;print,'slope',1-p9_mean,1-result[9],'  +/-',p9_plus,p9_min,p9_var
print,'A_V',result[10],'  +/-',p10_plus,p10_min,p10_var
print, 'chi2', b, '      dof',dof
nhpfmpl = 1-chisqr_pdf(b,dof)
print, 'NHP', nhpfmpl

end  

