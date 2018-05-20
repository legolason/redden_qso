function alambda,x,p
c1 = p(0)
c2 = p(1)
c3 = p(2)
c4 = p(3)
c5 = p(4)
gamma = p(5)
R_V = p(6)
x0    = p(7)
curve = x*0.

 xcutuv = 10000.0/2700.0
 xspluv = 10000.0/[2700.0,2600.0]
 iuv = where(x ge xcutuv, N_UV, complement = iopir, Ncomp = Nopir)
 IF (N_UV GT 0) THEN xuv = [xspluv,x[iuv]] ELSE  xuv = xspluv

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
          curve[iopir] = CSPLINE([xsplopir,xspluv],[ysplopir,yspluv],x[iopir])

return,curve  
 end

;================
pro getextfm

lam = findgen(500)/50+0.01
result = [-0.71,0.94,0.0,0.5,5.90,1.0,2.74,4.6]
;resultn = [-5.23,2.08,0.0,0.81,5.9,1.00,2.73,4.6]
;resultp = [-4.89,2.34,0.0,1.01,5.9,1.00,2.35,4.6]

curve = alambda(lam,result[0:7])/result(6)
;curven = alambda(lam,resultn[0:7])/resultn(6)
;curvep = alambda(lam,resultp[0:7])/resultp(6)

plot, lam,curve,yr=[0,8],xr=[0,8]
readcol,'SMC.pei.txt',x,y
oplot,x,y,color=250
;oplot,lam,curvep,color=250
;oplot,lam,curven,color=250

;wrtab,Transpose([[lam],[curve],[curven],[curvep]]),'qso.txt'

end
