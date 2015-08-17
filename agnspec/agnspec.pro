pro agnspec,mass=mass,mdot=mdot,angm=angm,alpha=alpha,mu=mu,    $
            rout=rout,tmin=tmin,deltar=deltar,nre=nre,          $
            rcut=rcut,tlim=tlim,dirprg=dirprg,                  $
            restore=restore,file=file,nocomp=nocomp,            $
            frmax=frmax,frmin=frmin,nfobs=nfobs,                $
            plmod=plmod,plspec=plspec,opl=opl,_extra=e, imu=imu, ispin=ispin, savename=savename                                
;
common model,tl,ml,ql,wls,int,pol,iconv
;
; initialization; default values
;
if n_elements(file) eq 0 then file='agnspec.save'
if n_elements(restore) gt 0 then restore,file
if n_elements(dirprg) eq 0 then dirprg='./'
;
if n_elements(mass) eq 0 then mass=1.e9
if n_elements(mdot) eq 0 then mdot=1.
if n_elements(angm) eq 0 then angm=0.998
if n_elements(alpha) eq 0 then alpha=0.01
if n_elements(mu) eq 0 then mu=0.6
if n_elements(rout) eq 0 then rout=0.
if n_elements(tmin) eq 0 then tmin=9000.
if n_elements(deltar) eq 0 then deltar=0.1
if n_elements(nre) eq 0 then nre=0
if n_elements(tlim) eq 0 then tlim=1000.
if n_elements(frmax) eq 0 then frmax=1.e17
if n_elements(plspec) eq 0 then plspec=1
if n_elements(imu) eq 0 then imu = 0
if n_elements(ispin) eq 0 then ispin = 0
;print, save
;
nfgrid=300
frgrid0=1.e17
frgrid1=1.e13
nretype=4
gr=1

readgrid
;
; -----------------------------------------------
; 1. parameters and spectra of individual annuli:
; -----------------------------------------------
;
; -----------------------------------------------
; 1a. set up input data for dispar3.f
; -----------------------------------------------
;
get_lun,lun1
openw,lun1,'tmp.5'
printf,lun1,mass,mdot,angm,alpha
printf,lun1,rout,tmin,deltar,nre
close,lun1
free_lun,lun1
;
; -----------------------------------------------
; 1b. run dispar3 and read results
; -----------------------------------------------
;
;a=dirprg+'dispar3 <tmp.5 >! tmp.6'
a=dirprg+'dispar3 <tmp.5 >tmp.6'
spawn,a
;
get_lun,lun1
openr,lun1,'tmp.6'
ii=0
r=fltarr(200) & tr=r & mr=r & qr=r
while not eof(lun1) do begin
   readf,lun1,r0,t0,m0,q0,teff,dm,qgrav
   r(ii)=r0 & tr(ii)=t0 & mr(ii)=m0 & qr(ii)=q0
   print,ii,r(ii),tr(ii),mr(ii),qr(ii),teff,dm,qgrav,$
         format='(i3,4f10.2,f10.1,2e10.2)'
   ii=ii+1
endwhile
close,lun1
free_lun,lun1
r=r(0:ii-1)
tr=tr(0:ii-1)
mr=mr(0:ii-1)
qr=qr(0:ii-1)
nr=n_elements(r)
if n_elements(rcut) eq 0 then rcut=r(nr-1)*(tlim/10.^tr(nr-1))^(-1.333)
;
; -----------------------------------------------
; 1c spectra for individual annuli
; -----------------------------------------------
;
tmpin=fltarr(20)
totr=fltarr(n_elements(wls))
get_lun,lun1
openw,lun1,'emrad.in'
for i=0,nr-1 do begin
   lnn=i-(i/5)*5
   ringspec,tr(i),mr(i),qr(i),int0
   printf,lun1,r(i),n_elements(wls)
   for j=0,n_elements(wls)-1 do begin
     for k=0,9 do tmpin(2*k)=exp(int0(0,0,0,j,k))
     for k=1,9 do tmpin(2*k+1)=0.
     printf,lun1,exp(wls(j)),tmpin(0)
     printf,lun1,tmpin
   endfor
   if n_elements(plmod) gt 0 then begin
      if plmod gt 1 then cc=r(i)^2*(10.^(deltar)-1.) else cc=1.
      totr=totr+exp(int0(0,0,0,*,0))*cc
      if plmod le 2 then begin
         if i eq 0 and n_elements(opl) eq 0 then $
         plot,3.e18/exp(wls),exp(int0(0,0,0,*,0))*cc,$
              /xlog,/ylog,_extra=e else $
         oplot,3.e18/exp(wls),exp(int0(0,0,0,*,0))*cc,line=lnn,_extra=e
      endif
      if plmod eq 3 then begin
         if i eq 0 and n_elements(opl) eq 0 then $
         plot,3.e18/exp(wls),totr,/xlog,/ylog,_extra=e else $
         oplot,3.e18/exp(wls),totr,line=lnn,_extra=e
      endif
   endif
endfor
if n_elements(plmod) gt 0 then oplot,3.e18/exp(wls),totr,thick=2,$
                                     _extra=e
close,lun1
free_lun,lun1
;
; -----------------------------------------------
; 2. integration over annuli to get the total spectrum
; -----------------------------------------------
;
; -----------------------------------------------
; 2a. set up input data for kerrtrans9: std input
; -----------------------------------------------
;
if n_elements(nocomp) gt 0 then return
;
get_lun,lun1
openw,lun1,'tmp.in'
printf,lun1,-fix(nfgrid)
printf,lun1,frgrid1,frgrid0
;
if frmax gt 0 then  begin
  if n_elements(frmin) eq 0 then frmin=1.e14  
  if n_elements(nfobs) eq 0 then nfobs=300
  printf,lun1,-fix(nfobs)
  printf,lun1,frmin,frmax
endif else begin
  get_lun,lun2
;  a='wc -l freq.in >! tmpf'
  a='wc -l freq.in >tmpf'
  spawn,a
  get_lun,lun2
  openr,lun2,'tmpf'
  readf,lun2,nko
  close,lun2
  free_lun,lun2
  printf,lun1,fix(nko)
endelse

printf,lun1,mu
printf,lun1,mass,mdot,angm
printf,lun1,rcut,tmin,nr,nretype,gr
close,lun1
free_lun,lun1
;
; -----------------------------------------------
; 2b. run kerrtrans9
; -----------------------------------------------
;
;a=dirprg+'kerrtrans9  <tmp.in >! sp.out'
;outname = 'sp' + STRTRIM(ispin,1)+'_mu' + STRTRIM(imu,1) + '.out'
outname=savename
print, outname
a=dirprg+'kerrtrans9  <tmp.in >'+outname
print, a
spawn,a
; comment out plot
;if n_elements(plspec) gt 0 then plt,'sp.out',nc=5,$
;                           /xlog,opl=opl,_extra=e
;
return
end


