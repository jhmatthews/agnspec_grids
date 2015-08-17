pro readgrid,save=save,name=name,lte=lte,$
             wlo=wlo,whi=whi,nws=nws 
;
common model,tl,ml,ql,wls,int,pol,iconv
;
if n_elements(name) eq 0 then name='agngrid.save'
;
if n_elements(save) gt 0 then begin
;
if keyword_set(lte) then llab='lt' else llab='nc'
if n_elements(wlo) eq 0 then wlo=30.
if n_elements(whi) eq 0 then whi=3.e5
if n_elements(nws) eq 0 then nws=300
wls=alog(wlo)+alog(whi/wlo)/(nws-1)*findgen(nws)
na=10
x0=fltarr(2*na+2)
na1=na-1
na2=na+2
nw0=500
wl=fltarr(nw0)
int0=fltarr(nw0,na)
pol0=int0
;
tl=4.+findgen(15)*0.1
ml=3.+findgen(5)
ql=-12.+findgen(9)
;
tlab=['40','41','42','43','44','45','46','47','48','49',$
      '50','51','52','53','54']
mlab=['30','40','50','60','70']
qlab=['120','110','100','90','80','70','60','50','40']
;
int=fltarr(n_elements(tl),n_elements(ml),n_elements(ql),nws,na)
pol=int
iconv=fltarr(n_elements(tl),n_elements(ml),n_elements(ql))
;
llabb=llab
if llabb eq 'nc' then llabb='[0-9]nc'
spawn,'ls *'+llabb+'.9 add/*'+llabb+'.9 >! tmpd'
get_lun,lun1
openr,lun1,'tmpd'
iex=0
ii=0
a=''
modl=strarr(1000)
while not eof(lun1) do begin
  readf,lun1,a
  modl(ii)=strtrim(a,2)
  ii=ii+1
endwhile
modl=modl(0:ii-1)
;
for i=0,n_elements(tl)-1 do begin
  td=tl(i)*10.
  if td eq fix(td)/2*2. then dire='' else dire='add/
  if td gt 51. then dire='add/'
  for j=0,n_elements(ml)-1 do begin
    for k=0,n_elements(ql)-1 do begin
      a=dire+'t'+tlab(i)+'m'+mlab(j)+'q-'+qlab(k)+llab
;      print,a
      iex=max(a+'.9' eq modl)
      print,a,iex
      if iex eq 0 then begin
         iconv(i,j,k)=0
       endif else begin
         reltot1,file=a,niter,chmax
         if chmax lt 0 and niter gt 1 then iconv(i,j,k)=1 else $
                                           iconv(i,j,k)=0
      endelse
      print,a,iconv(i,j,k)
      if iconv(i,j,k) eq 1 then begin
      get_lun,lun1
      openr,lun1,a+'.10'
      n=0L
      w0=0.
      wl=fltarr(nw0)
      int0=fltarr(nw0,na)
      pol0=int0
      while not eof(lun1) do begin
        readf,lun1,x0
        if x0(0) lt w0 then goto,endread 
        wl(n)=x0(0)
        for m=0,na1 do begin
          int0(n,m)=x0(2*m+2)
          pol0(n,m)=x0(2*m+3)
        endfor
        n=n+1
        w0=x0(0)
      endwhile
      endread: nw=n
      close,lun1
      free_lun,lun1
      wl=alog(wl(0:nw-1))
      int0=alog(int0(0:nw-1,0:na1)>1.e-36)
      pol0=pol0(0:nw-1,0:na1)
      for m=0,na1 do begin
        in0=int0(*,m)
        ins=interpol(in0,wl,wls)
        pl0=pol0(*,m)
        pls=interpol(pl0,wl,wls)
        for n=0,nws-1 do begin
          int(i,j,k,n,m)=ins(n)
          pol(i,j,k,n,m)=pls(n)
        endfor
      endfor
      endif
    endfor
  endfor
endfor
;
;help
save,tl,ml,ql,tlab,mlab,qlab,wls,int,pol,iconv,file=name
;
endif else restore,file=name
;
return
end


