pro ringspec,teff,dm0,qgrav,int0
;
common model,tl,ml,ql,wls,int,pol,iconv
;
nt=n_elements(tl)
it=findgen(nt)
ii=interpol(it,tl,teff) 
i=fix(ii)+1
i=i>1<(nt-1)
i1=i-1
a=(teff-tl(i1))/(tl(i)-tl(i1))
a1=1.-a
;
nm=n_elements(ml)
im=findgen(nm)
ii=interpol(im,ml,dm0) 
j=fix(ii)+1
j=j>1<(nm-1)
j1=j-1
b=(dm0-ml(j1))/(ml(j)-ml(j1))
b1=1.-b
;
nq=n_elements(ql)
iq=findgen(nq)
ii=interpol(iq,ql,qgrav) 
k=fix(ii)+1
k=k>1<(nq-1)
k1=k-1
c=(qgrav-ql(k1))/(ql(k)-ql(k1))
c1=1.-c
;
print,i,j,k,a,b,c,format='(3i4,3f10.3)'
in0=b*(a*int(i,j,k,*,*)+a1*int(i1,j,k,*,*)) + $
     b1*(a*int(i,j1,k,*,*)+a1*int(i1,j1,k,*,*)) 
in1=b*(a*int(i,j,k1,*,*)+a1*int(i1,j,k1,*,*)) + $
     b1*(a*int(i,j1,k1,*,*)+a1*int(i1,j1,k1,*,*)) 
int0=c*in0+c1*in1
;
return
end


