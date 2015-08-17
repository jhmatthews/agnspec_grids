pro make_agn_grid


spin = [0,0.998]
mass = 1e+9
mdot = 5.0
alpha = 0.1
; direction cosines in 10 degree steps from 0 to 90
mus = [1.00000000e+00,   9.84807753e-01,   9.39692622e-01,8.66025406e-01,   7.66044447e-01,   6.42787615e-01,5.00000008e-01,   3.42020153e-01,   1.73648189e-01,1.33491245e-08]
nmu = n_elements(mus)
nspin = 2
for i=0,nmu-1 do begin
  for j=0,nspin-1 do begin
    ; run agnspec
    agnspec, mdot=mdot, angm=spin(j), mu=mus(i), imu=i, ispin=j
    print, i, j
  endfor
endfor





; end of program
end

