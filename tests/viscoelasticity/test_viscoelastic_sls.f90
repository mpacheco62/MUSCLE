program test_sls
  use iso_fortran_env, only: real64
  use muscle_viscoelastic_sls
  use muscle_elasticity_base
  use muscle_tensors
  implicit none
  type(Viscoelastic_sls), target :: material
  class(Base_elasticity), pointer :: model
  type(Sls_internal) :: state,saved,independent
  type(ten_3D2Osym) :: eps,de,sig,again,pert,ds,predicted,eq
  type(ten_3D4O3sym) :: tangent
  real(real64) :: dt,exact,oldq(6),h
  integer :: j
  call material%set_parameters(1e7_real64,0.3_real64)
  model=>material
  dt=0.05_real64; h=1e-7_real64
  eps%vals=[0.01_real64,-0.003_real64,-0.003_real64,0.002_real64,0.001_real64,0.003_real64]
  de=eps; saved=state
  sig=model%stress(eps,de,dt,state)
  again=model%stress(eps,de,dt,state)
  if(any(state%q/=saved%q).or.any(sig%vals/=again%vals)) error stop 'impure trial'
  exact=1e5_real64*(1+sls_factor(dt,material%tau))
  if(abs(sig%vals(1)-exact)>1e-8_real64) error stop 'uniaxial analytic'
  tangent=model%dstress_dstrain(eps,de,dt,state)
  do j=1,6
    pert=eps; pert%vals(j)=pert%vals(j)+h
    ds=de; ds%vals(j)=ds%vals(j)+h
    again=model%stress(pert,ds,dt,state)
    ds%vals=0; ds%vals(j)=h
    predicted=tangent .ddot. ds
    if(maxval(abs(again%vals-sig%vals-predicted%vals))>1e-8_real64) error stop 'tangent'
  end do
  call model%stress_and_update(again,eps,de,dt,state)
  if(maxval(abs(again%vals-sig%vals))>1e-10_real64) error stop 'commit stress'
  if(any(independent%q/=0)) error stop 'shared Gauss history'
  oldq=state%q; de%vals=0
  call model%stress_and_update(sig,eps,de,dt,state)
  if(maxval(abs(state%q-exp(-dt/material%tau)*oldq))>1e-8_real64) error stop 'hold'
  if(abs(sls_factor(1e-12_real64,1.0_real64)-1)>1e-12_real64) error stop 'small dt'
  material%ratio=0; state%q=0
  sig=model%stress(eps,de,dt,state); eq=material%stress(eps)
  if(maxval(abs(sig%vals-eq%vals))>1e-8_real64) error stop 'elastic limit'
  print *, 'PASS MUSCLE SLS: polymorphism, trial/commit, tangent, relaxation, independent history'
end program
