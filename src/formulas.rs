
use cordic::sqrt;
use frunk::hlist;
use frunk::HList;
use physics_basic::stats::*;
use wacky_bag_fixed::utils::normal_cdf::normal_cdf;
use wacky_bag_fixed::utils::normal_pdf::normal_pdf;
use wacky_bag_fixed::vec_fix::VecFix;
use crate::stats::*;
//use wacky_bag::structures::delta::Delta;

type Delta<T>=T;

use crate::{constants::*, matters::Matters, num::Num};
#[inline]
pub fn mass_momentum_2_kenetic<const DIM:usize>(momentum:VecFix<DIM>,mass:Num)->Num {
    if mass.is_zero(){return Num::ZERO;}

    (momentum/mass).dot(&momentum)*NUMINV2

    //momentum/mass*momentum*NUMINV2
}

#[inline]
pub fn mass_kinetic_2_momentum<const DIM:usize>(kenetic:Num,mass:Num,dir_vec:VecFix<DIM>)->VecFix<DIM> {
    dir_vec*Num::sqrt(kenetic*mass*2)
}



/// for time dt, for ['MattersState'] with `volume`, to calculate how much [`Matters`] will cross the edge with `edge_len` and `edge_dir_vec` 
pub fn gas_cell_spread_to_side<const DIM:usize>(a:HList!(&Mass,&Vel<DIM>,&VelVarSq1Dir,&VelVar1Dir),volume:Num,edge_dir_vec:VecFix<DIM>,edge_len:Num,dt:Num)->Delta<Matters<DIM>>{
    let (mass_,v_mean_,v_var_sq_1dir_,v_var_1dir_)=a.into();
    let mass=mass_.0;
    let v_mean=v_mean_.0;
    let v_var_sq_1dir=v_var_sq_1dir_.0;
    let v_var_1dir=v_var_1dir_.0;
    let v_on_dir_mean=v_mean.dot(&edge_dir_vec);
    
    let mass_edge_volume_dt=mass*edge_len*dt/volume;

    let (cdf,pdf)=if v_var_1dir.is_zero() {
        ( if v_on_dir_mean>0 {Num::ONE} else {Num::ZERO} , Num::ZERO)
    }else{
        let frac=v_on_dir_mean/v_var_1dir;
        (normal_cdf(frac),normal_pdf(frac))
    };
    let v_on_dir_mean_p2=v_on_dir_mean*v_on_dir_mean;


    // E[max(v_on_dir,0)]
    let e_0=v_on_dir_mean*cdf+v_var_1dir*pdf;
    let e_1=(v_on_dir_mean_p2+v_var_sq_1dir)*cdf+v_on_dir_mean*v_var_1dir*pdf;
    let e_2=v_on_dir_mean*(v_on_dir_mean_p2+3*v_var_sq_1dir)*cdf+v_var_1dir*(v_on_dir_mean_p2+2*v_var_sq_1dir)*pdf;

    let pass_mass=mass_edge_volume_dt*e_0;

    let v_on_dir_vec=edge_dir_vec*v_on_dir_mean;
    let v_v_dir_vec=v_mean-v_on_dir_vec;

    let pass_momentum=(
        edge_dir_vec*e_1
        + v_v_dir_vec*e_0
    )*mass_edge_volume_dt;

    let pass_energy=(
        e_2+v_var_1dir*e_0
    )* Num::FRAC_1_SQRT_2 *mass_edge_volume_dt>>1;


    return //Delta
    hlist![Mass(pass_mass),Momentum(pass_momentum),Energy(pass_energy)];
    //Matters { mass: pass_mass, momentum: pass_momentum, energy: pass_energy };

}
/// for time dt, for ['MattersState'] with `volume`, to calculate how much [`Matters`] will cross the edge with `edge_len` and `edge_dir_vec` 

/*
pub fn gas_cell_spread(a:&GasCell,dt:Num,c:&mut Change<Matters>,n:&mut Change<Matters>,e:&mut Change<Matters>,s:&mut Change<Matters>,w:&mut Change<Matters>)->(){
    //let mut res=GasCellSpreadResult::default();
    let dme=gas_cell_spread_to_side(&a.matters,Num::ONE, Vec2Fix::new(Num::ONE, Num::ZERO),Num::ONE,dt);
    dme.transfer(c, e);
    let dmn=gas_cell_spread_to_side(&a.matters,Num::ONE, Vec2Fix::new(Num::ZERO, Num::ONE),Num::ONE,dt);
    dmn.transfer(c, n);
    let dmw=gas_cell_spread_to_side(&a.matters,Num::ONE, Vec2Fix::new(Num::NEG_ONE, Num::ZERO),Num::ONE,dt);
    dmw.transfer(c, w);
    let dms=gas_cell_spread_to_side(&a.matters,Num::ONE, Vec2Fix::new(Num::ZERO, Num::NEG_ONE),Num::ONE,dt);
    dms.transfer(c, s);
} */
 
//pub const interact_gas_cell_body_momentum_transfer:Num = ;

/// just a exist way
pub fn interact_gas_cell_body<const DIM:usize>(gc_m:(&Mass,&Momentum<DIM>,&Internal),b_m:(&Mass,&Momentum<DIM>,&Internal),b_radius:Num,half_life_period_factor_over_2_over_len:Num)->Delta<Matters<DIM>>{
    let gc_m_mass=gc_m.0.0;
    let gc_m_momentum=gc_m.1.0;
    let gc_m_internal=gc_m.2.0;
    let b_m_mass=b_m.0.0;
    let b_m_momentum=b_m.1.0;
    let b_m_internal=b_m.2.0;
    let factor=-b_radius*half_life_period_factor_over_2_over_len;
    let dmomentum=b_m_momentum-gc_m_momentum;
    let dmomentum_fac=dmomentum*factor;
    let dmomentum_fac_kenetic=
    mass_momentum_2_kenetic(dmomentum_fac, b_m_mass)-mass_momentum_2_kenetic(-dmomentum_fac, gc_m_mass);
    let dinternal=b_m_internal-gc_m_internal;
    let dinternal_fac=dinternal*factor;
    //Delta
    hlist![Mass(Num::ZERO),Momentum(dmomentum_fac),Energy(dmomentum_fac_kenetic+dinternal_fac)]
    /* 
    Matters{
        mass:Num::ZERO,
        momentum:dmomentum_fac,
        energy:dmomentum_fac_kenetic+dinternal_fac
    }*/
}

pub fn push_matters_by_work<const DIM:usize>(gc:(&Vel<DIM>,&Mass),work:(&Kinetic,&DirVec<DIM>,Vel<DIM>))->Delta<Matters<DIM>> {
    let (work_kinetic,dir_vec,worker_speed)=work;
    //let v1=(worker_speed-gc.v_mean())*dir_vec;
    let v1=(gc.0.0+worker_speed.0).dot( &dir_vec.0);

    let mass=gc.1.0;

    //v^2 = pmt

    // offset_v is init v, get final v when did `work` work
    // v2^2-v1^2=mass*work
    // mass v2^2=work+mass v1^2

    let v1_sq=v1*v1;
    let v2_sq=work_kinetic.0/mass+v1_sq;
    let v2=sqrt(v2_sq);
    let dv=v2-v1;


    // f=2p/v
    // I = 2 work / vel

    //Delta
    hlist![Mass(Num::ZERO),Momentum(dir_vec.0*dv*mass),Energy(work_kinetic.0)]
    // Matters { 
    //     mass: Mass(Num::ZERO), 
    //     momentum: Momentum(dir_vec.0*dv*mass),
    //     energy: Energy(work_kinetic.0), }
}

/*
    
    mass:&Mass,
    momentum:&Momentum,
    energy:&Energy,
    vel:&mut Vel,
    kinetic:&mut Kinetic,
    internal:&mut Internal,
    vel_var_sq:&mut VelVarSq,
    vel_var:&mut VelVar,
    vel_var_sq_1dir:&mut VelVarSq1Dir,
    vel_var_1dir:&mut VelVar1Dir,
*/

pub fn calculate_matters_state<const DIM:usize>(matters:HList!(&Mass,
    &Momentum<DIM>,
    &Energy,
    &mut Vel<DIM>,
    &mut Kinetic,
    &mut Internal,
    &mut VelVarSq,
    &mut VelVar,
    &mut VelVarSq1Dir,
    &mut VelVar1Dir)){

    let (mass,momentum,energy,vel,kinetic,internal,vel_var_sq,vel_var,vel_var_sq_1dir,vel_var_1dir)=matters.into();
    vel.0=momentum.0/mass.0;
    kinetic.0=  mass_momentum_2_kenetic(momentum.0,mass.0);//(momentum.0*momentum.0/mass.0) *NUMINV2;
    internal.0=energy.0-kinetic.0;
    vel_var_sq.0=2*internal.0/mass.0;
    vel_var.0=Num::sqrt(vel_var_sq.0);
    vel_var_sq_1dir.0=vel_var_sq.0/(DIM as i64);
    vel_var_1dir.0=vel_var.0/(SQRTS[DIM]);
    //calculate_matters_state::<DIM>(matters.into())
}
