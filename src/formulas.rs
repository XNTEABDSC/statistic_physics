
use cordic::sqrt;
//use wacky_bag::structures::delta::Delta;

type Delta<T>=T;

use crate::{constants::*, matters::{Matters, MattersState}, num::Num, vec2_fix::Vec2Fix};
#[inline]
pub fn mass_momentum_2_kenetic(momentum:Vec2Fix,mass:Num)->Num {
    if mass.is_zero(){return Num::ZERO;}
    momentum/mass*momentum*NUMINV2
}

#[inline]
pub fn mass_kinetic_2_momentum(kenetic:Num,mass:Num,dir_vec:Vec2Fix)->Vec2Fix {
    dir_vec*Num::sqrt(kenetic*mass*2)
}

pub fn normal_pdf(x:Num)->Num{
    return Num::FRAC_1_SQRT_2PI * cordic::exp(-x/2*x)
}

pub fn normal_cdf(x:Num)->Num {
    const A1: Result<Num, fixed::ParseFixedError> = Num::from_str("0.254829592");
    const A2: Result<Num, fixed::ParseFixedError> = Num::from_str("-0.284496736");
    const A3: Result<Num, fixed::ParseFixedError> = Num::from_str("1.421413741");
    const A4: Result<Num, fixed::ParseFixedError> = Num::from_str("-1.453152027");
    const A5: Result<Num, fixed::ParseFixedError> = Num::from_str("1.061405429");
    const P: Result<Num, fixed::ParseFixedError> = Num::from_str("0.3275911");
    let sign = Num::signum(x);
    let x2=x.abs()*Num::FRAC_1_SQRT_2;
    let t= Num::ONE/(Num::ONE+P.unwrap()*x2);
    let y=Num::ONE-(((((A5.unwrap()*t + A4.unwrap())*t) + A3.unwrap())*t + A2.unwrap())*t + A1.unwrap())*t*cordic::exp(-x*x);
    return (Num::ONE+sign*y)/Into::<Num>::into(2);
}



/// for time dt, for ['MattersState'] with `volume`, to calculate how much [`Matters`] will cross the edge with `edge_len` and `edge_dir_vec` 
pub fn gas_cell_spread_to_side(a:&MattersState,volume:Num,edge_dir_vec:Vec2Fix,edge_len:Num,dt:Num)->Delta<Matters>{
    let mass=a.mass;
    let v_mean=a.v_mean;
    let v_var_sq_1dir=a.v_var_sq_1dir;
    let v_var_1dir=a.v_var_1dir;
    let v_on_dir_mean=v_mean*edge_dir_vec;
    
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
    (Matters { mass: pass_mass, momentum: pass_momentum, energy: pass_energy });

}
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
pub fn interact_gas_cell_body(gc_m:MattersState,b_m:MattersState,b_radius:Num,half_life_period_factor_over_2_over_len:Num)->Delta<Matters>{
    let factor=-b_radius*half_life_period_factor_over_2_over_len;
    let dmomentum=b_m.momentum-gc_m.momentum;
    let dmomentum_fac=dmomentum*factor;
    let dmomentum_fac_kenetic=
    mass_momentum_2_kenetic(dmomentum_fac, b_m.mass)-mass_momentum_2_kenetic(-dmomentum_fac, gc_m.mass);
    let dinternal=b_m.internal-gc_m.internal;
    let dinternal_fac=dinternal*factor;
    //Delta
    (Matters{
        mass:Num::ZERO,
        momentum:dmomentum_fac,
        energy:dmomentum_fac_kenetic+dinternal_fac
    })
}

pub fn push_matters_by_work(gc:&MattersState,work:Num,dir_vec:Vec2Fix,worker_speed:Vec2Fix)->Delta<Matters> {
    //let v1=(worker_speed-gc.v_mean())*dir_vec;
    let v1=(gc.v_mean+worker_speed)*dir_vec;

    let mass=gc.mass;

    //v^2 = pmt

    // offset_v is init v, get final v when did `work` work
    // v2^2-v1^2=mass*work
    // mass v2^2=work+mass v1^2

    let v1_sq=v1*v1;
    let v2_sq=work/mass+v1_sq;
    let v2=sqrt(v2_sq);
    let dv=v2-v1;


    // f=2p/v
    // I = 2 work / vel

    //Delta
    (Matters { 
        mass: Num::ZERO, 
        momentum: dir_vec*dv*mass,
        energy: work, })
}