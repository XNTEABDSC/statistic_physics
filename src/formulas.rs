
use cordic::sqrt;
use fixed::traits::Fixed;

use crate::{
    body::CircleObject, constants::NUMINV2, delta::{Change, Delta}, gas_cell::GasCell, matters::{Matters, MattersState}, num::Num, vec2_fix::{self, dist, Vec2Fix}
};
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


pub fn gas_cell_spread_to_side_another(a:&MattersState,volume:Num,dir_vec:Vec2Fix,edge:Num,dt:Num)->Delta<Matters>{
    let mass=a.mass();
    let v_mean=a.v_mean();
    let v_var_sq=a.v_var_sq();
    let v_var=a.v_var();

    let v_mean_on_dir=v_mean*dir_vec;

    let factor=v_mean_on_dir/v_var;
    
    let pass_normal_cdf=if v_var.is_zero() { Num::ONE }else{normal_cdf( factor)};
    let pass_normal_pdf=if v_var.is_zero() { Num::ZERO }else{normal_pdf(factor)};

    
    let e_pass_v0=v_mean_on_dir*pass_normal_cdf+v_var*pass_normal_pdf;

    let pass_mass=( e_pass_v0 )*mass*dt*edge/volume;

    let pass_vector=v_mean+dir_vec*v_var*pass_normal_pdf;

    let pass_momentum=pass_vector*pass_mass;

    return todo!();
}


pub fn gas_cell_spread_to_side_copy(a:&MattersState,volume:Num,dir_vec:Vec2Fix,edge:Num,dt:Num)->Delta<Matters>{
    let mass=a.mass();
    let v_mean=a.v_mean();
    let v_var_sq=a.v_var_sq();
    let v_var=a.v_var();



    let v_mean_on_dir=(v_mean*dir_vec);
    let v_mean_on_dir_vec=dir_vec*v_mean_on_dir;
    //let v_dir_vec=dir_vec.rotate90();
    let v_mean_v_dir=v_mean-v_mean_on_dir_vec;

    
    let pass_normal_cdf=if false { Num::ONE }else{normal_cdf( v_mean_on_dir)};
    let pass_normal_pdf=if false { Num::ZERO }else{normal_pdf(v_mean_on_dir)};

    let e_pass_v0=v_mean_on_dir*pass_normal_cdf+v_var*pass_normal_pdf;
    let e_pass_v1=v_mean_on_dir*v_mean_on_dir*pass_normal_cdf
    +2*v_mean_on_dir*v_var*pass_normal_pdf
    + v_var*v_var*(Num::ONE-pass_normal_cdf);
    let e_pass_v2=
    v_mean_on_dir*v_mean_on_dir*v_mean_on_dir*pass_normal_cdf
    + 3*v_mean_on_dir*v_mean_on_dir*v_var*pass_normal_pdf
    + 3*v_mean_on_dir*v_var*v_var* pass_normal_cdf
    + v_var*v_var*v_var*pass_normal_pdf
    ;


    println!("cdf: {}, pdf: {}",pass_normal_cdf,pass_normal_pdf);

    // E( max(pass,0) * p.m)
    let pass_mass=( e_pass_v0 )*mass*dt*edge/volume;

    // E( max(pass,0) * p.m*p.v*p.v*0.5)
    let pass_energy=( e_pass_v2 )*mass*dt*edge/volume*NUMINV2;

    let pass_momentum=
        dir_vec*mass*dt*edge/volume*(
            e_pass_v1
        )+
        v_mean_v_dir*( e_pass_v0 )*mass*dt*edge/volume;

    return  Delta(Matters { mass: pass_mass, momentum: pass_momentum, energy: pass_energy });

}



pub fn gas_cell_spread_to_side(a:&MattersState,volume:Num,dir_vec:Vec2Fix,edge:Num,dt:Num)->Delta<Matters>{
    let mass=a.mass();
    let v_mean=a.v_mean();
    let v_var_sq=a.v_var_sq();
    let v_var=a.v_var();
    let v_on_dir_mean=v_mean*dir_vec;
    //let v_on_dir_var=v_var;
    let mass_edge_volume_dt=mass*edge*dt/volume;

    let (cdf,pdf)=if v_var.is_zero() {
        (if v_on_dir_mean>0 {Num::ONE} else {Num::ZERO},Num::ZERO)
    }else{
        let frac=v_on_dir_mean/v_var;
        (normal_cdf(frac),normal_pdf(frac))
    };
    let v_on_dir_mean_p2=v_on_dir_mean*v_on_dir_mean;


    // E[max(v_on_dir)]
    let e_0=v_on_dir_mean*cdf+v_var*pdf;
    let e_1=(v_on_dir_mean_p2+v_var_sq)*cdf+v_on_dir_mean*v_var*pdf;
    let e_2=v_on_dir_mean*(v_on_dir_mean_p2+3*v_var_sq)*cdf+v_var*(v_on_dir_mean_p2+2*v_var_sq)*pdf;

    let pass_mass=mass_edge_volume_dt*e_0;

    let v_on_dir_vec=dir_vec*v_on_dir_mean;
    let v_v_dir_vec=v_mean-v_on_dir_vec;

    let pass_momentum=(
        dir_vec*e_1
        + v_v_dir_vec*e_0
    )*mass_edge_volume_dt;

    let pass_energy=(
        e_2+v_var*e_0
    )* Num::FRAC_1_SQRT_2 *mass_edge_volume_dt>>1;


    return Delta(Matters { mass: pass_mass, momentum: pass_momentum, energy: pass_energy });

}


pub fn gas_cell_spread_to_side_old_0(a:&MattersState,volume:Num,dir_vec:Vec2Fix,edge:Num,dt:Num)->Delta<Matters>{
    let mass=a.mass();
    let v_mean=a.v_mean();
    let v_var_sq=a.v_var_sq();
    let v_var=a.v_var();
    let pass_factor=dt/volume;
    let edge_vec=dir_vec*edge;
    let pass_mean=v_mean*edge_vec*pass_factor;
    let pass_mean_sq=pass_mean*pass_mean;
    let pass_mean_pow_3=pass_mean_sq*pass_mean;
    let pass_var_sq=v_var_sq*pass_factor*pass_factor;
    let pass_var=v_var*pass_factor;
    let frac_pass_mean_pass_var=if pass_var.is_zero() {Num::ZERO} else {pass_mean/pass_var};

    let pass_normal_cdf=if pass_var.is_zero() { Num::ONE }else{normal_cdf(frac_pass_mean_pass_var)};
    let pass_normal_pdf=if pass_var.is_zero() { Num::ZERO }else{normal_pdf(frac_pass_mean_pass_var)};

    println!("cdf: {}, pdf: {}",pass_normal_cdf,pass_normal_pdf);
    
    let e_pass=pass_mean*pass_normal_cdf+pass_var*pass_normal_pdf;

    // E( max(pass,0) * p.m)
    let pass_mass=( e_pass )*mass;

    // E( max(pass,0) * p.m*p.v*p.v*0.5)
    let pass_energy=(
        pass_mean_pow_3*pass_normal_cdf
        + 3*pass_mean_sq*pass_var*pass_normal_pdf
        + 3*pass_mean*pass_var_sq * pass_normal_cdf
        + pass_var_sq*pass_var*pass_normal_pdf
    )*mass*NUMINV2;

    let v_mean_on_dir=dir_vec*(v_mean*dir_vec);
    //let v_dir_vec=dir_vec.rotate90();
    let v_mean_v_dir=v_mean-v_mean_on_dir;

    println!("v_mean_on_dir: {:?} , v_mean_v_dir: {:?}",v_mean_on_dir,v_mean_v_dir);
    println!("add_test: {:?}",v_mean_on_dir+v_mean_v_dir);
    println!("e_pass_v1: {:?} , e_pass_v0: {:?}",(
         (pass_mean_sq)*pass_normal_cdf
        + 2*pass_var*pass_mean*pass_normal_pdf
        + pass_var_sq*(Num::ONE-pass_normal_cdf)
    ),e_pass);

    let pass_momentum=
        (v_mean_on_dir*(
            (pass_mean_sq)*pass_normal_cdf
            + 2*pass_var*pass_mean*pass_normal_pdf
            +pass_var_sq*(Num::ONE-pass_normal_cdf)
        )
        +(v_mean_v_dir*e_pass));

    return  Delta(Matters { mass: pass_mass, momentum: pass_momentum, energy: pass_energy });

}

pub fn gas_cell_spread_to_side_bad(a:&GasCell,dirvec:Vec2Fix,dt:Num)->Delta<Matters>{
    /* */
    if a.matters.mass().is_zero(){return Delta::default();}
    //if a.matters.mass.is_zero() || a.matters.mass.is_negative(){return Delta::default();}
    //let mut res=Matters::default();
    let mut per=(a.matters.v_mean()*dirvec)*dt/a.edge;
    if per<Num::ZERO {per=Num::ZERO;}
    if per>Num::ONE {per=Num::ONE;}
    //res+=a.matters.take_matters_percent(per);
    //res
    let res1=a.matters.take_matters_percent(per);
    let vel_var_sq=Num::sqrt(a.matters.v_var());
    let per2_mean_vel=vel_var_sq*2;
    let mut per2=per2_mean_vel*dt/a.edge;
    if per2<Num::ZERO {per2=Num::ZERO;}
    if per2>Num::ONE {per2=Num::ONE;}
    let per2_mass=a.matters.mass()*per2;
    let per2_momentum=dirvec*per2_mass*per2_mean_vel;
    let per2_kenetic=a.matters.energy()*per2;
    let mut res=res1;
    res.0+=&Matters{mass:per2_mass,momentum:per2_momentum,energy:per2_kenetic};
    res
}

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
}

//pub const interact_gas_cell_body_momentum_transfer:Num = ;

pub fn interact_gas_cell_body(gc:&GasCell,b:&CircleObject,half_life_period_factor_over_2_over_len:Num)->Delta<Matters>{
    let factor=-b.radius*half_life_period_factor_over_2_over_len;
    let dmomentum=b.matters.momentum()-gc.matters.momentum();
    let dmomentum_fac=dmomentum*factor;
    let dmomentum_fac_kenetic=
    mass_momentum_2_kenetic(dmomentum_fac, b.matters.mass())-mass_momentum_2_kenetic(-dmomentum_fac, gc.matters.mass());
    let dinternal=b.matters.internal()-gc.matters.internal();
    let dinternal_fac=dinternal*factor;
    Delta(Matters{
        mass:Num::ZERO,
        momentum:dmomentum_fac,
        energy:dmomentum_fac_kenetic+dinternal_fac
    })
}


pub fn push_matters_by_work(gc:&MattersState,work:Num,dir_vec:Vec2Fix,worker_speed:Vec2Fix)->Delta<Matters> {
    //let v1=(worker_speed-gc.v_mean())*dir_vec;
    let v1=(gc.v_mean()+worker_speed)*dir_vec;

    let mass=gc.mass();

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

    Delta(Matters { 
        mass: Num::ZERO, 
        momentum: dir_vec*dv*mass,
        energy: work, })
}