
use crate::{
    body::CircleObject, constants::{HALF_LIFE_PERIOD_FRAME_FACTOR_OVER2, NUM2, STD_LENGTH}, delta::{Change, Delta, State}, gas_cell::GasCell, matters::Matters, num::Num, vec2_fix::Vec2Fix
};
#[inline]
pub fn momentum_kenetic(momentum:Vec2Fix,mass:Num)->Num {
    momentum*momentum/mass/NUM2
    //kenetic-1/2momentum^2/mass
}

pub fn gas_cell_spread_to_side(a:&GasCell,dirvec:Vec2Fix,dt:Num)->Delta<Matters>{
    //let mut res=Matters::default();
    let vel_var_sq=Num::sqrt(a.matters.v_var());
    let mut per=(a.matters.vmean()*dirvec+vel_var_sq)*dt/a.edge;
    if per<Num::ZERO {per=Num::ZERO;}
    if per>Num::ONE {per=Num::ONE;}
    //res+=a.matters.take_matters_percent(per);
    //res
    a.matters.take_matters_percent(per)
}
pub fn gas_cell_spread(a:&GasCell,dt:Num,c:&mut Change<Matters>,n:&mut Change<Matters>,e:&mut Change<Matters>,s:&mut Change<Matters>,w:&mut Change<Matters>)->(){
    //let mut res=GasCellSpreadResult::default();
    let dme=gas_cell_spread_to_side(a, Vec2Fix::new(Num::ONE, Num::ZERO),dt);
    dme.transfer(c, e);
    let dmn=gas_cell_spread_to_side(a, Vec2Fix::new(Num::ZERO, Num::ONE),dt);
    dmn.transfer(c, n);
    let dmw=gas_cell_spread_to_side(a, Vec2Fix::new(Num::NEG_ONE, Num::ZERO),dt);
    dmw.transfer(c, w);
    let dms=gas_cell_spread_to_side(a, Vec2Fix::new(Num::ZERO, Num::NEG_ONE),dt);
    dms.transfer(c, s);
}

//pub const interact_gas_cell_body_momentum_transfer:Num = ;

pub fn interact_gas_cell_body(gc:&GasCell,b:&CircleObject,half_life_period_factor_over_2_over_len:Num)->Delta<Matters>{
    let factor=-b.radius*half_life_period_factor_over_2_over_len;
    let dmomentum=b.matters.momentum-gc.matters.momentum;
    let dmomentum_fac=dmomentum*factor;
    let dmomentum_fac_kenetic=
    momentum_kenetic(dmomentum_fac, b.matters.mass)-momentum_kenetic(-dmomentum_fac, gc.matters.mass);
    let dinternal=b.matters.internal()-gc.matters.internal();
    let dinternal_fac=dinternal*factor;
    Delta(Matters{
        mass:Num::ZERO,
        momentum:dmomentum_fac,
        kinetic:dmomentum_fac_kenetic+dinternal_fac
    })
}