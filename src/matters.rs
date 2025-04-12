
use std::ops::AddAssign;

use fixed::traits::Fixed;

use crate::formulas::mass_momentum_2_kenetic;
use crate::num::Num;
use crate::vec2_fix::Vec2Fix;


/// mass = sum(a.m)
/// momentum = sum(a.v*a.m)
/// vmean = momentum/mass
/// kenetic = 1/2 sum(a/m*a.v^2)
/// internal = 1/2 sum(a.m*(a.v-vmean)^2) = kenetic-1/2momentum^2/mass
#[derive(Default,Clone,Debug)]
pub struct MattersState{
    /// total mass of all matters
    /// sum(a.m)
    pub mass:Num,
    /// sum(a.v*a.m)
    pub momentum:Vec2Fix,
    /// 1/2 sum(a.m*a.v^2)
    pub energy:Num,
    /// 1/2 sum(a.m*v_mean^2) = momentum*momentum/mass/2;
    pub kinetic:Num,
    /// energy-kinetic
    pub internal:Num,
    /// sum(a.m*a.v)/mass = momentum/mass
    pub v_mean:Vec2Fix,
    /// sum( a.m*(a.v-v_mean)^2 )/mass
    pub v_var_sq:Num,
    /// sqrt(v_var_sq)
    pub v_var:Num,
    /// sum( a.m*(a.v-v_mean)^2 )/mass
    pub v_var_sq_1dir:Num,
    /// sqrt(v_var_sq)
    pub v_var_1dir:Num,
    //pub internal:Num
}
/*
impl MattersState {
    pub fn new(matters:Matters)->MattersState{
        let mass=matters.mass;
        let momentum=matters.momentum;
        let energy=matters.energy;
        let kinetic=mass_momentum_2_kenetic(momentum, mass);
        let internal=energy-kinetic;
        let v_mean=if mass.is_zero() {Vec2Fix::default()} else { momentum/mass};
        let v_var_sq=if mass.is_zero() {Num::ZERO}else{internal/mass*2};
        let v_var=  Num::sqrt(v_var_sq );
        MattersState{
            mass,momentum,energy,kinetic,internal,v_mean,v_var_sq,v_var
        }
    }
    
    pub fn mass(&self) -> Num {
        self.mass
    }
    
    pub fn momentum(&self) -> Vec2Fix {
        self.momentum
    }
    
    pub fn energy(&self) -> Num{
        self.energy
    }
    
    pub fn kinetic(&self) -> Num {
        self.kinetic
    }
    
    pub fn internal(&self) -> Num {
        self.internal
    }
    
    pub fn v_mean(&self) -> Vec2Fix {
        self.v_mean
    }
    
    pub fn v_var_sq(&self) -> Num {
        self.v_var_sq
    }
    
    pub fn v_var(&self) -> Num {
        self.v_var
    }

    #[inline]
    pub fn take_matters_percent(&self,per:Num)->Delta <Matters>{
        Delta(Matters{
            mass:self.mass*per,
            momentum:self.momentum*per,
            energy:self.energy*per
            //internal:self.internal*per
        })
    }
    #[inline]
    pub fn take_matters(&self,mass:Num)->Delta<Matters> {
        let per=mass/self.mass;
        Delta(Matters{
            mass,
            momentum:self.momentum*per,
            energy:self.energy*per,
            //internal:self.internal*per
        })
    }

}
 */
/* */
#[derive(Default,Clone, Copy,Debug)]
pub struct Matters{
    /// total mass of all matters
    /// sum(a.m)
    pub mass:Num,
    /// sum(a.v*a.m)
    pub momentum:Vec2Fix,
    /// 1/2 sum(a.m*a.v^2)
    pub energy:Num,
    // 1/2 sum(a.m*(a.v-vmean)^2)
    
    //pub internal:Num
}
