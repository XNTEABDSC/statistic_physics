
use std::ops::AddAssign;

use fixed::traits::Fixed;

use crate::delta::Delta;
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
    mass:Num,
    /// sum(a.v*a.m)
    momentum:Vec2Fix,
    /// 1/2 sum(a.m*a.v^2)
    energy:Num,
    /// 1/2 sum(a.m*v_mean^2) = momentum*momentum/mass/2;
    kinetic:Num,
    /// energy-kinetic
    internal:Num,
    /// sum(a.m*a.v)/mass = momentum/mass
    v_mean:Vec2Fix,
    /// sum( a.m*(a.v-v_mean)^2 )/mass
    v_var_sq:Num,
    /// sqrt(v_var_sq)
    v_var:Num,
    //pub internal:Num
}

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

impl AddAssign<&Matters> for MattersState {
    fn add_assign(&mut self, rhs: &Matters) {
        *self=MattersState::new(Matters { mass: self.mass+rhs.mass, momentum: self.momentum+rhs.momentum, energy: self.energy+rhs.energy })
    }
}

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


impl std::ops::AddAssign<&Matters> for Matters {
    fn add_assign(&mut self, rhs: &Matters) {
        self.mass+=rhs.mass;
        self.momentum+=rhs.momentum;
        self.energy+=rhs.energy;
        //self.internal+=rhs.0.internal;
    }
}
impl std::ops::SubAssign<&Matters> for Matters {
    fn sub_assign(&mut self, rhs: &Matters) {
        self.mass-=rhs.mass;
        self.momentum-=rhs.momentum;
        self.energy-=rhs.energy;
        //self.internal-=rhs.0.internal;
    }
}


impl Matters {
    #[inline]
    pub fn momentum_kenetic(&self)->Num {
        mass_momentum_2_kenetic(self.momentum, self.mass)
        //kenetic-1/2momentum^2/mass
    }
    #[inline]
    pub fn internal(&self)->Num {
        self.energy-self.momentum_kenetic()
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
    // matters got pushed, assume all atom added same speed 
    pub fn from_momentum(&self,momentum:Vec2Fix)->Delta<Matters>{
        Delta(Matters{
            mass:Num::ZERO,
            momentum:momentum,
            energy:mass_momentum_2_kenetic(momentum, self.mass)
            //internal: Num::ZERO//(momentum/self.mass)*(momentum/Num::from_num(2)+self.momentum)
        })
    }
    pub fn vmean(&self)->Vec2Fix {
        self.momentum/self.mass
    }
    pub fn v_var(&self)->Num{
        self.internal()/self.mass*2
    }
}