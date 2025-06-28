


use crate::num::Num;
use crate::vec2_fix::Vec2Fix;


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
