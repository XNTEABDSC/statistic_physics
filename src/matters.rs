
use crate::delta::Delta;
use crate::formulas::momentum_kenetic;
use crate::num::Num;
use crate::vec2_fix::Vec2Fix;
/// mass = sum(a.m)
/// momentum = sum(a.v*a.m)
/// vmean = momentum/mass
/// kenetic = 1/2 sum(a/m*a.v^2)
/// internal = 1/2 sum(a.m*(a.v-vmean)^2) = kenetic-1/2momentum^2/mass
#[derive(Default,Clone, Copy,Debug)]
pub struct Matters{
    /// total mass of all matters
    /// sum(a.m)
    pub mass:Num,
    /// sum(a.v*a.m)
    pub momentum:Vec2Fix,
    /// 1/2 sum(a.m*a.v^2)
    pub kinetic:Num,
    // 1/2 sum(a.m*(a.v-vmean)^2)
    
    //pub internal:Num
}


impl std::ops::AddAssign<&Matters> for Matters {
    fn add_assign(&mut self, rhs: &Matters) {
        self.mass+=rhs.mass;
        self.momentum+=rhs.momentum;
        self.kinetic+=rhs.kinetic;
        //self.internal+=rhs.0.internal;
    }
}
impl std::ops::SubAssign<&Matters> for Matters {
    fn sub_assign(&mut self, rhs: &Matters) {
        self.mass-=rhs.mass;
        self.momentum-=rhs.momentum;
        self.kinetic-=rhs.kinetic;
        //self.internal-=rhs.0.internal;
    }
}

impl Matters {
    #[inline]
    pub fn momentum_kenetic(&self)->Num {
        momentum_kenetic(self.momentum, self.mass)
        //kenetic-1/2momentum^2/mass
    }
    #[inline]
    pub fn internal(&self)->Num {
        self.kinetic-self.momentum_kenetic()
    }
    #[inline]
    pub fn take_matters_percent(&self,per:Num)->Delta <Matters>{
        Delta(Matters{
            mass:self.mass*per,
            momentum:self.momentum*per,
            kinetic:self.kinetic*per
            //internal:self.internal*per
        })
    }
    #[inline]
    pub fn take_matters(&self,mass:Num)->Delta<Matters> {
        let per=mass/self.mass;
        Delta(Matters{
            mass,
            momentum:self.momentum*per,
            kinetic:self.kinetic*per,
            //internal:self.internal*per
        })
    }
    // matters got pushed, assume all atom added same speed 
    pub fn from_momentum(&self,momentum:Vec2Fix)->Delta<Matters>{
        Delta(Matters{
            mass:Num::ZERO,
            momentum:momentum,
            kinetic:momentum_kenetic(momentum, self.mass)
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