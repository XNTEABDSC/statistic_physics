use fixed::traits::Fixed;

use crate::{formulas::{self, gas_cell_spread_to_side, mass_momentum_2_kenetic}, matters::{Matters, MattersState}, num::Num, vec2_fix::Vec2Fix};


#[test]
fn test(){
    let mass=Num::ONE;
    let mom=Vec2Fix::new(Num::from_str("3.375").unwrap(), Num::ONE);
    let a_matters=MattersState::new(Matters{mass:mass,
        momentum:mom,
        energy:Num::from_str("291.92926").unwrap()+mass_momentum_2_kenetic(mom, mass)
    });
    let res=gas_cell_spread_to_side(&a_matters, Num::ONE, Vec2Fix::new(Num::ONE, Num::ZERO), Num::ONE, Num::ONE>>5);

    println!("mass: {} , momentum: {:?}",res.0.mass*Num::ONE<<5 , res.0.momentum* (Num::ONE<<5));
    
}