//use wacky_bag::derive_add_traits;

// use crate::num::Num;


use derive_more::{Add,AddAssign,Sub,SubAssign,Neg};

#[derive(Default,Clone,Copy,Debug,Add,AddAssign,Sub,SubAssign,Neg)]
pub struct Internal<Num>(pub Num);

//derive_add_traits!(Internal);

#[derive(Default, Clone, Copy, Debug)]
pub struct VelVarSq<Num>(pub Num);

//derive_add_traits!(VelVarSq);

#[derive(Default, Clone, Copy, Debug)]
pub struct VelVar<Num>(pub Num);

//derive_add_traits!(VelVar);

#[derive(Default, Clone, Copy, Debug)]
pub struct VelVarSq1Dir<Num>(pub Num);

//derive_add_traits!(VelVarSq1Dir);

#[derive(Default, Clone, Copy, Debug)]
pub struct VelVar1Dir<Num>(pub Num);

//derive_add_traits!(VelVar1Dir);
