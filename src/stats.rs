use wacky_bag::derive_add_traits;

use crate::num::Num;

#[derive(Default, Clone, Copy, Debug)]
pub struct Internal(pub Num);

derive_add_traits!(Internal);

#[derive(Default, Clone, Copy, Debug)]
pub struct VelVarSq(pub Num);

derive_add_traits!(VelVarSq);

#[derive(Default, Clone, Copy, Debug)]
pub struct VelVar(pub Num);

derive_add_traits!(VelVar);

#[derive(Default, Clone, Copy, Debug)]
pub struct VelVarSq1Dir(pub Num);

derive_add_traits!(VelVarSq1Dir);

#[derive(Default, Clone, Copy, Debug)]
pub struct VelVar1Dir(pub Num);

derive_add_traits!(VelVar1Dir);
