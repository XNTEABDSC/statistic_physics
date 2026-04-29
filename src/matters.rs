
use frunk::HList;
use physics_basic::stats::*;
use crate::stats::*;


/// 
/// 
/// They will be used to record changes, and used to calculate [MattersFull].
pub type MattersBasic<Num,const DIM:usize>=HList!(Mass<Num>,Momentum<Num,DIM>,Energy<Num>,Volume<Num>);

/// The basic stats of matters: [Mass] [Momentum] [Energy]
/// 
/// They will be used to record changes, and used to calculate [MattersFull].
pub type MattersBasicStat<Num,const DIM:usize>=HList!(Mass<Num>,Momentum<Num,DIM>,Energy<Num>);

/// All useful stats of matters
pub type MattersFull<Num,const DIM: usize>=HList!(
	Mass<Num>,
	Momentum<Num,DIM>,
	Energy<Num>,
    Vel<Num,DIM>,
    Kinetic<Num>,
    Internal<Num>,
    VelVarSq<Num>,
    VelVar<Num>,
    VelVarSq1Dir<Num>,
    VelVar1Dir<Num>,
	Volume<Num>,
	Density<Num>
);