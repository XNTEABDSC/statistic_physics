

use frunk::{HList, hlist, hlist_pat};
use nalgebra::{Const, DefaultAllocator, DimName, SVector, allocator::Allocator};
use physics_basic::{body::calculate_position_state, rotation::{AngularInertia, AngularMomentum, AngularVel, DimNameToSoDimName, DimNameToSoDimNameType, angular_kinetic_from_inertia_agv, angular_momentum_from_omega}, stats::*};
use simba::scalar::RealField;
use wacky_bag::{math::{normal_cdf::{NormalCdfConsts, normal_cdf}, normal_pdf::normal_pdf}, utils::{h_list_helpers::{HToMut, HToRef}, num_extend::NumExtends}};

use crate::matters::MattersBasic;
use crate::matters::MattersFull;
use crate::stats::*;
//use wacky_bag::structures::delta::Delta;

type Delta<T>=T;

use crate::{matters::MattersBasicStat};
#[inline]
pub fn mass_momentum_2_kenetic<Num:Copy+RealField+num_traits::Num,const DIM:usize>(momentum:SVector<Num,DIM>,mass:Num)->Num {
    if mass.is_zero(){return Num::zero();}

    (momentum/mass).dot(&momentum)/Num::p2()

    //momentum/mass*momentum*NUMINV2
}

#[inline]
pub fn mass_kinetic_2_momentum<Num:Copy+RealField+num_traits::Num,const DIM:usize>(kenetic:Num,mass:Num,dir_vec:SVector<Num,DIM>)->SVector<Num,DIM> {
    dir_vec*Num::sqrt(kenetic*mass* Num::p2() )
}

pub type GasCellSpreadToSideType<Num,const DIM:usize>=HList!(Density<Num>,Vel<Num,DIM>,VelVarSq1Dir<Num>,VelVar1Dir<Num>);

/// for time dt, for ['MattersState'] with `volume`, to calculate how much [`Matters`] will cross the edge with `edge_len` and `edge_dir_vec` 
pub fn gas_cell_spread_to_side<Num:Copy+RealField+num_traits::Num+NormalCdfConsts<Marker>,const DIM:usize,Marker>(a:HToRef<GasCellSpreadToSideType<Num,DIM>>,edge_dir_vec:SVector<Num,DIM>,edge_len:Num,dt:Num)->Delta<MattersBasicStat<Num,DIM>>{
    let (density_,v_mean_,v_var_sq_1dir_,v_var_1dir_)=a.into();
    let density=density_.0;
    let v_mean=v_mean_.0;
    let v_var_sq_1dir=v_var_sq_1dir_.0;
    let v_var_1dir=v_var_1dir_.0;
    let v_on_dir_mean=v_mean.dot(&edge_dir_vec);
    
    let mass_edge_volume_dt=density*edge_len*dt;

    let (cdf,pdf)=if v_var_1dir.is_zero() {
        ( if v_on_dir_mean.is_positive() {Num::one()} else {Num::zero()} , Num::zero())
    }else{
        let frac=v_on_dir_mean/v_var_1dir;
        (normal_cdf(frac),normal_pdf(frac))
    };
    // let v_on_dir_mean_p2=v_on_dir_mean*v_on_dir_mean;

	// \mu = v_on_dir_mean
	// \sigma = v_var_1dir

    // E[v_on_dir*max(v_on_dir,0)]
    let g_1=v_on_dir_mean*cdf+v_var_1dir*pdf;
    // // E[v_on_dir^2*max(v_on_dir,0)]
    // let e_2=(v_on_dir_mean_p2+v_var_sq_1dir)*cdf+v_on_dir_mean*v_var_1dir*pdf;
    // // E[v_on_dir^3*max(v_on_dir,0)]
    // let e_3=v_on_dir_mean*(v_on_dir_mean_p2+3*v_var_sq_1dir)*cdf+v_var_1dir*(v_on_dir_mean_p2+2*v_var_sq_1dir)*pdf;


    // let v_on_dir_vec=edge_dir_vec*v_on_dir_mean;
    // let v_v_dir_vec=v_mean-v_on_dir_vec;

	
    let pass_mass=mass_edge_volume_dt*g_1;

    let pass_momentum=(
        // edge_dir_vec*e_2
        // + v_v_dir_vec*g_1
		v_mean*g_1
		+edge_dir_vec*(v_var_sq_1dir*cdf)
    )*mass_edge_volume_dt;

    let pass_energy=(
        // e_3+v_var_1dir*g_1
		v_var_sq_1dir*cdf + g_1*(v_var_sq_1dir* Num::from_isize(DIM as isize-1).unwrap()+v_mean.dot(&v_mean))
    )/ Num::sqrt(Num::p2()) *mass_edge_volume_dt/Num::p2();
	//Num::FRAC_1_SQRT_2

    return //Delta
    hlist![Mass(pass_mass),Momentum(pass_momentum),Energy(pass_energy)];
    //Matters { mass: pass_mass, momentum: pass_momentum, energy: pass_energy };

}

/*
/// for time dt, for ['MattersState'] with `volume`, to calculate how much [`Matters`] will cross the edge with `edge_len` and `edge_dir_vec` 
pub fn gas_cell_spread_to_side<const DIM:usize>(a:HList!(&Mass,&Vel<DIM>,&VelVarSq1Dir,&VelVar1Dir),volume:Num,edge_dir_vec:VecFix<DIM>,edge_len:Num,dt:Num)->Delta<Matters<DIM>>{
    let (mass_,v_mean_,v_var_sq_1dir_,v_var_1dir_)=a.into();
    let mass=mass_.0;
    let v_mean=v_mean_.0;
    let v_var_sq_1dir=v_var_sq_1dir_.0;
    let v_var_1dir=v_var_1dir_.0;
    let v_on_dir_mean=v_mean.dot(&edge_dir_vec);
    
    let mass_edge_volume_dt=mass*edge_len*dt/volume;

    let (cdf,pdf)=if v_var_1dir.is_zero() {
        ( if v_on_dir_mean>0 {Num::ONE} else {Num::ZERO} , Num::ZERO)
    }else{
        let frac=v_on_dir_mean/v_var_1dir;
        (normal_cdf(frac),normal_pdf(frac))
    };
    let v_on_dir_mean_p2=v_on_dir_mean*v_on_dir_mean;


    // E[max(v_on_dir,0)]
    let e_0=v_on_dir_mean*cdf+v_var_1dir*pdf;
    let e_1=(v_on_dir_mean_p2+v_var_sq_1dir)*cdf+v_on_dir_mean*v_var_1dir*pdf;
    let e_2=v_on_dir_mean*(v_on_dir_mean_p2+3*v_var_sq_1dir)*cdf+v_var_1dir*(v_on_dir_mean_p2+2*v_var_sq_1dir)*pdf;

    let pass_mass=mass_edge_volume_dt*e_0;

    let v_on_dir_vec=edge_dir_vec*v_on_dir_mean;
    let v_v_dir_vec=v_mean-v_on_dir_vec;

    let pass_momentum=(
        edge_dir_vec*e_1
        + v_v_dir_vec*e_0
    )*mass_edge_volume_dt;

    let pass_energy=(
        e_2+v_var_1dir*e_0
    )* Num::FRAC_1_SQRT_2 *mass_edge_volume_dt>>1;


    return //Delta
    hlist![Mass(pass_mass),Momentum(pass_momentum),Energy(pass_energy)];
    //Matters { mass: pass_mass, momentum: pass_momentum, energy: pass_energy };

}
*/

/// for time dt, for ['MattersState'] with `volume`, to calculate how much [`Matters`] will cross the edge with `edge_len` and `edge_dir_vec` 

/*
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
} */
 
//pub const interact_gas_cell_body_momentum_transfer:Num = ;

pub type InteractGasCellBodyGridCellMatters<Num,const DIM:usize>=HList!(
	Mass<Num>,Vel<Num,DIM>,
	VelVarSq<Num>,
	Density<Num>);
pub type InteractGasCellBodyBodyMatters<Num,const DIM:usize>=HList!(
	Mass<Num>,Vel<Num,DIM>,
	AngularInertia<Num,DIM>,AngularVel<Num,DIM>,
	VelVarSq<Num>,
	Density<Num>);
pub type InteractGasCellBodyGridCellChange<Num,const DIM:usize> = MattersBasicStat<Num,DIM>;
pub type InteractGasCellBodyBodyChange<Num,const DIM:usize> = HList!(Momentum<Num,DIM>,AngularMomentum<Num,DIM>,Energy<Num>);

/// just a way
pub fn interact_gas_cell_body_simple<Num:RealField+Copy,const DIM:usize>(
	hlist_pat![
		_gc_m_mass,gc_m_vel,
		gc_m_vel_var_sq,
		gc_m_density]:
	HToRef<InteractGasCellBodyGridCellMatters<Num,DIM>>,

	hlist_pat![
		b_m_mass,b_m_vel,
		b_m_agi,b_n_agv,
		b_m_vel_var_sq,
		_b_m_density]:
	HToRef<InteractGasCellBodyBodyMatters<Num,DIM>>,
	
	half_life_period_factor_over_len:Num
)->(InteractGasCellBodyGridCellChange<Num,DIM>,InteractGasCellBodyBodyChange<Num,DIM>)
	where 
    Const<DIM>: DimNameToSoDimName + DimName,
    DefaultAllocator: Allocator<DimNameToSoDimNameType<DIM>>
        + Allocator<DimNameToSoDimNameType<DIM>, DimNameToSoDimNameType<DIM>>,
{
	
    let factor=
		// *gc_m_mass / *gc_m_volume * *b_m_volume / *b_m_mass
		gc_m_density.0 / b_m_mass.0
		* half_life_period_factor_over_len
	;
	let d_vel = *b_m_vel-*gc_m_vel;
	let d_vel_fac = Vel(d_vel.0*factor);

	let hlist_pat![d_momentum] = mass_vel_2_momentum(hlist![b_m_mass,&d_vel_fac]);
	let hlist_pat![d_kinetic] = mass_vel_2_kinetic(hlist![b_m_mass,&d_vel_fac]);

	let d_agv_fac = AngularVel(b_n_agv.0*factor);
	let hlist_pat![d_ag_mom] = angular_momentum_from_omega(hlist![b_m_agi,&d_agv_fac]);
	let hlist_pat![d_ag_kin] = angular_kinetic_from_inertia_agv(hlist![b_m_agi,&d_agv_fac]);

	let d_vel_var_sq=*b_m_vel_var_sq-*gc_m_vel_var_sq;
	let d_vel_var_sq_fac=VelVarSq(d_vel_var_sq.0*factor);

	let hlist_pat![d_internal] = mass_vel_var_sq_2_internal(hlist![b_m_mass,&d_vel_var_sq_fac]);

	let d_energy=Energy(d_kinetic.0+d_ag_kin.0+d_internal.0);

	return (
		hlist![Mass(Num::zero()),-d_momentum,-d_energy],
		hlist![d_momentum,d_ag_mom,d_energy]
	);
}

pub fn push_matters_by_work<Num:RealField+Copy,const DIM:usize>(gc:(&Vel<Num,DIM>,&Mass<Num>),work:(&Kinetic<Num>,&DirVec<Num,DIM>,Vel<Num,DIM>))->Delta<MattersBasicStat<Num,DIM>> {
    let (work_kinetic,dir_vec,worker_speed)=work;
    //let v1=(worker_speed-gc.v_mean())*dir_vec;
    let v1=(gc.0.0+worker_speed.0).dot( &dir_vec.0);

    let mass=gc.1.0;

    //v^2 = pmt

    // offset_v is init v, get final v when did `work` work
    // v2^2-v1^2=mass*work
    // mass v2^2=work+mass v1^2

    let v1_sq=v1*v1;
    let v2_sq=work_kinetic.0/mass+v1_sq;
    let v2=v2_sq.sqrt();
    let dv=v2-v1;


    // f=2p/v
    // I = 2 work / vel

    //Delta
    hlist![Mass(Num::zero()),Momentum(dir_vec.0*dv*mass),Energy(work_kinetic.0)]
    // Matters { 
    //     mass: Mass(Num::ZERO), 
    //     momentum: Momentum(dir_vec.0*dv*mass),
    //     energy: Energy(work_kinetic.0), }
}

/*
    
    mass:&Mass,
    momentum:&Momentum,
    energy:&Energy,
    vel:&mut Vel,
    kinetic:&mut Kinetic,
    internal:&mut Internal,
    vel_var_sq:&mut VelVarSq,
    vel_var:&mut VelVar,
    vel_var_sq_1dir:&mut VelVarSq1Dir,
    vel_var_1dir:&mut VelVar1Dir,
*/

pub fn mass_vel_var_sq_2_internal<Num:RealField+Copy>(hlist_pat![mass,vel_var_sq]:HList!(&Mass<Num>,&VelVarSq<Num>))->HList!(Internal<Num>){
	hlist![Internal(vel_var_sq.0 * Num::frac_1_2() * mass.0)]
}

pub fn calculate_vel_var<Num:RealField+Copy,const DIM:usize>(hlist_pat![internal,mass]:HList!(&Internal<Num>,&Mass<Num>))->HList!(VelVarSq<Num>,VelVar<Num>,VelVarSq1Dir<Num>,VelVar1Dir<Num>){
	let vel_var_sq=VelVarSq(Num::p2()*internal.0/mass.0);
    let vel_var=VelVar(
		if vel_var_sq.0.is_negative(){
			-Num::sqrt(-vel_var_sq.0)
		}else{
			Num::sqrt(vel_var_sq.0)
		}
	);
    let vel_var_sq_1dir=VelVarSq1Dir(vel_var_sq.0/Num::from_isize(DIM as isize).unwrap());
    let vel_var_1dir=VelVar1Dir(vel_var.0/(Num::from_isize(DIM as isize).unwrap().sqrt()));
	hlist![vel_var_sq,vel_var,vel_var_sq_1dir,vel_var_1dir]
}

pub fn calculate_density<Num:RealField+Copy>(hlist_pat![mass,volume]:HList!(&Mass<Num>,&Volume<Num>))
->HList!(Density<Num>){
	hlist![Density(mass.0/volume.0)]
}

// pub fn calculate_internal<Num:RealField+Copy>(total_energy:Energy<Num>,other_energy:impl IntoIterator<Item :>)

pub fn calculate_matters_state<Num:RealField+Copy,const DIM:usize>(matters:MattersBasic<Num,DIM>)->MattersFull<Num,DIM>{
	let (mass,momentum,energy,volume)=matters.into();

	let hlist_pat![vel,kinetic] = calculate_position_state(hlist![&mass,&momentum]);
    let internal=Internal(energy.0-kinetic.0);

	let hlist_pat![vel_var_sq,vel_var,vel_var_sq_1dir,vel_var_1dir] = calculate_vel_var::<_,DIM>(hlist![&internal,&mass]);
	
	let hlist_pat![density] = calculate_density(hlist![&mass,&volume]);
	hlist![mass,momentum,energy,vel,kinetic,internal,vel_var_sq,vel_var,vel_var_sq_1dir,vel_var_1dir,volume,density]
}


pub fn calculate_matters_state_inplace_m<Num:RealField+Copy,const DIM:usize>(
	matters:
		// HList!(
		// &mut Mass<Num>,
		// &mut Momentum<Num,DIM>,
		// &mut Energy<Num>,
		// &mut Vel<Num,DIM>,
		// &mut Kinetic<Num>,
		// &mut Internal<Num>,
		// &mut VelVarSq<Num>,
		// &mut VelVar<Num>,
		// &mut VelVarSq1Dir<Num>,
		// &mut VelVar1Dir<Num>)
		HToMut<MattersFull<Num,DIM>>
){
		
    let (mass,momentum,energy,volume,vel,kinetic,internal,vel_var_sq,vel_var,vel_var_sq_1dir,vel_var_1dir,density)=matters.into();
    calculate_matters_state_inplace(hlist!(mass,momentum,energy,volume,vel,kinetic,internal,vel_var_sq,vel_var,vel_var_sq_1dir,vel_var_1dir,density));
	//vel.0=momentum.0/mass.0;
    //kinetic.0=  mass_momentum_2_kenetic(momentum.0,mass.0);//(momentum.0*momentum.0/mass.0) *NUMINV2;
    //internal.0=energy.0-kinetic.0;
    //vel_var_sq.0=2*internal.0/mass.0;
    //vel_var.0=Num::sqrt(vel_var_sq.0);
    //vel_var_sq_1dir.0=vel_var_sq.0/(DIM as i64);
    //vel_var_1dir.0=vel_var.0/(SQRTS[DIM]);
    //calculate_matters_state::<DIM>(matters.into())
}



pub fn calculate_matters_state_inplace<Num:RealField+Copy,const DIM:usize>(matters:HList!(
	&Mass<Num>,
    &Momentum<Num,DIM>,
    &Energy<Num>,
    &mut Vel<Num,DIM>,
    &mut Kinetic<Num>,
    &mut Internal<Num>,
    &mut VelVarSq<Num>,
    &mut VelVar<Num>,
    &mut VelVarSq1Dir<Num>,
    &mut VelVar1Dir<Num>,
	&Volume<Num>,
	&mut Density<Num>))
{
    let (mass,momentum,energy,vel,kinetic,internal,vel_var_sq,vel_var,vel_var_sq_1dir,vel_var_1dir,volume,density)=matters.into();
    vel.0=momentum.0/mass.0;
    kinetic.0=  mass_momentum_2_kenetic(momentum.0,mass.0);//(momentum.0*momentum.0/mass.0) *NUMINV2;
    internal.0=energy.0-kinetic.0;
    vel_var_sq.0=Num::p2()*internal.0/mass.0;
    vel_var.0=if vel_var_sq.0.is_negative(){
		-Num::sqrt(-vel_var_sq.0)
	}else{
		Num::sqrt(vel_var_sq.0)
	};
    vel_var_sq_1dir.0=vel_var_sq.0/Num::from_isize(DIM as isize).unwrap();
    vel_var_1dir.0=vel_var.0/(Num::sqrt(Num::from_isize(DIM as isize).unwrap()));
	density.0=mass.0/volume.0;
    //calculate_matters_state::<DIM>(matters.into())
}

#[cfg(test)]
mod tests{
	use num_traits::{One, zero};

	use super::*;

	type Num=simba::scalar::FixedI32F32;

	const DIM:usize=2;
	#[test]
	pub fn bug_test(){
		
		let one=Num::one();
		let f_1_2=Num::from_num(0.5);
		let input=hlist!(Mass(one),Momentum([f_1_2*f_1_2*f_1_2*f_1_2,zero()].into()),Energy(zero()),Volume(one));
		let output=calculate_matters_state(input);
		println!("{:?}",output);
	}

}