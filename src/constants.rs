use std::array;

use crate::{num::Num, vec2_fix::Vec2Fix};

#[test]
fn getbits(){
    
}

pub const NUMINV2:Num=Num::unwrapped_shr(Num::ONE, 1);

pub const NUM2:Num=Num::unwrapped_shl(Num::ONE, 1);
pub const NUM4:Num=Num::unwrapped_shl(Num::ONE, 2);
pub const NUM8:Num=Num::unwrapped_shl(Num::ONE, 3);
pub const NUM16:Num=Num::unwrapped_shl(Num::ONE, 4);
pub const NUM32:Num=Num::unwrapped_shl(Num::ONE, 5);
/// 32=2^5
pub const FRAME_PER_SECOND:Num=mut_frame_per_second(Num::ONE);

pub const HALF_LIFE_PERIOD_SECOND:Num=NUM16;
/// 64*32=2048=2^11
/// For any 2 value to go closer
pub const HALF_LIFE_PERIOD_FRAME:Num=Num::unwrapped_mul(HALF_LIFE_PERIOD_SECOND, FRAME_PER_SECOND);

pub const HALF_LIFE_PERIOD_FRAME_LOG2:u32=11;
/// 1 - (1/2)^(1/2048)~=1 - (1+1/2048*ln(1/2)) = ln(2) / 2048
/// mult this to delta per frame to make value drain
pub const HALF_LIFE_PERIOD_FRAME_FACTOR:Num=Num::unwrapped_div(Num::LN_2, HALF_LIFE_PERIOD_FRAME);

/// HALF_LIFE_PERIOD_FRAME_FACTOR/2
/// mult this to 2 value delta per frame to make 2 value approach eachother
pub const HALF_LIFE_PERIOD_FRAME_FACTOR_OVER2:Num=Num::unwrapped_shr(HALF_LIFE_PERIOD_FRAME_FACTOR, 1);

pub const SECOND_PER_FRAME:Num=mut_second_per_frame(Num::ONE);

pub const INV_HALF_LIFE_PERIOD:Num=Num::unwrapped_div(Num::ONE,HALF_LIFE_PERIOD_SECOND);
//length of grid
pub const STD_LENGTH:Num=Num::ONE;

pub const fn mut_frame_per_second(n:Num)->Num{
    Num::unwrapped_shl(n,5)
}
pub const fn mut_second_per_frame(n:Num)->Num{
    Num::unwrapped_shr(n,5)
}

pub const fn mut_second_per_frame_vec(v:Vec2Fix)->Vec2Fix{
    Vec2Fix::new(mut_second_per_frame(v.0[0]), mut_second_per_frame(v.0[1]))
}