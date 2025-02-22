use crate::{delta::State, matters::MattersState, num::Num};
#[derive(Debug)]
pub struct CircleObject{
    pub matters:State<MattersState>,
    pub radius:Num,
}