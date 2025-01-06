use crate::{delta::State, matters::Matters, num::Num};
#[derive(Debug)]
pub struct CircleObject{
    pub matters:State<Matters>,
    pub radius:Num,
}