use crate::{delta::State, matters::Matters, num::Num};

#[derive(Debug,Clone,Default)]
pub struct GasCell{
    pub matters:State<Matters>,
    pub edge:Num,
}

