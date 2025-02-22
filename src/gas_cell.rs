use crate::{delta::State, matters::{Matters, MattersState}, num::Num};

#[derive(Debug,Clone,Default)]
pub struct GasCell{
    pub matters:State<MattersState>,
    pub edge:Num,
}

