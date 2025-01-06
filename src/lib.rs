pub mod matters;
pub mod num;
pub mod vec2_fix;
pub mod gas_cell;
pub mod body;
pub mod formulas;
pub mod constants;
pub mod delta;

pub fn add(left: u64, right: u64) -> u64 {
    left + right
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }
}
