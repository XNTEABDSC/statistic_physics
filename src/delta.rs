use std::{mem, ops::{AddAssign, Deref}};

#[derive(Debug,Clone,Default)]
pub struct Delta<T>(pub T);
#[derive(Default,Debug,Clone)]
pub struct State<T>(T); 
#[derive(Default,Debug)]
pub struct Change<T>(T);

impl<T> State<T> 
{
    pub const fn spawn(v:T)->Self {
        Self(v)
    }
    pub const fn v(&self)->&T {
        &self.0
    }
    pub fn apply_change<'a,T2>(&mut self,change:&'a mut Change<T2>)
        
        where T:for<'b >std::ops::AddAssign<&'b T2>,
        T2:Default
    {
        self.0 += &mem::replace::<T2>(&mut change.0, T2::default());
    }
}
impl<T> Change<T> 
    where T:Default+for<'b >std::ops::AddAssign<&'b T>+for<'b >std::ops::SubAssign<&'b T>
{
    pub const fn spawn(v:T)->Self {
        Self(v)
    }
    pub fn merge_change(&mut self,b:Change<T>) {
        self.0+=&b.0;
    }
    pub fn merge_delta(&mut self,b:Delta<T>){
        self.0+=&b.0;

    }
    pub const fn v(&self)->&T {
        &self.0
    }
}
impl<T> Delta<T> 
    where T:Default+for<'b >std::ops::AddAssign<&'b T>+for<'b >std::ops::SubAssign<&'b T>
{
    pub fn transfer(self,from:&mut Change<T>,to:& mut Change<T>) {
        from.0-=&self.0;
        to.0+=&self.0;
        
    }
}
impl<T> Deref for State<T> {
    type Target=T;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
impl<T> Deref for Change<T> {
    type Target=T;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}